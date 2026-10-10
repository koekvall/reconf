# Solve A %*% x = b for symmetric A. Tries Cholesky first (fast); falls back
# to eigendecomposition-based pseudoinverse when A is not positive definite.
.solve_sym_eigen <- function(A, b, tol = 1e-10) {
  ch <- tryCatch(chol(A), error = function(e) NULL)
  if (!is.null(ch)) return(backsolve(ch, backsolve(ch, b, transpose = TRUE)))
  ed <- eigen(A, symmetric = TRUE)
  threshold <- tol * max(abs(ed$values))
  inv_vals <- ifelse(abs(ed$values) > threshold, 1 / ed$values, 0)
  ed$vectors %*% (inv_vals * (t(ed$vectors) %*% b))
}

# Resolve parameters given as names or as indices into the covariance
# parameters to indices; NULL stays NULL. nms are the parameter names and
# arg names the argument in error messages.
.parm_index <- function(parm, nms, arg) {
  if (is.null(parm)) return(NULL)
  if (is.character(parm)) {
    idx <- match(parm, nms)
    if (anyNA(idx)) {
      stop("unknown ", arg, ": ", paste(parm[is.na(idx)], collapse = ", "),
           ". The parameters are ", paste(nms, collapse = ", "))
    }
    parm <- idx
  }
  .check_idx(parm, arg, length(nms), unique = TRUE)
  as.integer(parm)
}

# Feasibility check: Sigma = psi_r (I_n + Z Psi_r Z') is positive definite
# iff I + F Psi_r F' is positive definite for any F with F'F = Z'Z, because
# the nonzero eigenvalues of F Psi_r F' and Psi_r Z'Z coincide. F is the
# precomputed q x q factor R with R'R = Z'Z when available and Z otherwise.
# Decided by attempting a sparse Cholesky factorization, which fails iff the
# matrix is not positive definite.
#
# When a gate cache is supplied (see get_precomp), the symbolic analysis of
# the factor is reused across evaluations and only the numeric factorization
# is redone. gate$pat0 is an all-zero matrix carrying the pattern of
# R Psi(1) R', which contains the pattern of R Psi_r R' for every psi;
# adding it pads the parent to the analyzed pattern.
.sigma_pd <- function(Psi_r, R = NULL, Z = NULL, gate = NULL) {
  if (!is.null(gate)) {
    Mx <- Matrix::forceSymmetric(Matrix::tcrossprod(R %*% Psi_r, R)) + gate$pat0
    return(!inherits(tryCatch(suppressWarnings(update(gate$ch, Mx, mult = 1)),
                              error = identity),
                     "condition"))
  }
  M <- if (!is.null(R)) Matrix::tcrossprod(R %*% Psi_r, R)
       else Matrix::tcrossprod(Z %*% Psi_r, Z)
  M <- Matrix::forceSymmetric(M + Matrix::Diagonal(nrow(M)))
  # LDL = FALSE forces a true Cholesky, which fails iff M is not positive
  # definite; the LDL' default would complete for indefinite M.
  !inherits(tryCatch(suppressWarnings(Matrix::Cholesky(M, perm = TRUE,
                                                       LDL = FALSE)),
                     error = identity),
            "condition")
}

# Quantities reusable across likelihood evaluations. H concatenates the
# structure matrices; R is a q x q factor with R'R = ZtZ used by the
# feasibility check, from .ztz_factor (NULL when no factor is available, in
# which case the check uses Z itself).
#
# With method = "n_side", instead precompute for the dense n-by-n likelihood
# path (loglik_n, loglik_res_n): the concatenated K = [K_1 ... K_{r - 1}] with
# K_j = Z H_j Z'. The K_j do not depend on psi, so every likelihood
# evaluation reuses them and the cost is independent of q.
#
# With method = "spectral" (requires r = 2, a single structure matrix),
# precompute for the O(n)-per-evaluation path (loglik_spectral,
# loglik_res_spectral): the eigendecomposition K = Z H_1 Z' = U diag(d) U'
# and the rotated data Yt = U'Y, Xt = U'X. Because Sigma = psi_1 K + psi_2 I
# shares the eigenvectors of K for every psi, the one-time O(n^3)
# decomposition replaces the per-evaluation factorization entirely.
#
# method = "auto" picks a dense path iff q >= n and Z is dense, and among
# the dense paths the spectral one when r = 2. The density condition matters:
# with sparse Z (e.g., crossed intercepts), Z'Z has O(n) off-diagonals
# however large q is, and the sparse q-side is as fast as or faster than the
# dense n-side even for q >> n (benchmarked 2026-07-12: at n = 1000,
# q = 6000 crossed, q-side ~4 ms vs n-side ~150 ms). With dense Z the q-side
# degenerates to dense q x q algebra and the n-side is 50-140x faster (same
# benchmark). The density threshold is a heuristic; method forces a path.
get_precomp <- function(Y, X, Z, Hlist = NULL,
                        method = c("auto", "q_side", "n_side", "spectral")) {
  method <- match.arg(method)
  if (method == "auto") {
    dens <- Matrix::nnzero(Z) / (as.double(nrow(Z)) * ncol(Z))
    if (ncol(Z) >= nrow(Z) && dens > 0.1) {
      method <- if (length(Hlist) == 1) "spectral" else "n_side"
    } else {
      method <- "q_side"
    }
  }
  if (method == "spectral") {
    if (!(is.list(Hlist) && length(Hlist) == 1)) {
      stop("method = 'spectral' requires exactly one structure matrix (r = 2)",
           call. = FALSE)
    }
    K <- as.matrix(Matrix::tcrossprod(Z %*% Hlist[[1]], Z))
    ed <- eigen(K, symmetric = TRUE)
    return(list("d" = ed$values,
                "Yt" = as.vector(crossprod(ed$vectors, Y)),
                "Xt" = as.matrix(crossprod(ed$vectors, X)),
                "method" = "spectral"))
  }
  if (method == "n_side") {
    if (!is.list(Hlist)) {
      stop("Hlist is required when method = 'n_side'", call. = FALSE)
    }
    n <- nrow(Z)
    K <- do.call(cbind, lapply(Hlist, function(Hj) {
      as.matrix(Matrix::tcrossprod(Z %*% Hj, Z))
    }))
    if (is.null(K)) K <- matrix(0, n, 0) # r = 1: only the error variance
    return(list("K" = K, "method" = "n_side"))
  }
  ZtZ <- methods::as(crossprod(Z), "generalMatrix")
  R <- .ztz_factor(Z, ZtZ)
  precomp <- list("XtX" = as.matrix(crossprod(X)),
                  "XtZ" = as.matrix(crossprod(X, Z)),
                  "ZtZ" = ZtZ,
                  "XtY" = as.vector(crossprod(X, Y)),
                  "ZtY" = as.vector(crossprod(Z, Y)),
                  "R" = R,
                  "method" = "q_side")
  if (!is.null(Hlist)) {
    precomp$H <- methods::as(do.call(cbind, Hlist), "generalMatrix")
    if (!is.null(R)) {
      # Cache the symbolic Cholesky analysis for the feasibility check. The
      # pattern of R Psi_r R' is contained in that of R Psi(1) R', where
      # Psi(1) = sum of the H_j has a one in every structurally nonzero
      # position and is positive semidefinite for the block structures
      # produced by lme4, so the analysis of I + R Psi(1) R' covers every
      # psi. If the initial factorization fails, .sigma_pd falls back to
      # factorizing from scratch at each evaluation.
      M1 <- Matrix::forceSymmetric(Matrix::tcrossprod(R %*% Reduce(`+`, Hlist), R))
      ch <- tryCatch(suppressWarnings(
        Matrix::Cholesky(M1 + Matrix::Diagonal(nrow(M1)), perm = TRUE,
                         LDL = FALSE)),
        error = function(e) NULL)
      if (!is.null(ch)) precomp$gate <- list(ch = ch, pat0 = M1 * 0)
    }
  }
  precomp
}

# A q x q factor R with R'R = ZtZ for the feasibility check, or NULL if none
# is available: the sparse Cholesky of ZtZ or, if that fails because ZtZ is
# not positive definite to working precision, the R factor of a sparse QR of
# Z from .qr_factor. On a rank-deficient ZtZ, such as crossed random
# intercepts, the Cholesky usually completes with a pivot near zero; R'R =
# ZtZ then still holds to rounding, which is all the check needs.
.ztz_factor <- function(Z, ZtZ) {
  R <- tryCatch(suppressWarnings(Matrix::chol(Matrix::forceSymmetric(ZtZ))),
                error = function(e) NULL)
  if (is.null(R)) R <- .qr_factor(Z)
  R
}

# The R factor of a sparse QR of Z with the column permutation undone, so
# that crossprod(R) equals Z'Z; NULL if the factorization fails. Matrix's
# sparse QR requires nrow(Z) >= ncol(Z), so for a wider Z this is NULL.
.qr_factor <- function(Z) {
  tryCatch(suppressWarnings({
    qrZ <- Matrix::qr(Z)
    Matrix::qr.R(qrZ)[, order(qrZ@q + 1L), drop = FALSE]
  }), error = function(e) NULL)
}

# Where each covariance parameter appears in Psi = sum_j psi_j H_j, from the
# structure matrices alone: on_diag[i, j] is TRUE if psi_j is the i-th
# diagonal entry of Psi, and n_off[i, j] counts the off-diagonal entries of
# row i that equal psi_j. is_var (length r) marks the variances: the
# parameters on the diagonal of Psi, and the error variance.
.psi_structure <- function(Hlist) {
  if (length(Hlist) == 0) {
    return(list(on_diag = matrix(FALSE, 0, 0), n_off = matrix(0, 0, 0),
                is_var = TRUE))
  }
  q <- ncol(Hlist[[1]])
  on_diag <- matrix(vapply(Hlist, function(H) as.vector(Matrix::diag(H) != 0),
                           logical(q)), q)
  n_off <- matrix(vapply(Hlist, function(H) {
    nz <- H != 0
    as.vector(Matrix::rowSums(nz)) - as.vector(Matrix::diag(nz))
  }, numeric(q)), q)
  list(on_diag = on_diag, n_off = n_off, is_var = c(colSums(on_diag) > 0, TRUE))
}

# The arguments vc_ci and vc_test accept in ... and pass to the optimizer:
# the trust::trust settings a caller can change, and the maximize_loglik
# flag warn_nonconv. The other formals of trust and maximize_loglik are set
# by the functions themselves. Any other name is rejected by .check_dots,
# because the CI search takes an error inside a profile evaluation as an
# infeasible point.
.optimizer_args <- c("rinit", "rmax", "iterlim", "fterm", "mterm",
                     "warn_nonconv")

# Validate the list of arguments in ... against .optimizer_args
.check_dots <- function(dots) {
  if (length(dots) > 0 && (is.null(names(dots)) || any(names(dots) == ""))) {
    stop("arguments in ... must be named", call. = FALSE)
  }
  bad <- setdiff(names(dots), .optimizer_args)
  if (length(bad) > 0) {
    stop("unused argument(s): ", paste(bad, collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

# Validate the model arguments of make_loglik. Conditions
# are joined by && so later ones are evaluated only when earlier ones hold.
.check_model <- function(Y, X, Z, Hlist, precomp = NULL) {
  if (!(is.vector(Y, mode = "numeric") && length(Y) > 0)) {
    stop("Y should be a numeric vector of positive length", call. = FALSE)
  }
  if (!(is.matrix(X) && nrow(X) == length(Y))) {
    stop("X should be a matrix with nrow(X) == length(Y)", call. = FALSE)
  }
  if (!(is(Z, "sparseMatrix") && nrow(Z) == length(Y) && ncol(Z) > 0)) {
    stop("Z should be a sparse matrix with nrow(Z) == length(Y) and ",
         "ncol(Z) > 0", call. = FALSE)
  }
  if (!(is.list(Hlist) &&
        all(vapply(Hlist, methods::is, logical(1), "sparseMatrix")) &&
        all(vapply(Hlist, function(H) all(dim(H) == ncol(Z)), logical(1))))) {
    stop("Hlist should be a list of q x q sparse matrices, q = ncol(Z)",
         call. = FALSE)
  }
  if (!(is.null(precomp) || is.list(precomp))) {
    stop("precomp should be NULL or a list", call. = FALSE)
  }
  invisible(TRUE)
}

# Validate a vector of r covariance parameters
.check_psi <- function(psi, r) {
  if (!(is.vector(psi, mode = "numeric") && length(psi) == r &&
        all(is.finite(psi)))) {
    stop("psi should be a numeric vector of finite values with length r = ",
         r, call. = FALSE)
  }
  invisible(TRUE)
}

# Validate a named list of single logical arguments
.check_flags <- function(flags) {
  for (nm in names(flags)) {
    if (!(is.logical(flags[[nm]]) && length(flags[[nm]]) == 1 &&
          !is.na(flags[[nm]]))) {
      stop(nm, " should be a single logical value", call. = FALSE)
    }
  }
  invisible(TRUE)
}

# Validate a vector of indices into a parameter vector of length n_par
.check_idx <- function(idx, idx_name, n_par, unique = FALSE) {
  if (!(is.vector(idx, mode = "numeric") && length(idx) > 0 &&
        !anyNA(idx) && all(idx == floor(idx)) && all(idx > 0))) {
    stop(idx_name, " should be a vector of positive integers", call. = FALSE)
  }
  if (max(idx) > n_par) {
    stop(idx_name, " values must not exceed ", n_par, call. = FALSE)
  }
  if (unique && anyDuplicated(idx) > 0) {
    stop(idx_name, " should not contain duplicate values", call. = FALSE)
  }
  invisible(TRUE)
}
