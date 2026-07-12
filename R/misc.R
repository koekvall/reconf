# get the row and column index from the vectorization of the
# lower triangular part of an n x n matrix
get_row_col_ltri <- function(idx, n)
{
  last_idx <- as.integer(n * (n + 1) / 2)
  stopifnot(idx <= last_idx)
  # Count backwards from the end to our idx
  back_idx <- last_idx - idx + 1
  # Column counting from the right
  back_col <- ceiling(0.5 * (-1 + sqrt(1 + 8 * back_idx)))
  # Column counting from the left
  col <- n - back_col + 1

  # Get the row
  row <- 0.5 * back_col * (back_col + 1)
  row <- (n - back_col) + (row - back_idx) + 1
  c(row, col)
}

# Get the index in vectorization of lower triangular
get_idx_ltri <- function(row, col, n)
{
  stopifnot(row >= col && row <= n && col <= n)
  back_col <- n - col + 1
  # Indexing counting backwards from the end
  back_idx <- 0.5 * back_col * (back_col + 1) + col - row
  0.5 * n * (n + 1) - back_idx + 1
}

# Human-readable names for covariance parameters, from the data frame
# as.data.frame(VarCorr(fit), order = "lower.tri"). The residual row has
# var1 = var2 = NA and is named by grp alone; variance rows have var2 = NA;
# covariance rows have both variables.
.param_names <- function(vc, idx = seq_len(nrow(vc))) {
  vapply(idx, function(i) {
    if (is.na(vc$var1[i])) return(vc$grp[i])
    paste0(vc$var1[i],
           if (is.na(vc$var2[i])) "" else paste0(":", vc$var2[i]),
           " | ", vc$grp[i])
  }, character(1))
}

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

# Feasibility check: Sigma = psi_r (I_n + Z Psi_r Z') is positive definite
# iff I + F Psi_r F' is positive definite for any F with F'F = Z'Z, because
# the nonzero eigenvalues of F Psi_r F' and Psi_r Z'Z coincide. Uses the
# precomputed q x q factor R = chol(ZtZ) when available and falls back to
# F = Z (an n x n check) when ZtZ is singular. Decided by attempting a sparse
# Cholesky factorization, which fails iff the matrix is not positive definite.
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
# structure matrices; R is the Cholesky factor of ZtZ used by the feasibility
# check (NULL when ZtZ is singular, e.g., crossed random intercepts).
#
# With method = "n_side", instead precompute for the dense n-by-n likelihood
# path (loglik_n, loglik_res_n): the concatenated K = [K_1 ... K_{r - 1}] with
# K_j = Z H_j Z'. The K_j do not depend on psi, so every likelihood
# evaluation reuses them; this is what makes the n-side path independent of q.
#
# method = "auto" picks n_side iff q >= n and Z is dense-ish. Dimensions
# alone are not the right criterion: with sparse Z (e.g., crossed
# intercepts), Z'Z has O(n) off-diagonals however large q is, and the sparse
# q-side stays as fast or faster than the dense n-side even for q >> n
# (benchmarked 2026-07-12: at n = 1000, q = 6000 crossed, q-side ~4 ms vs
# n-side ~150 ms). The n-side wins when Z is dense -- kernel-, kinship-, or
# spline-type designs -- where the q-side degenerates to dense q x q algebra
# (same benchmark: 50-140x in favor of n-side). The density threshold is a
# heuristic; callers can always force a path via method.
get_precomp <- function(Y, X, Z, REML = TRUE, Hlist = NULL,
                        method = c("auto", "q_side", "n_side")) {
  method <- match.arg(method)
  if (method == "auto") {
    dens <- Matrix::nnzero(Z) / (as.double(nrow(Z)) * ncol(Z))
    method <- if (ncol(Z) >= nrow(Z) && dens > 0.1) "n_side" else "q_side"
  }
  if (method == "n_side") {
    assertthat::assert_that(is.list(Hlist),
                            msg = "Hlist is required when method = 'n_side'")
    n <- nrow(Z)
    K <- do.call(cbind, lapply(Hlist, function(Hj) {
      as.matrix(Matrix::tcrossprod(Z %*% Hj, Z))
    }))
    if (is.null(K)) K <- matrix(0, n, 0) # r = 1: only the error variance
    return(list("K" = K, "method" = "n_side"))
  }
  ZtZ <- methods::as(crossprod(Z), "generalMatrix")
  R <- tryCatch(suppressWarnings(Matrix::chol(Matrix::forceSymmetric(ZtZ))),
                error = function(e) NULL)
  if (is.null(R)) {
    # ZtZ singular (e.g., a group with a single observation, or crossed
    # random intercepts). Use the R factor of a sparse QR of Z instead:
    # after undoing the column permutation, crossprod(R) equals ZtZ, which
    # is all the feasibility check needs.
    R <- tryCatch(suppressWarnings({
      qrZ <- Matrix::qr(Z)
      Matrix::qr.R(qrZ)[, order(qrZ@q + 1L), drop = FALSE]
    }), error = function(e) NULL)
  }
  precomp <- list("XtX" = as.matrix(crossprod(X)),
                  "XtZ" = as.matrix(crossprod(X, Z)),
                  "ZtZ" = ZtZ,
                  "R" = R,
                  "method" = "q_side")
  if (REML) {
    precomp$XtY <- as.vector(crossprod(X, Y))
    precomp$ZtY <- as.vector(crossprod(Z, Y))
  }
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
