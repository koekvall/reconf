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
get_precomp <- function(Y, X, Z, REML = TRUE, Hlist = NULL) {
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
                  "R" = R)
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
