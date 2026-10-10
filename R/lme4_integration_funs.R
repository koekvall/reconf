#' Structure matrices of an lme4 fit
#'
#' The matrices \eqn{H_j} with \eqn{\Psi = \sum_j \psi_j H_j} as in
#' \code{make_loglik}: \eqn{H_j} has ones at the positions of \eqn{\psi_j}
#' in \eqn{\Psi} and zeros elsewhere.
#'
#' @param lmerfit An \code{lmerMod} object from \code{lme4::lmer}.
#'
#' @return A list of \eqn{r - 1} sparse symmetric matrices, one per
#'   random-effect covariance parameter, in the order of
#'   \code{as.data.frame(VarCorr(lmerfit), order = "lower.tri")}.
#'
#' @noRd
get_Hlist_lmer <- function(lmerfit)
{
  # Psi, and hence H, has the same structure as Lambdat
  H <- lme4::getME(lmerfit, "Lambdat")
  param_idx <- lme4::getME(lmerfit, "Lind")
  q <- nrow(H)
  m <- lme4::getME(lmerfit, "m")

  # Build each indicator matrix from the nonzero positions of Lambdat whose
  # Lind entry is i. Expand the column pointers to a column index per entry
  col_idx <- rep(seq_len(q), diff(H@p))
  row_idx <- H@i + 1L  # 0-based to 1-based

  lapply(seq_len(m), function(i) {
    keep <- which(param_idx == i)
    M <- Matrix::sparseMatrix(i = row_idx[keep], j = col_idx[keep],
                              x = 1, dims = c(q, q))
    Matrix::forceSymmetric(M, uplo = "U")
  })
}

# Extract Y, X, Z from an lmer fit with the offset and prior weights applied.
# Premultiplying by W^(1/2) turns the weighted model, Var(E) = psi_r W^{-1},
# into the unit-variance model the likelihood code implements, with identical
# parameters (beta, Psi, psi_r); the offset enters the mean only. All
# psi-dependent quantities (scores, information, statistics, intervals) are
# invariant under the transformation; only the log-likelihood value
# changes, by the constant 0.5 * sum(log(w)). The diagonal scaling leaves
# Z's sparsity pattern unchanged.
.lmer_matrices <- function(lmerfit) {
  Y <- lme4::getME(lmerfit, "y")
  X <- lme4::getME(lmerfit, "X")
  Z <- lme4::getME(lmerfit, "Z")
  off <- lme4::getME(lmerfit, "offset")
  if (length(off) > 0 && any(off != 0)) Y <- Y - off
  w <- stats::weights(lmerfit)
  if (length(w) > 0 && any(w != 1)) {
    if (any(w <= 0)) stop("prior weights must be positive")
    sw <- sqrt(w)
    Y <- sw * Y
    X <- sw * X
    Z <- Matrix::Diagonal(x = sw) %*% Z
  }
  list(Y = Y, X = X, Z = Z)
}

# The vc_model of an lme4 fit, with lme4's estimates and parameter names;
# the only place the inference functions depend on lme4.
#' @rdname vc_model
#' @export
vc_model.lmerMod <- function(x, method = c("auto", "q_side", "n_side",
                                           "spectral"), ...) {
  m  <- .lmer_matrices(x)  # offset and prior weights applied
  vc <- as.data.frame(lme4::VarCorr(x), order = "lower.tri")
  .make_setup(Y = m$Y, X = m$X, Z = m$Z, Hlist = get_Hlist_lmer(x),
              REML = lme4::getME(x, "REML") != 0,
              psi_hat = get_psi_hat_lmer(x), names = .param_names(vc),
              method = match.arg(method))
}

# Names for covariance parameters, from the data frame
# as.data.frame(VarCorr(fit), order = "lower.tri"). They follow lme4's names
# for standard deviations and correlations (sd_x|g, cor_y.x|g) on the
# variance scale: var_x|g, cov_y.x|g, and var_Residual for the error
# variance. The residual row has var1 = var2 = NA, variance rows have
# var2 = NA, and covariance rows have both variables.
.param_names <- function(vc, idx = seq_len(nrow(vc))) {
  vapply(idx, function(i) {
    if (is.na(vc$var1[i])) return(paste0("var_", vc$grp[i]))
    if (is.na(vc$var2[i])) return(paste0("var_", vc$var1[i], "|", vc$grp[i]))
    paste0("cov_", vc$var2[i], ".", vc$var1[i], "|", vc$grp[i])
  }, character(1))
}

#' Estimated covariance parameters of an lme4 fit
#'
#' @param lmerfit An \code{lmerMod} object from \code{lme4::lmer}.
#'
#' @return The \code{vcov} column of
#'   \code{as.data.frame(VarCorr(lmerfit), order = "lower.tri")}: the
#'   \eqn{r} covariance parameters, the error variance last.
#'
#' @noRd
get_psi_hat_lmer <- function(lmerfit)
{
  as.data.frame(lme4::VarCorr(lmerfit), order = "lower.tri")$vcov
}
