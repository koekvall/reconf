#' Get structure matrices for covariance parameterization
#'
#' Extracts the list of structure matrices (H matrices) from an lme4 fit that
#' determine how covariance parameters map to the covariance matrix structure.
#' These matrices are used in likelihood computations where the covariance matrix
#' is expressed as a linear combination: Psi = sum(psi\[i\] * H\[\[i\]\]).
#'
#' @param lmerfit An `lmerMod` object from fitting a linear mixed model using
#'   `lme4::lmer`.
#'
#' @return A list of sparse symmetric matrices (dsCMatrix), one for each
#'   covariance parameter (excluding error variance). The length of the list
#'   equals `getME(lmerfit, "m")`, which is r - 1 where r is the total number
#'   of covariance parameters including error variance.
#'
#' @details
#' Each matrix in the returned list is an indicator matrix showing which elements
#' of the random effects covariance matrix are associated with each parameter.
#' The i-th matrix has 1s in positions determined by the i-th covariance parameter
#' and 0s elsewhere.
#'
#' @noRd
get_Hlist_lmer <- function(lmerfit)
{
  .check_lmerfit(lmerfit)

  # Psi, and hence H, has the same structure as Lambdat
  H <- lme4::getME(lmerfit, "Lambdat")
  param_idx <- lme4::getME(lmerfit, "Lind")
  q <- nrow(H)
  m <- lme4::getME(lmerfit, "m")

  # Build each indicator matrix directly from matching nonzero positions,
  # avoiding r-1 full copies of H followed by drop0
  # Expand column-pointer format to per-element column indices
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

#' Extract estimated covariance parameters from lme4 fit
#'
#' Extracts all estimated variance and covariance parameters from a fitted
#' linear mixed model, including the error variance.
#'
#' @param lmerfit An `lmerMod` object from fitting a linear mixed model using
#'   `lme4::lmer`.
#'
#' @return A numeric vector containing all estimated covariance parameters,
#'   ordered as in the `vcov` column of 
#'   `as.data.frame(VarCorr(lmerfit), order = "lower.tri")`. The last element
#'   is the error variance. The vector has length r, where r is the total
#'   number of covariance parameters.
#'
#' @noRd
get_psi_hat_lmer <- function(lmerfit)
{
  .check_lmerfit(lmerfit)
  as.data.frame(lme4::VarCorr(lmerfit), order = "lower.tri")$vcov
}
