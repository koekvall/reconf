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
#' @keywords internal
get_Hlist_lmer <- function(lmerfit)
{
  # Validate input
  if (!inherits(lmerfit, "lmerMod")) {
    stop("lmerfit must be an lmerMod object from lme4::lmer")
  }
  
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
# exactly invariant under the transformation; only the log-likelihood value
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

#' Get precomputed quantities from lme4 fit
#'
#' Extracts and computes quantities from an lme4 fit that can be reused in
#' likelihood calculations to avoid redundant computations. The quantities
#' computed depend on whether REML or ML estimation is used.
#'
#' @param lmerfit An `lmerMod` object from fitting a linear mixed model using
#'   `lme4::lmer`.
#' @param REML Logical indicating whether to compute quantities for REML
#'   (\code{TRUE}) or ML (\code{FALSE}). If \code{NULL} (default), uses the
#'   estimation method from the fitted model.
#' @param Hlist Optional list of structure matrices (see
#'   \code{\link{get_Hlist_lmer}}); when supplied, their concatenation is
#'   stored in the result so likelihood evaluations do not rebuild it.
#' @param method Computational path to precompute for; see
#'   \code{?loglikelihood}. The default \code{"auto"} picks the dense n-by-n
#'   path iff \eqn{q \ge n} and \eqn{Z} is dense.
#'
#' @return A list containing precomputed cross-products and related quantities:
#'   \code{XtX}, \code{XtZ}, \code{ZtZ}, \code{R} (sparse Cholesky factor of
#'   \code{ZtZ}, or \code{NULL} if singular), for REML also \code{XtY} and
#'   \code{ZtY}, and \code{H} if \code{Hlist} was supplied. For the n-side
#'   path, instead the dense matrix \code{K} (see \code{?loglik_n}). Either
#'   way the list is tagged with \code{method}.
#'
#' @keywords internal
get_precomp_lmer <- function(lmerfit, REML = NULL, Hlist = NULL,
                             method = c("auto", "q_side", "n_side")){
  # Validate input
  if (!inherits(lmerfit, "lmerMod")) {
    stop("lmerfit must be an lmerMod object from lme4::lmer")
  }

  if(is.null(REML)){
    # 0 indicates ML, non-zero indicates REML
    REML <- lme4::getME(lmerfit, "REML") != 0
  } else {
    if (!is.logical(REML) || length(REML) != 1) {
      stop("REML must be a single logical value")
    }
  }

  m <- .lmer_matrices(lmerfit)

  get_precomp(Y = m$Y, X = m$X, Z = m$Z, REML = REML, Hlist = Hlist,
              method = match.arg(method))
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
#' @details
#' This function extracts the full vector of covariance parameter estimates
#' (often denoted psi or theta in the package), including:
#' \itemize{
#'   \item Variances and covariances of random effects
#'   \item Error variance (last element)
#' }
#' 
#' The ordering follows lme4's internal parameterization with "lower.tri" ordering.
#'
#' @keywords internal
get_psi_hat_lmer <- function(lmerfit)
{
  # Validate input
  if (!inherits(lmerfit, "lmerMod")) {
    stop("lmerfit must be an lmerMod object from lme4::lmer")
  }
  
  # Extract variance components
  vcov_vec <- as.data.frame(lme4::VarCorr(lmerfit), order = "lower.tri")$vcov
  
  if (!is.numeric(vcov_vec) || length(vcov_vec) == 0) {
    stop("Failed to extract variance components from lmerfit")
  }
  
  vcov_vec
}

#' Score test for linear mixed model fitted with lme4
#'
#' Computes a score test statistic for a linear mixed model fitted using lme4::lmer.
#' This is a convenience wrapper around \code{\link{score_stat}} that extracts the
#' necessary components from an lmerMod object.
#'
#' @param lmerfit An `lmerMod` object from fitting a linear mixed model using
#'   `lme4::lmer`.
#' @param theta_null Numeric vector of parameter values under the null hypothesis.
#'   For REML fits, this should be a vector of length r (covariance parameters only).
#'   For ML fits, this should be a vector of length p + r (fixed effects followed by
#'   covariance parameters). If \code{NULL}, defaults to zero random effects and
#'   unit error variance.
#' @param test_idx Integer vector specifying which elements of \code{theta_null} to
#'   test. If \code{NULL}, tests all covariance parameters except error variance
#'   (i.e., tests for zero random effects).
#' @param efficient Logical. If \code{TRUE} (default), use efficient information
#'   that accounts for estimation of nuisance parameters.
#' @param expected Logical. If \code{TRUE} (default), use expected Fisher information;
#'   otherwise use observed information.
#' @param profile Logical. If \code{TRUE} (default), optimize nuisance parameters
#'   under the null hypothesis before computing the test statistic.
#' @param known_idx Integer vector or \code{NULL} specifying which elements of
#'   \code{theta_null} (other than \code{test_idx}) have known values and should
#'   not be optimized when \code{profile = TRUE}. If \code{NULL}, all parameters
#'   except \code{test_idx} are treated as nuisance parameters.
#' @param ... Additional arguments passed to the trust-region optimizer used
#'   when \code{profile = TRUE}, such as \code{iterlim} (default 1000 here:
#'   profiling at a distant null can involve many cheap iterations along
#'   nearly flat directions).
#'
#' @return A named numeric vector with elements:
#'   \item{stat}{The score test statistic}
#'   \item{p_val}{P-value from chi-squared distribution}
#'   \item{df}{Degrees of freedom (length of \code{test_idx})}
#'
#' @details
#' This function provides a convenient interface for score testing in lme4 fits.
#' It automatically extracts the model matrices, variance structure, and other
#' components needed by \code{\link{score_stat}}.
#'
#' When \code{profile = TRUE}, the function uses \code{\link{maximize_loglik}} to
#' optimize nuisance parameters under the null hypothesis before computing the test
#' statistic. This typically yields more powerful tests.
#'
#' @examples
#' library(lme4)
#' fit <- lmer(Reaction ~ Days + (1 | Subject), data = sleepstudy)
#'
#' # Default: test H0 that the random intercept variance is zero
#' score_test_lmer(fit)
#'
#' # Test only the random intercept variance (index 1), against null value 0
#' score_test_lmer(fit, test_idx = 1L)
#'
#' @export
score_test_lmer <- function(lmerfit,
                            theta_null = NULL,
                            test_idx = NULL,
                            efficient = TRUE,
                            expected = TRUE,
                            profile = TRUE,
                            known_idx = NULL,
                            ...)
{
  # Validate input
  if (!inherits(lmerfit, "lmerMod")) {
    stop("lmerfit must be an lmerMod object from lme4::lmer")
  }

  # Extract model components (offset and prior weights applied)
  m <- .lmer_matrices(lmerfit)
  Y <- m$Y
  X <- m$X
  Z <- m$Z
  Hlist <- get_Hlist_lmer(lmerfit)
  REML <- lme4::getME(lmerfit, "REML") != 0
  precomp <- get_precomp_lmer(lmerfit, REML = REML, Hlist = Hlist)
  
  p <- ncol(X)
  r <- length(Hlist) + 1
  
  # Set up theta_null
  if(is.null(theta_null)){
    # Default: zero random effects, unit error variance
    if(REML){
      theta_null <- c(rep(0, r - 1), 1)
    } else {
      theta_null <- c(lme4::fixef(lmerfit), rep(0, r - 1), 1)
    }
  }
  
  # Validate theta_null
  expected_length <- if(REML) r else p + r
  if(length(theta_null) != expected_length){
    stop("theta_null should have length ", expected_length,
         " (", if(REML) "r" else "p + r", " for REML = ", REML, ")")
  }
  
  # Check error variance is positive
  psi_r_idx <- if(REML) r else p + r
  if(theta_null[psi_r_idx] <= 0){
    stop("Error variance (last element of theta_null) must be positive")
  }
  
  # Set up test_idx
  if(is.null(test_idx)){
    # Default: test all random effect parameters (not error variance)
    if(REML){
      test_idx <- seq_len(r - 1)
    } else {
      test_idx <- (p + 1):(p + r - 1)
    }
  }
  
  k <- length(test_idx)
  if(k == 0){
    stop("test_idx must have length > 0")
  }
  
  # Profile nuisance parameters if requested
  if(profile){
    # Determine which parameters to optimize (exclude test_idx and known_idx)
    exclude_idx <- c(test_idx, known_idx)
    opt_idx <- seq_along(theta_null)[-exclude_idx]
    
    if(length(opt_idx) > 0){
      # Optimize nuisance parameters. Profiling at a null far from the
      # estimates can require many cheap iterations along nearly flat
      # directions (Fisher scoring converges linearly there), so the
      # default iteration limit is generous.
      dots <- list(...)
      if (is.null(dots$iterlim)) dots$iterlim <- 1000L
      theta_null <- do.call(maximize_loglik,
                            c(list(start_val = theta_null,
                                   opt_idx = opt_idx,
                                   Y = Y,
                                   X = X,
                                   Z = Z,
                                   Hlist = Hlist,
                                   expected = expected,
                                   REML = REML,
                                   precomp = precomp),
                              dots))$arg
    }
  }
  
  # Compute score test statistic
  test_stat <- score_stat(theta = theta_null,
                          test_idx = test_idx,
                          Y = Y,
                          X = X,
                          Z = Z,
                          Hlist = Hlist,
                          REML = REML,
                          expected = expected,
                          efficient = efficient,
                          signed = FALSE,
                          known_idx = known_idx,
                          precomp = precomp)
  
  # Return results
  c("stat" = as.numeric(test_stat),
    "p_val" = stats::pchisq(as.numeric(test_stat), df = k, lower.tail = FALSE),
    "df" = k)
}

#' Score tests for all covariance parameters individually
#'
#' Performs individual score tests for each covariance parameter in a linear mixed
#' model fitted with lme4. This function tests each parameter separately while
#' profiling over all other parameters. The default null value is zero for
#' every parameter, so variance rows test whether the corresponding random
#' effect is needed and covariance rows test zero covariance.
#'
#' @param lmerfit An `lmerMod` object from fitting a linear mixed model using
#'   `lme4::lmer`.
#' @param theta_null Numeric vector of parameter values under the null hypothesis.
#'   For REML fits, this should be a vector of length r (covariance parameters only).
#'   For ML fits, this should be a vector of length p + r (fixed effects followed by
#'   covariance parameters). If \code{NULL} (default), each parameter is tested
#'   at zero, with the remaining parameters started at their estimates when
#'   profiling.
#' @param test_idx Integer vector specifying which covariance parameters to test.
#'   If \code{NULL}, tests all covariance parameters except error variance.
#' @param efficient Logical. If \code{TRUE} (default), use efficient information
#'   that accounts for estimation of nuisance parameters.
#' @param expected Logical. If \code{TRUE} (default), use expected Fisher information;
#'   otherwise use observed information.
#' @param ... Additional arguments passed to \code{\link{score_test_lmer}}.
#'
#' @return A matrix of class \code{reconf_test} with one row per tested
#'   parameter, named as in \code{\link{ci_all_lmer}}, and columns
#'   \code{statistic} (the score test statistic), \code{df} (degrees of
#'   freedom, 1 for individual tests), and \code{p.value} (from the
#'   chi-squared distribution). Supports \code{print} and
#'   \code{\link[generics]{tidy}}.
#'
#' @details
#' For each parameter specified in \code{test_idx}, this function:
#' \enumerate{
#'   \item Sets up a null hypothesis with that parameter at its null value
#'   \item Profiles over all other parameters to maximize the likelihood under the null
#'   \item Computes the score test statistic
#' }
#'
#' When testing covariance parameters (off-diagonal elements), the function ensures
#' that the starting values for optimization yield a positive semi-definite covariance
#' matrix by appropriately adjusting the corresponding variance parameters.
#'
#' @examples
#' \donttest{
#' library(lme4)
#' fit <- lmer(Reaction ~ Days + (Days | Subject), data = sleepstudy)
#'
#' # Test each variance/covariance parameter against zero
#' score_test_all_lmer(fit)
#' }
#' @seealso \code{\link{score_test_lmer}} for joint tests and user-supplied
#'   null values, \code{\link{ci_all_lmer}} for confidence intervals.
#' @export
score_test_all_lmer <- function(lmerfit,
                          theta_null = NULL,
                          test_idx = NULL,
                          efficient = TRUE,
                          expected = TRUE,
                          ...)
{
  # Validate input
  if (!inherits(lmerfit, "lmerMod")) {
    stop("lmerfit must be an lmerMod object from lme4::lmer")
  }
  
  # Get model dimensions
  REML <- lme4::getME(lmerfit, "REML") != 0
  r_i <- lme4::getME(lmerfit, "m_i")
  r <- sum(r_i) + 1
  p <- ncol(lme4::getME(lmerfit, "X"))
  
  # Set up theta_null. By default each parameter is tested at zero while the
  # nuisance parameters start from their estimates, which keeps the profiling
  # well conditioned regardless of the response scale.
  default_null <- is.null(theta_null)
  if(default_null){
    psi_hat <- get_psi_hat_lmer(lmerfit)
    theta_null <- if(REML) psi_hat else c(as.vector(lme4::fixef(lmerfit)), psi_hat)
  }

  # Validate theta_null length
  expected_length <- if(REML) r else p + r
  if(length(theta_null) != expected_length){
    stop("theta_null should have length ", expected_length,
         " (", if(REML) "r" else "p + r", " for REML = ", REML, ")")
  }

  # Set up test_idx (indices in theta_null space)
  if(is.null(test_idx)){
    # Default: test all covariance parameters except error variance
    if(REML){
      test_idx <- seq_len(r - 1)
    } else {
      test_idx <- (p + 1):(p + r - 1)
    }
  }
  
  k <- length(test_idx)
  if(k == 0){
    stop("test_idx must have length > 0")
  }
  
  # Get variance/covariance indicator and parameter names
  vc <- as.data.frame(lme4::VarCorr(lmerfit), order = "lower.tri")
  is_var_param <- is.na(vc$var2)

  # Used to determine which term a parameter belongs to
  last_idx <- lme4::getME(lmerfit, "Tp")

  # Prepare output matrix; test_idx indexes theta, psi_idx the psi vector
  psi_idx <- if (REML) test_idx else test_idx - p
  out <- matrix(NA, nrow = k, ncol = 3,
                dimnames = list(.param_names(vc, psi_idx),
                                c("statistic", "df", "p.value")))
  
  # Loop over parameters to test
  for(i in seq_along(test_idx)){
    param_idx <- test_idx[i]
    
    # Adjust for REML vs ML indexing
    psi_param_idx <- if(REML) param_idx else param_idx - p
    
    # Create starting point for this test; under the default, the tested
    # parameter's null value is zero
    theta_start <- theta_null
    if (default_null) theta_start[param_idx] <- 0

    if (psi_param_idx <= (r - 1)) {
      # Term bookkeeping. Tp has a leading zero: term t covers parameters
      # (Tp[t] + 1):Tp[t + 1].
      term_idx <- max(which(last_idx < psi_param_idx))
      first_param <- last_idx[term_idx] + 1L
      dim_i <- as.integer(0.5 * (-1 + sqrt(1 + 8 * r_i[term_idx])))
      jj <- psi_param_idx - last_idx[term_idx]
      row_col <- get_row_col_ltri(jj, n = dim_i)
      shift <- if (REML) 0L else p

      if (!is_var_param[psi_param_idx]) {
        # Tested parameter is a covariance: ensure the two corresponding
        # variances give a positive semidefinite starting value
        var_idx1 <- get_idx_ltri(row = row_col[2], col = row_col[2], n = dim_i)
        var_idx2 <- get_idx_ltri(row = row_col[1], col = row_col[1], n = dim_i)
        adj <- shift + first_param + c(var_idx1, var_idx2) - 1
        theta_start[adj] <- pmax(theta_start[adj], abs(theta_start[param_idx]))
      } else if (default_null && dim_i > 1) {
        # Tested parameter is a variance set to zero: also start the
        # covariances involving the same variable at zero. They are still
        # profiled; the modified start avoids slow optimization along the
        # weakly identified covariance direction at zero variance.
        v <- row_col[1]
        cov_idx <- vapply(setdiff(seq_len(dim_i), v), function(i)
          get_idx_ltri(row = max(i, v), col = min(i, v), n = dim_i), numeric(1))
        theta_start[shift + first_param + cov_idx - 1] <- 0
      }
    }

    # Perform score test for this parameter
    res <- score_test_lmer(lmerfit = lmerfit,
                           theta_null = theta_start,
                           test_idx = param_idx,
                           efficient = efficient,
                           expected = expected,
                           profile = TRUE, ...)
    out[i, ] <- res[c("stat", "df", "p_val")]
  }
  structure(out, class = c("reconf_test", class(out)),
            method = if (REML) "REML" else "ML")
}



