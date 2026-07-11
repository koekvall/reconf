#' Get the covariance matrix of random effects
#'
#' Returns the covariance matrix of the random effects vector in a linear mixed model,
#' either as estimated by lme4::lmer or evaluated at particular parameter values
#' supplied as an argument.
#'
#' @param lmerfit An `lmerMod` object from fitting a linear mixed model using `lme4::lmer`.
#' @param psi_mr Optional numeric vector with covariance parameters, not including the 
#'   error variance; i.e., `psi` minus its r:th element, where r is the total number
#'   of covariance parameters (see details).
#'
#' @return A sparse symmetric covariance matrix (dsCMatrix) of the random effects vector.
#'
#' @details
#' If `psi_mr` is not supplied, the estimated covariance matrix from the fitted model 
#' is returned. If `psi_mr` is supplied, the covariance matrix is constructed using 
#' those parameter values instead.
#' 
#' When provided, `psi_mr` should be a numeric vector of length `getME(lmerfit, "m")`,
#' which equals r - 1 where r is the total number of covariance parameters (including
#' error variance). The elements should be ordered as in the `vcov` column of
#' `as.data.frame(VarCorr(lmerfit), order = "lower.tri")`, with the last element 
#' (error variance) omitted.
#'
#' @keywords internal
get_Psi_lmer <- function(lmerfit, psi_mr = NULL){
  # Validate input
  if (!inherits(lmerfit, "lmerMod")) {
    stop("lmerfit must be an lmerMod object from lme4::lmer")
  }
  
  rm1 <- lme4::getME(lmerfit, "m")
  
  if(is.null(psi_mr)){
    # Extract variances and covariances of random effects ordered as in the covmat
    # lower.tri despite creating upper triangular Psi since it appears consistent
    # with the indexing in getME(, "Lind")
    psi_mr <- as.data.frame(lme4::VarCorr(lmerfit), order = "lower.tri")$vcov[1:rm1]
  } else {
    # Validate psi_mr
    if (!is.numeric(psi_mr)) {
      stop("psi_mr must be a numeric vector")
    }
    if (length(psi_mr) != rm1) {
      stop("psi_mr must have length ", rm1, " (getME(lmerfit, 'm'))")
    }
  }

  # Lambdat has the right structure, but not the same entries as Psi
  Psi_half <- lme4::getME(lmerfit, "Lambdat")

  # Fill in the upper triangular part of Psi with the extracted elements
  param_idx <- lme4::getME(lmerfit, "Lind")
  Psi_half@x <- psi_mr[param_idx]

  # Return symmetric matrix
  Matrix::forceSymmetric(Psi_half, uplo = "U")
}

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
#'
#' @return A list containing precomputed cross-products and related quantities:
#'   \code{XtX}, \code{XtZ}, \code{ZtZ}, \code{R} (sparse Cholesky factor of
#'   \code{ZtZ}, or \code{NULL} if singular), for REML also \code{XtY} and
#'   \code{ZtY}, and \code{H} if \code{Hlist} was supplied.
#'
#' @keywords internal
get_precomp_lmer <- function(lmerfit, REML = NULL, Hlist = NULL){
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

  X <- lme4::getME(lmerfit, "X")
  Z <- lme4::getME(lmerfit, "Z")
  Y <- lme4::getME(lmerfit, "y")

  get_precomp(Y = Y, X = X, Z = Z, REML = REML, Hlist = Hlist)
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
                            known_idx = NULL)
{
  # Validate input
  if (!inherits(lmerfit, "lmerMod")) {
    stop("lmerfit must be an lmerMod object from lme4::lmer")
  }

  # Extract model components
  Y <- lme4::getME(lmerfit, "y")
  X <- lme4::getME(lmerfit, "X")
  Z <- lme4::getME(lmerfit, "Z")
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
      # Optimize nuisance parameters
      theta_null <- maximize_loglik(start_val = theta_null,
                                    opt_idx = opt_idx,
                                    Y = Y,
                                    X = X,
                                    Z = Z,
                                    Hlist = Hlist,
                                    expected = expected,
                                    REML = REML,
                                    precomp = precomp)$arg
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

#' Test all covariance parameters individually
#'
#' Performs individual score tests for each covariance parameter in a linear mixed
#' model fitted with lme4. This function tests each parameter separately while
#' profiling over all other parameters.
#'
#' @param lmerfit An `lmerMod` object from fitting a linear mixed model using
#'   `lme4::lmer`.
#' @param theta_null Numeric vector of parameter values under the null hypothesis.
#'   For REML fits, this should be a vector of length r (covariance parameters only).
#'   For ML fits, this should be a vector of length p + r (fixed effects followed by
#'   covariance parameters). If \code{NULL}, defaults to zero random effects and
#'   unit error variance for each test.
#' @param test_idx Integer vector specifying which covariance parameters to test.
#'   If \code{NULL}, tests all covariance parameters except error variance.
#' @param efficient Logical. If \code{TRUE} (default), use efficient information
#'   that accounts for estimation of nuisance parameters.
#' @param expected Logical. If \code{TRUE} (default), use expected Fisher information;
#'   otherwise use observed information.
#'
#' @return A matrix with one row per tested parameter and three columns:
#'   \item{stat}{The score test statistic}
#'   \item{p_val}{P-value from chi-squared distribution with 1 degree of freedom}
#'   \item{df}{Degrees of freedom (always 1 for individual tests)}
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
#' @keywords internal
test_all_lmer <- function(lmerfit,
                          theta_null = NULL,
                          test_idx = NULL,
                          efficient = TRUE,
                          expected = TRUE)
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
  
  # Set up theta_null
  if(is.null(theta_null)){
    # Default: zero random effects, unit error variance
    if(REML){
      theta_null <- c(rep(0, r - 1), 1)
    } else {
      theta_null <- c(lme4::fixef(lmerfit), rep(0, r - 1), 1)
    }
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
  
  # Get variance/covariance indicator
  is_var_param <- is.na(as.data.frame(lme4::VarCorr(lmerfit), order = "lower.tri")$var2)
  
  # Used to determine which term a parameter belongs to
  last_idx <- lme4::getME(lmerfit, "Tp")
  
  # Prepare output matrix
  out <- matrix(NA, nrow = k, ncol = 3)
  colnames(out) <- c("stat", "p_val", "df")
  
  # Loop over parameters to test
  for(i in seq_along(test_idx)){
    param_idx <- test_idx[i]
    
    # Adjust for REML vs ML indexing
    psi_param_idx <- if(REML) param_idx else param_idx - p
    
    # Create starting point for this test
    theta_start <- theta_null
    
    # For covariance parameters, ensure positive semi-definite starting values
    if(psi_param_idx <= (r - 1) && !is_var_param[psi_param_idx]){
      # Find which term this parameter belongs to
      term_idx <- which(psi_param_idx <= last_idx)[1]
      term_param_idxs <- if(term_idx == 1) {
        seq_len(last_idx[1])
      } else {
        (last_idx[term_idx - 1] + 1):last_idx[term_idx]
      }
      
      # Dimension of covariance matrix for this term
      dim_i <- as.integer(0.5 * (-1 + sqrt(1 + 8 * r_i[term_idx])))
      
      # Position within term
      jj <- psi_param_idx - ifelse(term_idx == 1, 0, last_idx[term_idx - 1])
      
      # Get row and column indices for this covariance
      row_col <- get_row_col_ltri(jj, n = dim_i)
      
      # Find corresponding variance parameters
      var_idx1 <- get_idx_ltri(row = row_col[2], col = row_col[2], n = dim_i)
      var_idx2 <- get_idx_ltri(row = row_col[1], col = row_col[1], n = dim_i)
      
      # Adjust to full parameter vector
      var_idx1_full <- term_param_idxs[1] + var_idx1 - 1
      var_idx2_full <- term_param_idxs[1] + var_idx2 - 1
      
      if(REML){
        theta_start[c(var_idx1_full, var_idx2_full)] <- 
          pmax(theta_start[c(var_idx1_full, var_idx2_full)], abs(theta_start[psi_param_idx]))
      } else {
        theta_start[p + c(var_idx1_full, var_idx2_full)] <- 
          pmax(theta_start[p + c(var_idx1_full, var_idx2_full)], abs(theta_start[param_idx]))
      }
    }
    
    # Perform score test for this parameter
    out[i, ] <- score_test_lmer(lmerfit = lmerfit,
                                theta_null = theta_start,
                                test_idx = param_idx,
                                efficient = efficient,
                                expected = expected,
                                profile = TRUE)
  }
  out
}



