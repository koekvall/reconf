#' Score test for covariance parameters
#'
#' Computes the score test statistic for a hypothesis about the covariance
#' parameters of a linear mixed model.
#'
#' @inheritParams vc_ci
#' @param parm the tested parameters, as names or indices as in
#'   \code{\link{vc_ci}}. If \code{NULL}, all random-effect covariance
#'   parameters not in \code{known}.
#' @param null_value values of the tested parameters under the null
#'   hypothesis, recycled to the length of \code{parm}.
#'
#' @return An object of class \code{"htest"} with components
#'   \code{statistic}; \code{parameter}, the degrees of freedom, which is the
#'   number of tested parameters; \code{p.value}, from the chi-squared
#'   distribution; \code{null.value}, the null values named as in
#'   \code{\link{vc_ci}}; \code{alternative}; \code{method}; and
#'   \code{data.name}.
#'
#' @details The statistic is the profile score statistic: the score for the
#'   tested parameters, evaluated at the maximizer of the likelihood with the
#'   tested parameters at their null values and the parameters in
#'   \code{known} at their estimates, and standardized by the efficient
#'   information. The maximization is over the other covariance parameters
#'   and, for a maximum likelihood fit, the fixed effects. It starts from the
#'   estimates in the model, with the covariances of a random effect
#'   whose variance is tested at zero set to zero. If the covariance matrix
#'   of the response is not positive definite at this start, free entries
#'   of \eqn{\Psi} are changed to make its blocks weakly diagonally
#'   dominant, and the error variance is doubled until it is.
#'
#' @references Shedden, M. and Ekvall, K. O. (2026). Score-based confidence
#'   intervals for variance-covariance parameters in linear mixed models.
#'   arXiv:2610.04181.
#'
#' @examples
#' library(lme4)
#' fit <- lmer(Reaction ~ Days + (Days | Subject), data = sleepstudy)
#'
#' # All random-effect parameters zero
#' vc_test(fit)
#'
#' # Random slope variance zero
#' vc_test(fit, parm = "var_Days|Subject")
#'
#' @export
vc_test <- function(object, parm = NULL, null_value = 0, expected = TRUE,
                    known = NULL,
                    method = c("auto", "q_side", "n_side", "spectral"), ...)
{
  setup <- .as_vc_model(object, match.arg(method), !missing(method))
  r     <- setup$r
  st    <- setup$structure
  .check_flags(list(expected = expected))
  .check_dots(list(...))

  # Indices into the covariance parameters, as in vc_ci
  known_idx <- .parm_index(known, setup$names, "known")
  test_idx  <- .parm_index(parm, setup$names, "parm")
  if (is.null(test_idx)) {
    test_idx <- setdiff(seq_len(r - 1), known_idx)
    if (length(test_idx) == 0) {
      stop("no parameters to test: every random-effect parameter is in known")
    }
  } else if (length(intersect(test_idx, known_idx)) > 0) {
    stop("parm and known should not overlap")
  }
  k <- length(test_idx)
  if (!is.numeric(null_value) || !(length(null_value) %in% c(1L, k))) {
    stop("null_value should be numeric of length 1 or length(parm)")
  }
  null_value <- rep_len(null_value, k)
  if (r %in% test_idx && null_value[test_idx == r] <= 0) {
    stop("The null value of the error variance must be positive")
  }

  # Start: the estimates, with the tested parameters at their null values.
  # Covariances of a random effect whose variance is tested at zero start at
  # zero; this avoids slow optimization along the weakly identified
  # covariance direction at zero variance.
  psi <- setup$psi_hat
  psi[test_idx] <- null_value
  fixed <- c(test_idx, known_idx)
  for (j in test_idx[null_value == 0 & st$is_var[test_idx] & test_idx < r]) {
    in_rows <- colSums(st$n_off[st$on_diag[, j], , drop = FALSE]) > 0
    cov_j <- which(in_rows & !st$is_var[-r])
    psi[setdiff(cov_j, fixed)] <- 0
  }
  # Construct a feasible start if this one is not; see .feasible_start
  psi <- .feasible_start(psi, fixed, setup)
  if (is.null(psi)) {
    stop("No feasible starting point: the error variance is fixed and a ",
         "fixed variance is below the absolute fixed covariances in its row ",
         "of Psi")
  }

  # Maximize over the nuisance parameters; the default iteration limit is
  # that of maximize_loglik
  opt_idx <- seq_len(r)[-fixed]
  if (length(opt_idx) > 0) {
    psi <- maximize_loglik(start_val = psi, opt_idx = opt_idx, ll = setup$ll,
                           expected = expected, ...)$arg
  }

  test_stat <- as.numeric(score_stat(psi = psi, test_idx = test_idx,
                                     ll = setup$ll, expected = expected,
                                     signed = FALSE, known_idx = known_idx))

  structure(list(statistic = c(S = test_stat),
                 parameter = c(df = k),
                 p.value = stats::pchisq(test_stat, df = k,
                                         lower.tail = FALSE),
                 null.value = stats::setNames(null_value,
                                              setup$names[test_idx]),
                 alternative = "two.sided",
                 method = paste0("Score test (",
                                 if (setup$REML) "REML" else "ML", ")"),
                 data.name = paste(deparse(substitute(object)),
                                   collapse = " ")),
            class = "htest")
}
