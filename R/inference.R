#' Maximize log-likelihood with respect to specified parameters
#'
#' Optimizes a subset of the covariance parameters of a linear mixed effects
#' model while holding the others fixed. Uses the trust region algorithm for
#' optimization.
#'
#' @param start_val Numeric vector of starting values for the covariance
#'   parameters \eqn{\psi}, of length \eqn{r}.
#' @param opt_idx Integer vector specifying which elements of \code{start_val}
#'   to optimize. All other parameters are held fixed at their starting values.
#'   Must have length > 0 and contain unique positive integers not exceeding
#'   \code{length(start_val)}.
#' @param ll The log-likelihood as a function of \eqn{\psi}, from
#'   \code{make_loglik}.
#' @param expected Logical. If \code{TRUE}, use expected Fisher information
#'   matrix; otherwise use observed information. Default is \code{TRUE}.
#' @param rinit Initial trust-region radius passed to \code{\link[trust]{trust}}.
#'   Default is 1. With the default \code{parscale}, radii are lengths in the
#'   diagonal information metric at \code{start_val}, so 1 is roughly one
#'   standard error.
#' @param rmax Maximum trust-region radius passed to \code{\link[trust]{trust}}.
#'   Default is 100, in the same units as \code{rinit}.
#' @param parscale Numeric vector of positive scaling factors for the free
#'   parameters, passed to \code{\link[trust]{trust}}; the trust region is
#'   \eqn{\|diag(parscale) s\| \le r} in the step \eqn{s}. If \code{NULL}
#'   (default), the square root of the diagonal of the information matrix at
#'   \code{start_val} is used, so the region approximates a ball in the
#'   information metric; if that diagonal is not finite and positive, no
#'   rescaling is done.
#' @param warn_nonconv Logical. If \code{TRUE} (default), a warning is issued
#'   when the trust-region optimizer does not converge. Set to \code{FALSE}
#'   when non-convergence is expected by design (e.g., when \code{iterlim = 1L}
#'   is used for a one-step update).
#' @param check If \code{TRUE} (default), validate arguments. Internal callers
#'   in loops set \code{FALSE} to skip redundant validation.
#' @param ... Additional arguments passed to \code{\link[trust]{trust}} optimizer,
#'   such as tolerance settings or iteration limits.
#'
#' @return A list with components:
#'   \item{arg}{Numeric vector of optimized parameter values. Parameters not in
#'     \code{opt_idx} retain their starting values from \code{start_val}.}
#'   \item{value}{The maximized log-likelihood value.}
#'   \item{conv}{Logical indicating whether the optimization converged.}
#'   \item{iter}{Integer giving the number of iterations performed.}
#'
#' @noRd
maximize_loglik <- function(start_val, opt_idx, ll, expected = TRUE,
                            rinit = 1, rmax = 100, parscale = NULL,
                            warn_nonconv = TRUE, check = TRUE, ...) {
  if (check) {
    .check_psi(start_val, r = attr(ll, "r"))
    .check_idx(opt_idx, "opt_idx", length(start_val), unique = TRUE)
    .check_flags(list(expected = expected))
  }

  # Objective for trust(): negative log-likelihood in the free parameters
  obj_fun <- function(x) {
    psi <- start_val
    psi[opt_idx] <- x
    ll_things <- ll(psi, expected = expected)
    list("value" = -ll_things$value, "gradient" = -ll_things$score[opt_idx],
         "hessian" = as.matrix(ll_things$inf_mat[opt_idx, opt_idx]))
  }

  if (is.null(parscale)) {
    # trust's region is ||diag(parscale) s|| <= r; the square root of the
    # information diagonal at the start makes it the diagonal information
    # metric there.
    ps <- sqrt(diag(as.matrix(obj_fun(start_val[opt_idx])$hessian)))
    if (all(is.finite(ps)) && all(ps > 0) && all(is.finite(1 / ps))) {
      parscale <- ps
    }
  }

  trust_args <- list(objfun = obj_fun, parinit = start_val[opt_idx],
                     rinit = rinit, rmax = rmax, ...)
  if (!is.null(parscale)) trust_args$parscale <- parscale
  fit <- do.call(trust::trust, trust_args)

  if (!fit$converged && warn_nonconv) {
    warning("The optimizer did not converge in ", fit$iterations,
            " iterations")
  }

  start_val[opt_idx] <- fit$argument
  names(start_val) <- paste0("psi", seq_along(start_val))
  list("arg" = start_val, "value" = -fit$value, "conv" = fit$converged,
       "iter" = fit$iterations)
}

#' Score test statistic
#'
#' Computes the score test statistic for testing hypotheses about the
#' covariance parameters of a linear mixed effects model, standardized by the
#' efficient information, as the quadratic form (unsigned) or the signed root
#' statistic.
#'
#' @param psi Numeric vector of covariance parameters, of length \eqn{r}, at
#'   which to evaluate the statistic.
#' @param test_idx Integer vector specifying which elements of \code{psi}
#'   are being tested. These are the parameters constrained by the null hypothesis.
#' @param ll The log-likelihood as a function of \eqn{\psi}, from
#'   \code{make_loglik}.
#' @param expected Logical. If \code{TRUE}, use expected Fisher information
#'   matrix; otherwise use observed information. Default is \code{TRUE}.
#' @param signed Logical. If \code{TRUE}, return the signed root statistic
#'   (vector). If \code{FALSE}, return the quadratic form statistic (scalar
#'   for single parameter tests). Default is \code{FALSE}.
#' @param known_idx Integer vector or \code{NULL} specifying which elements of
#'   \code{psi} (other than those in \code{test_idx}) have known values and
#'   should be held fixed (not treated as nuisance parameters to be profiled over).
#'   If \code{NULL} (default), all parameters not in \code{test_idx} are treated
#'   as nuisance parameters. Must not overlap with \code{test_idx}.
#' @param check If \code{TRUE} (default), validate arguments and warn when the
#'   information matrix is poorly conditioned. Internal callers in loops set
#'   \code{FALSE}.
#'
#' @return Numeric value or vector containing the score test statistic with
#'   attributes \code{"score"} and \code{"info"}. If \code{signed = FALSE},
#'   returns a scalar (the quadratic form). If \code{signed = TRUE}, returns a
#'   vector (the signed root statistic). The \code{"score"} attribute contains
#'   the score vector for the test parameter(s), and \code{"info"} contains
#'   the used information for the test parameter(s). Under the
#'   null hypothesis, the squared statistic asymptotically follows a chi-squared
#'   distribution with degrees of freedom equal to \code{length(test_idx)}.
#'
#' @details
#' The score test statistic is computed as:
#' \deqn{T = S_t^T I_{tt}^{-1} S_t}
#' where \eqn{S_t} is the score vector for the test parameters and \eqn{I_{tt}}
#' is the efficient information matrix. When there are nuisance parameters
#' (parameters not in \code{test_idx} or \code{known_idx}), it is
#' \deqn{I_{tt}^{eff} = I_{tt} - I_{tn} I_{nn}^{-1} I_{nt}}
#' where subscripts \eqn{t} denote test parameters and \eqn{n} denote nuisance
#' parameters. For a maximum likelihood fit the fixed effects are profiled out
#' of \code{ll}, so they are nuisance parameters as well.
#'
#' When \code{known_idx} is specified, those parameters are treated as fixed and
#' known (not as nuisance parameters). This is useful when some parameters have
#' been estimated separately or are constrained to specific values.
#'
#' @noRd
score_stat <- function(psi, test_idx, ll, expected = TRUE, signed = FALSE,
                       known_idx = NULL, check = TRUE)
{
  if (check) {
    .check_psi(psi, r = attr(ll, "r"))
    .check_flags(list(expected = expected, signed = signed))
    .check_idx(test_idx, "test_idx", length(psi))
    if (!is.null(known_idx) && length(known_idx) > 0) {
      .check_idx(known_idx, "known_idx", length(psi))
      if (length(intersect(test_idx, known_idx)) > 0) {
        stop("test_idx and known_idx should not overlap", call. = FALSE)
      }
    }
  }

  test_idx <- unique(test_idx)
  if (!is.null(known_idx) && length(known_idx) > 0) {
    known_idx <- unique(known_idx)
  } else {
    known_idx <- NULL  # Treat empty vector as NULL
  }

  ll_things <- ll(psi, get_val = FALSE, expected = expected)
  # Check condition of information matrix (skipped in hot loops via check).
  # The matrix is small (r-dimensional), so the exact spectral
  # condition number is affordable; the threshold says three quarters of the
  # double-precision digits are gone.
  if (check) {
    cond <- tryCatch(
      kappa(ll_things$inf_mat, exact = TRUE),
      error = function(e) Inf
    )
    if (cond > .Machine$double.eps^-0.75) {
      warning("The information matrix has condition number ",
              format(cond, digits = 3, scientific = TRUE))
    }
  }

  inf_mat <- ll_things$inf_mat[test_idx, test_idx, drop = FALSE]

  # Efficient information when there are nuisance parameters; known
  # parameters are excluded entirely (fixed, not profiled)
  exclude_idx <- if (is.null(known_idx)) test_idx else c(test_idx, known_idx)
  nuis_idx <- setdiff(seq_along(ll_things$score), exclude_idx)

  if (length(nuis_idx) > 0) {
    # Efficient information from one Cholesky of the joint (nuisance, test)
    # block: with R'R = I[(n,t),(n,t)], the Schur complement
    # I_tt - I_tn I_nn^{-1} I_nt equals R22'R22, positive semidefinite by
    # construction. The explicit subtraction suffers catastrophic
    # cancellation when the information is nearly singular, which is the
    # near-boundary regime this method targets. It remains only as the
    # fallback when the joint block is not positive definite (possible for
    # observed information).
    joint_idx <- c(nuis_idx, test_idx)
    ch <- tryCatch(chol(ll_things$inf_mat[joint_idx, joint_idx]),
                   error = function(e) NULL)
    if (!is.null(ch)) {
      tt <- length(nuis_idx) + seq_along(test_idx)
      inf_mat <- crossprod(ch[tt, tt, drop = FALSE])
    } else {
      A_nt <- ll_things$inf_mat[nuis_idx, test_idx, drop = FALSE]
      I_nn <- ll_things$inf_mat[nuis_idx, nuis_idx, drop = FALSE]
      inf_mat <- inf_mat - crossprod(A_nt, .solve_sym_eigen(I_nn, A_nt))
    }
  }

  if (signed) {
    ed <- eigen(inf_mat, symmetric = TRUE)
    # Floor eigenvalues at a tolerance relative to the largest before taking
    # 1/sqrt, and say so rather than clamp silently (suppressed in hot loops
    # via check, like the condition diagnostic above)
    tol <- max(ed$values, .Machine$double.eps) * .Machine$double.eps
    if (check && any(ed$values < tol)) {
      warning("The efficient information for the tested parameters is ",
              "singular to working precision; eigenvalues below ",
              format(tol, digits = 3), " were set to that value")
    }
    ev <- pmax(ed$values, tol)
    inf_root <- ed$vectors %*% (sqrt(ev) * t(ed$vectors))
    test_stat <- solve(inf_root, ll_things$score[test_idx])
  } else {
    test_stat <- crossprod(ll_things$score[test_idx], solve(inf_mat,
                                    ll_things$score[test_idx]))
  }
  test_stat <- as.vector(test_stat)
  attr(test_stat, "score") <- ll_things$score[test_idx]
  attr(test_stat, "info") <- inf_mat
  test_stat
}
