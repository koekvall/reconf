#' Maximize the log-likelihood over a subset of the covariance parameters
#'
#' Maximizes the log-likelihood over the covariance parameters in
#' \code{opt_idx}, holding the others fixed, by the trust region method of
#' \code{\link[trust]{trust}}.
#'
#' @param start_val Numeric vector of starting values for the covariance
#'   parameters \eqn{\psi}, of length \eqn{r}.
#' @param opt_idx Indices of the elements of \code{start_val} to maximize
#'   over; the others are held at their starting values. Unique positive
#'   integers, at most \code{length(start_val)}.
#' @param ll The log-likelihood as a function of \eqn{\psi}, from
#'   \code{make_loglik}.
#' @param expected Logical. If \code{TRUE}, the expected information is the
#'   Hessian of the objective; otherwise the observed information.
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
#'   \code{start_val} when that is finite and positive, so the region
#'   approximates a ball in the information metric, and otherwise the
#'   default of \code{trust}.
#' @param iterlim Maximum number of iterations passed to
#'   \code{\link[trust]{trust}}. Default is 1000: profiling far from the
#'   estimates can require many cheap iterations along nearly flat
#'   directions, where the trust-region method converges linearly.
#' @param warn_nonconv Logical. If \code{TRUE} (default), a warning is issued
#'   when the optimizer does not converge; \code{FALSE} for the one-step
#'   update with \code{iterlim = 1}.
#' @param check If \code{TRUE} (default), validate arguments. Internal callers
#'   in loops set \code{FALSE} to skip redundant validation.
#' @param ... Additional arguments passed to \code{\link[trust]{trust}}, such
#'   as the tolerances \code{fterm} and \code{mterm}.
#'
#' @return A list with components:
#'   \item{arg}{The parameters, with those in \code{opt_idx} at the maximizer
#'     and the others at their values in \code{start_val}.}
#'   \item{value}{The maximized log-likelihood.}
#'   \item{conv}{Logical, whether the optimizer converged.}
#'   \item{iter}{The number of iterations.}
#'
#' @noRd
maximize_loglik <- function(start_val, opt_idx, ll, expected = TRUE,
                            rinit = 1, rmax = 100, parscale = NULL,
                            iterlim = 1000L, warn_nonconv = TRUE,
                            check = TRUE, ...) {
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
    if (all(is.finite(ps)) && all(ps > 0)) {
      parscale <- ps
    }
  }

  trust_args <- list(objfun = obj_fun, parinit = start_val[opt_idx],
                     rinit = rinit, rmax = rmax, iterlim = iterlim, ...)
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

#' Score statistic
#'
#' Computes the score statistic for the covariance parameters in
#' \code{test_idx}, standardized by their efficient information, as a
#' quadratic form or as a signed root.
#'
#' @param psi Numeric vector of covariance parameters, of length \eqn{r}, at
#'   which to evaluate the statistic.
#' @param test_idx Indices of the tested parameters.
#' @param ll The log-likelihood as a function of \eqn{\psi}, from
#'   \code{make_loglik}.
#' @param expected Logical. If \code{TRUE}, the expected information;
#'   otherwise the observed information.
#' @param signed Logical. If \code{TRUE}, the signed root statistic, a vector
#'   of length \code{length(test_idx)}; otherwise the quadratic form, a
#'   scalar.
#' @param known_idx Indices of parameters held fixed at their values in
#'   \code{psi}, or \code{NULL}. The parameters in neither \code{test_idx}
#'   nor \code{known_idx} are the nuisance parameters. Must not overlap with
#'   \code{test_idx}.
#' @param check If \code{TRUE} (default), validate arguments and warn when the
#'   information matrix is poorly conditioned. Internal callers in loops set
#'   \code{FALSE}.
#'
#' @return The statistic, with attributes \code{"score"}, the score for the
#'   tested parameters, and \code{"info"}, their efficient information.
#'
#' @details
#' The quadratic form is \eqn{S_t' I_{tt}^{-1} S_t}, with \eqn{S_t} the score
#' for the tested parameters and \eqn{I_{tt}} their efficient information
#' \eqn{I_{tt} - I_{tn} I_{nn}^{-1} I_{nt}}, where \eqn{n} indexes the
#' nuisance parameters. The signed root is \eqn{I_{tt}^{-1/2} S_t} with the
#' symmetric square root. For a maximum likelihood fit the fixed effects are
#' profiled out of \code{ll}.
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
