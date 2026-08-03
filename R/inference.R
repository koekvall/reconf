#' Maximize log-likelihood with respect to specified parameters
#'
#' Optimizes a subset of parameters in a linear mixed effects model while
#' holding others fixed. Uses the trust region algorithm for optimization.
#'
#' @param start_val Numeric vector of starting parameter values. If
#'   \code{REML = FALSE}, this should be of length \eqn{p + r} with
#'   \code{start_val = c(beta, psi)}, where \eqn{p} is the number of fixed
#'   effects and \eqn{r} is the number of variance parameters. If
#'   \code{REML = TRUE}, this should be of length \eqn{r} with
#'   \code{start_val = psi}.
#' @param opt_idx Integer vector specifying which elements of \code{start_val}
#'   to optimize. All other parameters are held fixed at their starting values.
#'   Must have length > 0 and contain unique positive integers not exceeding
#'   \code{length(start_val)}.
#' @param Y Numeric vector of length \eqn{n} containing the response values.
#' @param X Numeric matrix of size \eqn{n \times p} containing fixed effect
#'   predictors.
#' @param Z Sparse matrix of size \eqn{n \times q} containing the random effect
#'   design matrix.
#' @param Hlist List of sparse matrices determining how \eqn{\psi} is mapped to
#'   the covariance matrix \eqn{\Psi} (see \code{?loglikelihood}).
#' @param expected Logical. If \code{TRUE}, use expected Fisher information
#'   matrix; otherwise use observed information. Default is \code{TRUE}.
#' @param REML Logical. If \code{TRUE}, use restricted maximum likelihood;
#'   otherwise use maximum likelihood. Default is \code{TRUE}.
#' @param precomp List or \code{NULL} containing precomputed quantities to speed
#'   up computation (see \code{?get_precomp}). If \code{NULL}, quantities are
#'   computed internally.
#' @param rinit Initial trust-region radius passed to \code{\link[trust]{trust}}.
#'   Default is 1.
#' @param rmax Maximum trust-region radius passed to \code{\link[trust]{trust}}.
#'   Default is 100.
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
#' @seealso \code{\link{score_stat}},
#'   \code{\link{loglikelihood}}, \code{\link[trust]{trust}}
#'
#' @keywords internal
maximize_loglik <- function(start_val, opt_idx, Y, X, Z, Hlist, expected = TRUE,
                             REML = TRUE, precomp = NULL,
                             rinit = 1, rmax = 100, warn_nonconv = TRUE,
                             check = TRUE, ...) {
  r <- length(Hlist) + 1
  p <- ncol(X)

  if (check) {
    .check_lmm_args(start_val, "start_val", Y = Y, X = X, Z = Z,
                    Hlist = Hlist, REML = REML, precomp = precomp,
                    flags = list(REML = REML, expected = expected))
    .check_idx(opt_idx, "opt_idx", length(start_val), unique = TRUE)
  }

  if(is.null(precomp)) {
    precomp <- get_precomp(Y = Y, X = X, Z = Z, REML = REML, Hlist = Hlist)
  }

  # Objective for trust(): negative log-likelihood in the free parameters
  obj_fun <- function(x) {
    theta <- start_val
    theta[opt_idx] <- x
    psi <- if (REML || p == 0) theta else theta[(p + 1):(p + r)]
    b <- if (!REML && p > 0) theta[seq_len(p)] else NULL
    ll_things <- loglikelihood(psi = psi, b = b, Y = Y, X = X, Z = Z,
                               Hlist = Hlist, REML = REML,
                               get_val = TRUE, get_score = TRUE, get_inf = TRUE,
                               get_beta = !REML, expected = expected,
                               precomp = precomp, check = FALSE)
    list("value" = -ll_things$value, "gradient" = -ll_things$score[opt_idx],
         "hessian" = as.matrix(ll_things$inf_mat[opt_idx, opt_idx]))
  }

  fit <- trust::trust(objfun = obj_fun, parinit = start_val[opt_idx],
                      rinit = rinit, rmax = rmax, ...)

  if (!fit$converged && warn_nonconv) {
    warning("Optimization did not converge. Results may be unreliable. ",
            "Iterations: ", fit$iterations)
  }

  start_val[opt_idx] <- fit$argument
  names(start_val) <- if(REML || p == 0) paste0("psi", 1:r) else c(paste0("b", 1:p), paste0("psi", 1:r))
  list("arg" = start_val, "value" = -fit$value, "conv" = fit$converged,
       "iter" = fit$iterations)
}

#' Score test statistic
#'
#' Computes the score test statistic for testing hypotheses about parameters
#' in a linear mixed effects model. The test statistic can be computed with or
#' without efficient information (accounting for nuisance parameters), and can
#' return either the quadratic form (unsigned) or the signed root statistic.
#'
#' @param theta Numeric vector of parameter values at which to evaluate the
#'   score test statistic. If \code{REML = FALSE}, this should be of length
#'   \eqn{p + r} with \code{theta = c(beta, psi)}, where \eqn{p} is the number
#'   of fixed effects and \eqn{r} is the number of variance parameters. If
#'   \code{REML = TRUE}, this should be of length \eqn{r} with \code{theta = psi}.
#' @param test_idx Integer vector specifying which elements of \code{theta}
#'   are being tested. These are the parameters constrained by the null hypothesis.
#' @param Y Numeric vector of length \eqn{n} containing the response values.
#' @param X Numeric matrix of size \eqn{n \times p} containing fixed effect
#'   predictors.
#' @param Z Sparse matrix of size \eqn{n \times q} containing the random effect
#'   design matrix.
#' @param Hlist List of sparse matrices determining how \eqn{\psi} is mapped to
#'   the covariance matrix \eqn{\Psi} (see \code{?loglikelihood}).
#' @param REML Logical. If \code{TRUE}, use restricted maximum likelihood;
#'   otherwise use maximum likelihood. Default is \code{TRUE}.
#' @param expected Logical. If \code{TRUE}, use expected Fisher information
#'   matrix; otherwise use observed information. Default is \code{TRUE}.
#'   Note: Observed information is not available for REML.
#' @param efficient Logical. If \code{TRUE}, use efficient information that
#'   accounts for estimation of nuisance parameters. If \code{FALSE}, use
#'   the full information matrix without adjustment. Default is \code{TRUE}.
#' @param signed Logical. If \code{TRUE}, return the signed root statistic
#'   (vector). If \code{FALSE}, return the quadratic form statistic (scalar
#'   for single parameter tests). Default is \code{FALSE}.
#' @param known_idx Integer vector or \code{NULL} specifying which elements of
#'   \code{theta} (other than those in \code{test_idx}) have known values and
#'   should be held fixed (not treated as nuisance parameters to be profiled over).
#'   If \code{NULL} (default), all parameters not in \code{test_idx} are treated
#'   as nuisance parameters. Must not overlap with \code{test_idx}.
#' @param precomp List or \code{NULL} containing precomputed quantities to speed
#'   up computation (see \code{?get_precomp}). If \code{NULL}, all quantities
#'   are computed from scratch.
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
#' is the information matrix (possibly efficient if \code{efficient = TRUE}).
#'
#' When \code{efficient = TRUE} and there are nuisance parameters (parameters
#' not in \code{test_idx} or \code{known_idx}), the efficient information is:
#' \deqn{I_{tt}^{eff} = I_{tt} - I_{tn} I_{nn}^{-1} I_{nt}}
#' where subscripts \eqn{t} denote test parameters and \eqn{n} denote nuisance
#' parameters.
#'
#' When \code{known_idx} is specified, those parameters are treated as fixed and
#' known (not as nuisance parameters). This is useful when some parameters have
#' been estimated separately or are constrained to specific values.
#'
#' @seealso \code{\link{loglikelihood}}
#'
#' @keywords internal
score_stat <- function(theta, test_idx, Y, X, Z, Hlist, REML = TRUE,
                       expected = TRUE, efficient = TRUE, signed = FALSE,
                       known_idx = NULL,
                       precomp = NULL, check = TRUE)
{
  p <- ncol(X)
  r <- length(Hlist) + 1

  if (check) {
    .check_lmm_args(theta, "theta", Y = Y, X = X, Z = Z, Hlist = Hlist,
                    REML = REML, precomp = precomp,
                    flags = list(REML = REML, expected = expected,
                                 efficient = efficient, signed = signed))
    .check_idx(test_idx, "test_idx", length(theta))
    if (!is.null(known_idx) && length(known_idx) > 0) {
      .check_idx(known_idx, "known_idx", length(theta))
      assertthat::assert_that(length(intersect(test_idx, known_idx)) == 0,
                              msg = "test_idx and known_idx should not overlap")
    }
    if(!expected && REML){
      warning("Observed information not available for restricted likelihood; using
              expected.")
    }
  }

  test_idx <- unique(test_idx)
  if (!is.null(known_idx) && length(known_idx) > 0) {
    known_idx <- unique(known_idx)
  } else {
    known_idx <- NULL  # Treat empty vector as NULL
  }

  psi <- if(REML || p < 1) theta else theta[-seq_len(p)]
  b <- if(!REML && p >= 1) theta[seq_len(p)] else NULL

  ll_things <- loglikelihood(psi = psi,
                             b = b,
                             Y = Y,
                             X = X,
                             Z = Z,
                             Hlist = Hlist,
                             REML = REML,
                             get_val = FALSE,
                             get_score = TRUE,
                             get_inf = TRUE,
                             get_beta = (!REML && p>=1),
                             expected = expected,
                             precomp = precomp,
                             check = FALSE)
  # Check condition of information matrix (skipped in hot loops via check).
  # The matrix is small (at most (p + r)-dimensional), so the exact spectral
  # condition number is affordable; the threshold says three quarters of the
  # double-precision digits are gone.
  if (check) {
    cond <- tryCatch(
      kappa(ll_things$inf_mat, exact = TRUE),
      error = function(e) Inf
    )
    if (cond > .Machine$double.eps^-0.75) {
      warning("Information matrix is poorly conditioned (condition number: ",
              format(cond, scientific = TRUE), "). Results may be unreliable.")
    }
  }

  inf_mat <- ll_things$inf_mat[test_idx, test_idx, drop = FALSE]

  # Use efficient information only if there are nuisance parameters; known
  # parameters are excluded entirely (fixed, not profiled)
  exclude_idx <- if (is.null(known_idx)) test_idx else c(test_idx, known_idx)
  nuis_idx <- setdiff(seq_along(ll_things$score), exclude_idx)

  if (efficient && length(nuis_idx) > 0) {
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
    # 1/sqrt; a floored eigenvalue means the statistic is not reliable in
    # that direction, so flag it rather than clamp silently (suppressed in
    # hot loops via check, like the condition diagnostic above)
    tol <- max(ed$values, .Machine$double.eps) * .Machine$double.eps
    if (check && any(ed$values < tol)) {
      warning("Efficient information for the test parameter(s) is nearly ",
              "singular; the signed statistic may be unreliable.")
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
