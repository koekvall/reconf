#' Confidence intervals for covariance parameters
#'
#' Computes confidence intervals for covariance parameters of a linear mixed
#' model by inverting a signed profile score or likelihood ratio statistic.
#' See the references for the theory.
#'
#' @param object a fitted linear mixed model of class \code{lmerMod} from
#'   \code{lme4::lmer}, or a \code{\link{vc_model}}.
#' @param parm parameters, as names (see Value) or as indices in the order
#'   of \code{\link{vc_model}}, which for an lme4 fit is that of
#'   \code{as.data.frame(VarCorr(fit), order = "lower.tri")}, with the error
#'   variance last. If \code{NULL}, all covariance parameters not in
#'   \code{known}.
#' @param level confidence level.
#' @param statistic \code{"score"} for the signed profile score statistic,
#'   standardized by the efficient information, or \code{"rlrt"} for the
#'   signed root of the profile likelihood ratio statistic, which for a
#'   restricted likelihood fit is the restricted likelihood ratio. See
#'   Details.
#' @param expected if \code{TRUE} use the expected information, otherwise the
#'   observed information, both in the statistic and in the maximization over
#'   nuisance parameters.
#' @param known parameters, as for \code{parm}, held fixed at their
#'   estimates.
#' @param nonneg if \code{TRUE}, the interval for a variance, the error
#'   variance included, is intersected with \eqn{[0, \infty)}. An empty
#'   intersection gives \code{NA} bounds and a warning.
#' @param onestep if \code{TRUE}, the nuisance parameters at each step are
#'   one iteration of the optimizer, with unbounded trust region, from those
#'   at the previous step. Requires \code{statistic = "score"}.
#' @param step_size initial step length of the search (see Details), on the
#'   scale of the parameter. If \code{NULL}, one fortieth of the standard
#'   error from the expected information at the estimate.
#' @param num_points maximum number of steps in each direction. The search
#'   also stops at distance \code{num_points * step_size} from its origin.
#' @param method computational method for the likelihood: \code{"auto"},
#'   \code{"q_side"}, \code{"n_side"}, or \code{"spectral"}. See Details.
#'   A \code{vc_model} carries the method it was built with.
#' @param ... arguments passed to the optimizer \code{\link[trust]{trust}}:
#'   \code{rinit}, \code{rmax}, \code{iterlim}, \code{fterm}, and
#'   \code{mterm}; and \code{warn_nonconv}, which if \code{FALSE} suppresses
#'   the warning when the optimizer does not converge. The default
#'   \code{iterlim} is 1000.
#'
#' @return A matrix of class \code{vc_ci} with one row per parameter and
#'   columns \code{estimate}, \code{lower}, and \code{upper}. For an lme4
#'   fit the estimates are those of lme4, and rows are named
#'   \code{"var_x|g"} for the variance of the random effect \code{x} with
#'   grouping factor \code{g}, \code{"cov_y.x|g"} for the covariance of
#'   \code{x} and \code{y}, and \code{"var_Residual"} for the error
#'   variance, following the names \code{lme4} uses for standard deviations
#'   and correlations.
#'
#' @details
#' The interval at level \eqn{1 - \alpha} for a parameter \eqn{\psi_j} is
#' \deqn{\{\psi_j : |T(\psi_j)| \leq z_{1-\alpha/2}\},}
#' where \eqn{T} is the signed statistic with the other parameters profiled
#' out and \eqn{z_{1-\alpha/2}} is the standard normal quantile. The
#' parameter set is the set of covariance parameters for which the covariance
#' matrix of the response is positive definite. It includes negative
#' variances, and the nuisance parameters are maximized over it.
#'
#' The bounds are found by a search outward from the estimate in each
#' direction, with the nuisance parameters at each step computed from those
#' at the previous step. Step lengths are chosen by secant extrapolation of
#' the statistic, each at most twice the previous one. If the covariance
#' matrix of the response is not positive definite at the start of a step,
#' the step is halved, up to 20 times, after which a starting point is
#' constructed as in \code{\link{vc_test}}. A bound is where the
#' statistic crosses \eqn{\pm z_{1-\alpha/2}}, located to within
#' \code{step_size}. If there is no crossing within the search distance, the
#' bound is \code{-Inf} or \code{Inf} and a warning is issued, reporting
#' the statistic at the end of the search. As a variance grows, the
#' statistic tends to a finite limit, which for a grouping factor with few
#' levels can be above \eqn{-z_{1-\alpha/2}}; the upper bound is then
#' infinite for any search distance. If the statistic at the estimate
#' exceeds the critical value, the search starts instead from the maximizer
#' over the parameter set, with a warning.
#'
#' With \code{statistic = "rlrt"}, the search starts from the maximizer over
#' the parameter set. A warning is issued if that maximizer has a negative
#' variance; the interval for that variance can then lie entirely below zero.
#'
#' The likelihood is computed by one of three methods. \code{"q_side"} uses
#' sparse \eqn{q \times q} matrices, where \eqn{q} is the number of random
#' effects. \code{"n_side"} uses dense \eqn{n \times n} matrices. For models
#' with one random-effect variance, \code{"spectral"} uses one
#' eigendecomposition and \eqn{O(n)} operations per evaluation.
#' \code{"auto"} uses \code{"spectral"} or \code{"n_side"} if
#' \eqn{q \geq n} and more than 10 percent of the entries of the
#' random-effects design matrix are nonzero, and \code{"q_side"} otherwise.
#'
#' @references Shedden, M. and Ekvall, K. O. (2026). Score-based confidence
#'   intervals for variance-covariance parameters in linear mixed models.
#'   arXiv:2610.04181.
#'
#' @seealso \code{\link{vc_model}}, \code{\link{vc_test}}
#'
#' @examples
#' library(lme4)
#' fit <- lmer(Reaction ~ Days + (Days | Subject), data = sleepstudy)
#'
#' # 95% intervals for all covariance parameters
#' vc_ci(fit)
#'
#' # The random slope variance, by name or by index
#' vc_ci(fit, parm = "var_Days|Subject")
#' vc_ci(fit, parm = 3)
#' @export
vc_ci <- function(object, parm = NULL, level = 0.95,
                  statistic = c("score", "rlrt"), expected = TRUE,
                  known = NULL, nonneg = TRUE, onestep = FALSE,
                  step_size = NULL, num_points = 500L,
                  method = c("auto", "q_side", "n_side", "spectral"), ...) {
  statistic <- match.arg(statistic)
  setup <- .as_vc_model(object, match.arg(method), !missing(method))
  known <- .parm_index(known, setup$names, "known")
  parm  <- .parm_index(parm, setup$names, "parm")
  if (is.null(parm)) {
    parm <- setdiff(seq_len(setup$r), known)
    if (length(parm) == 0) {
      stop("no parameters to compute intervals for: every parameter is in ",
           "known")
    }
  } else if (length(intersect(parm, known)) > 0) {
    stop("parm and known should not overlap")
  }

  .check_flags(list(expected = expected, nonneg = nonneg, onestep = onestep))
  if (!(is.numeric(level) && length(level) == 1L && level > 0 && level < 1)) {
    stop("level must be a single number in (0, 1)")
  }
  if (!(is.null(step_size) ||
        (is.numeric(step_size) && length(step_size) == 1L &&
         is.finite(step_size) && step_size > 0))) {
    stop("step_size must be NULL or a single positive number")
  }
  if (!(is.numeric(num_points) && length(num_points) == 1L &&
        num_points >= 2)) {
    stop("num_points must be a single integer >= 2")
  }
  if (statistic == "rlrt" && onestep) {
    stop("onestep = TRUE is not available with statistic = 'rlrt': ",
         "the likelihood ratio requires fully profiled nuisance parameters")
  }
  .check_dots(list(...))

  # The information at the estimate sets the default step size for every
  # parameter, so it is computed once here
  inf_hat <- if (is.null(step_size)) {
    setup$ll(setup$psi_hat, get_val = FALSE, get_score = FALSE)$inf_mat
  }

  # The rlrt reference maximum does not depend on the parameter, so compute
  # it once and share it across parameters
  rlrt_ref <- if (statistic == "rlrt") {
    .rlrt_reference(setup, known, expected, ...)
  }

  ci <- do.call(rbind, lapply(parm, function(j) {
    .ci_lmer_core(setup, test_idx = j, level = level, step_size = step_size,
                  num_points = num_points, expected = expected,
                  known_idx = known, onestep = onestep, nonneg = nonneg,
                  statistic = statistic, rlrt_ref = rlrt_ref,
                  inf_hat = inf_hat, ...)
  }))
  .as_vc_ci(ci, level = level, REML = setup$REML, statistic = statistic)
}


# Internal helpers ---------------------------------------------------------

# Compute the CI for one parameter given a prebuilt setup; see vc_ci,
# which validates the arguments. test_idx and known_idx index the covariance
# parameters. accelerate = FALSE gives the fixed-step search, every step of
# length step_size, which the tests compare against. inf_hat is the
# information at the estimate, which vc_ci computes once for all
# parameters; if NULL and step_size is NULL, it is computed here.
.ci_lmer_core <- function(setup, test_idx, level = 0.95, step_size = NULL,
                          num_points = 500L, expected = TRUE,
                          known_idx = NULL,
                          onestep = FALSE, nonneg = TRUE, accelerate = TRUE,
                          statistic = c("score", "rlrt"), rlrt_ref = NULL,
                          inf_hat = NULL, ...) {

  statistic <- match.arg(statistic)

  # The rlrt reference maximum is taken over the same extended set on which
  # the nuisance parameters are maximized; a smaller reference would let
  # the profile exceed it. The search starts at the maximizer, where the
  # statistic is zero. The reference does not depend on the tested
  # parameter; vc_ci precomputes it and passes rlrt_ref.
  theta_origin <- setup$psi_hat
  ll_max <- NULL
  if (statistic == "rlrt") {
    if (is.null(rlrt_ref)) {
      rlrt_ref <- .rlrt_reference(setup, known_idx, expected, ...)
    }
    theta_origin <- rlrt_ref$arg
    ll_max <- rlrt_ref$value
  }

  # Determine step size from expected information if not provided.
  # Use SE/40 so roughly 40 steps cover one Wald CI half-width on each side.
  if (is.null(step_size)) {
    if (is.null(inf_hat)) {
      inf_hat <- setup$ll(setup$psi_hat, get_val = FALSE,
                          get_score = FALSE)$inf_mat
    }
    se_approx <- tryCatch(
      sqrt(solve(inf_hat)[test_idx, test_idx]),
      error = function(e) sqrt(1 / inf_hat[test_idx, test_idx])
    )
    step_size <- se_approx / 40
    if (!is.finite(step_size) || step_size <= 0) {
      stop("step_size could not be set from the information at the ",
           "estimate; supply step_size")
    }
  }

  z_crit <- stats::qnorm((1 + level) / 2)

  # Optimizer arguments shared by every profile evaluation. The one-step
  # Newton update is a trust-region solve with iterlim = 1 and a radius
  # large enough that the Newton step is unconstrained; see vc_ci.
  dots <- list(...)
  if (onestep) {
    dots[c("iterlim", "rinit", "rmax", "warn_nonconv")] <- NULL
    opt_args <- c(list(iterlim = 1L, rinit = 1e10, rmax = 1e10,
                       warn_nonconv = FALSE), dots)
  } else {
    opt_args <- dots
  }
  opt_idx <- seq_along(setup$psi_hat)[-unique(c(test_idx, known_idx))]

  # Everything the profile evaluations and the outward search need that is
  # fixed for this parameter
  ctx <- list(setup = setup, test_idx = test_idx, known_idx = known_idx,
              opt_idx = opt_idx, statistic = statistic, ll_max = ll_max,
              expected = expected, opt_args = opt_args, z_crit = z_crit,
              step_size = step_size, max_steps = as.integer(num_points),
              accelerate = accelerate)

  # Evaluate the signed statistic at the search origin. It vanishes at an
  # interior maximizer but not at a boundary estimate, and an assumed zero
  # can produce a spurious sign change in the first step. If the evaluation
  # fails, fall back to zero. Both search directions share the result.
  .origin_eval <- function() {
    ev <- .profile_stat(theta_origin, theta_origin[test_idx], ctx)
    if (is.null(ev)) ev <- list(theta = theta_origin, stat = 0)
    ev
  }
  origin_eval <- .origin_eval()

  # At a boundary estimate the statistic can exceed the critical value at
  # the origin itself. The confidence set then does not contain the
  # estimate, and an outward search from the estimate finds no bounds.
  # Restart from the extended-set maximizer, where the profile score
  # vanishes. The rlrt origin is that maximizer already.
  if (statistic == "score" && abs(origin_eval$stat) > z_crit) {
    warning("The signed score statistic at the estimate is ",
            format(origin_eval$stat, digits = 3), ", beyond the critical ",
            "value; searching from the extended-set maximizer instead.")
    if (is.null(rlrt_ref)) {
      rlrt_ref <- .rlrt_reference(setup, known_idx, expected, ...)
    }
    theta_origin <- rlrt_ref$arg
    origin_eval <- .origin_eval()
  }

  # Whether the tested parameter is a variance, from the structure matrices
  # (see .psi_structure); the nonneg clamp applies to variances only
  is_variance <- setup$structure$is_var[test_idx]
  clamp0 <- nonneg && is_variance
  origin_val <- theta_origin[test_idx]

  # The nonneg clamp is the intersection of the search interval with
  # [0, Inf). When the origin is nonnegative, the lower search takes the
  # intersection by stopping at 0. A negative origin needs no lower search:
  # the interval contains the origin, so the intersected lower bound is 0
  # if the interval reaches past 0, and the intersection is empty
  # otherwise, reported as NA below. Clamping a search that starts below 0
  # would return a lower bound above the upper bound.
  lower_clamp <- if (clamp0 && origin_val >= 0) 0 else -Inf

  upper <- .outward_bound(ctx, origin_eval, direction = 1L)
  lower <- if (clamp0 && origin_val < 0) {
    if (!is.na(upper) && upper >= 0) 0 else NA_real_
  } else {
    .outward_bound(ctx, origin_eval, direction = -1L,
                   lower_clamp = lower_clamp)
  }

  # Diagnostics for a variance whose extended-set estimate is negative
  if (is_variance && is.finite(upper) && upper < 0) {
    warning("The ", statistic, " interval for '", setup$names[test_idx],
            "' contains no nonnegative values (upper bound ",
            format(upper, digits = 3), ").",
            if (clamp0) " Reporting NA bounds." else "")
    if (clamp0) upper <- NA_real_
  } else if (is_variance && origin_val < 0) {
    warning("The extended-set estimate of '", setup$names[test_idx],
            "' is negative (", format(origin_val, digits = 3), ").",
            if (clamp0) " The lower bound is truncated at 0." else "")
  }

  if (is.finite(lower) && is.finite(upper) && lower > upper) {
    warning("Lower bound is not less than upper bound. ",
            "Consider decreasing step_size or increasing num_points.")
  }

  matrix(c(setup$psi_hat[test_idx], lower, upper), nrow = 1L,
         dimnames = list(setup$names[test_idx],
                         c("estimate", "lower", "upper")))
}

# Attach the class and display attributes of vc_ci output
.as_vc_ci <- function(ci, level, REML, statistic) {
  structure(ci, class = c("vc_ci", class(ci)), level = level,
            statistic = statistic, method = if (REML) "REML" else "ML")
}

# Maximize the (restricted) log-likelihood over the extended parameter set,
# free in every parameter except those in known_idx. Returns the
# maximize_loglik list; arg and value are the reference maximizer and
# maximum for the rlrt statistic. An optimizer error propagates, because
# without a reference no interval exists.
.rlrt_reference <- function(setup, known_idx, expected = TRUE, ...) {
  free_idx <- seq_len(setup$r)
  if (length(known_idx) > 0L) free_idx <- setdiff(free_idx, known_idx)
  do.call(maximize_loglik,
          c(list(start_val = setup$psi_hat, opt_idx = free_idx,
                 ll = setup$ll, expected = expected, check = FALSE),
            list(...)))
}

# Log-likelihood value at psi; -Inf if the evaluation fails (infeasible psi).
# ll is the likelihood function of make_loglik.
.ll_value <- function(psi, ll) {
  tryCatch(ll(psi, get_score = FALSE, get_inf = FALSE)$value,
           error = function(e) -Inf)
}

# A starting point for the maximization over the nuisance parameters at
# which the covariance matrix of the response, Sigma, is positive definite.
# fixed indexes the covariance parameters that keep their values. If psi is
# infeasible, each row of Psi is made weakly diagonally dominant as far as
# the free parameters allow, using setup$structure: in a row whose fixed
# diagonal entry is below the row's sum of absolute off-diagonal entries,
# the free covariances are set to zero, and each free variance is raised to
# the largest such sum over its rows. A weakly diagonally dominant Psi with
# nonnegative diagonal is positive semidefinite, so Sigma is then positive
# definite unless a fixed diagonal entry is below the fixed off-diagonal
# entries in its row. A free error variance is then doubled until Sigma is
# positive definite, which holds for a large enough error variance. Returns
# NULL only if the error variance is fixed and no feasible point was found.
.feasible_start <- function(psi, fixed, setup) {
  if (is.finite(.ll_value(psi, setup$ll))) return(psi)
  r       <- setup$r
  st      <- setup$structure
  is_free <- !(seq_len(r) %in% fixed)
  if (r > 1) {
    re  <- seq_len(r - 1)
    cov <- re[!st$is_var[re]]                  # off the diagonal of Psi only
    off_sum <- function(psi) as.vector(st$n_off %*% abs(psi[re]))
    # The parameter on each diagonal entry of Psi, 0 if that entry is zero
    d_idx  <- as.vector(st$on_diag %*% re)
    d_val  <- ifelse(d_idx > 0, psi[pmax(d_idx, 1L)], 0)
    d_free <- d_idx > 0 & is_free[pmax(d_idx, 1L)]
    # Rows whose fixed diagonal entry is below the row's absolute
    # off-diagonal sum: set the free covariances in them to zero
    bad <- which(!d_free & d_val < off_sum(psi))
    if (length(bad) > 0) {
      in_bad <- colSums(st$n_off[bad, , drop = FALSE]) > 0
      psi[intersect(cov[in_bad[cov]], which(is_free))] <- 0
    }
    # Raise each free variance to the largest off-diagonal sum in its rows
    s_off <- off_sum(psi)
    for (v in re[st$is_var[re] & is_free[re]]) {
      psi[v] <- max(psi[v], s_off[st$on_diag[, v]])
    }
  }
  if (is_free[r] && psi[r] <= 0) psi[r] <- setup$psi_hat[r]
  for (k in 0:200) {
    if (is.finite(.ll_value(psi, setup$ll))) return(psi)
    if (!is_free[r]) return(NULL)
    psi[r] <- 2 * psi[r]
  }
  NULL
}

# Profile the nuisance parameters from a warm start and evaluate the signed
# statistic at the result; NULL if either step fails. origin_val fixes the
# sign of the rlrt root. The rlrt statistic reuses the maximized
# log-likelihood from maximize_loglik, so profiling and evaluation are one
# call; its signed root is on the scale of the signed score. ctx is the
# search context built in .ci_lmer_core; ctx$opt_args holds the one-step
# settings and user optimizer arguments.
.profile_stat <- function(theta_start, origin_val, ctx) {
  s <- ctx$setup
  ll_prof <- NA_real_
  if (length(ctx$opt_idx) > 0L) {
    opt <- tryCatch(
      do.call(maximize_loglik,
              c(list(start_val = theta_start, opt_idx = ctx$opt_idx,
                     ll = s$ll, expected = ctx$expected, check = FALSE),
                ctx$opt_args)),
      error = function(e) NULL
    )
    if (is.null(opt)) return(NULL)
    theta_start <- opt$arg
    ll_prof <- opt$value
  } else if (ctx$statistic == "rlrt") {
    ll_prof <- .ll_value(theta_start, s$ll)
  }
  stat <- if (ctx$statistic == "rlrt") {
    if (!is.finite(ll_prof)) return(NULL)
    lr <- 2 * (ctx$ll_max - ll_prof)
    if (!is.finite(lr)) return(NULL)
    sign(origin_val - theta_start[ctx$test_idx]) * sqrt(max(lr, 0))
  } else {
    tryCatch(
      as.numeric(score_stat(psi = theta_start, test_idx = ctx$test_idx,
                            ll = s$ll, expected = ctx$expected, signed = TRUE,
                            known_idx = ctx$known_idx, check = FALSE)),
      error = function(e) NA_real_
    )
  }
  if (is.na(stat)) return(NULL)
  list(theta = theta_start, stat = stat)
}


# Search outward from the origin in one direction until the signed statistic
# crosses the critical value, then interpolate to find the CI bound.
#
# direction: -1L to search left (lower bound), +1L to search right (upper bound)
# target:  +z_crit for lower bound, -z_crit for upper bound
#
# The search steps monotonically outward from the origin. Each accepted
# point warm-starts the next, so the nuisance optimum is tracked
# continuously. With accelerate = TRUE the step is a secant prediction of
# the crossing, overshot by 10% so the target is bracketed, and capped at
# twice the last accepted step. The nuisance warm start is extrapolated
# from the last two accepted points, with the plain warm start as fallback
# when the extrapolation is infeasible. A bracketed bound is refined by
# regula falsi until the bracket is no wider than step_size, the
# resolution of the fixed-step search. With accelerate = FALSE every step
# is step_size, which reproduces the fixed-step search.
#
# An infeasible warm start halves the step, up to 20 times. The growth cap
# then keeps later steps small near the feasibility boundary. The search
# radius is max_steps * step_size in both modes.
#
# ctx is the search context built in .ci_lmer_core and origin_eval the
# profiled evaluation at the search origin, shared by both directions.
.outward_bound <- function(ctx, origin_eval, direction, lower_clamp = -Inf) {

  test_idx   <- ctx$test_idx
  step_size  <- ctx$step_size
  max_steps  <- ctx$max_steps
  accelerate <- ctx$accelerate
  target     <- if (direction == -1L) ctx$z_crit else -ctx$z_crit
  growth     <- if (accelerate) 2 else 1
  max_radius <- max_steps * step_size
  # Profiling leaves the tested component unchanged, so this is the origin
  origin_val <- origin_eval$theta[test_idx]
  # Covariance parameters that keep their values
  fixed_psi <- c(test_idx, ctx$known_idx)

  .ll_val  <- function(theta) .ll_value(theta, ctx$setup$ll)
  .eval_at <- function(theta_start) .profile_stat(theta_start, origin_val, ctx)

  # Linear interpolation of the crossing between two bracketing points
  .interp <- function(x1, y1, x2, y2) x1 - (y1 - target) * (x2 - x1) / (y2 - y1)

  # Accepted-path state: current point and the one before it
  theta_cur <- origin_eval$theta
  val_cur   <- origin_val
  stat_cur  <- origin_eval$stat
  theta_prev <- NULL;   val_prev <- NA_real_;         stat_prev <- NA_real_
  step <- step_size

  for (i in seq_len(max_steps)) {

    # Propose next test value; clamp the lower search at a hard boundary
    # (e.g. 0 for a variance parameter). If the profile has not crossed
    # z_crit by the boundary, the boundary is returned as the CI bound.
    prop_val <- val_cur + direction * step
    hit_boundary <- direction == -1L && prop_val < lower_clamp
    if (hit_boundary) prop_val <- lower_clamp

    # Warm start: linear extrapolation of the accepted path (continuation),
    # falling back to the previous solution if the extrapolation is
    # infeasible. Known parameters are identical along the path, so the
    # extrapolation leaves them unchanged.
    theta_warm <- theta_cur
    theta_warm[test_idx] <- prop_val
    if (accelerate && !is.null(theta_prev) && val_cur != val_prev) {
      theta_pred <- theta_cur +
        ((prop_val - val_cur) / (val_cur - val_prev)) * (theta_cur - theta_prev)
      theta_pred[test_idx] <- prop_val
      if (is.finite(.ll_val(theta_pred))) theta_warm <- theta_pred
    }

    # If the warm start is infeasible, halve the step up to 20 times
    prop0 <- prop_val; step0 <- step; hit0 <- hit_boundary
    n_halve <- 0L
    feasible <- is.finite(.ll_val(theta_warm))
    while (!feasible && n_halve < 20L) {
      step <- step / 2
      prop_val <- val_cur + direction * step
      if (direction == -1L && prop_val < lower_clamp) {
        prop_val <- lower_clamp
        hit_boundary <- TRUE
      }
      theta_warm <- theta_cur
      theta_warm[test_idx] <- prop_val
      n_halve <- n_halve + 1L
      feasible <- is.finite(.ll_val(theta_warm))
    }
    if (!feasible) {
      # Construct a feasible start at the original proposal instead
      step <- step0; prop_val <- prop0; hit_boundary <- hit0
      theta_warm <- theta_cur
      theta_warm[test_idx] <- prop_val
      theta_warm <- .feasible_start(theta_warm, fixed_psi, ctx$setup)
      if (is.null(theta_warm)) {
        # Only possible with the error variance fixed, as at a lower clamp
        # of zero for the error variance
        if (direction == -1L && is.finite(lower_clamp)) return(lower_clamp)
        side <- if (direction == -1L) "lower" else "upper"
        warning("No feasible starting point was found for the next step of ",
                "the ", side, " search; returning ",
                if (direction == -1L) "-Inf." else "Inf.")
        return(if (direction == -1L) -Inf else Inf)
      }
    }

    res <- .eval_at(theta_warm)
    if (isTRUE(getOption("reconf.trace"))) {
      message(sprintf("propose %.4f (step %.4f, halve %d): stat %s | cur (%.4f, %.4f)",
                      prop_val, step, n_halve,
                      if (is.null(res)) "FAIL" else sprintf("%.4f", res$stat),
                      val_cur, stat_cur))
    }
    if (is.null(res)) {
      step <- step / 2
      next
    }

    # Crossing of the target: refine the bracket by regula falsi down to the
    # resolution of the fixed-step search, then interpolate.
    if ((stat_cur - target) * (res$stat - target) <= 0) {
      lo_val <- val_cur;  lo_stat <- stat_cur;  lo_theta <- theta_cur
      hi_val <- prop_val; hi_stat <- res$stat;  hi_theta <- res$theta
      for (k in seq_len(8L)) {
        if (abs(hi_val - lo_val) <= step_size) break
        # Alternate regula falsi with bisection: with a convex profile,
        # regula falsi alone can stagnate against a pinned endpoint.
        mid_val <- if (k %% 2L == 0L) 0.5 * (lo_val + hi_val)
                   else .interp(lo_val, lo_stat, hi_val, hi_stat)
        theta_mid <- if (abs(mid_val - lo_val) < abs(hi_val - mid_val))
          lo_theta else hi_theta
        theta_mid[test_idx] <- mid_val
        res_mid <- if (is.finite(.ll_val(theta_mid))) .eval_at(theta_mid) else NULL
        if (isTRUE(getOption("reconf.trace"))) {
          message(sprintf("  refine %d: mid %.4f stat %s | bracket [%.4f (%.4f), %.4f (%.4f)]",
                          k, mid_val,
                          if (is.null(res_mid)) "FAIL" else sprintf("%.4f", res_mid$stat),
                          lo_val, lo_stat, hi_val, hi_stat))
        }
        if (is.null(res_mid)) break
        if ((lo_stat - target) * (res_mid$stat - target) <= 0) {
          hi_val <- mid_val; hi_stat <- res_mid$stat; hi_theta <- res_mid$theta
        } else {
          lo_val <- mid_val; lo_stat <- res_mid$stat; lo_theta <- res_mid$theta
        }
      }
      return(.interp(lo_val, lo_stat, hi_val, hi_stat))
    }

    # Reached the lower clamp without crossing z_crit: return the clamp.
    if (hit_boundary) return(lower_clamp)

    # Advance the accepted path
    theta_prev <- theta_cur; val_prev <- val_cur; stat_prev <- stat_cur
    theta_cur <- res$theta;  val_cur <- prop_val; stat_cur <- res$stat

    # Choose the next step. Accelerated: secant prediction of the remaining
    # distance to the crossing, overshot by 10% to force a bracket, capped
    # at growth times the step just accepted (so halvings near the
    # feasibility boundary keep subsequent steps small). Falls back to pure
    # growth when the local slope is uninformative. Fixed: step_size always.
    # Chosen here, after acceptance, so failure-driven halvings above are
    # not overwritten.
    if (accelerate) {
      last_step <- abs(val_cur - val_prev)
      # Slope of the statistic per unit distance walked in `direction`;
      # the signed difference matters (direction * (val_cur - val_prev) is
      # the positive walked distance in both directions)
      slope <- (stat_cur - stat_prev) / (direction * (val_cur - val_prev))
      d_pred <- (target - stat_cur) / slope
      step <- if (is.finite(d_pred) && d_pred > 0) {
        min(1.1 * d_pred, growth * last_step)
      } else {
        growth * last_step
      }
      step <- max(step, 1e-3 * step_size)
    } else {
      step <- step_size
    }

    # Same search radius as the fixed-step search with max_steps steps
    if (direction * (val_cur - origin_val) > max_radius) break
  }

  side <- if (direction == -1L) "lower" else "upper"
  warning("CI ", side, " bound not found within the search radius (",
          format(max_radius, digits = 4), " = num_points * step_size from ",
          "the search origin): the statistic is ",
          format(stat_cur, digits = 3), " at ", format(val_cur, digits = 4),
          ", critical value ", format(target, digits = 3), ". The bound is ",
          "infinite if the statistic levels off short of the critical ",
          "value; otherwise increase num_points or step_size.")
  if (direction == -1L) -Inf else Inf
}
