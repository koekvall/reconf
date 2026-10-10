library(lme4)

fit_ri <- lmer(Reaction ~ Days + (1 | Subject), data = sleepstudy, REML = TRUE)
fit_rs <- lmer(Reaction ~ Days + (Days | Subject), data = sleepstudy, REML = TRUE)

# Fixed-step search, every step of length step_size, through the internal
# core; vc_ci always uses the accelerated search
ci_fixed_step <- function(fit, statistic = "score") {
  setup <- vc_model(fit)
  do.call(rbind, lapply(seq_len(setup$r), function(j)
    reconf:::.ci_lmer_core(setup, test_idx = j, accelerate = FALSE,
                           statistic = statistic)))
}

# ── vc_ci ────────────────────────────────────────────────────────────────────

test_that("vc_ci returns a 1-row vc_ci matrix", {
  skip_on_cran()
  ci <- vc_ci(fit_ri, parm = 1L)
  expect_true(is.matrix(ci))
  expect_s3_class(ci, "vc_ci")
  expect_equal(nrow(ci), 1L)
  expect_equal(colnames(ci), c("estimate", "lower", "upper"))
  expect_true(is.numeric(ci))
})

test_that("vc_ci: lower < MLE < upper, and estimate column matches VarCorr", {
  skip_on_cran()
  ci <- vc_ci(fit_ri, parm = 1L)
  mle <- as.data.frame(VarCorr(fit_ri), order = "lower.tri")$vcov[1]
  expect_equal(as.numeric(ci[1, "estimate"]), mle)
  expect_lt(ci[1, "lower"], mle)
  expect_gt(ci[1, "upper"], mle)
})

test_that("vc_ci: lower < upper", {
  skip_on_cran()
  ci <- vc_ci(fit_ri, parm = 1L)
  expect_lt(ci[1, "lower"], ci[1, "upper"])
})

test_that("vc_ci: lower bound is non-negative for variance parameter", {
  skip_on_cran()
  ci <- vc_ci(fit_ri, parm = 1L)
  expect_gte(ci[1, "lower"], 0)
})

test_that("print and tidy methods work for confidence intervals", {
  skip_on_cran()
  ci <- vc_ci(fit_ri)
  expect_output(print(ci), paste("95% confidence intervals for",
                                 "variance-covariance parameters \\(score, REML\\)"))
  ci_r <- vc_ci(fit_ri, parm = 1L, statistic = "rlrt")
  expect_output(print(ci_r), "\\(rlrt, REML\\)")
  td <- tidy(ci)
  expect_s3_class(td, "data.frame")
  expect_equal(names(td), c("term", "estimate", "conf.low", "conf.high"))
  expect_equal(nrow(td), nrow(ci))
})

test_that("single-parameter score tests agree with the intervals", {
  skip_on_cran()
  # The test of psi_j = 0 uses the statistic the interval inverts, so it
  # rejects at level 0.05 iff 0 is outside the 95% interval
  ci <- vc_ci(fit_rs, nonneg = FALSE)
  for (j in 1:3) {
    st <- vc_test(fit_rs, parm = j)
    expect_identical(names(st$null.value), rownames(ci)[j])
    expect_identical(st$p.value < 0.05,
                     ci[j, "lower"] > 0 || ci[j, "upper"] < 0)
  }
})

test_that("accelerated and fixed-step searches agree", {
  skip_on_cran()
  # Both searches resolve the crossing to within one step_size (SE/40), so
  # bounds must agree to that resolution. Guards the secant step and the
  # regula-falsi refinement against sign errors in the slope, which degrade
  # one search direction only.
  ci_a <- vc_ci(fit_ri)
  ci_f <- ci_fixed_step(fit_ri)
  expect_equal(as.numeric(ci_a), as.numeric(ci_f), tolerance = 1e-3)
})

# ── several parameters ───────────────────────────────────────────────────────

test_that("a subset of parameters, by name or index, matches the full call", {
  skip_on_cran()
  ci_all <- vc_ci(fit_rs)
  ci_sub <- vc_ci(fit_rs, parm = c("var_Days|Subject", "var_Residual"))
  expect_identical(rownames(ci_sub), c("var_Days|Subject", "var_Residual"))
  expect_equal(unclass(ci_sub)[, ], unclass(ci_all)[3:4, ],
               ignore_attr = TRUE)
  expect_equal(unclass(vc_ci(fit_rs, parm = c(3L, 4L)))[, ],
               unclass(ci_sub)[, ], ignore_attr = TRUE)
  expect_error(vc_ci(fit_rs, parm = "var_Hours|Subject"), "unknown parm")
})

test_that("vc_ci returns matrix with correct dimensions", {
  skip_on_cran()
  ci <- vc_ci(fit_ri)
  expect_true(is.matrix(ci))
  expect_equal(ncol(ci), 3L)
  expect_equal(colnames(ci), c("estimate", "lower", "upper"))
  # random intercept model: 1 RE variance + error variance
  expect_equal(nrow(ci), 2L)
})

test_that("vc_ci: all lower < upper", {
  skip_on_cran()
  ci <- vc_ci(fit_rs)
  expect_true(all(ci[, "lower"] < ci[, "upper"]))
})

test_that("vc_ci: all MLEs inside CIs", {
  skip_on_cran()
  ci <- vc_ci(fit_rs)
  vc <- as.data.frame(VarCorr(fit_rs), order = "lower.tri")
  mles <- vc$vcov
  expect_true(all(ci[, "lower"] <= mles & mles <= ci[, "upper"]))
})

test_that("vc_ci: correct number of rows for random slope model", {
  skip_on_cran()
  ci <- vc_ci(fit_rs)
  # intercept var, covariance, slope var, error var = 4 parameters
  expect_equal(nrow(ci), 4L)
  expect_identical(rownames(ci),
                   c("var_(Intercept)|Subject", "cov_Days.(Intercept)|Subject",
                     "var_Days|Subject", "var_Residual"))
})

test_that("vc_ci respects the parm argument", {
  skip_on_cran()
  ci <- vc_ci(fit_rs, parm = 1L)
  expect_equal(nrow(ci), 1L)
})

test_that("onestep CIs agree with full profiling on sleepstudy", {
  skip_on_cran()
  ci_full <- vc_ci(fit_rs)
  ci_one  <- vc_ci(fit_rs, onestep = TRUE)
  expect_equal(ci_one, ci_full, tolerance = 0.05)
})

# ── nonneg clamping ──────────────────────────────────────────────────────────
# Build a small model whose random-effect variance MLE is near zero, so that
# the unclamped score CI crosses zero on the lower side. The default nonneg
# behavior should truncate the lower bound at 0; disabling nonneg should
# recover a negative lower bound.

make_tiny_var_fit <- function(seed = 11L) {
  set.seed(seed)
  n_grp <- 6L
  n_obs <- 4L
  grp   <- factor(rep(seq_len(n_grp), each = n_obs))
  # Effectively no random effect: small signal relative to residual
  b     <- rnorm(n_grp, sd = 0.05)
  y     <- b[grp] + rnorm(n_grp * n_obs, sd = 1)
  dat   <- data.frame(y = y, grp = grp)
  suppressMessages(suppressWarnings(
    lmer(y ~ 1 + (1 | grp), data = dat, REML = TRUE)
  ))
}

test_that("nonneg = TRUE clamps the variance lower bound at 0", {
  skip_on_cran()
  fit <- make_tiny_var_fit()
  ci  <- suppressWarnings(vc_ci(fit, parm = 1L))
  expect_gte(ci[1, "lower"], 0)
})

test_that("nonneg = FALSE allows negative lower bound for a variance", {
  skip_on_cran()
  fit     <- make_tiny_var_fit()
  ci_on   <- suppressWarnings(vc_ci(fit, parm = 1L, nonneg = TRUE))
  ci_off  <- suppressWarnings(vc_ci(fit, parm = 1L, nonneg = FALSE))
  # The raw lower bound need not be negative for this seed; the test
  # enforces only that turning off the clamp does not increase it.
  expect_lte(ci_off[1, "lower"], ci_on[1, "lower"] + 1e-8)
})

test_that("nonneg does not clamp covariance parameters", {
  skip_on_cran()
  # fit_rs has a random-intercept/slope covariance at index 2.
  ci_on  <- vc_ci(fit_rs, nonneg = TRUE)
  ci_off <- vc_ci(fit_rs, nonneg = FALSE)
  expect_equal(ci_on[2, ], ci_off[2, ])
})

test_that("all variance lower bounds are nonneg under default", {
  skip_on_cran()
  ci <- vc_ci(fit_rs)
  vc <- as.data.frame(VarCorr(fit_rs), order = "lower.tri")
  is_var <- is.na(vc$var2)
  expect_true(all(ci[is_var, "lower"] >= 0))
})

# ── statistic = "rlrt" ───────────────────────────────────────────────────────

test_that("rlrt intervals are finite and bracket the estimates", {
  skip_on_cran()
  ci <- vc_ci(fit_rs, statistic = "rlrt")
  expect_true(all(is.finite(ci)))
  expect_true(all(ci[, "lower"] <= ci[, "estimate"]))
  expect_true(all(ci[, "estimate"] <= ci[, "upper"]))
})

test_that("rlrt and score intervals roughly agree on sleepstudy", {
  skip_on_cran()
  # The statistics are first-order equivalent and all estimates in this fit
  # are interior. The tolerance is half the score interval's width because
  # rlrt upper bounds for variances are tighter (1562 vs 2331 for the
  # intercept variance); sign errors or a wrong reference maximum shift
  # bounds by far more.
  ci_r <- vc_ci(fit_rs, statistic = "rlrt")
  ci_s <- vc_ci(fit_rs, statistic = "score")
  width <- ci_s[, "upper"] - ci_s[, "lower"]
  expect_true(all(abs(ci_r - ci_s) / width < 0.5))
})

test_that("rlrt intervals do not depend on which other parameters are requested", {
  skip_on_cran()
  # One reference maximum is shared across the requested parameters
  ci_all <- vc_ci(fit_ri, statistic = "rlrt")
  ci_one <- vc_ci(fit_ri, parm = 2L, statistic = "rlrt")
  expect_equal(as.numeric(ci_all[2, ]), as.numeric(ci_one[1, ]),
               tolerance = 1e-6)
})

test_that("rlrt: accelerated and fixed-step searches agree", {
  skip_on_cran()
  ci_a <- vc_ci(fit_ri, statistic = "rlrt")
  ci_f <- ci_fixed_step(fit_ri, statistic = "rlrt")
  expect_equal(as.numeric(ci_a), as.numeric(ci_f), tolerance = 1e-3)
})

test_that("rlrt respects the nonneg clamp", {
  skip_on_cran()
  fit    <- make_tiny_var_fit()
  ci_on  <- suppressWarnings(vc_ci(fit, parm = 1L, statistic = "rlrt"))
  ci_off <- suppressWarnings(vc_ci(fit, parm = 1L, statistic = "rlrt",
                                     nonneg = FALSE))
  expect_gte(ci_on[1, "lower"], 0)
  expect_lte(ci_off[1, "lower"], ci_on[1, "lower"] + 1e-8)
})

test_that("rlrt rejects onestep = TRUE", {
  skip_on_cran()
  expect_error(vc_ci(fit_ri, parm = 1L, statistic = "rlrt",
                       onestep = TRUE),
               "onestep")
})

# ── boundary behavior: negative extended-set estimate ────────────────────────
# Partial group-centering (c = 0.5) induces negative within-group correlation
# with implied psi1/psi2 around -0.19, inside the feasibility bound
# psi1 > -psi2/n_obs, so the extended-set estimate of the group variance is
# negative while lme4 reports 0. The rlrt interval then lies
# entirely below zero and its intersection with [0, Inf) is empty.

make_neg_icc_fit <- function(seed = 3L, n_grp = 40L, n_obs = 4L) {
  set.seed(seed)
  grp <- factor(rep(seq_len(n_grp), each = n_obs))
  e   <- rnorm(n_grp * n_obs)
  y   <- e - 0.5 * ave(e, grp)
  suppressMessages(suppressWarnings(
    lmer(y ~ 1 + (1 | grp), data = data.frame(y = y, grp = grp), REML = TRUE)
  ))
}

test_that("rlrt: empty nonneg intersection gives NA bounds with a warning", {
  skip_on_cran()
  fit <- make_neg_icc_fit()
  expect_warning(
    ci <- vc_ci(fit, parm = 1L, statistic = "rlrt"),
    "no nonnegative values")
  expect_true(is.na(ci[1, "lower"]))
  expect_true(is.na(ci[1, "upper"]))
})

test_that("rlrt: unclamped boundary interval lies entirely below zero", {
  skip_on_cran()
  fit <- make_neg_icc_fit()
  ci <- suppressWarnings(
    vc_ci(fit, parm = 1L, statistic = "rlrt", nonneg = FALSE))
  expect_lt(ci[1, "lower"], ci[1, "upper"])
  expect_lt(ci[1, "upper"], 0)
})

test_that("score: origin outside the confidence set restarts and reports NA", {
  skip_on_cran()
  # The signed score at the estimate exceeds the critical value, so the
  # search restarts from the extended-set maximizer and the nonneg
  # intersection is empty.
  fit <- make_neg_icc_fit()
  w <- character()
  ci <- withCallingHandlers(
    vc_ci(fit, parm = 1L, statistic = "score"),
    warning = function(cnd) {
      w <<- c(w, conditionMessage(cnd))
      invokeRestart("muffleWarning")
    })
  expect_true(any(grepl("beyond the critical value", w)))
  expect_true(any(grepl("no nonnegative values", w)))
  expect_true(is.na(ci[1, "lower"]) && is.na(ci[1, "upper"]))
})

test_that("score: unclamped boundary interval lies entirely below zero", {
  skip_on_cran()
  fit <- make_neg_icc_fit()
  ci <- suppressWarnings(
    vc_ci(fit, parm = 1L, statistic = "score", nonneg = FALSE))
  expect_lt(ci[1, "lower"], ci[1, "upper"])
  expect_lt(ci[1, "upper"], 0)
})

# ── no lme4 objects past the setup ───────────────────────────────────────────

test_that("intervals from a model built from matrices match the lme4 path", {
  skip_on_cran()
  m <- reconf:::.lmer_matrices(fit_rs)
  setup <- reconf:::.make_setup(Y = m$Y, X = m$X, Z = m$Z,
                                Hlist = reconf:::get_Hlist_lmer(fit_rs),
                                REML = TRUE,
                                psi_hat = reconf:::get_psi_hat_lmer(fit_rs))
  # The variances, found from the structure matrices alone
  vc <- as.data.frame(VarCorr(fit_rs), order = "lower.tri")
  expect_identical(setup$structure$is_var, is.na(vc$var2))
  ci_m <- do.call(rbind, lapply(seq_len(setup$r), function(j)
    reconf:::.ci_lmer_core(setup, test_idx = j)))
  expect_identical(rownames(ci_m), paste0("psi", 1:4))
  expect_identical(as.numeric(ci_m), as.numeric(vc_ci(fit_rs)))
})

# ── weights and offsets ──────────────────────────────────────────────────────

test_that("prior weights are handled exactly via the W^(1/2) transformation", {
  skip_on_cran()
  set.seed(7)
  w <- runif(nrow(sleepstudy), 0.2, 5)
  for (reml in c(TRUE, FALSE)) {
    fitw <- lmer(Reaction ~ Days + (Days | Subject), data = sleepstudy,
                 REML = reml, weights = w)
    m <- reconf:::.lmer_matrices(fitw)
    psi <- reconf:::get_psi_hat_lmer(fitw)
    ll <- reconf:::make_loglik(Y = m$Y, X = m$X, Z = m$Z,
                               Hlist = reconf:::get_Hlist_lmer(fitw),
                               REML = reml)(psi, get_inf = FALSE)
    # Value matches logLik() up to the weight-transformation Jacobian
    expect_equal(ll$value + 0.5 * sum(log(w)), as.numeric(logLik(fitw)),
                 tolerance = 1e-6)
    # Score vanishes at the weighted estimates
    expect_lt(max(abs(ll$score)), 1e-2)
  }
  # CI machinery runs on a weighted fit and brackets the estimate
  fitw <- lmer(Reaction ~ Days + (1 | Subject), data = sleepstudy, weights = w)
  ci <- vc_ci(fitw, parm = 1L)
  mle <- as.data.frame(VarCorr(fitw), order = "lower.tri")$vcov[1]
  expect_lt(ci[1, "lower"], mle)
  expect_gt(ci[1, "upper"], mle)
})

# ── argument checks ──────────────────────────────────────────────────────────

test_that("vc_ci rejects arguments the search would otherwise swallow", {
  skip_on_cran()
  # Inside the search an error in a profile evaluation counts as an
  # infeasible point, so a misnamed or invalid argument must be rejected
  # up front rather than turn into an infinite interval
  expect_error(vc_ci(fit_ri, parm = 1L, check = FALSE), "unused argument")
  expect_error(vc_ci(fit_ri, parm = 1L, minimize = FALSE), "unused argument")
  expect_error(vc_ci(fit_ri, parm = 1L, expected = "yes"), "single logical")
  expect_error(vc_ci(fit_ri, parm = 1L, nonneg = NA), "single logical")
  expect_error(vc_ci(fit_ri, parm = 1L, step_size = 0), "step_size")
  expect_error(vc_ci(fit_ri, parm = 1L, step_size = c(1, 2)), "step_size")
  expect_error(vc_ci(fit_ri, parm = 1L, known = 1L), "should not overlap")
  expect_error(vc_ci(fit_ri, known = 1:2), "every parameter is in known")
  # Allowed optimizer settings pass through
  ci <- vc_ci(fit_ri, parm = 1L, iterlim = 50L, warn_nonconv = FALSE)
  expect_true(all(is.finite(ci)))
})

test_that("the default parm excludes known parameters", {
  skip_on_cran()
  ci <- vc_ci(fit_rs, known = "cov_Days.(Intercept)|Subject")
  expect_identical(rownames(ci), c("var_(Intercept)|Subject",
                                   "var_Days|Subject", "var_Residual"))
  expect_equal(unclass(ci)[1, ],
               unclass(vc_ci(fit_rs, parm = 1L, known = 2L))[1, ])
})

test_that("a statistic that levels off gives an infinite bound and a warning", {
  skip_on_cran()
  # The sample factor has 6 levels: as its variance grows, the signed score
  # statistic tends to -sqrt((6 - 1) / 2) = -1.58, short of -1.96, so the
  # 95% upper bound is infinite whatever the search distance
  fit <- lmer(diameter ~ (1 | plate) + (1 | sample), Penicillin)
  expect_warning(ci <- vc_ci(fit, parm = "var_(Intercept)|sample"),
                 "levels off")
  expect_identical(ci[1, "upper"], Inf)
  expect_true(is.finite(ci[1, "lower"]))
})

# ── offsets ──────────────────────────────────────────────────────────────────

test_that("offsets are subtracted before the analysis", {
  skip_on_cran()
  set.seed(8)
  off <- runif(nrow(sleepstudy), -10, 10)
  fito <- lmer(Reaction ~ Days + (1 | Subject), data = sleepstudy,
               REML = TRUE, offset = off)
  m <- reconf:::.lmer_matrices(fito)
  psi <- reconf:::get_psi_hat_lmer(fito)
  ll <- reconf:::make_loglik(Y = m$Y, X = m$X, Z = m$Z,
                             Hlist = reconf:::get_Hlist_lmer(fito),
                             REML = TRUE)(psi, get_inf = FALSE)
  expect_equal(ll$value, as.numeric(logLik(fito)), tolerance = 1e-6)
  expect_lt(max(abs(ll$score)), 1e-2)
})
