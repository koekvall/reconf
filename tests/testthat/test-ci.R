library(lme4)

fit_ri <- lmer(Reaction ~ Days + (1 | Subject), data = sleepstudy, REML = TRUE)
fit_rs <- lmer(Reaction ~ Days + (Days | Subject), data = sleepstudy, REML = TRUE)

# ── ci_lmer ──────────────────────────────────────────────────────────────────

test_that("ci_lmer returns a 1-row reconf_ci matrix", {
  skip_on_cran()
  ci <- ci_lmer(fit_ri, test_idx = 1L)
  expect_true(is.matrix(ci))
  expect_s3_class(ci, "reconf_ci")
  expect_equal(nrow(ci), 1L)
  expect_equal(colnames(ci), c("estimate", "lower", "upper"))
  expect_true(is.numeric(ci))
})

test_that("ci_lmer: lower < MLE < upper, and estimate column matches VarCorr", {
  skip_on_cran()
  ci <- ci_lmer(fit_ri, test_idx = 1L)
  mle <- as.data.frame(VarCorr(fit_ri), order = "lower.tri")$vcov[1]
  expect_equal(as.numeric(ci[1, "estimate"]), mle)
  expect_lt(ci[1, "lower"], mle)
  expect_gt(ci[1, "upper"], mle)
})

test_that("ci_lmer: lower < upper", {
  skip_on_cran()
  ci <- ci_lmer(fit_ri, test_idx = 1L)
  expect_lt(ci[1, "lower"], ci[1, "upper"])
})

test_that("ci_lmer: lower bound is non-negative for variance parameter", {
  skip_on_cran()
  ci <- ci_lmer(fit_ri, test_idx = 1L)
  expect_gte(ci[1, "lower"], 0)
})

test_that("print and tidy methods work for confidence intervals", {
  skip_on_cran()
  ci <- ci_all_lmer(fit_ri)
  expect_output(print(ci), "95% score-based confidence intervals \\(REML\\)")
  td <- tidy(ci)
  expect_s3_class(td, "data.frame")
  expect_equal(names(td), c("term", "estimate", "conf.low", "conf.high"))
  expect_equal(nrow(td), nrow(ci))
})

test_that("score_test_all_lmer returns a named test table", {
  skip_on_cran()
  st <- score_test_all_lmer(fit_rs)
  expect_s3_class(st, "reconf_test")
  expect_equal(colnames(st), c("statistic", "df", "p.value"))
  expect_equal(nrow(st), 3L)  # two variances and one covariance
  expect_true(all(rownames(st) != ""))
  expect_true(all(st[, "p.value"] >= 0 & st[, "p.value"] <= 1))
  expect_true(all(is.finite(st[, "statistic"])))
  expect_output(print(st), "Score tests of zero covariance parameters")
  td <- tidy(st)
  expect_equal(names(td), c("term", "statistic", "df", "p.value"))
})

test_that("accelerated and fixed-step searches agree", {
  skip_on_cran()
  # Both searches resolve the crossing to within one step_size (SE/40), so
  # bounds must agree to that resolution. Guards the secant step logic and
  # the regula-falsi refinement against regressions such as sign errors in
  # the slope (which silently degrade one search direction only).
  ci_a <- ci_all_lmer(fit_ri, accelerate = TRUE)
  ci_f <- ci_all_lmer(fit_ri, accelerate = FALSE)
  expect_equal(ci_a, ci_f, tolerance = 1e-3)
})

# ── ci_all_lmer ──────────────────────────────────────────────────────────────

test_that("ci_all_lmer returns matrix with correct dimensions", {
  skip_on_cran()
  ci <- ci_all_lmer(fit_ri)
  expect_true(is.matrix(ci))
  expect_equal(ncol(ci), 3L)
  expect_equal(colnames(ci), c("estimate", "lower", "upper"))
  # random intercept model: 1 RE variance + error variance
  expect_equal(nrow(ci), 2L)
})

test_that("ci_all_lmer: all lower < upper", {
  skip_on_cran()
  ci <- ci_all_lmer(fit_rs)
  expect_true(all(ci[, "lower"] < ci[, "upper"]))
})

test_that("ci_all_lmer: all MLEs inside CIs", {
  skip_on_cran()
  ci <- ci_all_lmer(fit_rs)
  vc <- as.data.frame(VarCorr(fit_rs), order = "lower.tri")
  mles <- vc$vcov
  expect_true(all(ci[, "lower"] <= mles & mles <= ci[, "upper"]))
})

test_that("ci_all_lmer: correct number of rows for random slope model", {
  skip_on_cran()
  ci <- ci_all_lmer(fit_rs)
  # intercept var, covariance, slope var, error var = 4 parameters
  expect_equal(nrow(ci), 4L)
})

test_that("ci_all_lmer respects test_idx argument", {
  skip_on_cran()
  ci <- ci_all_lmer(fit_rs, test_idx = 1L)
  expect_equal(nrow(ci), 1L)
})

test_that("onestep CIs agree with full profiling on sleepstudy", {
  skip_on_cran()
  ci_full <- ci_all_lmer(fit_rs)
  ci_one  <- ci_all_lmer(fit_rs, onestep = TRUE)
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
  ci  <- suppressWarnings(ci_lmer(fit, test_idx = 1L))
  expect_gte(ci[1, "lower"], 0)
})

test_that("nonneg = FALSE allows negative lower bound for a variance", {
  skip_on_cran()
  fit     <- make_tiny_var_fit()
  ci_on   <- suppressWarnings(ci_lmer(fit, test_idx = 1L, nonneg = TRUE))
  ci_off  <- suppressWarnings(ci_lmer(fit, test_idx = 1L, nonneg = FALSE))
  # Either the raw bound is negative (so clamp changes the answer) or both
  # agree when the raw bound happens to be nonneg. Both outcomes are valid;
  # the test enforces that turning off the clamp does not increase the lower.
  expect_lte(ci_off[1, "lower"], ci_on[1, "lower"] + 1e-8)
})

test_that("nonneg does not clamp covariance parameters", {
  skip_on_cran()
  # fit_rs has a random-intercept/slope covariance at index 2.
  ci_on  <- ci_all_lmer(fit_rs, nonneg = TRUE)
  ci_off <- ci_all_lmer(fit_rs, nonneg = FALSE)
  expect_equal(ci_on[2, ], ci_off[2, ])
})

test_that("all variance lower bounds are nonneg under default", {
  skip_on_cran()
  ci <- ci_all_lmer(fit_rs)
  vc <- as.data.frame(VarCorr(fit_rs), order = "lower.tri")
  is_var <- is.na(vc$var2)
  expect_true(all(ci[is_var, "lower"] >= 0))
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
    b <- if (reml) NULL else as.vector(fixef(fitw))
    ll <- reconf:::loglikelihood(psi = psi, b = b, Y = m$Y, X = m$X, Z = m$Z,
                                 Hlist = reconf:::get_Hlist_lmer(fitw),
                                 REML = reml, get_inf = FALSE)
    # Value matches logLik() up to the weight-transformation Jacobian
    expect_equal(ll$value + 0.5 * sum(log(w)), as.numeric(logLik(fitw)),
                 tolerance = 1e-6)
    # Score vanishes at the weighted estimates
    expect_lt(max(abs(ll$score)), 1e-2)
  }
  # CI machinery runs on a weighted fit and brackets the estimate
  fitw <- lmer(Reaction ~ Days + (1 | Subject), data = sleepstudy, weights = w)
  ci <- ci_lmer(fitw, test_idx = 1L)
  mle <- as.data.frame(VarCorr(fitw), order = "lower.tri")$vcov[1]
  expect_lt(ci[1, "lower"], mle)
  expect_gt(ci[1, "upper"], mle)
})

test_that("offsets are subtracted before the analysis", {
  skip_on_cran()
  set.seed(8)
  off <- runif(nrow(sleepstudy), -10, 10)
  fito <- lmer(Reaction ~ Days + (1 | Subject), data = sleepstudy,
               REML = TRUE, offset = off)
  m <- reconf:::.lmer_matrices(fito)
  psi <- reconf:::get_psi_hat_lmer(fito)
  ll <- reconf:::loglikelihood(psi = psi, Y = m$Y, X = m$X, Z = m$Z,
                               Hlist = reconf:::get_Hlist_lmer(fito),
                               REML = TRUE, get_inf = FALSE)
  expect_equal(ll$value, as.numeric(logLik(fito)), tolerance = 1e-6)
  expect_lt(max(abs(ll$score)), 1e-2)
})
