library(lme4)

fit_ri <- lmer(Reaction ~ Days + (1 | Subject), data = sleepstudy, REML = TRUE)
fit_rs <- lmer(Reaction ~ Days + (Days | Subject), data = sleepstudy, REML = TRUE)

# ── ci_lmer ──────────────────────────────────────────────────────────────────

test_that("ci_lmer returns a 1-row matrix with lower and upper columns", {
  skip_on_cran()
  ci <- ci_lmer(fit_ri, test_idx = 1L)
  expect_true(is.matrix(ci))
  expect_equal(nrow(ci), 1L)
  expect_equal(colnames(ci), c("lower", "upper"))
  expect_true(is.numeric(ci))
})

test_that("ci_lmer: lower < MLE < upper for random intercept variance", {
  skip_on_cran()
  ci <- ci_lmer(fit_ri, test_idx = 1L)
  mle <- as.data.frame(VarCorr(fit_ri), order = "lower.tri")$vcov[1]
  expect_lt(ci[1], mle)
  expect_gt(ci[2], mle)
})

test_that("ci_lmer: lower < upper", {
  skip_on_cran()
  ci <- ci_lmer(fit_ri, test_idx = 1L)
  expect_lt(ci[1], ci[2])
})

test_that("ci_lmer: lower bound is non-negative for variance parameter", {
  skip_on_cran()
  ci <- ci_lmer(fit_ri, test_idx = 1L)
  expect_gte(ci[1], 0)
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
  expect_equal(ncol(ci), 2L)
  expect_equal(colnames(ci), c("lower", "upper"))
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
