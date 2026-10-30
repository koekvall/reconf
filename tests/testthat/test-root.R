library(lme4)

fit_rs <- lmer(Reaction ~ Days + (Days | Subject), data = sleepstudy, REML = TRUE)

# Correlation 0.999, 20 groups of 10: at the true covariance the nuisance
# maximization stops away from a root
set.seed(34)
g <- factor(rep(1:20, each = 10))
x <- runif(200, -1, 1)
u <- matrix(rnorm(40), 20) %*% chol(matrix(c(2.5, 2.4975, 2.4975, 2.5), 2))
y <- u[g, 1] + u[g, 2] * x + rnorm(200)
fit_miss <- suppressMessages(lmer(y ~ x + (x | g)))

# ── root check ──────────────────────────────────────────────────────────────

# The misses below depend on the optimizer path, which can differ across
# platforms, so they are not run on CRAN
test_that("maximize_loglik flags whether it found a root", {
  m <- vc_model(fit_rs)
  expect_true(reconf:::maximize_loglik(m$psi_hat, 1:4, m$ll)$root)

  skip_on_cran()
  m <- vc_model(fit_miss)
  psi <- m$psi_hat
  psi[2] <- 2.4975
  psi <- reconf:::.feasible_start(psi, 2L, m)
  o <- reconf:::maximize_loglik(psi, c(1L, 3L, 4L), m$ll, warn_nonconv = FALSE)
  expect_false(o$root)
  expect_gt(o$score_max, 1)
})

test_that("vc_test warns when no root is found", {
  w <- capture_warnings(vc_test(fit_rs))
  expect_false(any(grepl("root of their score", w)))

  skip_on_cran()
  expect_warning(vc_test(fit_miss, parm = 2L, null_value = 2.4975),
                 "root of their score")
})

test_that("vc_ci warns once per interval when no root is found", {
  skip_on_cran()
  w <- capture_warnings(vc_ci(fit_miss, parm = 2L, num_points = 20))
  expect_equal(sum(grepl("root of their score", w)), 1L)
  w <- capture_warnings(vc_ci(fit_miss, parm = 2L, num_points = 20,
                              onestep = TRUE))
  expect_false(any(grepl("root of their score", w)))
})

test_that("the extended-set maximization warns when no root is found", {
  w <- capture_warnings(reconf:::.rlrt_reference(vc_model(fit_rs), NULL))
  expect_false(any(grepl("extended parameter set", w)))

  skip_on_cran()
  expect_warning(reconf:::.rlrt_reference(vc_model(fit_miss), NULL),
                 "extended parameter set")
})
