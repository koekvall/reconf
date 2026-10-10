library(lme4)

fit_ri <- lmer(Reaction ~ Days + (1 | Subject), data = sleepstudy, REML = TRUE)
fit_rs <- lmer(Reaction ~ Days + (Days | Subject), data = sleepstudy, REML = TRUE)

# ── vc_test return structure ────────────────────────────────────────────────

test_that("vc_test returns an htest object", {
  res <- vc_test(fit_ri)
  expect_s3_class(res, "htest")
  expect_named(res$statistic, "S")
  expect_named(res$parameter, "df")
  expect_equal(res$null.value, c("var_(Intercept)|Subject" = 0))
  expect_output(print(res), "Score test \\(REML\\)")
})

test_that("p-value is in [0, 1]", {
  res <- vc_test(fit_ri)
  expect_gte(res$p.value, 0)
  expect_lte(res$p.value, 1)
})

test_that("df equals the number of tested parameters", {
  res <- vc_test(fit_ri, parm = 1L)
  expect_equal(res$parameter[["df"]], 1)
})

# ── score at MLE should be near zero ────────────────────────────────────────

test_that("score stat near zero when the null value equals the estimate", {
  psi_hat <- reconf:::get_psi_hat_lmer(fit_ri)
  res <- vc_test(fit_ri, parm = 1L, null_value = psi_hat[1])
  expect_lt(abs(res$statistic[["S"]]), 0.01)
})

test_that("indices are into the covariance parameters for ML fits too", {
  fit_ml <- lmer(Reaction ~ Days + (Days | Subject), data = sleepstudy,
                 REML = FALSE)
  res <- vc_test(fit_ml, parm = 3L)
  expect_identical(names(res$null.value), "var_Days|Subject")
  expect_output(print(res), "Score test \\(ML\\)")
})

test_that("a covariance null outside the estimated variances is feasible", {
  # null_value^2 exceeds the product of the estimated variances, so the
  # start is infeasible until the nuisance variances are scaled up
  psi_hat <- reconf:::get_psi_hat_lmer(fit_rs)
  c_null <- 3 * sqrt(psi_hat[1] * psi_hat[3])
  res <- vc_test(fit_rs, parm = 2L, null_value = c_null)
  expect_true(is.finite(res$statistic[["S"]]))
  expect_lt(res$p.value, 0.05)
})

test_that("known and argument checks", {
  res <- vc_test(fit_rs, parm = 3L, known = 1L)
  expect_true(is.finite(res$statistic[["S"]]))
  expect_error(vc_test(fit_rs, parm = 3L, known = 3L),
               "should not overlap")
  expect_error(vc_test(fit_rs, parm = 1:2, null_value = 1:3),
               "null_value")
  expect_error(vc_test(fit_rs, parm = 4L),
               "error variance must be positive")
  # The default parm excludes known parameters
  res <- vc_test(fit_rs, known = 2L)
  expect_identical(names(res$null.value),
                   c("var_(Intercept)|Subject", "var_Days|Subject"))
  expect_error(vc_test(fit_rs, known = 1:3),
               "every random-effect parameter is in known")
  # Misnamed or invalid arguments are rejected, not passed to the optimizer
  expect_error(vc_test(fit_rs, foo = 1), "unused argument")
  expect_error(vc_test(fit_rs, expected = "yes"), "single logical")
})

test_that("score test agrees across computational methods", {
  res_q <- vc_test(fit_ri, method = "q_side")
  res_n <- vc_test(fit_ri, method = "n_side")
  expect_equal(res_n$statistic, res_q$statistic, tolerance = 1e-6)
})

# ── rejects zero variance when between-subject variability is strong ─────────

test_that("random intercept is significant in sleepstudy", {
  res <- vc_test(fit_ri)
  expect_lt(res$p.value, 0.001)
})

# ── multiple random effects ──────────────────────────────────────────────────

test_that("works with random slope model, testing all RE params", {
  res <- vc_test(fit_rs)
  expect_s3_class(res, "htest")
  expect_equal(res$parameter[["df"]], 3)  # intercept var, covariance, slope var
})

# ── Efficient information: factorized Schur complement ───────────────────────

test_that("efficient information matches the explicit Schur subtraction", {
  # score_stat computes I_tt - I_tn I_nn^{-1} I_nt from one Cholesky of the
  # joint information; the explicit subtraction is the oracle here
  eff_sub <- function(inf, test_idx, known_idx = NULL) {
    nuis <- setdiff(seq_len(nrow(inf)), c(test_idx, known_idx))
    inf[test_idx, test_idx, drop = FALSE] -
      inf[test_idx, nuis, drop = FALSE] %*%
      solve(inf[nuis, nuis], inf[nuis, test_idx, drop = FALSE])
  }

  psi_hat <- reconf:::get_psi_hat_lmer(fit_rs)
  Y <- getME(fit_rs, "y"); X <- getME(fit_rs, "X"); Z <- getME(fit_rs, "Z")
  Hlist <- reconf:::get_Hlist_lmer(fit_rs)
  ll <- reconf:::make_loglik(Y = Y, X = X, Z = Z, Hlist = Hlist, REML = TRUE)

  # REML at the estimate and at a non-stationary point; single and joint
  # tests, and with a known parameter excluded from the nuisances
  for (psi in list(psi_hat, psi_hat * c(1.3, 0.8, 1.1, 0.9))) {
    inf_full <- ll(psi, get_val = FALSE)$inf_mat
    for (cfg in list(list(t = 1L, k = NULL), list(t = c(1L, 3L), k = NULL),
                     list(t = 1L, k = 2L))) {
      st <- reconf:::score_stat(psi = psi, test_idx = cfg$t,
                                known_idx = cfg$k, ll = ll, signed = TRUE)
      expect_equal(attr(st, "info"), eff_sub(inf_full, cfg$t, cfg$k),
                   tolerance = 1e-8, ignore_attr = TRUE)
    }
  }

  # ML, with the fixed effects profiled out of the likelihood, with expected
  # and observed information
  fit_ml <- lmer(Reaction ~ Days + (Days | Subject), data = sleepstudy,
                 REML = FALSE)
  ll_ml <- reconf:::make_loglik(Y = Y, X = X, Z = Z, Hlist = Hlist,
                                REML = FALSE)
  psi_ml <- reconf:::get_psi_hat_lmer(fit_ml) * c(1.3, 0.8, 1.1, 0.9)
  for (expected in c(TRUE, FALSE)) {
    inf_full <- ll_ml(psi_ml, get_val = FALSE, expected = expected)$inf_mat
    st <- reconf:::score_stat(psi = psi_ml, test_idx = 1L, ll = ll_ml,
                              expected = expected, signed = TRUE)
    expect_equal(attr(st, "info"), eff_sub(inf_full, 1L),
                 tolerance = 1e-8, ignore_attr = TRUE)
  }
})

test_that("singular efficient information warns instead of clamping silently", {
  # Two variance parameters with the same structure matrix enter the
  # likelihood only through their sum, so the information for them is
  # exactly singular: the joint Cholesky fails, the subtraction fallback
  # produces a rank-deficient 2 x 2 block, and the relative eigenvalue floor
  # must flag it
  set.seed(9)
  n <- 40
  x <- rnorm(n)
  X <- cbind(1, x)
  g <- factor(rep(1:8, each = 5))
  Z <- Matrix::t(Matrix::fac2sparse(g))
  H <- methods::as(Matrix::Diagonal(8), "generalMatrix")
  Y <- as.vector(2 + x + rnorm(8)[g] + rnorm(n))
  ll <- reconf:::make_loglik(Y = Y, X = X, Z = Z, Hlist = list(H, H),
                             REML = TRUE)

  # The singular design also trips the condition-number diagnostic, so
  # collect all warnings and look for the eigenvalue-floor one
  w <- capture_warnings(
    reconf:::score_stat(psi = c(0.5, 0.5, 1), test_idx = c(1L, 2L), ll = ll,
                        signed = TRUE)
  )
  expect_true(any(grepl("singular to working precision", w)))
})

test_that("feasible starts are constructed from infeasible ones", {
  setup <- vc_model(fit_rs)
  psi <- setup$psi_hat
  feasible <- function(th) is.finite(reconf:::.ll_value(th, setup$ll))

  # A covariance fixed far beyond its variances: the free variances are
  # raised until the block is diagonally dominant
  th <- psi; th[2] <- 50 * sqrt(psi[1] * psi[3])
  expect_false(feasible(th))
  out <- reconf:::.feasible_start(th, fixed = 2L, setup = setup)
  expect_true(feasible(out))
  expect_identical(out[2], th[2])
  expect_gte(min(out[1], out[3]), abs(th[2]))

  # A negative fixed variance: only a larger error variance helps
  th <- psi; th[1] <- -50 * psi[1]
  out <- reconf:::.feasible_start(th, fixed = 1L, setup = setup)
  expect_true(feasible(out))
  expect_identical(out[1:3], c(th[1], 0, th[3]))
  expect_gt(out[4], psi[4])

  # The same with the error variance fixed has no feasible start this way
  expect_null(reconf:::.feasible_start(th, fixed = c(1L, 4L), setup = setup))
})
