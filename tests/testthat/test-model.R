library(lme4)

fit_rs <- lmer(Reaction ~ Days + (Days | Subject), data = sleepstudy, REML = TRUE)

# ── vc_model ─────────────────────────────────────────────────────────────────

test_that("vc_model holds the parameters and is accepted in place of the fit", {
  skip_on_cran()
  m <- vc_model(fit_rs)
  expect_s3_class(m, "vc_model")
  expect_identical(m$names, c("var_(Intercept)|Subject",
                              "cov_Days.(Intercept)|Subject",
                              "var_Days|Subject", "var_Residual"))
  expect_output(print(m), "var_Days\\|Subject")
  expect_identical(vc_ci(m), vc_ci(fit_rs))
  expect_identical(vc_test(m, parm = 3L)$statistic,
                   vc_test(fit_rs, parm = 3L)$statistic)
  # A vc_model fixes the computational method
  expect_error(vc_ci(m, method = "n_side"), "method is set when the model")
  expect_error(vc_model(lm(Reaction ~ Days, sleepstudy)), "no method")
})

# ── vc_loglik ────────────────────────────────────────────────────────────────

test_that("vc_loglik returns named value, score, and information", {
  ll <- vc_loglik(fit_rs)
  psi_hat <- vc_model(fit_rs)$psi_hat
  out <- ll(psi_hat)
  expect_equal(out$value, as.numeric(logLik(fit_rs)), tolerance = 1e-8)
  expect_named(out$score, vc_model(fit_rs)$names)
  expect_identical(dimnames(out$information),
                   list(names(out$score), names(out$score)))
  expect_lt(max(abs(out$score)), 1e-2)
  # Components not requested are NULL
  out_v <- ll(psi_hat, get_score = FALSE, get_inf = FALSE)
  expect_null(out_v$score)
  expect_null(out_v$information)
  expect_identical(out_v$value, out$value)
  # Outside the parameter set
  expect_identical(ll(c(-1e6, 0, 35, 650))$value, -Inf)
  expect_error(ll(psi_hat[-1]), "length r = 4")
  expect_error(ll(c(NA, 10, 35, 650)), "finite")
})
