library(lme4)

# ── Our optimizer recovers lme4's MLE ────────────────────────────────────────

test_that("maximize_loglik (REML) recovers lme4 estimates", {
  fit <- lmer(Reaction ~ Days + (Days | Subject), data = sleepstudy, REML = TRUE)

  psi_hat <- reconf:::get_psi_hat_lmer(fit)
  r <- length(psi_hat)

  ll <- reconf:::make_loglik(Y = getME(fit, "y"), X = getME(fit, "X"),
                             Z = getME(fit, "Z"),
                             Hlist = reconf:::get_Hlist_lmer(fit), REML = TRUE)
  fit_our <- reconf:::maximize_loglik(start_val = psi_hat, opt_idx = seq_len(r),
                                      ll = ll)

  expect_equal(unname(fit_our$arg), unname(psi_hat), tolerance = 1e-4)
})

test_that("maximize_loglik (ML) recovers lme4 estimates", {
  fit <- lmer(Reaction ~ Days + (Days | Subject), data = sleepstudy, REML = FALSE)

  psi_hat <- reconf:::get_psi_hat_lmer(fit)
  r       <- length(psi_hat)

  # The fixed effects are profiled out of the ML likelihood, so only psi is
  # optimized; start away from the estimate
  ll <- reconf:::make_loglik(Y = getME(fit, "y"), X = getME(fit, "X"),
                             Z = getME(fit, "Z"),
                             Hlist = reconf:::get_Hlist_lmer(fit), REML = FALSE)
  fit_our <- reconf:::maximize_loglik(start_val = psi_hat * c(1.5, 0.5, 1.5, 0.8),
                                      opt_idx = seq_len(r), ll = ll)

  expect_equal(unname(fit_our$arg), unname(psi_hat), tolerance = 1e-4)
  expect_equal(fit_our$value, as.numeric(logLik(fit)), tolerance = 1e-8)
})

# ── get_psi_hat_lmer returns the VarCorr vcov vector ─────────────────────────

test_that("get_psi_hat_lmer matches as.data.frame(VarCorr(...))", {
  fit <- lmer(Reaction ~ Days + (Days | Subject), data = sleepstudy)
  psi_hat <- reconf:::get_psi_hat_lmer(fit)
  vc <- as.data.frame(VarCorr(fit), order = "lower.tri")$vcov
  expect_equal(psi_hat, vc, tolerance = 1e-12)
})
