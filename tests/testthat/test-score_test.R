library(lme4)

fit_ri <- lmer(Reaction ~ Days + (1 | Subject), data = sleepstudy, REML = TRUE)
fit_rs <- lmer(Reaction ~ Days + (Days | Subject), data = sleepstudy, REML = TRUE)

# ── score_test_lmer return structure ────────────────────────────────────────

test_that("score_test_lmer returns named numeric vector of length 3", {
  res <- score_test_lmer(fit_ri)
  expect_named(res, c("stat", "p_val", "df"))
  expect_length(res, 3L)
  expect_true(is.numeric(res))
})

test_that("p-value is in [0, 1]", {
  res <- score_test_lmer(fit_ri)
  expect_gte(res[["p_val"]], 0)
  expect_lte(res[["p_val"]], 1)
})

test_that("df equals length of test_idx", {
  res <- score_test_lmer(fit_ri, test_idx = 1L)
  expect_equal(res[["df"]], 1)
})

# ── score at MLE should be near zero ────────────────────────────────────────

test_that("score stat near zero when theta_null equals MLE", {
  psi_hat <- reconf:::get_psi_hat_lmer(fit_ri)
  # Test with null = MLE: stat should be ~0, p-value ~1
  res <- score_test_lmer(fit_ri,
                         theta_null = psi_hat,
                         test_idx = 1L,
                         profile = FALSE)
  expect_lt(abs(res[["stat"]]), 0.01)
})

# ── rejects zero variance when between-subject variability is strong ─────────

test_that("random intercept is significant in sleepstudy", {
  res <- score_test_lmer(fit_ri)
  expect_lt(res[["p_val"]], 0.001)
})

# ── multiple random effects ──────────────────────────────────────────────────

test_that("works with random slope model, testing all RE params", {
  res <- score_test_lmer(fit_rs)
  expect_named(res, c("stat", "p_val", "df"))
  expect_equal(res[["df"]], 3)  # intercept var, covariance, slope var
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
  b_hat <- as.vector(fixef(fit_rs))
  Y <- getME(fit_rs, "y"); X <- getME(fit_rs, "X"); Z <- getME(fit_rs, "Z")
  Hlist <- reconf:::get_Hlist_lmer(fit_rs)

  # REML at the estimate and at a non-stationary point; single and joint
  # tests, and with a known parameter excluded from the nuisances
  for (psi in list(psi_hat, psi_hat * c(1.3, 0.8, 1.1, 0.9))) {
    inf_full <- reconf:::loglikelihood(psi = psi, Y = Y, X = X, Z = Z,
                                       Hlist = Hlist, REML = TRUE,
                                       get_val = FALSE)$inf_mat
    for (cfg in list(list(t = 1L, k = NULL), list(t = c(1L, 3L), k = NULL),
                     list(t = 1L, k = 2L))) {
      st <- reconf:::score_stat(theta = psi, test_idx = cfg$t,
                                known_idx = cfg$k, Y = Y, X = X, Z = Z,
                                Hlist = Hlist, REML = TRUE, signed = TRUE)
      expect_equal(attr(st, "info"), eff_sub(inf_full, cfg$t, cfg$k),
                   tolerance = 1e-10, ignore_attr = TRUE)
    }
  }

  # ML with fixed effects in theta
  fit_ml <- lmer(Reaction ~ Days + (Days | Subject), data = sleepstudy,
                 REML = FALSE)
  theta <- c(as.vector(fixef(fit_ml)), reconf:::get_psi_hat_lmer(fit_ml))
  inf_full <- reconf:::loglikelihood(psi = theta[-(1:2)], b = theta[1:2],
                                     Y = getME(fit_ml, "y"),
                                     X = getME(fit_ml, "X"),
                                     Z = getME(fit_ml, "Z"),
                                     Hlist = Hlist, REML = FALSE,
                                     get_val = FALSE, get_beta = TRUE)$inf_mat
  st <- reconf:::score_stat(theta = theta, test_idx = 3L,
                            Y = getME(fit_ml, "y"), X = getME(fit_ml, "X"),
                            Z = getME(fit_ml, "Z"), Hlist = Hlist,
                            REML = FALSE, signed = TRUE)
  expect_equal(attr(st, "info"), eff_sub(inf_full, 3L),
               tolerance = 1e-10, ignore_attr = TRUE)
})

test_that("singular efficient information warns instead of clamping silently", {
  # Duplicated regressor: the ML information for the two collinear
  # coefficients is exactly singular, so the joint Cholesky fails, the
  # subtraction fallback produces a rank-deficient 2 x 2 block, and the
  # relative eigenvalue floor must flag it
  set.seed(9)
  n <- 40
  x <- rnorm(n)
  X <- cbind(1, x, x)
  g <- factor(rep(1:8, each = 5))
  Z <- Matrix::t(Matrix::fac2sparse(g))
  Hlist <- list(methods::as(Matrix::Diagonal(8), "generalMatrix"))
  Y <- as.vector(2 + x + rnorm(8)[g] + rnorm(n))
  theta <- c(2, 0.5, 0.5, 1, 1)  # (beta, psi)

  # The singular design also trips the condition-number diagnostic, so
  # collect all warnings and look for the eigenvalue-floor one
  w <- capture_warnings(
    reconf:::score_stat(theta = theta, test_idx = c(2L, 3L), Y = Y, X = X,
                        Z = Z, Hlist = Hlist, REML = FALSE, signed = TRUE)
  )
  expect_true(any(grepl("nearly singular", w)))
})
