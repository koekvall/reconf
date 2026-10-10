library(lme4)
library(Matrix)

# All tests use the sleepstudy random-slope model. Fits and extracted
# components are shared across tests; Y, X, Z, and Hlist are identical for
# the ML and REML fits.
fit_ml   <- lmer(Reaction ~ Days + (Days | Subject), data = sleepstudy,
                 REML = FALSE)
fit_reml <- lmer(Reaction ~ Days + (Days | Subject), data = sleepstudy,
                 REML = TRUE)
Y <- getME(fit_ml, "y"); X <- getME(fit_ml, "X"); Z <- getME(fit_ml, "Z")
Hlist    <- reconf:::get_Hlist_lmer(fit_ml)
H        <- do.call(cbind, Hlist)
psi_ml   <- reconf:::get_psi_hat_lmer(fit_ml)
psi_reml <- reconf:::get_psi_hat_lmer(fit_reml)
r        <- length(psi_ml)

# The ML likelihood is profiled over beta (generalized least squares)
ll_ml   <- reconf:::make_loglik(Y = Y, X = X, Z = Z, Hlist = Hlist, REML = FALSE)
ll_reml <- reconf:::make_loglik(Y = Y, X = X, Z = Z, Hlist = Hlist, REML = TRUE)
val_ml   <- function(psi) ll_ml(psi, get_score = FALSE, get_inf = FALSE)$value
val_reml <- function(psi) ll_reml(psi, get_score = FALSE, get_inf = FALSE)$value

# Generalized least squares estimate of beta at psi, by dense algebra
gls_beta <- function(psi) {
  Psi <- as.matrix(reconf:::Psi_from_H_cpp(psi_mr = psi[-r], H = H))
  Sigma <- as.matrix(Z %*% Psi %*% t(Z)) + psi[r] * diag(length(Y))
  as.vector(solve(crossprod(X, solve(Sigma, X)), crossprod(X, solve(Sigma, Y))))
}

# ── Numerical derivative checks ───────────────────────────────────────────────

test_that("analytical ML score for psi agrees with numerical gradient", {
  skip_on_cran()
  skip_if_not_installed("numDeriv")
  # The profile likelihood's gradient is the score at the profiled beta;
  # checked at the estimate and away from it
  for (psi in list(psi_ml, psi_ml * c(1.4, 0.7, 1.2, 0.8))) {
    num_score <- numDeriv::grad(val_ml, psi)
    ana_score <- ll_ml(psi, get_val = FALSE, get_inf = FALSE)$score
    expect_equal(num_score, ana_score, tolerance = 1e-4)
  }
})

test_that("analytical ML observed information agrees with negative numerical Hessian", {
  skip_on_cran()
  skip_if_not_installed("numDeriv")
  # The negative Hessian of the profile likelihood is the Schur complement
  # of the beta block in the observed information for (beta, psi)
  for (psi in list(psi_ml, psi_ml * c(1.4, 0.7, 1.2, 0.8))) {
    num_hess <- -numDeriv::hessian(val_ml, psi)
    ana_hess <- ll_ml(psi, get_val = FALSE, get_score = FALSE,
                      expected = FALSE)$inf_mat
    expect_equal(num_hess, ana_hess, tolerance = 1e-4)
  }
})

test_that("analytical REML score agrees with numerical gradient", {
  skip_on_cran()
  skip_if_not_installed("numDeriv")
  num_score <- numDeriv::grad(val_reml, psi_reml)
  ana_score <- ll_reml(psi_reml, get_val = FALSE, get_inf = FALSE)$score
  expect_equal(num_score, ana_score, tolerance = 1e-4)
})

test_that("analytical REML observed information agrees with negative numerical Hessian", {
  skip_on_cran()
  skip_if_not_installed("numDeriv")
  # Checked at the REML estimate and away from it, where observed and
  # expected differ by about 30 percent. At the estimate they nearly
  # coincide, but only because sleepstudy is balanced; do not compare the
  # numerical Hessian to the expected information on that basis.
  for (psi in list(psi_reml, psi_reml * c(1.4, 0.7, 1.2, 0.8))) {
    num_hess <- -numDeriv::hessian(val_reml, psi)
    ana_hess <- ll_reml(psi, get_val = FALSE, get_score = FALSE,
                        expected = FALSE)$inf_mat
    expect_equal(num_hess, ana_hess, tolerance = 1e-4)
  }
})

# ── R and C++ implementations agree ──────────────────────────────────────────

test_that("Psi_from_H_cpp and Psi_from_Hlist agree", {
  Psi_cpp <- reconf:::Psi_from_H_cpp(psi_mr = psi_ml[-r], H = H)
  Psi_R   <- Psi_from_Hlist(psi_mr = psi_ml[-r], Hlist = Hlist)

  expect_equal(as.matrix(Psi_cpp), as.matrix(Psi_R), tolerance = 1e-12)
})

test_that("C++ and R log-likelihood values agree at MLE (ML)", {
  ll_cpp <- val_ml(psi_ml)

  # R implementation at the generalized least squares estimate of beta
  e     <- Y - X %*% gls_beta(psi_ml)
  Psi_r <- reconf:::Psi_from_H_cpp(psi_mr = psi_ml[-r], H = H) / psi_ml[r]

  ll_R <- loglik_psi(
    Z = Z,
    ZtZXe = crossprod(Z, cbind(Z, X, e)),
    e = e, H = H, Psi_r = Psi_r, psi_r = psi_ml[r],
    get_val = TRUE, get_score = FALSE, get_inf = FALSE
  )$value

  expect_equal(ll_cpp, ll_R, tolerance = 1e-6)
})

test_that("C++ and R restricted log-likelihood values agree at MLE (REML)", {
  Psi_r <- reconf:::Psi_from_H_cpp(psi_mr = psi_reml[-r], H = H) / psi_reml[r]

  ll_cpp <- val_reml(psi_reml)

  ll_R <- res_ll(
    XtX = crossprod(X), XtY = crossprod(X, Y),
    XtZ = crossprod(X, Z), ZtZ = crossprod(Z),
    ZtY = crossprod(Z, Y),
    Y = Y, X = X, Z = Z, H = H,
    Psi_r = Psi_r, psi_r = psi_reml[r],
    get_val = TRUE, get_score = FALSE, get_inf = FALSE
  )$value

  expect_equal(ll_cpp, ll_R, tolerance = 1e-6)
})

test_that("log-likelihood is -Inf when Sigma is not positive definite", {
  # sleepstudy is balanced, so all subjects' blocks cross the feasibility
  # boundary together: the determinant of I + Psi_r Z'Z stays positive at
  # infeasible parameters and a determinant sign check alone cannot detect
  # them. Guards the Cholesky-based gate in make_loglik.
  psi_bad <- psi_ml
  psi_bad[1] <- -1e6  # Sigma indefinite; det(I + Psi_r Z'Z) still positive
  expect_identical(val_ml(psi_bad), -Inf)
  expect_identical(val_reml(psi_bad), -Inf)

  # Nonpositive error variance is infeasible
  psi_bad2 <- psi_ml
  psi_bad2[r] <- 0
  expect_identical(val_reml(psi_bad2), -Inf)

  # A mildly negative variance parameter can still give Sigma PD; the gate
  # must not reject it (the parameter space is Sigma PD, not Psi PSD)
  psi_ok <- psi_ml
  psi_ok[1] <- -20
  expect_true(is.finite(val_reml(psi_ok)))
})

test_that("C++ and R score vectors agree at MLE (REML)", {
  Psi_r <- reconf:::Psi_from_H_cpp(psi_mr = psi_reml[-r], H = H) / psi_reml[r]

  score_cpp <- ll_reml(psi_reml, get_val = FALSE, get_inf = FALSE)$score

  score_R <- res_ll(
    XtX = crossprod(X), XtY = crossprod(X, Y),
    XtZ = crossprod(X, Z), ZtZ = crossprod(Z),
    ZtY = crossprod(Z, Y),
    Y = Y, X = X, Z = Z, H = H,
    Psi_r = Psi_r, psi_r = psi_reml[r],
    get_val = FALSE, get_score = TRUE, get_inf = FALSE
  )$score

  expect_equal(score_cpp, score_R, tolerance = 1e-6)
})
