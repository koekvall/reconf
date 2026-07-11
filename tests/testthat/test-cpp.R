library(lme4)
library(Matrix)

# ── Numerical derivative checks ───────────────────────────────────────────────

test_that("analytical ML score for psi agrees with numerical gradient", {
  skip_on_cran()
  skip_if_not_installed("numDeriv")
  fit <- lmer(Reaction ~ Days + (Days | Subject), data = sleepstudy, REML = FALSE)
  psi_hat <- reconf:::get_psi_hat_lmer(fit)
  b_hat   <- getME(fit, "beta")
  Y <- getME(fit, "y"); X <- getME(fit, "X"); Z <- getME(fit, "Z")
  Hlist <- reconf:::get_Hlist_lmer(fit)

  ll_fun <- function(psi) {
    reconf:::loglikelihood(psi = psi, b = b_hat, Y = Y, X = X, Z = Z,
                           Hlist = Hlist, REML = FALSE,
                           get_val = TRUE, get_score = FALSE, get_inf = FALSE)$value
  }
  num_score <- numDeriv::grad(ll_fun, psi_hat)
  ana_score <- reconf:::loglikelihood(psi = psi_hat, b = b_hat, Y = Y, X = X,
                                      Z = Z, Hlist = Hlist, REML = FALSE,
                                      get_val = FALSE, get_score = TRUE,
                                      get_inf = FALSE)$score
  expect_equal(num_score, ana_score, tolerance = 1e-4)
})

test_that("analytical ML observed information agrees with negative numerical Hessian", {
  skip_on_cran()
  skip_if_not_installed("numDeriv")
  fit <- lmer(Reaction ~ Days + (Days | Subject), data = sleepstudy, REML = FALSE)
  psi_hat <- reconf:::get_psi_hat_lmer(fit)
  b_hat   <- getME(fit, "beta")
  Y <- getME(fit, "y"); X <- getME(fit, "X"); Z <- getME(fit, "Z")
  Hlist <- reconf:::get_Hlist_lmer(fit)

  ll_fun <- function(psi) {
    reconf:::loglikelihood(psi = psi, b = b_hat, Y = Y, X = X, Z = Z,
                           Hlist = Hlist, REML = FALSE,
                           get_val = TRUE, get_score = FALSE, get_inf = FALSE)$value
  }
  num_hess <- -numDeriv::hessian(ll_fun, psi_hat)
  ana_hess <- reconf:::loglikelihood(psi = psi_hat, b = b_hat, Y = Y, X = X,
                                     Z = Z, Hlist = Hlist, REML = FALSE,
                                     get_val = FALSE, get_score = FALSE,
                                     get_inf = TRUE, expected = FALSE)$inf_mat
  expect_equal(num_hess, ana_hess, tolerance = 1e-4)
})

test_that("analytical REML score agrees with numerical gradient", {
  skip_on_cran()
  skip_if_not_installed("numDeriv")
  fit <- lmer(Reaction ~ Days + (Days | Subject), data = sleepstudy, REML = TRUE)
  psi_hat <- reconf:::get_psi_hat_lmer(fit)
  Y <- getME(fit, "y"); X <- getME(fit, "X"); Z <- getME(fit, "Z")
  Hlist <- reconf:::get_Hlist_lmer(fit)

  ll_fun <- function(psi) {
    reconf:::loglikelihood(psi = psi, Y = Y, X = X, Z = Z, Hlist = Hlist,
                           REML = TRUE, get_val = TRUE, get_score = FALSE,
                           get_inf = FALSE)$value
  }
  num_score <- numDeriv::grad(ll_fun, psi_hat)
  ana_score <- reconf:::loglikelihood(psi = psi_hat, Y = Y, X = X, Z = Z,
                                      Hlist = Hlist, REML = TRUE,
                                      get_val = FALSE, get_score = TRUE,
                                      get_inf = FALSE)$score
  expect_equal(num_score, ana_score, tolerance = 1e-4)
})

test_that("analytical REML expected information agrees with negative numerical Hessian", {
  skip_on_cran()
  skip_if_not_installed("numDeriv")
  fit <- lmer(Reaction ~ Days + (Days | Subject), data = sleepstudy, REML = TRUE)
  psi_hat <- reconf:::get_psi_hat_lmer(fit)
  Y <- getME(fit, "y"); X <- getME(fit, "X"); Z <- getME(fit, "Z")
  Hlist <- reconf:::get_Hlist_lmer(fit)

  ll_fun <- function(psi) {
    reconf:::loglikelihood(psi = psi, Y = Y, X = X, Z = Z, Hlist = Hlist,
                           REML = TRUE, get_val = TRUE, get_score = FALSE,
                           get_inf = FALSE)$value
  }
  num_hess <- -numDeriv::hessian(ll_fun, psi_hat)
  ana_hess <- reconf:::loglikelihood(psi = psi_hat, Y = Y, X = X, Z = Z,
                                     Hlist = Hlist, REML = TRUE,
                                     get_val = FALSE, get_score = FALSE,
                                     get_inf = TRUE, expected = TRUE)$inf_mat
  expect_equal(num_hess, ana_hess, tolerance = 1e-4)
})

# ── R and C++ implementations agree ──────────────────────────────────────────

test_that("Psi_from_H_cpp and Psi_from_Hlist agree", {
  fit <- lmer(Reaction ~ Days + (Days | Subject), data = sleepstudy, REML = FALSE)
  psi_hat <- reconf:::get_psi_hat_lmer(fit)
  Hlist   <- reconf:::get_Hlist_lmer(fit)
  H       <- do.call(cbind, Hlist)
  r       <- length(psi_hat)

  Psi_cpp <- reconf:::Psi_from_H_cpp(psi_mr = psi_hat[-r], H = H)
  Psi_R   <- reconf:::Psi_from_Hlist(psi_mr = psi_hat[-r], Hlist = Hlist)

  expect_equal(as.matrix(Psi_cpp), as.matrix(Psi_R), tolerance = 1e-12)
})

test_that("C++ and R log-likelihood values agree at MLE (ML)", {
  fit <- lmer(Reaction ~ Days + (Days | Subject), data = sleepstudy, REML = FALSE)
  psi_hat <- reconf:::get_psi_hat_lmer(fit)
  b_hat   <- getME(fit, "beta")
  Hlist   <- reconf:::get_Hlist_lmer(fit)

  ll_cpp <- reconf:::loglikelihood(
    psi = psi_hat, b = b_hat,
    Y = getME(fit, "y"), X = getME(fit, "X"), Z = getME(fit, "Z"),
    Hlist = Hlist, REML = FALSE,
    get_val = TRUE, get_score = FALSE, get_inf = FALSE
  )$value

  # R implementation
  H     <- do.call(cbind, Hlist)
  r     <- length(psi_hat)
  e     <- getME(fit, "y") - getME(fit, "X") %*% b_hat
  Z     <- getME(fit, "Z")
  Psi_r <- reconf:::Psi_from_H_cpp(psi_mr = psi_hat[-r], H = H) / psi_hat[r]

  ll_R <- reconf:::loglik_psi(
    Z = Z,
    ZtZXe = crossprod(Z, cbind(Z, getME(fit, "X"), e)),
    e = e, H = H, Psi_r = Psi_r, psi_r = psi_hat[r],
    get_val = TRUE, get_score = FALSE, get_inf = FALSE
  )$value

  expect_equal(ll_cpp, ll_R, tolerance = 1e-6)
})

test_that("C++ and R restricted log-likelihood values agree at MLE (REML)", {
  fit <- lmer(Reaction ~ Days + (Days | Subject), data = sleepstudy, REML = TRUE)
  psi_hat <- reconf:::get_psi_hat_lmer(fit)
  Hlist   <- reconf:::get_Hlist_lmer(fit)
  H       <- do.call(cbind, Hlist)
  r       <- length(psi_hat)
  Y       <- getME(fit, "y")
  X       <- getME(fit, "X")
  Z       <- getME(fit, "Z")
  Psi_r   <- reconf:::Psi_from_H_cpp(psi_mr = psi_hat[-r], H = H) / psi_hat[r]

  ll_cpp <- reconf:::loglikelihood(
    psi = psi_hat,
    Y = Y, X = X, Z = Z, Hlist = Hlist, REML = TRUE,
    get_val = TRUE, get_score = FALSE, get_inf = FALSE
  )$value

  ll_R <- reconf:::res_ll(
    XtX = crossprod(X), XtY = crossprod(X, Y),
    XtZ = crossprod(X, Z), ZtZ = crossprod(Z),
    ZtY = crossprod(Z, Y),
    Y = Y, X = X, Z = Z, H = H,
    Psi_r = Psi_r, psi_r = psi_hat[r],
    get_val = TRUE, get_score = FALSE, get_inf = FALSE
  )$value

  expect_equal(ll_cpp, ll_R, tolerance = 1e-6)
})

test_that("log-likelihood is -Inf when Sigma is not positive definite", {
  # sleepstudy is balanced, so all subjects' blocks cross the feasibility
  # boundary together: the determinant of I + Psi_r Z'Z stays positive at
  # infeasible parameters and a determinant sign check alone cannot detect
  # them. Guards the Cholesky-based gate in loglikelihood.
  fit <- lmer(Reaction ~ Days + (Days | Subject), data = sleepstudy, REML = FALSE)
  psi_hat <- reconf:::get_psi_hat_lmer(fit)
  b_hat   <- as.vector(getME(fit, "beta"))
  Y <- getME(fit, "y"); X <- getME(fit, "X"); Z <- getME(fit, "Z")
  Hlist <- reconf:::get_Hlist_lmer(fit)

  psi_bad <- psi_hat
  psi_bad[1] <- -1e6  # Sigma indefinite; det(I + Psi_r Z'Z) still positive

  ll_ml <- reconf:::loglikelihood(psi = psi_bad, b = b_hat, Y = Y, X = X, Z = Z,
                                  Hlist = Hlist, REML = FALSE,
                                  get_val = TRUE, get_score = FALSE,
                                  get_inf = FALSE)
  expect_identical(ll_ml$value, -Inf)

  ll_reml <- reconf:::loglikelihood(psi = psi_bad, Y = Y, X = X, Z = Z,
                                    Hlist = Hlist, REML = TRUE,
                                    get_val = TRUE, get_score = FALSE,
                                    get_inf = FALSE)
  expect_identical(ll_reml$value, -Inf)

  # Nonpositive error variance is infeasible
  psi_bad2 <- psi_hat
  psi_bad2[4] <- 0
  ll0 <- reconf:::loglikelihood(psi = psi_bad2, Y = Y, X = X, Z = Z,
                                Hlist = Hlist, REML = TRUE,
                                get_val = TRUE, get_score = FALSE,
                                get_inf = FALSE)
  expect_identical(ll0$value, -Inf)

  # A mildly negative variance parameter can still give Sigma PD; the gate
  # must not reject it (the parameter space is Sigma PD, not Psi PSD)
  psi_ok <- psi_hat
  psi_ok[1] <- -20
  ll_ok <- reconf:::loglikelihood(psi = psi_ok, Y = Y, X = X, Z = Z,
                                  Hlist = Hlist, REML = TRUE,
                                  get_val = TRUE, get_score = FALSE,
                                  get_inf = FALSE)
  expect_true(is.finite(ll_ok$value))
})

test_that("C++ and R score vectors agree at MLE (REML)", {
  fit <- lmer(Reaction ~ Days + (Days | Subject), data = sleepstudy, REML = TRUE)
  psi_hat <- reconf:::get_psi_hat_lmer(fit)
  Hlist   <- reconf:::get_Hlist_lmer(fit)
  H       <- do.call(cbind, Hlist)
  r       <- length(psi_hat)
  Y       <- getME(fit, "y")
  X       <- getME(fit, "X")
  Z       <- getME(fit, "Z")
  Psi_r   <- reconf:::Psi_from_H_cpp(psi_mr = psi_hat[-r], H = H) / psi_hat[r]

  score_cpp <- reconf:::loglikelihood(
    psi = psi_hat,
    Y = Y, X = X, Z = Z, Hlist = Hlist, REML = TRUE,
    get_val = FALSE, get_score = TRUE, get_inf = FALSE
  )$score

  score_R <- reconf:::res_ll(
    XtX = crossprod(X), XtY = crossprod(X, Y),
    XtZ = crossprod(X, Z), ZtZ = crossprod(Z),
    ZtY = crossprod(Z, Y),
    Y = Y, X = X, Z = Z, H = H,
    Psi_r = Psi_r, psi_r = psi_hat[r],
    get_val = FALSE, get_score = TRUE, get_inf = FALSE
  )$score

  expect_equal(score_cpp, score_R, tolerance = 1e-6)
})
