library(lme4)
library(Matrix)

# The spectral kernels (loglik_spectral, loglik_res_spectral) apply when
# r = 2, where Sigma = psi_1 K + psi_2 I_n shares the eigenvectors of
# K = Z H_1 Z' for every psi. They are validated against the q-side path
# through loglikelihood(), which is itself checked against lme4 and
# numerical derivatives in test-cpp.R; the dense-Z design below also has
# q > n, the case where method = "auto" selects the spectral path.

# Dense kernel-type design with q > n and a single random-effect variance;
# Z Z' is positive definite, so psi_1 alone can make Sigma "look" positive
# definite and the explicit psi_2 > 0 check matters.
make_dense_r2 <- function(n = 25, q = 30, seed = 1) {
  set.seed(seed)
  Z <- methods::as(Matrix::Matrix(matrix(rnorm(n * q), n, q) / sqrt(q),
                                  sparse = TRUE), "generalMatrix")
  Hlist <- list(methods::as(Matrix::Diagonal(q), "generalMatrix"))
  X <- cbind(1, rnorm(n))
  psi_true <- c(0.8, 0.5)
  u <- rnorm(q, sd = sqrt(psi_true[1]))
  Y <- as.vector(X %*% c(2, -1) + Z %*% u + rnorm(n, sd = sqrt(psi_true[2])))
  list(Y = Y, X = X, Z = Z, Hlist = Hlist, psi = psi_true, n = n, q = q)
}

# Spectral evaluation of the same likelihood; ll_q (helper-paths.R) is the oracle
ll_s <- function(psi, Y, X, Z, Hlist, REML, b = NULL, expected = TRUE) {
  pc <- reconf:::get_precomp(Y = Y, X = X, Z = Z, REML = REML, Hlist = Hlist,
                             method = "spectral")
  if (REML) {
    reconf:::loglik_res_spectral(d = pc$d, Yt = pc$Yt, Xt = pc$Xt, psi = psi)
  } else {
    reconf:::loglik_spectral(d = pc$d, Yt = pc$Yt, Xt = pc$Xt, psi = psi,
                             b = b, expected = expected)
  }
}

# ── Path agreement, grouped design (sleepstudy random intercept, q < n) ──────

test_that("spectral agrees with q-side on sleepstudy (ML, REML)", {
  fit <- lmer(Reaction ~ Days + (1 | Subject), data = sleepstudy, REML = FALSE)
  psi_hat <- reconf:::get_psi_hat_lmer(fit)
  b_hat   <- as.vector(getME(fit, "beta"))
  Y <- getME(fit, "y"); X <- getME(fit, "X"); Z <- getME(fit, "Z")
  Hlist <- reconf:::get_Hlist_lmer(fit)

  # At the MLE and at a nearby non-stationary point, where the score is nonzero
  for (psi in list(psi_hat, psi_hat * c(1.3, 0.8))) {
    for (expected in c(TRUE, FALSE)) {
      out_q <- ll_q(psi, Y, X, Z, Hlist, REML = FALSE, b = b_hat,
                    expected = expected)
      out_s <- ll_s(psi, Y, X, Z, Hlist, REML = FALSE, b = b_hat,
                    expected = expected)
      expect_equal(out_s$value, out_q$value, tolerance = 1e-8)
      expect_equal(out_s$score, out_q$score, tolerance = 1e-6)
      expect_equal(out_s$inf_mat, out_q$inf_mat, tolerance = 1e-6)
    }
    out_q <- ll_q(psi, Y, X, Z, Hlist, REML = TRUE)
    out_s <- ll_s(psi, Y, X, Z, Hlist, REML = TRUE)
    expect_equal(out_s$value, out_q$value, tolerance = 1e-8)
    expect_equal(out_s$score, out_q$score, tolerance = 1e-6)
    expect_equal(out_s$inf_mat, out_q$inf_mat, tolerance = 1e-6)
  }
})

# ── Path agreement, dense design with q > n ──────────────────────────────────

test_that("spectral agrees with q-side in a dense design with q > n", {
  d <- make_dense_r2()
  expect_true(d$q > d$n)
  b <- c(2, -1)

  for (psi in list(d$psi, c(0.4, 1.2))) {
    for (expected in c(TRUE, FALSE)) {
      out_q <- ll_q(psi, d$Y, d$X, d$Z, d$Hlist, REML = FALSE, b = b,
                    expected = expected)
      out_s <- ll_s(psi, d$Y, d$X, d$Z, d$Hlist, REML = FALSE, b = b,
                    expected = expected)
      expect_equal(out_s$value, out_q$value, tolerance = 1e-8)
      expect_equal(out_s$score, out_q$score, tolerance = 1e-6)
      expect_equal(out_s$inf_mat, out_q$inf_mat, tolerance = 1e-6)
    }
    out_q <- ll_q(psi, d$Y, d$X, d$Z, d$Hlist, REML = TRUE)
    out_s <- ll_s(psi, d$Y, d$X, d$Z, d$Hlist, REML = TRUE)
    expect_equal(out_s$value, out_q$value, tolerance = 1e-8)
    expect_equal(out_s$score, out_q$score, tolerance = 1e-6)
    expect_equal(out_s$inf_mat, out_q$inf_mat, tolerance = 1e-6)
  }
})

test_that("spectral REML beta and I_b_chol match direct GLS computations", {
  d <- make_dense_r2()

  out_s <- ll_s(d$psi, d$Y, d$X, d$Z, d$Hlist, REML = TRUE)

  # Direct dense computation of the GLS quantities
  K <- as.matrix(Matrix::tcrossprod(d$Z))
  Sigma <- d$psi[1] * K + d$psi[2] * diag(d$n)
  XtSiX <- crossprod(d$X, solve(Sigma, d$X))
  beta_gls <- solve(XtSiX, crossprod(d$X, solve(Sigma, d$Y)))
  expect_equal(out_s$beta, as.vector(beta_gls), tolerance = 1e-6)
  expect_equal(crossprod(out_s$I_b_chol), XtSiX, tolerance = 1e-6,
               ignore_attr = TRUE)
})

# ── Numerical derivative spot-checks on the spectral value itself ────────────

test_that("spectral derivatives agree with numerical ones (q > n)", {
  skip_on_cran()
  skip_if_not_installed("numDeriv")
  d <- make_dense_r2()
  b <- c(2, -1)
  pc <- reconf:::get_precomp(Y = d$Y, X = d$X, Z = d$Z, Hlist = d$Hlist,
                             method = "spectral")

  # ML: score and observed information against the spectral value
  val_ml <- function(psi) {
    reconf:::loglik_spectral(d = pc$d, Yt = pc$Yt, Xt = pc$Xt, psi = psi,
                             b = b, get_score = FALSE, get_inf = FALSE)$value
  }
  out <- reconf:::loglik_spectral(d = pc$d, Yt = pc$Yt, Xt = pc$Xt,
                                  psi = d$psi, b = b, expected = FALSE)
  p <- ncol(d$X)
  expect_equal(numDeriv::grad(val_ml, d$psi),
               out$score[(p + 1):(p + 2)], tolerance = 1e-4)
  expect_equal(-numDeriv::hessian(val_ml, d$psi),
               out$inf_mat[(p + 1):(p + 2), (p + 1):(p + 2)],
               tolerance = 1e-4)

  # REML: score against the spectral value
  val_reml <- function(psi) {
    reconf:::loglik_res_spectral(d = pc$d, Yt = pc$Yt, Xt = pc$Xt, psi = psi,
                                 get_score = FALSE, get_inf = FALSE)$value
  }
  out_reml <- reconf:::loglik_res_spectral(d = pc$d, Yt = pc$Yt, Xt = pc$Xt,
                                           psi = d$psi)
  expect_equal(numDeriv::grad(val_reml, d$psi), out_reml$score,
               tolerance = 1e-4)
})

# ── Feasibility gate agreement ───────────────────────────────────────────────

test_that("spectral gate matches q-side feasibility decisions", {
  d <- make_dense_r2()

  # Infeasible: Sigma indefinite
  psi_bad <- c(-100, 0.5)
  out_q <- ll_q(psi_bad, d$Y, d$X, d$Z, d$Hlist, REML = TRUE)
  out_s <- ll_s(psi_bad, d$Y, d$X, d$Z, d$Hlist, REML = TRUE)
  expect_identical(out_q$value, -Inf)
  expect_identical(out_s$value, -Inf)
  expect_true(all(out_s$score == 0) && all(out_s$inf_mat == 0))
  expect_true(all(is.na(out_s$beta)) && all(is.na(out_s$I_b_chol)))

  # Mildly negative variance with Sigma still positive definite: feasible on
  # both paths and equal (the parameter space is Sigma PD, not Psi PSD)
  psi_ok <- c(-0.05, 0.7)
  out_q <- ll_q(psi_ok, d$Y, d$X, d$Z, d$Hlist, REML = TRUE)
  out_s <- ll_s(psi_ok, d$Y, d$X, d$Z, d$Hlist, REML = TRUE)
  expect_true(is.finite(out_q$value))
  expect_equal(out_s$value, out_q$value, tolerance = 1e-8)

  # Nonpositive error variance is infeasible even though K = Z Z' is
  # positive definite here (q > n); caught by the explicit psi_2 check
  for (psi_r_bad in c(0, -0.1)) {
    out_s <- ll_s(c(1, psi_r_bad), d$Y, d$X, d$Z, d$Hlist, REML = TRUE)
    expect_identical(out_s$value, -Inf)
    out_s_ml <- ll_s(c(1, psi_r_bad), d$Y, d$X, d$Z, d$Hlist, REML = FALSE,
                     b = c(0, 0))
    expect_identical(out_s_ml$value, -Inf)
  }
})

# ── Precompute layer ─────────────────────────────────────────────────────────

test_that("get_precomp spectral stores eigenvalues and a rotation", {
  d <- make_dense_r2()
  pc <- reconf:::get_precomp(Y = d$Y, X = d$X, Z = d$Z, Hlist = d$Hlist,
                             method = "spectral")
  # Eigenvalues of K = Z H_1 Z'
  K <- as.matrix(Matrix::tcrossprod(d$Z))
  expect_equal(pc$d, eigen(K, symmetric = TRUE, only.values = TRUE)$values,
               tolerance = 1e-10)
  # Yt and Xt are the data in a common orthonormal basis: all inner products
  # are preserved
  expect_equal(sum(pc$Yt^2), sum(d$Y^2), tolerance = 1e-10)
  expect_equal(crossprod(pc$Xt), crossprod(d$X), tolerance = 1e-10,
               ignore_attr = TRUE)
  expect_equal(as.vector(crossprod(pc$Xt, pc$Yt)),
               as.vector(crossprod(d$X, d$Y)), tolerance = 1e-10)

  # r > 2 is not representable on this path
  expect_error(reconf:::get_precomp(Y = d$Y, X = d$X, Z = d$Z,
                                    Hlist = c(d$Hlist, d$Hlist),
                                    method = "spectral"))
})

# ── Dispatch through loglikelihood() ─────────────────────────────────────────

test_that("loglikelihood dispatches on method and precomp tag", {
  d <- make_dense_r2()
  b <- c(2, -1)

  # auto picks the spectral path iff the dense regime applies (q >= n, dense
  # Z) and r = 2; sparse designs keep the q-side even when r = 2
  expect_identical(
    reconf:::get_precomp(d$Y, d$X, d$Z, Hlist = d$Hlist,
                         method = "auto")$method,
    "spectral")
  fit <- lmer(Reaction ~ Days + (1 | Subject), data = sleepstudy)
  expect_identical(
    reconf:::get_precomp_lmer(fit, Hlist = reconf:::get_Hlist_lmer(fit))$method,
    "q_side")

  # Explicit method: paths agree through the wrapper, REML and ML
  for (psi in list(d$psi, c(0.4, 1.2))) {
    out_q <- reconf:::loglikelihood(psi = psi, Y = d$Y, X = d$X, Z = d$Z,
                                    Hlist = d$Hlist, REML = TRUE,
                                    method = "q_side")
    out_s <- reconf:::loglikelihood(psi = psi, Y = d$Y, X = d$X, Z = d$Z,
                                    Hlist = d$Hlist, REML = TRUE,
                                    method = "spectral")
    expect_equal(out_s$value, out_q$value, tolerance = 1e-8)
    expect_equal(out_s$score, out_q$score, tolerance = 1e-6)
    expect_equal(out_s$inf_mat, out_q$inf_mat, tolerance = 1e-6)

    # ML without get_beta exercises the wrapper's stripping of the beta block
    out_q <- reconf:::loglikelihood(psi = psi, b = b, Y = d$Y, X = d$X,
                                    Z = d$Z, Hlist = d$Hlist, REML = FALSE,
                                    get_beta = FALSE, method = "q_side")
    out_s <- reconf:::loglikelihood(psi = psi, b = b, Y = d$Y, X = d$X,
                                    Z = d$Z, Hlist = d$Hlist, REML = FALSE,
                                    get_beta = FALSE, method = "spectral")
    expect_equal(length(out_s$score), length(psi))
    expect_equal(out_s$score, out_q$score, tolerance = 1e-6)
    expect_equal(out_s$inf_mat, out_q$inf_mat, tolerance = 1e-6)
  }

  # A supplied precomp determines the path regardless of method
  pc_s <- reconf:::get_precomp(d$Y, d$X, d$Z, Hlist = d$Hlist,
                               method = "spectral")
  out_pc <- reconf:::loglikelihood(psi = d$psi, Y = d$Y, X = d$X, Z = d$Z,
                                   Hlist = d$Hlist, REML = TRUE,
                                   precomp = pc_s, method = "q_side")
  out_s <- reconf:::loglikelihood(psi = d$psi, Y = d$Y, X = d$X, Z = d$Z,
                                  Hlist = d$Hlist, REML = TRUE,
                                  method = "spectral")
  expect_identical(out_pc, out_s)

  # Infeasible parameters keep the q-side return shapes
  out_bad <- reconf:::loglikelihood(psi = c(-100, 0.5), Y = d$Y, X = d$X,
                                    Z = d$Z, Hlist = d$Hlist, REML = TRUE,
                                    method = "spectral")
  expect_identical(out_bad$value, -Inf)
  expect_identical(length(out_bad$score), 2L)
})

# ── End to end through the lme4 front door ───────────────────────────────────

test_that("ci_all_lmer with method = 'spectral' matches the q-side", {
  skip_on_cran()
  fit <- lmer(Reaction ~ Days + (1 | Subject), data = sleepstudy, REML = TRUE)

  ci_q <- ci_all_lmer(fit, method = "q_side")
  ci_s <- ci_all_lmer(fit, method = "spectral")
  expect_true(all(is.finite(ci_q[, "lower"])) &&
                all(is.finite(ci_q[, "upper"])))
  expect_equal(unclass(ci_s), unclass(ci_q), tolerance = 1e-4)
})
