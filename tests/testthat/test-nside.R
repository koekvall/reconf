library(lme4)
library(Matrix)

# The n-side kernels (loglik_n, loglik_res_n) are validated against the
# q-side path through loglikelihood(), which is itself checked against lme4
# and numerical derivatives in test-cpp.R. The q-side algebra does not
# require q < n, so the q-side path serves as the oracle also in the
# crossed q > n design below (kept small so cost is irrelevant).

# Crossed random intercepts with q = l1 + l2 > n; Z'Z is singular, and for
# suitable psi, Z Psi Z' alone is positive definite.
make_crossed <- function(n = 30, l1 = 12, l2 = 20, seed = 1) {
  set.seed(seed)
  # Unused levels give zero columns in Z, which the likelihood handles
  f1 <- factor(sample.int(l1, n, replace = TRUE), levels = seq_len(l1))
  f2 <- factor(sample.int(l2, n, replace = TRUE), levels = seq_len(l2))
  Z <- cbind(Matrix::t(Matrix::fac2sparse(f1, drop.unused.levels = FALSE)),
             Matrix::t(Matrix::fac2sparse(f2, drop.unused.levels = FALSE)))
  q <- l1 + l2
  Hlist <- list(
    methods::as(Matrix::Diagonal(q, c(rep(1, l1), rep(0, l2))), "generalMatrix"),
    methods::as(Matrix::Diagonal(q, c(rep(0, l1), rep(1, l2))), "generalMatrix")
  )
  X <- cbind(1, rnorm(n))
  psi_true <- c(1, 0.5, 0.7)
  u <- c(rnorm(l1, sd = sqrt(psi_true[1])), rnorm(l2, sd = sqrt(psi_true[2])))
  Y <- as.vector(X %*% c(2, -1) + Z %*% u + rnorm(n, sd = sqrt(psi_true[3])))
  list(Y = Y, X = X, Z = Z, Hlist = Hlist, psi = psi_true, n = n, q = q)
}

# q-side and n-side evaluations of the same likelihood, for comparison
ll_q <- function(psi, Y, X, Z, Hlist, REML, b = NULL, expected = TRUE) {
  reconf:::loglikelihood(psi = psi, b = b, Y = Y, X = X, Z = Z, Hlist = Hlist,
                         REML = REML, get_val = TRUE, get_score = TRUE,
                         get_inf = TRUE, get_beta = TRUE, expected = expected,
                         check = FALSE)
}

ll_n <- function(psi, Y, X, Z, Hlist, REML, b = NULL, expected = TRUE) {
  K <- reconf:::get_precomp(Y = Y, X = X, Z = Z, REML = REML, Hlist = Hlist,
                            method = "n_side")$K
  if (REML) {
    reconf:::loglik_res_n(K = K, psi = psi, Y = Y, X = X)
  } else {
    e <- as.vector(Y - X %*% b)
    reconf:::loglik_n(K = K, psi = psi, e = e, X = X, expected = expected)
  }
}

# ── Path agreement, grouped design (sleepstudy, q < n) ───────────────────────

test_that("n-side agrees with q-side on sleepstudy (ML, expected and observed)", {
  fit <- lmer(Reaction ~ Days + (Days | Subject), data = sleepstudy, REML = FALSE)
  psi_hat <- reconf:::get_psi_hat_lmer(fit)
  b_hat   <- as.vector(getME(fit, "beta"))
  Y <- getME(fit, "y"); X <- getME(fit, "X"); Z <- getME(fit, "Z")
  Hlist <- reconf:::get_Hlist_lmer(fit)

  # At the MLE and at a nearby non-stationary point, where the score is nonzero
  for (psi in list(psi_hat, psi_hat * c(1.3, 0.8, 1.1, 0.9))) {
    for (expected in c(TRUE, FALSE)) {
      out_q <- ll_q(psi, Y, X, Z, Hlist, REML = FALSE, b = b_hat, expected = expected)
      out_n <- ll_n(psi, Y, X, Z, Hlist, REML = FALSE, b = b_hat, expected = expected)
      expect_equal(out_n$value, out_q$value, tolerance = 1e-8)
      expect_equal(out_n$score, out_q$score, tolerance = 1e-6)
      expect_equal(out_n$inf_mat, out_q$inf_mat, tolerance = 1e-6)
    }
  }
})

test_that("n-side agrees with q-side on sleepstudy (REML)", {
  fit <- lmer(Reaction ~ Days + (Days | Subject), data = sleepstudy, REML = TRUE)
  psi_hat <- reconf:::get_psi_hat_lmer(fit)
  Y <- getME(fit, "y"); X <- getME(fit, "X"); Z <- getME(fit, "Z")
  Hlist <- reconf:::get_Hlist_lmer(fit)

  for (psi in list(psi_hat, psi_hat * c(1.3, 0.8, 1.1, 0.9))) {
    out_q <- ll_q(psi, Y, X, Z, Hlist, REML = TRUE)
    out_n <- ll_n(psi, Y, X, Z, Hlist, REML = TRUE)
    expect_equal(out_n$value, out_q$value, tolerance = 1e-8)
    expect_equal(out_n$score, out_q$score, tolerance = 1e-6)
    expect_equal(out_n$inf_mat, out_q$inf_mat, tolerance = 1e-6)
  }
})

test_that("n-side REML beta and I_b_chol match direct GLS computations", {
  fit <- lmer(Reaction ~ Days + (Days | Subject), data = sleepstudy, REML = TRUE)
  psi_hat <- reconf:::get_psi_hat_lmer(fit)
  Y <- getME(fit, "y"); X <- getME(fit, "X"); Z <- getME(fit, "Z")
  Hlist <- reconf:::get_Hlist_lmer(fit)
  r <- length(psi_hat)

  out_n <- ll_n(psi_hat, Y, X, Z, Hlist, REML = TRUE)

  # Direct dense computation of the GLS quantities
  Psi <- as.matrix(reconf:::Psi_from_H_cpp(psi_hat[-r], do.call(cbind, Hlist)))
  Sigma <- as.matrix(Z %*% Psi %*% t(Z)) + psi_hat[r] * diag(length(Y))
  XtSiX <- crossprod(X, solve(Sigma, X))
  beta_gls <- solve(XtSiX, crossprod(X, solve(Sigma, Y)))
  expect_equal(out_n$beta, as.vector(beta_gls), tolerance = 1e-6)
  expect_equal(crossprod(out_n$I_b_chol), XtSiX, tolerance = 1e-6,
               ignore_attr = TRUE)
})

# ── Path agreement, crossed design with q > n ────────────────────────────────

test_that("n-side agrees with q-side in a crossed design with q > n", {
  d <- make_crossed()
  expect_true(d$q > d$n)
  b <- c(2, -1)

  for (psi in list(d$psi, c(0.4, 1.2, 0.3))) {
    for (expected in c(TRUE, FALSE)) {
      out_q <- ll_q(psi, d$Y, d$X, d$Z, d$Hlist, REML = FALSE, b = b,
                    expected = expected)
      out_n <- ll_n(psi, d$Y, d$X, d$Z, d$Hlist, REML = FALSE, b = b,
                    expected = expected)
      expect_equal(out_n$value, out_q$value, tolerance = 1e-8)
      expect_equal(out_n$score, out_q$score, tolerance = 1e-6)
      expect_equal(out_n$inf_mat, out_q$inf_mat, tolerance = 1e-6)
    }
    out_q <- ll_q(psi, d$Y, d$X, d$Z, d$Hlist, REML = TRUE)
    out_n <- ll_n(psi, d$Y, d$X, d$Z, d$Hlist, REML = TRUE)
    expect_equal(out_n$value, out_q$value, tolerance = 1e-8)
    expect_equal(out_n$score, out_q$score, tolerance = 1e-6)
    expect_equal(out_n$inf_mat, out_q$inf_mat, tolerance = 1e-6)
  }
})

# ── Numerical derivative spot-checks on the n-side value itself ──────────────

test_that("n-side derivatives agree with numerical ones (q > n)", {
  skip_on_cran()
  skip_if_not_installed("numDeriv")
  d <- make_crossed()
  b <- c(2, -1)
  K <- reconf:::get_precomp(Y = d$Y, X = d$X, Z = d$Z, Hlist = d$Hlist,
                            method = "n_side")$K
  e <- as.vector(d$Y - d$X %*% b)

  # ML: score and observed information against the n-side value
  val_ml <- function(psi) {
    reconf:::loglik_n(K = K, psi = psi, e = e, X = d$X, get_score = FALSE,
                      get_inf = FALSE)$value
  }
  out <- reconf:::loglik_n(K = K, psi = d$psi, e = e, X = d$X, expected = FALSE)
  p <- ncol(d$X)
  r <- length(d$psi)
  expect_equal(numDeriv::grad(val_ml, d$psi),
               out$score[(p + 1):(p + r)], tolerance = 1e-4)
  expect_equal(-numDeriv::hessian(val_ml, d$psi),
               out$inf_mat[(p + 1):(p + r), (p + 1):(p + r)], tolerance = 1e-4)

  # REML: score against the n-side value
  val_reml <- function(psi) {
    reconf:::loglik_res_n(K = K, psi = psi, Y = d$Y, X = d$X,
                          get_score = FALSE, get_inf = FALSE)$value
  }
  out_reml <- reconf:::loglik_res_n(K = K, psi = d$psi, Y = d$Y, X = d$X)
  expect_equal(numDeriv::grad(val_reml, d$psi), out_reml$score,
               tolerance = 1e-4)
})

# ── Feasibility gate agreement ───────────────────────────────────────────────

test_that("n-side gate matches q-side feasibility decisions", {
  d <- make_crossed()
  b <- c(2, -1)

  # Clearly infeasible: Sigma indefinite
  psi_bad <- c(-100, 0.5, 0.7)
  out_q <- ll_q(psi_bad, d$Y, d$X, d$Z, d$Hlist, REML = TRUE)
  out_n <- ll_n(psi_bad, d$Y, d$X, d$Z, d$Hlist, REML = TRUE)
  expect_identical(out_q$value, -Inf)
  expect_identical(out_n$value, -Inf)
  expect_true(all(out_n$score == 0) && all(out_n$inf_mat == 0))

  # Mildly negative variance with Sigma still positive definite: feasible on
  # both paths and equal (the parameter space is Sigma PD, not Psi PSD)
  psi_ok <- c(1, -0.05, 0.7)
  out_q <- ll_q(psi_ok, d$Y, d$X, d$Z, d$Hlist, REML = TRUE)
  out_n <- ll_n(psi_ok, d$Y, d$X, d$Z, d$Hlist, REML = TRUE)
  expect_true(is.finite(out_q$value))
  expect_equal(out_n$value, out_q$value, tolerance = 1e-8)

  # Nonpositive error variance is infeasible even when Z Psi Z' alone is
  # positive definite (possible because q >= n); this must be caught by the
  # explicit psi_r check, not the Cholesky. Constructed so that
  # Z Psi Z' = psi_1 (I_n + 1_n 1_n') is positive definite by design.
  n2 <- 10
  Z2 <- methods::as(cbind(Matrix::Diagonal(n2), Matrix::Matrix(1, n2, 1)),
                    "generalMatrix")
  Hlist2 <- list(methods::as(Matrix::Diagonal(n2 + 1), "generalMatrix"))
  set.seed(2)
  X2 <- cbind(1, rnorm(n2))
  Y2 <- rnorm(n2)
  for (psi_r_bad in c(0, -0.1)) {
    out_n <- ll_n(c(1, psi_r_bad), Y2, X2, Z2, Hlist2, REML = TRUE)
    expect_identical(out_n$value, -Inf)
    out_n_ml <- ll_n(c(1, psi_r_bad), Y2, X2, Z2, Hlist2, REML = FALSE,
                     b = c(0, 0))
    expect_identical(out_n_ml$value, -Inf)
  }
})

# ── Precompute layer ─────────────────────────────────────────────────────────

test_that("get_precomp n_side builds K = [Z H_1 Z' ... Z H_{r-1} Z']", {
  d <- make_crossed()
  pc <- reconf:::get_precomp(Y = d$Y, X = d$X, Z = d$Z, Hlist = d$Hlist,
                             method = "n_side")
  expect_equal(dim(pc$K), c(d$n, d$n * length(d$Hlist)))
  for (j in seq_along(d$Hlist)) {
    Kj <- as.matrix(d$Z %*% d$Hlist[[j]] %*% Matrix::t(d$Z))
    expect_equal(pc$K[, (j - 1) * d$n + seq_len(d$n)], Kj, tolerance = 1e-12,
                 ignore_attr = TRUE)
  }
})

# ── Dispatch through loglikelihood() ─────────────────────────────────────────

test_that("loglikelihood dispatches on method and precomp tag", {
  d <- make_crossed()
  b <- c(2, -1)

  # auto keeps the q-side for sparse Z even when q > n (crossed intercepts:
  # Z'Z has O(n) off-diagonals however large q is) and picks the n-side only
  # for dense Z with q >= n, where the sparse path degenerates
  expect_identical(
    reconf:::get_precomp(d$Y, d$X, d$Z, Hlist = d$Hlist, method = "auto")$method,
    "q_side")
  fit <- lmer(Reaction ~ Days + (Days | Subject), data = sleepstudy)
  expect_identical(
    reconf:::get_precomp_lmer(fit, Hlist = reconf:::get_Hlist_lmer(fit))$method,
    "q_side")
  set.seed(4)
  Zd <- methods::as(Matrix::Matrix(matrix(rnorm(20 * 25), 20, 25),
                                   sparse = TRUE), "generalMatrix")
  # Two structure matrices (r = 3): the dense regime picks the n-side; with
  # a single structure matrix (r = 2) it picks the spectral path instead
  # (tested in test-spectral.R)
  Hd <- list(
    methods::as(Matrix::Diagonal(25, rep(c(1, 0), c(10, 15))), "generalMatrix"),
    methods::as(Matrix::Diagonal(25, rep(c(0, 1), c(10, 15))), "generalMatrix")
  )
  expect_identical(
    reconf:::get_precomp(rnorm(20), cbind(1, rnorm(20)), Zd, Hlist = Hd,
                         method = "auto")$method,
    "n_side")

  # Explicit method: paths agree through the wrapper, REML and ML
  for (psi in list(d$psi, c(0.4, 1.2, 0.3))) {
    out_q <- reconf:::loglikelihood(psi = psi, Y = d$Y, X = d$X, Z = d$Z,
                                    Hlist = d$Hlist, REML = TRUE,
                                    method = "q_side")
    out_n <- reconf:::loglikelihood(psi = psi, Y = d$Y, X = d$X, Z = d$Z,
                                    Hlist = d$Hlist, REML = TRUE,
                                    method = "n_side")
    expect_equal(out_n$value, out_q$value, tolerance = 1e-8)
    expect_equal(out_n$score, out_q$score, tolerance = 1e-6)
    expect_equal(out_n$inf_mat, out_q$inf_mat, tolerance = 1e-6)

    # ML without get_beta exercises the wrapper's stripping of the beta block
    out_q <- reconf:::loglikelihood(psi = psi, b = b, Y = d$Y, X = d$X,
                                    Z = d$Z, Hlist = d$Hlist, REML = FALSE,
                                    get_beta = FALSE, method = "q_side")
    out_n <- reconf:::loglikelihood(psi = psi, b = b, Y = d$Y, X = d$X,
                                    Z = d$Z, Hlist = d$Hlist, REML = FALSE,
                                    get_beta = FALSE, method = "n_side")
    expect_equal(length(out_n$score), length(psi))
    expect_equal(out_n$score, out_q$score, tolerance = 1e-6)
    expect_equal(out_n$inf_mat, out_q$inf_mat, tolerance = 1e-6)
  }

  # A supplied precomp determines the path regardless of method
  pc_n <- reconf:::get_precomp(d$Y, d$X, d$Z, Hlist = d$Hlist,
                               method = "n_side")
  out_pc <- reconf:::loglikelihood(psi = d$psi, Y = d$Y, X = d$X, Z = d$Z,
                                   Hlist = d$Hlist, REML = TRUE,
                                   precomp = pc_n, method = "q_side")
  out_n <- reconf:::loglikelihood(psi = d$psi, Y = d$Y, X = d$X, Z = d$Z,
                                  Hlist = d$Hlist, REML = TRUE,
                                  method = "n_side")
  expect_identical(out_pc, out_n)

  # Infeasible parameters keep the q-side return shapes
  out_bad <- reconf:::loglikelihood(psi = c(-100, 0.5, 0.7), Y = d$Y, X = d$X,
                                    Z = d$Z, Hlist = d$Hlist, REML = TRUE,
                                    method = "n_side")
  expect_identical(out_bad$value, -Inf)
  expect_identical(length(out_bad$score), 3L)
})

# ── End to end with q > n through the lme4 front door ────────────────────────

test_that("score_stat and ci_all_lmer work with q > n on both paths", {
  skip_on_cran()
  # Deterministic crossed design where every level occurs, so lmer keeps all
  # q = 15 + 16 = 31 > n = 30 random effects
  n <- 30
  set.seed(3)
  df <- data.frame(x = rnorm(n),
                   f1 = factor(rep(seq_len(15), each = 2)),
                   f2 = factor(rep(seq_len(16), length.out = n)))
  df$y <- 2 - df$x + rnorm(15, sd = 1)[df$f1] + rnorm(16, sd = 0.8)[df$f2] +
    rnorm(n, sd = 0.7)

  fit <- suppressMessages(lmer(
    y ~ x + (1 | f1) + (1 | f2), data = df, REML = TRUE,
    control = lmerControl(check.nobs.vs.nlev = "ignore",
                          check.nobs.vs.nRE = "ignore",
                          calc.derivs = FALSE)))

  Hlist <- reconf:::get_Hlist_lmer(fit)
  psi_hat <- reconf:::get_psi_hat_lmer(fit)
  Y <- getME(fit, "y"); X <- getME(fit, "X"); Z <- getME(fit, "Z")
  expect_gt(ncol(Z), nrow(Z))

  # Same score statistic whichever path computes it
  pc_q <- reconf:::get_precomp(Y, X, Z, REML = TRUE, Hlist = Hlist,
                               method = "q_side")
  pc_n <- reconf:::get_precomp(Y, X, Z, REML = TRUE, Hlist = Hlist,
                               method = "n_side")
  psi0 <- pmax(psi_hat, 0.1)  # away from the boundary
  st_q <- reconf:::score_stat(theta = psi0, test_idx = 1, Y = Y, X = X, Z = Z,
                              Hlist = Hlist, REML = TRUE, precomp = pc_q)
  st_n <- reconf:::score_stat(theta = psi0, test_idx = 1, Y = Y, X = X, Z = Z,
                              Hlist = Hlist, REML = TRUE, precomp = pc_n)
  expect_equal(as.vector(st_n), as.vector(st_q), tolerance = 1e-6)

  # Full CI pipeline gives the same interval on both paths
  ci_q <- ci_all_lmer(fit, test_idx = 1, method = "q_side")
  ci_n <- ci_all_lmer(fit, test_idx = 1, method = "n_side")
  expect_true(is.finite(ci_q[1, "lower"]) && is.finite(ci_q[1, "upper"]))
  expect_lte(ci_q[1, "lower"], ci_q[1, "estimate"])
  expect_gte(ci_q[1, "upper"], ci_q[1, "estimate"])
  expect_equal(unclass(ci_n), unclass(ci_q), tolerance = 1e-4)
})
