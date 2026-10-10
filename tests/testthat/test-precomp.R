library(Matrix)

# The feasibility check of the q-side path uses a q x q factor R with
# R'R = Z'Z when one is available and Z itself otherwise (get_precomp,
# .sigma_pd). Three cases: the sparse Cholesky of Z'Z; the QR fallback when
# the Cholesky fails; and no factor when Z is wider than tall, where
# Matrix's sparse QR does not apply. In each case the feasibility decisions
# and the likelihood must agree with the dense n-side path, which decides
# feasibility by a Cholesky of Sigma.

# Random intercepts for 8 groups of 5 with the first column of Z
# duplicated, so that Z'Z (9 x 9) is exactly singular while n = 40 > q = 9
make_dup_col <- function(seed = 2) {
  set.seed(seed)
  n <- 40
  g <- factor(rep(1:8, each = 5))
  Z0 <- Matrix::t(Matrix::fac2sparse(g))
  Z <- methods::as(Matrix::cbind2(Z0, Z0[, 1]), "generalMatrix")
  H <- methods::as(Matrix::Diagonal(ncol(Z)), "generalMatrix")
  X <- cbind(1, rnorm(n))
  Y <- as.vector(X %*% c(2, -1) + Z %*% rnorm(ncol(Z)) + rnorm(n))
  # Z H Z' has largest eigenvalue 10 (group 1, whose column is doubled), so
  # Sigma = psi_1 Z H Z' + psi_2 I is positive definite iff psi_1 > -psi_2/10
  list(Y = Y, X = X, Z = Z, Hlist = list(H), psi_bound = -1 / 10)
}

# Crossed random intercepts with q = 12 + 20 > n = 30, as in test-nside.R
make_wide <- function(n = 30, l1 = 12, l2 = 20, seed = 1) {
  set.seed(seed)
  f1 <- factor(sample.int(l1, n, replace = TRUE), levels = seq_len(l1))
  f2 <- factor(sample.int(l2, n, replace = TRUE), levels = seq_len(l2))
  Z <- cbind(Matrix::t(Matrix::fac2sparse(f1, drop.unused.levels = FALSE)),
             Matrix::t(Matrix::fac2sparse(f2, drop.unused.levels = FALSE)))
  Z <- methods::as(Z, "generalMatrix")
  q <- l1 + l2
  Hlist <- list(
    methods::as(Matrix::Diagonal(q, c(rep(1, l1), rep(0, l2))),
                "generalMatrix"),
    methods::as(Matrix::Diagonal(q, c(rep(0, l1), rep(1, l2))),
                "generalMatrix")
  )
  X <- cbind(1, rnorm(n))
  u <- c(rnorm(l1), rnorm(l2, sd = 0.7))
  Y <- as.vector(X %*% c(2, -1) + Z %*% u + rnorm(n, sd = 0.8))
  # Z H_1 Z' has eigenvalues equal to the level counts of f1 and Z H_2 Z' is
  # positive semidefinite, so with psi_3 = 0.7, Sigma is positive definite
  # whenever psi_1 > -0.7 / max(table(f1))
  list(Y = Y, X = X, Z = Z, Hlist = Hlist, psi_bound = -0.7 / max(table(f1)))
}

# The feasibility check of .sigma_pd at psi, with the factor R or, if NULL,
# with Z
gate_at <- function(d, psi, R) {
  r <- length(psi)
  H <- methods::as(do.call(cbind, d$Hlist), "generalMatrix")
  Psi_r <- reconf:::Psi_from_H_cpp(psi_mr = psi[-r] / psi[r], H = H)
  reconf:::.sigma_pd(Psi_r, R = R, Z = d$Z)
}

# Value, score, and information on both paths at psi; the feasibility
# decisions must match, and the results at a feasible psi
expect_paths_agree <- function(d, psi, feasible) {
  out_q <- reconf:::make_loglik(d$Y, d$X, d$Z, d$Hlist, REML = TRUE,
                                method = "q_side")(psi)
  out_n <- reconf:::make_loglik(d$Y, d$X, d$Z, d$Hlist, REML = TRUE,
                                method = "n_side")(psi)
  expect_identical(is.finite(out_q$value), feasible)
  expect_identical(is.finite(out_n$value), feasible)
  if (feasible) {
    expect_equal(out_q$value, out_n$value, tolerance = 1e-8)
    expect_equal(out_q$score, out_n$score, tolerance = 1e-6)
    expect_equal(out_q$inf_mat, out_n$inf_mat, tolerance = 1e-6)
  }
}

test_that("the QR fallback gives R'R = Z'Z when the Cholesky of Z'Z fails", {
  d <- make_dup_col()
  ZtZ <- as.matrix(Matrix::crossprod(d$Z))
  R <- reconf:::.qr_factor(d$Z)
  expect_equal(as.matrix(Matrix::crossprod(R)), ZtZ, tolerance = 1e-12)
  # get_precomp supplies a factor by one route or the other
  pc <- reconf:::get_precomp(d$Y, d$X, d$Z, Hlist = d$Hlist, method = "q_side")
  expect_false(is.null(pc$R))
  expect_equal(as.matrix(Matrix::crossprod(pc$R)), ZtZ, tolerance = 1e-12)
  # With the QR factor, the check decides as the dense Cholesky of Sigma:
  # Psi indefinite with Sigma positive definite, and Sigma indefinite
  psi_in  <- c(0.5 * d$psi_bound, 1)
  psi_out <- c(2 * d$psi_bound, 1)
  expect_true(gate_at(d, c(1, 1), R))
  expect_true(gate_at(d, psi_in, R))
  expect_false(gate_at(d, psi_out, R))
  expect_paths_agree(d, c(1, 1), feasible = TRUE)
  expect_paths_agree(d, psi_in, feasible = TRUE)
  expect_paths_agree(d, psi_out, feasible = FALSE)
})

test_that("wide Z (q > n) works without a factor and matches the n-side path", {
  d <- make_wide()
  expect_gt(ncol(d$Z), nrow(d$Z))
  expect_null(reconf:::.qr_factor(d$Z))
  # The check with Z itself, as used when no factor is available
  psi_in  <- c(0.5 * d$psi_bound, 0.5, 0.7)
  psi_out <- c(-5, 0.5, 0.7)
  expect_true(gate_at(d, c(1, 0.5, 0.7), NULL))
  expect_true(gate_at(d, psi_in, NULL))
  expect_false(gate_at(d, psi_out, NULL))
  expect_paths_agree(d, c(1, 0.5, 0.7), feasible = TRUE)
  expect_paths_agree(d, psi_in, feasible = TRUE)
  expect_paths_agree(d, psi_out, feasible = FALSE)
})
