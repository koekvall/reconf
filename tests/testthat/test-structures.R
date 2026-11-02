library(lme4)

# Covariance structures need lme4 2.0
has_structures <- "getReCovs" %in% getNamespaceExports("lme4")

# H as reconf 0.2.0 built it, from Lambdat and Lind; correct for
# unstructured terms
Hlist_from_Lind <- function(fit) {
  H <- getME(fit, "Lambdat")
  param_idx <- getME(fit, "Lind")
  q <- nrow(H)
  col_idx <- rep(seq_len(q), diff(H@p))
  row_idx <- H@i + 1L
  lapply(seq_len(getME(fit, "m")), function(i) {
    keep <- which(param_idx == i)
    M <- Matrix::sparseMatrix(i = row_idx[keep], j = col_idx[keep],
                              x = 1, dims = c(q, q))
    Matrix::forceSymmetric(M, uplo = "U")
  })
}

set.seed(11)
n_grp <- 60
dat <- data.frame(g = gl(n_grp, 8), h = factor(rep(1:12, 40)),
                  x1 = rnorm(n_grp * 8), x2 = rnorm(n_grp * 8))
b <- matrix(rnorm(3 * n_grp), n_grp) %*% chol(matrix(c(1, 0.4, 0.4,
                                                        0.4, 0.8, 0.4,
                                                        0.4, 0.4, 0.6), 3))
dat$y <- with(dat, 1 + x1 + b[g, 1] + b[g, 2] * x1 + b[g, 3] * x2 +
                rnorm(12)[h] + rnorm(n_grp * 8))
fit_s <- function(re) {
  suppressMessages(suppressWarnings(
    lmer(stats::as.formula(paste("y ~ x1 + x2 +", re)), data = dat)))
}

# ── Unstructured terms: unchanged from reconf 0.2.0 ─────────────────────────

test_that("unstructured terms give the structure matrices of reconf 0.2.0", {
  fits <- list(
    lmer(Reaction ~ Days + (Days | Subject), data = sleepstudy),
    lmer(Reaction ~ Days + (Days || Subject), data = sleepstudy),
    fit_s("(1 + x1 + x2 | g) + (1 | h)"))
  for (fit in fits) {
    s <- reconf:::.lmer_structure(fit)
    expect_identical(s$Hlist, Hlist_from_Lind(fit))
    vc <- as.data.frame(VarCorr(fit), order = "lower.tri")
    expect_equal(s$psi_hat, vc$vcov, tolerance = 1e-12)
  }
  expect_identical(
    vc_model(fits[[2]])$names,
    c("var_(Intercept)|Subject", "var_Days|Subject.1", "var_Residual"))
})

# ── Linear structures ────────────────────────────────────────────────────────

test_that("diag() terms have no covariances and match the || fit", {
  skip_if_not(has_structures)
  fit_d <- lmer(Reaction ~ Days + diag(Days | Subject), data = sleepstudy)
  fit_u <- lmer(Reaction ~ Days + (Days || Subject), data = sleepstudy)
  m <- vc_model(fit_d)
  expect_identical(m$names, c("var_(Intercept)|Subject", "var_Days|Subject",
                              "var_Residual"))
  expect_equal(m$psi_hat, vc_model(fit_u)$psi_hat, tolerance = 1e-6)
  # Same model as the || fit; reconf 0.2.0 gave infinite intervals
  ci_d <- vc_ci(fit_d)
  expect_true(all(is.finite(as.matrix(ci_d[, c("lower", "upper")]))))
  expect_equal(unname(as.matrix(ci_d[, c("lower", "upper")])),
               unname(as.matrix(vc_ci(fit_u)[, c("lower", "upper")])),
               tolerance = 1e-4)
})

test_that("diag(hom = TRUE) and cs(hom = TRUE) share parameters", {
  skip_if_not(has_structures)
  fit_dh <- fit_s("diag(1 + x1 + x2 | g, hom = TRUE)")
  m <- vc_model(fit_dh)
  expect_identical(m$names, c("var_*|g", "var_Residual"))
  expect_equal(m$psi_hat[1], as.numeric(VarCorr(fit_dh)$g[1, 1]),
               tolerance = 1e-10)
  expect_equal(as.matrix(m$precomp$H), as.matrix(Matrix::Diagonal(180)))

  fit_ch <- fit_s("cs(1 + x1 + x2 | g, hom = TRUE)")
  m <- vc_model(fit_ch)
  expect_identical(m$names, c("var_*|g", "cov_*|g", "var_Residual"))
  S <- VarCorr(fit_ch)$g
  expect_equal(m$psi_hat[1:2], c(S[1, 1], S[2, 1]), tolerance = 1e-10)
  expect_identical(m$structure$is_var, c(TRUE, FALSE, TRUE))
  # Score near zero at the estimates
  expect_lt(max(abs(vc_loglik(m)(m$psi_hat)$score)), 1e-2)
  ci <- vc_ci(m)
  expect_true(all(ci[, "lower"] < m$psi_hat & m$psi_hat < ci[, "upper"]))
})

test_that("two-effect cs() and ar1() terms reduce to linear structures", {
  skip_if_not(has_structures)
  # With two effects, cs() is us
  fit_cs <- fit_s("cs(1 + x1 | g)")
  fit_us <- fit_s("(1 + x1 | g)")
  expect_identical(vc_model(fit_cs)$names, vc_model(fit_us)$names)
  expect_equal(unname(as.matrix(vc_ci(fit_cs)[, c("lower", "upper")])),
               unname(as.matrix(vc_ci(fit_us)[, c("lower", "upper")])),
               tolerance = 1e-4)
  # With two effects, ar1() is cs(hom = TRUE); its one covariance keeps its
  # own name
  fit_ar <- fit_s("ar1(1 + x1 | g)")
  fit_ch <- fit_s("cs(1 + x1 | g, hom = TRUE)")
  expect_identical(vc_model(fit_ar)$names,
                   c("var_*|g", "cov_x1.(Intercept)|g", "var_Residual"))
  expect_equal(vc_model(fit_ar)$psi_hat, vc_model(fit_ch)$psi_hat,
               tolerance = 1e-4)
})

test_that("terms with different structures combine", {
  skip_if_not(has_structures)
  fit <- fit_s("diag(1 + x1 | g, hom = TRUE) + cs(x1 + x2 | h, hom = TRUE)")
  m <- vc_model(fit)
  expect_identical(m$names, c("var_*|g", "var_*|h", "cov_*|h",
                              "var_Residual"))
  expect_true(all(is.finite(m$psi_hat)))
})

# ── Unsupported structures and layout checks ────────────────────────────────

test_that("terms not linear in their variances and covariances are an error", {
  skip_if_not(has_structures)
  expect_error(vc_model(fit_s("cs(1 + x1 + x2 | g)")),
               "cs\\(hom = FALSE\\) with 3 effects, which is not linear")
  expect_error(vc_model(fit_s("ar1(1 + x1 + x2 | g)")),
               "ar1\\(hom = TRUE\\) with 3 effects, which is not linear")
})

test_that(".term_pattern refuses unknown structures", {
  expect_null(reconf:::.term_pattern(2L, NA_character_, NA))
  expect_null(reconf:::.term_pattern(3L, "diag", NA))
  expect_null(reconf:::.term_pattern(3L, "cs", FALSE))
  expect_identical(reconf:::.term_pattern(1L, "ar1", TRUE),
                   matrix(1L, 1, 1))
  expect_identical(reconf:::.term_pattern(3L, "us", NA),
                   matrix(c(1L, 2L, 3L, 2L, 4L, 5L, 3L, 5L, 6L), 3))
})

test_that("a structure that does not reproduce the fit is an error", {
  fit <- lmer(Reaction ~ Days + (Days | Subject), data = sleepstudy)
  s <- reconf:::.lmer_structure(fit)
  Psi_fit <- sigma(fit)^2 * Matrix::tcrossprod(getME(fit, "Lambda"))
  psi <- s$psi_hat[-4]
  expect_true(reconf:::.check_reproduces(Psi_fit, psi, s$Hlist))
  # Covariance read as zero, as by a diag() pattern
  expect_error(reconf:::.check_reproduces(Psi_fit, psi * c(1, 0, 1),
                                          s$Hlist),
               "laid out differently")
  # Parameters in the wrong order
  expect_error(reconf:::.check_reproduces(Psi_fit, psi[3:1], s$Hlist),
               "laid out differently")
})

test_that("the check compares each entry on its own scale", {
  # Variance 1e4 next to a two-level cs(hom = TRUE) term with variance 1e-4
  # and covariance 5e-5
  sym <- function(M) Matrix::forceSymmetric(as(M, "CsparseMatrix"))
  P   <- reconf:::.term_pattern(3L, "cs", TRUE)
  Hk  <- function(k) sym(Matrix::bdiag(0, (P == k) * 1, (P == k) * 1))
  Hlist <- list(sym(Matrix::bdiag(1, matrix(0, 6, 6))), Hk(1), Hk(2))
  psi <- c(1e4, 1e-4, 5e-5)
  blk <- matrix(5e-5, 3, 3)
  diag(blk) <- 1e-4
  expect_true(reconf:::.check_reproduces(sym(Matrix::bdiag(1e4, blk, blk)),
                                         psi, Hlist))
  # One covariance off by less than 1e-8 of the largest variance
  bad <- blk
  bad[3, 1] <- bad[1, 3] <- 2e-5
  expect_error(reconf:::.check_reproduces(sym(Matrix::bdiag(1e4, blk, bad)),
                                          psi, Hlist),
               "laid out differently")
  Psi_nan <- sym(Matrix::bdiag(1e4, blk, blk))
  Psi_nan[1, 1] <- NaN
  expect_error(reconf:::.check_reproduces(Psi_nan, psi, Hlist), "not finite")
})

test_that("a fit with residual standard deviation NaN is an error", {
  # As in some lme4 2.0-6 fits
  fit <- lmer(Reaction ~ Days + (Days | Subject), data = sleepstudy)
  fit@devcomp$cmp[["sigmaREML"]] <- NaN
  expect_error(vc_model(fit), "residual standard deviation of NaN")
})

test_that("vc_model setup checks the lengths of psi_hat and names", {
  fit <- lmer(Reaction ~ Days + (Days | Subject), data = sleepstudy)
  s <- reconf:::.lmer_structure(fit)
  expect_error(reconf:::.make_setup(Y = getME(fit, "y"), X = getME(fit, "X"),
                                    Z = getME(fit, "Z"), Hlist = s$Hlist,
                                    REML = TRUE, psi_hat = s$psi_hat[-1]),
               "length r = 4")
})
