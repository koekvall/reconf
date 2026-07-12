# Coverage simulation for the README figure.
#
# Model: y_ij = 1 + b_i + e_ij with b_i ~ N(0, psi1), e_ij ~ N(0, 1),
# i = 1, ..., n_grp groups, j = 1, ..., m observations per group.
#
# Coverage is estimated by test-interval duality: psi1 is inside a
# confidence interval obtained by inverting a test if and only if that test
# fails to reject at the true value, so no interval needs to be
# constructed. For each replicate, with the nuisance error variance
# profiled at psi1 fixed to the truth:
#   - Wald:    |psi1_hat - psi1| <= z * SE, SE from the expected information
#              at the (refined) REML estimates;
#   - profile: 2 * (l_max - l_profile(psi1)) <= qchisq(0.95, 1), the
#              likelihood-ratio inversion underlying profile intervals;
#   - score:   the efficient score statistic at psi1 <= qchisq(0.95, 1),
#              the test inverted by reconf::ci_lmer.
#
# Run from the package root:
#   PILOT=1 Rscript scripts/coverage_simulation.R   # small run, sanity check
#   Rscript scripts/coverage_simulation.R           # full run + figure
#
# The full run writes man/figures/coverage.png and prints the coverage table.

suppressMessages({
  library(lme4)
  library(reconf)
  library(parallel)
})

pilot <- nzchar(Sys.getenv("PILOT"))
n_grp <- 30L
m     <- 5L
level <- 0.95
z     <- qnorm((1 + level) / 2)
crit  <- qchisq(level, df = 1)
psi1_grid <- if (pilot) c(0.05, 0.5) else
  c(0.01, 0.02, 0.05, 0.1, 0.15, 0.2, 0.35, 0.5, 0.75, 1, 1.5)
N     <- if (pilot) 100L else 2000L
cores <- max(1L, detectCores() - 2L)

grp <- factor(rep(seq_len(n_grp), each = m))

one_rep <- function(psi1, rep_id) {
  set.seed(100000L * round(1000 * psi1) + rep_id)
  y <- 1 + rnorm(n_grp, sd = sqrt(psi1))[grp] + rnorm(n_grp * m)
  tryCatch({
    fit <- suppressMessages(suppressWarnings(
      lmer(y ~ 1 + (1 | grp), REML = TRUE,
           control = lmerControl(check.conv.singular = "ignore",
                                 calc.derivs = FALSE))
    ))
    Y <- getME(fit, "y"); X <- getME(fit, "X"); Z <- getME(fit, "Z")
    Hlist <- reconf:::get_Hlist_lmer(fit)
    pc <- reconf:::get_precomp(Y, X, Z, REML = TRUE, Hlist = Hlist)

    # lme4's constrained REML estimates: the quantities practitioners use
    psi_hat <- reconf:::get_psi_hat_lmer(fit)

    # Value and expected information at the estimates
    ll_hat <- reconf:::loglikelihood(
      psi = psi_hat, Y = Y, X = X, Z = Z, Hlist = Hlist, REML = TRUE,
      get_val = TRUE, get_score = FALSE, get_inf = TRUE,
      precomp = pc, check = FALSE)

    # Wald on the variance scale
    wald_cover <- abs(psi_hat[1] - psi1) <=
      z * sqrt(solve(ll_hat$inf_mat)[1, 1])

    # Nuisance error variance profiled at psi1 fixed to the truth
    prof <- reconf:::maximize_loglik(
      start_val = c(psi1, psi_hat[2]), opt_idx = 2L,
      Y = Y, X = X, Z = Z, Hlist = Hlist, REML = TRUE, precomp = pc,
      check = FALSE, warn_nonconv = FALSE, iterlim = 1000L)

    # Profile-likelihood interval coverage via the LRT at the truth
    # (clamped at zero: optimizer tolerance can make the difference
    # marginally negative when psi1 is near the estimate)
    prof_cover <- max(0, 2 * (ll_hat$value - prof$value)) <= crit

    # Score interval coverage via the efficient score test at the truth
    st <- reconf:::score_stat(
      theta = prof$arg, test_idx = 1L,
      Y = Y, X = X, Z = Z, Hlist = Hlist, REML = TRUE,
      precomp = pc, check = FALSE)
    score_cover <- as.numeric(st) <= crit

    c(wald_cover, prof_cover, score_cover)
  }, error = function(e) c(NA, NA, NA))
}

t0 <- Sys.time()
results <- lapply(psi1_grid, function(psi1) {
  out <- mclapply(seq_len(N), function(r) one_rep(psi1, r), mc.cores = cores)
  do.call(rbind, out)
})
cat(sprintf("simulation: %.1f min on %d cores\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins")), cores))

tab <- do.call(rbind, lapply(seq_along(psi1_grid), function(i) {
  b <- results[[i]]
  n_ok <- sum(!is.na(b[, 1]))
  cov <- colMeans(b, na.rm = TRUE)
  data.frame(psi1    = psi1_grid[i],
             wald    = cov[1], profile = cov[2], score = cov[3],
             se_wald    = sqrt(cov[1] * (1 - cov[1]) / n_ok),
             se_profile = sqrt(cov[2] * (1 - cov[2]) / n_ok),
             se_score   = sqrt(cov[3] * (1 - cov[3]) / n_ok),
             n_fail  = sum(is.na(b[, 1])))
}))
print(tab[, c("psi1", "wald", "profile", "score", "n_fail")], digits = 3)

if (pilot) quit(save = "no")

saveRDS(list(tab = tab, n_grp = n_grp, m = m, N = N),
        "scripts/coverage_results.rds")

## Figure ---------------------------------------------------------------------
# Palette: validated categorical slots (blue, aqua, yellow) on white; the
# two lower-contrast hues are relieved by direct labels on every line.
col_score <- "#2a78d6"
col_wald  <- "#1baf7a"
col_prof  <- "#eda100"
ink       <- "#1a1a19"
ink2      <- "#6f6e66"
grid_col  <- "#e8e7e2"

dir.create("man/figures", showWarnings = FALSE, recursive = TRUE)
png("man/figures/coverage.png", width = 2100, height = 1350, res = 300)
op <- par(mar = c(4.2, 4.2, 3.2, 7.5), mgp = c(2.6, 0.7, 0), xpd = FALSE)

ymin <- floor(min(tab[, c("wald", "profile", "score")]) * 40) / 40
ylim <- c(ymin, 1)
plot(NA, xlim = range(psi1_grid), ylim = ylim, log = "x", axes = FALSE,
     xlab = expression("True random-intercept variance " * psi[1]),
     ylab = "Empirical coverage")
title(main = "Coverage of nominal 95% confidence intervals for a variance",
      adj = 0, cex.main = 1.05, col.main = ink)
mtext(sprintf(
  "%d groups of %d observations, error variance 1; %d replicates and ±2 Monte Carlo SE bands per point",
  n_grp, m, N),
  side = 3, line = 0.2, adj = 0, cex = 0.75, col = ink2)

abline(h = level, col = ink2, lwd = 1, lty = 2)
xlabs <- c(0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1)  # thinned to avoid collisions
axis(1, at = xlabs, labels = xlabs, col = NA, col.ticks = ink2,
     col.axis = ink2, cex.axis = 0.8, tcl = -0.25)
axis(2, at = c(seq(ymin, 1, by = 0.025)), las = 1, col = NA,
     col.ticks = ink2, col.axis = ink2, cex.axis = 0.8, tcl = -0.25)

# Semi-transparent versions of the line colors for the ribbons
faint <- function(col) adjustcolor(col, alpha.f = 0.12)
draw <- function(y, se, col) {
  polygon(c(psi1_grid, rev(psi1_grid)),
          c(pmin(y + 2 * se, 1), rev(y - 2 * se)),
          col = faint(col), border = NA)
  lines(psi1_grid, y, col = col, lwd = 2)
}
draw(tab$wald,    tab$se_wald,    col_wald)
draw(tab$profile, tab$se_profile, col_prof)
draw(tab$score,   tab$se_score,   col_score)

par(xpd = NA)
lab <- function(y, txt, col) {
  text(max(psi1_grid) * 1.12, y[length(y)], txt, col = col, adj = 0,
       cex = 0.85, font = 2)
}
lab(tab$score,   "Score (reconf)", col_score)
lab(tab$wald,    "Wald",           col_wald)
lab(tab$profile, "Profile",        col_prof)
par(op)
dev.off()
cat("wrote man/figures/coverage.png\n")
