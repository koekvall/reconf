# Benchmark the log-likelihood with different get_* flag combinations, for
# the restricted likelihood and the likelihood maximized over the fixed
# effects. Uses the unvalidated function inside vc_model(); vc_loglik() adds
# argument checks on each call. Install a release build first (R CMD
# INSTALL, or pkgbuild::compile_dll(debug = FALSE)); debug builds of the
# Eigen code are about ten times slower.
#
# Run from the package root: Rscript scripts/benchmark_loglik.R

library(reconf)
library(lme4)

data(fev1)
n_rep <- 500

for (reml in c(FALSE, TRUE)) {
  fit <- lmer(logfev1 ~ age + ht + baseage + baseht + (age | id), data = fev1,
              REML = reml)
  model <- vc_model(fit)
  ll <- model$ll
  psi <- model$psi_hat
  cat("===", if (reml) "REML" else "ML", "(", n_rep, "reps, ms per call) ===\n")
  cases <- list(
    "value only"          = function() ll(psi, get_score = FALSE, get_inf = FALSE),
    "score only"          = function() ll(psi, get_val = FALSE, get_inf = FALSE),
    "value and score"     = function() ll(psi, get_inf = FALSE),
    "expected information" = function() ll(psi),
    "observed information" = function() ll(psi, expected = FALSE))
  for (nm in names(cases)) {
    f <- cases[[nm]]
    t <- system.time(for (i in seq_len(n_rep)) f())[["elapsed"]]
    cat(sprintf("  %-22s %.3f\n", nm, 1e3 * t / n_rep))
  }
}
