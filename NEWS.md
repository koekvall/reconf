# reconf 0.1.1

* Fixed: the log-likelihood could be finite at parameters where
  `Sigma = Z Psi Z' + psi_r I` is not positive definite, corrupting CI
  searches. Feasibility is now decided by an exact Cholesky test.
* Fixed: models with prior weights or an offset silently used the unweighted
  likelihood. They are now handled exactly by transforming Y, X, and Z with
  the square-root weights.
* The log-likelihood value now includes all constants and equals `logLik()`
  from lme4 (weighted fits: up to half the sum of log weights).
* Large speedups; `ci_all_lmer()` on the FEV1 example runs in 0.6 s instead
  of 43 s, with identical intervals. Main changes: sparsity-aware solves via
  Matrix, no dense q x q matrices in the restricted information, cached
  symbolic factorizations, reuse of precomputed quantities, and a
  secant-accelerated outward CI search (`accelerate = FALSE` restores the
  fixed-step search; both resolve bounds to within SE/40).
* Internal `loglik()` and `loglik_res()` now take the matrix A and its
  log-determinant as arguments; consistent return names in `res_ll()`;
  documentation fixes.

# reconf 0.1

* Initial release.
* `ci_lmer()`: score-based confidence interval for a single covariance parameter
  in a linear mixed model fitted with lme4.
* `ci_all_lmer()`: score-based confidence intervals for all (or a subset of)
  covariance parameters.
* `score_test_lmer()`: score test for covariance parameters against a user-supplied
  null hypothesis.
