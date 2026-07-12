# reconf 0.2.0

## New features

* New exported `score_test_all_lmer()`: per-parameter score tests against
  zero, with named rows, p-values, and `print`/`tidy` methods.
* `ci_lmer()` and `ci_all_lmer()` gain prior-weight and offset support,
  handled exactly by transforming Y, X, and Z with the square-root weights.
* New `accelerate` argument (default `TRUE`) for a secant-accelerated
  outward search; set `FALSE` for the previous fixed-step search.
* `tidy()` methods (broom-style) and `print` methods for the confidence
  interval and score-test objects.

## Changes that may affect existing code

* `ci_lmer()` and `ci_all_lmer()` now return an object of class
  `reconf_ci` with an added `estimate` column (columns `estimate`,
  `lower`, `upper`). Code that indexed the previous two-column output by
  position should select columns by name.
* The log-likelihood value now includes all normalizing constants and
  equals `logLik()` from an equivalent lme4 fit (weighted fits: up to
  half the sum of log weights).

## Bug fixes

* The log-likelihood could be finite at parameters where
  `Sigma = Z Psi Z' + psi_r I` is not positive definite, corrupting CI
  searches. Feasibility is now decided by an exact Cholesky test.
* Score tests of covariance parameters in correlated-random-effect models
  used a wrong term lookup; profiling now also starts nuisance parameters
  from their estimates, fixing non-convergence on large response scales.

## Performance

* Large speedups; `ci_all_lmer()` on the FEV1 example runs in about 0.6 s
  instead of 43 s, with identical intervals. Sparsity-aware solves via
  Matrix, no dense q-by-q matrices in the restricted information, cached
  symbolic factorizations, and reuse of precomputed quantities.

## Infrastructure

* GitHub Actions R CMD check (Linux, macOS, Windows).
* Internal `loglik()` and `loglik_res()` take the matrix A and its
  log-determinant as arguments. Dead code removed and pure-R reference
  implementations moved to the test suite.

# reconf 0.1

* Initial release.
* `ci_lmer()`: score-based confidence interval for a single covariance parameter
  in a linear mixed model fitted with lme4.
* `ci_all_lmer()`: score-based confidence intervals for all (or a subset of)
  covariance parameters.
* `score_test_lmer()`: score test for covariance parameters against a user-supplied
  null hypothesis.
