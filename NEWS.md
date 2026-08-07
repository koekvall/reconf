# reconf 0.2.0

## New features

* New `statistic` argument for `ci_lmer()` and `ci_all_lmer()`. The default
  `"score"` is the previous behavior; `"rlrt"` instead inverts the profile
  restricted likelihood ratio statistic on the extended parameter set, with
  nuisance parameters maximized over that set as for the score statistic and
  the reference maximum taken over the same set. For a variance whose
  extended-set estimate is negative the interval is centered below zero and
  can exclude the estimate or contain no nonnegative values; warnings flag
  both cases, and with `nonneg = TRUE` an interval with no nonnegative
  values is reported as `NA` bounds (possible with probability about
  (1 - level)/2 under a true zero variance; suggestive of an inadequate
  covariance structure if it recurs). Experimental: these intervals have no
  supporting theory at present. Requires `onestep = FALSE`; `ci_all_lmer()`
  computes the reference maximum once and shares it across parameters.
* The outward CI search now evaluates the statistic at the search origin
  instead of assuming it vanishes there. If the estimate is outside the
  confidence set, which can happen for estimates on the boundary of the
  parameter space, the search warns and restarts from the extended-set
  maximizer, where the statistic vanishes. The `nonneg` clamp is applied
  as a set intersection, so the lower bound cannot exceed the upper bound;
  an empty intersection is reported as `NA` bounds for either statistic.

* New exported `score_test_all_lmer()`: per-parameter score tests against
  zero, with named rows, p-values, and `print`/`tidy` methods.
* `ci_lmer()` and `ci_all_lmer()` gain prior-weight and offset support,
  handled exactly by transforming Y, X, and Z with the square-root weights.
* New `accelerate` argument (default `TRUE`) for a secant-accelerated
  outward search; set `FALSE` for the previous fixed-step search.
* New `method` argument selecting the computational path. The default
  `"auto"` uses the sparse q-by-q path, a dense n-by-n path when `q >= n`
  and `Z` is dense, or an `O(n)`-per-evaluation spectral path when `r = 2`.
  On a dense genomic model (n = 742, q = 1484) one evaluation drops from
  6.8 s to 0.07 s.
* `expected = FALSE` now works with the restricted likelihood: all three
  computational paths return the observed information, where previously a
  warning was issued and the expected information used. The added cost is
  of the same order as the score.
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
