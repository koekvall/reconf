# reconf 0.2.0

## New interface

* The functions are `vc_ci()` for confidence intervals, `vc_test()` for
  score tests, and `vc_loglik()` for the log-likelihood as a function of
  the covariance parameters. Each accepts an lme4 fit or a model object
  from `vc_model()`, which does the precomputation once and can be reused
  across calls. `vc_model()` is generic, so other model sources can be
  added as methods.
* `vc_ci()` replaces `ci_lmer()` and `ci_all_lmer()`: `vc_ci(fit)` computes
  all intervals, and `parm` selects parameters by name or index. Only the
  requested intervals are computed.
* `vc_test()` replaces `score_test_lmer()`.
* `vc_loglik()` returns the value, score, and information, named by the
  parameters, at any covariance parameters in the parameter set.

## New features

* New `statistic` argument for `vc_ci()`. The default
  `"score"` is the previous behavior; `"rlrt"` inverts the profile
  likelihood ratio statistic, restricted for restricted likelihood fits, with
  nuisance parameters and the reference maximum taken over the extended
  parameter set. If the maximizer has a negative variance, the interval for
  that variance can lie below zero, and a warning is issued. Requires
  `onestep = FALSE`.
* The outward CI search now evaluates the statistic at the search origin
  instead of assuming it vanishes there. If the estimate is outside the
  confidence set, the search warns and restarts from the extended-set
  maximizer. The `nonneg` clamp is applied as a set intersection; an empty
  intersection gives `NA` bounds and a warning.
* Prior weights and offsets in lme4 fits are supported, handled exactly by
  transforming Y, X, and Z with the square-root weights.
* The outward search chooses step lengths by secant extrapolation, each
  at most twice the previous one. If the nuisance parameters from the
  previous step are infeasible at a proposed value, a feasible starting
  point is constructed instead of ending the search.
* New `method` argument selecting the computational path. The default
  `"auto"` uses the sparse q-by-q path, a dense n-by-n path when `q >= n`
  and `Z` is dense, or an `O(n)`-per-evaluation spectral path when `r = 2`.
  On a dense genomic model (n = 742, q = 1484) one evaluation drops from
  6.8 s to 0.07 s.
* `expected = FALSE` now works with the restricted likelihood: all three
  computational paths return the observed information, where previously a
  warning was issued and the expected information used. The added cost is
  of the same order as the score.
* `print` and `tidy()` (broom-style) methods for the confidence intervals.

## Changes that may affect existing code

* Intervals are returned as an object of class `vc_ci` with an added
  `estimate` column (columns `estimate`, `lower`, `upper`). Code that
  indexed the previous two-column output by position should select columns
  by name.
* Interval rows are named `var_x|g` for a variance, `cov_y.x|g` for a
  covariance, and `var_Residual` for the error variance, following lme4's
  names for standard deviations and correlations.
* `test_idx` and `known_idx` are renamed `parm` and `known`; both index the
  covariance parameters, also for maximum likelihood fits, and parameters
  in `known` are held at their estimates. The `REML` argument is removed;
  the likelihood is the one the fit maximized.
* `vc_test()` returns an object of class `htest` with components
  `statistic`, `parameter`, `p.value`, and `null.value`, instead of a named
  vector with elements `stat`, `p_val`, and `df`. It takes `null_value`,
  the values of the tested parameters (default 0), in place of the full
  parameter vector `theta_null`. The arguments `efficient` and `profile` are
  removed; the statistic is the profile score statistic, standardized by the
  efficient information.
* For maximum likelihood fits the fixed effects are profiled out of the
  likelihood by generalized least squares, which leaves the statistics
  unchanged.
* The log-likelihood value now includes all normalizing constants and
  equals `logLik()` from an equivalent lme4 fit (weighted fits: up to
  half the sum of log weights).

## Bug fixes

* The log-likelihood could be finite at parameters where
  `Sigma = Z Psi Z' + psi_r I` is not positive definite, corrupting CI
  searches. Feasibility is now decided by an exact Cholesky test.

## Performance

* Large speedups; `vc_ci()` on the FEV1 example runs in about 0.6 s
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
