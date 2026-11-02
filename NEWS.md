# reconf 0.2.1

* Builds with Eigen 5 (#4, @eddelbuettel).
* Supports the lme4 2.0 covariance structures `diag()`, `cs(hom = TRUE)`,
  and `cs()` and `ar1()` with two effects. A variance or covariance shared
  by a term's effects is one parameter, named `var_*|g` or `cov_*|g`. With
  three or more effects, `cs()` with separate variances and `ar1()` are an
  error. Version 0.2.0 misread all of these terms and could return wrong
  intervals without an error.
* `vc_model()` stops if the covariance structure does not reproduce the
  fitted covariance matrix, or if lme4 reports a residual standard
  deviation that is not finite, as lme4 2.0-6 can.
* `vc_test()` and `vc_ci()` warn when an internal maximization stops away
  from a root of the score.
* No spurious "NaNs produced" warnings in nearly singular fits.

# reconf 0.2.0

* First CRAN release.

## New interface

* `vc_ci()` computes confidence intervals, `vc_test()` score tests, and
  `vc_loglik()` the log-likelihood, score, and information as functions of
  the covariance parameters. Each takes an lme4 fit or a model from the
  generic `vc_model()`, which does the precomputation once.
* `vc_ci()` replaces `ci_lmer()` and `ci_all_lmer()`, and `vc_test()`
  replaces `score_test_lmer()`. `parm` selects parameters by name or index.

## Changes that may affect existing code

* `vc_ci()` returns an object of class `vc_ci` with columns `estimate`,
  `lower`, and `upper`. Code that indexed the old two-column output by
  position should select columns by name.
* Parameters are named `var_x|g` for a variance, `cov_y.x|g` for a
  covariance, and `var_Residual` for the error variance.
* `test_idx` and `known_idx` are renamed `parm` and `known`. Both index the
  covariance parameters, also for maximum likelihood fits. Parameters in
  `known` are held at their estimates.
* The `REML` argument is removed. The likelihood is the one the fit
  maximized.
* `vc_test()` returns an `htest` object instead of a named vector, and takes
  `null_value`, the values of the tested parameters (default 0), in place of
  `theta_null`. The arguments `efficient` and `profile` are removed. The
  statistic is the profile score statistic, standardized by the efficient
  information.
* Log-likelihood values include all normalizing constants and equal
  `logLik()` of the lme4 fit, up to half the sum of the log weights in
  weighted fits.

## New features

* `vc_ci(statistic = "rlrt")` inverts the profile likelihood ratio
  statistic, restricted for restricted likelihood fits, with nuisance
  parameters and the reference maximum taken over the extended parameter
  set. If the maximizer has a negative variance, the interval for that
  variance can extend below zero, with a warning. Requires
  `onestep = FALSE`. The default is `statistic = "score"`.
* `method` selects the computational path. The default `"auto"` uses a
  sparse q-by-q path, a dense n-by-n path when `q >= n` and `Z` is dense,
  or a spectral path with cost O(n) per evaluation when `r = 2`. On a dense
  genomic model (n = 742, q = 1484) one evaluation takes 0.07 s instead of
  6.8 s.
* `expected = FALSE` gives the observed information also for the
  restricted likelihood. Previously the expected information was used, with
  a warning.
* Prior weights and offsets in lme4 fits are supported.
* `print()` and `tidy()` methods for `vc_ci` objects.
* The interval search evaluates the statistic at the estimate. If the
  estimate is outside the confidence set, the search warns and restarts
  from the maximizer over the extended parameter set.
* With `nonneg = TRUE` the interval for a variance is intersected with
  [0, Inf). An empty intersection gives `NA` bounds and a warning.
* If no bound is found within the search distance, the warning reports the
  statistic where the search ended. The bound is infinite if the statistic
  levels off below the critical value, as for a variance whose grouping
  factor has few levels.
* The search continues from a feasible starting point when the nuisance
  parameters from the previous step are infeasible at the next value.
* `vc_ci()` and `vc_test()` reject misnamed optimizer arguments and invalid
  flags. The default `parm` excludes the parameters in `known`.

## Bug fixes

* The log-likelihood is `-Inf` wherever `Z Psi Z' + psi_r I` is not
  positive definite. It could be finite there before, which corrupted
  interval searches.

## Performance

* `vc_ci()` on the `fev1` example takes about 0.6 s instead of 43 s, with
  identical intervals.

# reconf 0.1

* Initial release.
* `ci_lmer()`: score-based confidence interval for a covariance parameter of
  an lme4 fit.
* `ci_all_lmer()`: intervals for all or some covariance parameters.
* `score_test_lmer()`: score test for covariance parameters.
