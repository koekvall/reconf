# reconf 0.1.1

* The restricted-information computation no longer forms dense q x q
  matrices: all traces use the decomposition Z'PZ = Z'Sigma^{-1}Z - G'E1,
  whose second term has rank p, so the cost is linear in q for
  block-structured models. A full REML evaluation on the FEV1 example drops
  from 7.6 ms to 1.1 ms and `ci_all_lmer()` from 43 s to 6 s (26 s to 4 s
  with `onestep = TRUE`) relative to reconf 0.1, with identical results.
* The feasibility check caches the symbolic Cholesky analysis (pattern from
  `Psi(1)`, the sum of the structure matrices) and redoes only the numeric
  factorization at each evaluation.

* The log-likelihood is now `-Inf` whenever the marginal covariance matrix
  `Sigma = Z Psi Z' + psi_r I` is not positive definite. Previously the
  log-determinant was computed from its absolute value, so infeasible
  parameters could yield a finite value; in balanced designs this could not
  be detected from the determinant's sign because several eigenvalues cross
  zero together. Feasibility is now decided by an exact Cholesky-based test.
* The solve producing `A = (I + Psi_r Z'Z)^{-1} Psi_r` moved from Eigen's
  SparseLU, which treats each right-hand-side column as dense, to
  `Matrix::solve`, whose sparsity-aware triangular solve costs time
  proportional to the nonzeros in the solution. For block-structured random
  effects this makes the solve linear rather than quadratic in the number of
  random effects; `ci_all_lmer()` is roughly 1.7x faster on the FEV1 example
  and large grouped models gain much more.
* Likelihood evaluations in the confidence-interval search now reuse
  precomputed quantities (including the concatenated structure matrix `H`)
  and skip redundant argument validation; `ci_all_lmer()` extracts model
  components once instead of once per parameter.
* `loglik()` and `loglik_res()` (internal) now take `A` and the
  log-determinant as arguments instead of factorizing internally.
* Fixed inconsistent element names in the error path of `res_ll()` and
  renamed the returned Cholesky factor to `I_b_chol` to match what is
  computed (the factor of the information, not its inverse).
* Documentation fixes: the default `ci_lmer()` step size is `SE / 40`, and
  the unused residual entry was removed from the precompute list.

# reconf 0.1

* Initial release.
* `ci_lmer()`: score-based confidence interval for a single covariance parameter
  in a linear mixed model fitted with lme4.
* `ci_all_lmer()`: score-based confidence intervals for all (or a subset of)
  covariance parameters.
* `score_test_lmer()`: score test for covariance parameters against a user-supplied
  null hypothesis.
