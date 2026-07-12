# q-side evaluation through loglikelihood(), used as the oracle for the
# n-side and spectral paths (test-nside.R, test-spectral.R). The q-side is
# itself checked against lme4 and numerical derivatives in test-cpp.R.
ll_q <- function(psi, Y, X, Z, Hlist, REML, b = NULL, expected = TRUE) {
  reconf:::loglikelihood(psi = psi, b = b, Y = Y, X = X, Z = Z, Hlist = Hlist,
                         REML = REML, get_val = TRUE, get_score = TRUE,
                         get_inf = TRUE, get_beta = TRUE, expected = expected,
                         check = FALSE)
}
