# q-side evaluation through make_loglik(), used as the oracle for the
# n-side and spectral paths (test-nside.R, test-spectral.R). The q-side is
# itself checked against lme4 and numerical derivatives in test-cpp.R.
ll_q <- function(psi, Y, X, Z, Hlist, REML, expected = TRUE) {
  reconf:::make_loglik(Y = Y, X = X, Z = Z, Hlist = Hlist, REML = REML,
                       method = "q_side")(psi, expected = expected)
}
