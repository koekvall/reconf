#' Log-likelihood as a function of the covariance parameters
#'
#' Returns the log-likelihood, or the restricted log-likelihood, of a linear
#' mixed effects model as a function of its covariance parameters, with the
#' model and its precomputations fixed. The likelihood is maximized over the
#' fixed effects.
#'
#' @param Y Vector of length \eqn{n} of responses.
#' @param X Dense matrix of size \eqn{n\times p} of regressors, or
#'   \code{NULL} for no fixed effects.
#' @param Z Sparse design matrix for the random effects of size \eqn{n\times q}.
#' @param Hlist A list of matrices determining how \eqn{\psi} is mapped to
#'   \eqn{\Psi} (see details).
#' @param REML If \code{TRUE}, use the restricted likelihood; otherwise the
#'   likelihood.
#' @param method The computational path: \code{"q_side"} works with sparse
#'   \eqn{q \times q} matrices via the Woodbury identity; \code{"n_side"}
#'   works with dense \eqn{n \times n} matrices, at a cost independent of
#'   \eqn{q}; \code{"spectral"} requires \eqn{r = 2} and evaluates in
#'   \eqn{O(n)} time after a one-time eigendecomposition stored in the
#'   precomputations (see \code{?loglik_spectral}). The default
#'   \code{"auto"} is \code{"spectral"} or \code{"n_side"} if \eqn{q \ge n}
#'   and more than 10 percent of the entries of \eqn{Z} are nonzero, and
#'   \code{"q_side"} otherwise. A supplied \code{precomp} fixes the path.
#' @param precomp Optional list of precomputed quantities from
#'   \code{get_precomp}.
#'
#' @return A function of \code{psi}, the vector of \eqn{r} covariance
#'   parameters, with arguments \code{get_val}, \code{get_score},
#'   \code{get_inf} (whether to compute the value, the score, and the
#'   information) and \code{expected} (expected or observed information),
#'   all \code{TRUE} by default. It returns a list with components
#'  \item{value}{The value of the log-likelihood}
#'  \item{score}{The score, or gradient of the log-likelihood, for \eqn{\psi}}
#'  \item{inf_mat}{The information matrix for \eqn{\psi}}
#'   The function has attribute \code{"r"}, the number of covariance
#'   parameters. The caller validates \code{psi}.
#'
#' @details
#' The model is \deqn{Y = X\beta + Z U + E,} where \eqn{U \sim N_q(0, \Psi)}
#' and \eqn{E \sim N_n(0, \psi_r I_n)}. The last element \code{psi[r]} of
#' \eqn{\psi} is the error variance and the first \eqn{r - 1} elements are
#' variances and covariances of random effects.
#'
#' The covariance matrix of the random effects is
#' \eqn{\Psi = \sum_{j = 1}^{r - 1}\psi_j H_j}, where each \eqn{H_j} is a
#' \eqn{q\times q} matrix of zeros and ones; \code{Hlist} is the list of
#' \eqn{H_1, \dots, H_{r - 1}}. Each element of \eqn{\Psi} is then zero or
#' one of \eqn{\psi_1, \dots, \psi_{r - 1}}.
#'
#' For \code{REML = FALSE}, the likelihood is evaluated at the generalized
#' least squares estimate of \eqn{\beta} given \eqn{\psi}: this profile
#' log-likelihood has the score for \eqn{\psi} as its gradient, and its
#' information is the Schur complement of the \eqn{\beta} block in the
#' information for \eqn{(\beta, \psi)} (see \code{?loglik}). Without fixed
#' effects the restricted likelihood is the likelihood.
#'
#' The parameters are feasible if \code{psi[r] > 0} and
#' \eqn{\Sigma = Z \Psi Z' + \psi_r I_n} is positive definite; this allows
#' \eqn{\Psi} to be indefinite. At infeasible parameters the value is
#' \code{-Inf} and score and information are zero.
#'
#' The value includes all constants and matches \code{logLik()} of an
#' equivalent \code{lme4} fit (for fits with prior weights, up to the
#' constant \eqn{\sum_i \log(w_i)/2} from the weight transformation).
#'
#' @noRd
make_loglik <- function(Y, X, Z, Hlist, REML = TRUE,
                        method = c("auto", "q_side", "n_side", "spectral"),
                        precomp = NULL) {
  method <- match.arg(method)
  if (is.null(X)) X <- matrix(0, length(Y), 0)
  if (!is(Z, "generalMatrix")) Z <- as(Z, "generalMatrix")
  .check_model(Y = Y, X = X, Z = Z, Hlist = Hlist, precomp = precomp)
  .check_flags(list(REML = REML))
  if (ncol(X) == 0) REML <- FALSE
  if (is.null(precomp)) {
    precomp <- get_precomp(Y = Y, X = X, Z = Z, Hlist = Hlist, method = method)
  }
  r <- length(Hlist) + 1L
  # A supplied precomp determines the path; untagged lists are q-side ones
  path <- if (is.null(precomp$method)) "q_side" else precomp$method

  ll <- switch(path,
    # Dense n-by-n path: everything, including the Sigma feasibility gate,
    # happens inside the kernels; see ?loglik_n
    n_side = {
      kernel <- if (REML) loglik_res_n else loglik_n
      function(psi, get_val = TRUE, get_score = TRUE, get_inf = TRUE,
               expected = TRUE) {
        kernel(K = precomp$K, psi = as.numeric(psi), Y = Y, X = X,
               get_val = get_val, get_score = get_score, get_inf = get_inf,
               expected = expected)
      }
    },
    # Spectral r = 2 path: Sigma = psi_1 K + psi_2 I_n shares the
    # eigenvectors of K = Z H_1 Z' for every psi, so the one-time
    # decomposition stored in the precomp makes each evaluation O(n) up to
    # fixed-dimension factors; see ?loglik_spectral. Feasibility is decided
    # inside the kernels, as on the n side.
    spectral = {
      kernel <- if (REML) loglik_res_spectral else loglik_spectral
      function(psi, get_val = TRUE, get_score = TRUE, get_inf = TRUE,
               expected = TRUE) {
        kernel(d = precomp$d, Yt = precomp$Yt, Xt = precomp$Xt,
               psi = as.numeric(psi), get_val = get_val,
               get_score = get_score, get_inf = get_inf, expected = expected)
      }
    },
    q_side = {
      H <- precomp$H
      if (is.null(H)) H <- methods::as(do.call(cbind, Hlist), "generalMatrix")
      q <- ncol(Z)
      infeasible <- list("value" = -Inf, "score" = numeric(r),
                         "inf_mat" = matrix(0, r, r))
      function(psi, get_val = TRUE, get_score = TRUE, get_inf = TRUE,
               expected = TRUE) {
        # Feasibility gate. The parameters are feasible iff psi[r] > 0 and
        # Sigma = Z Psi Z' + psi_r I_n is positive definite, which holds iff
        # I_q + F Psi_r F' is positive definite for any F with F'F = Z'Z
        # (the nonzero eigenvalues of F Psi_r F' and Psi_r Z'Z coincide).
        # The determinant sign alone is not sufficient: in balanced designs
        # several eigenvalues cross zero together and the sign may not flip.
        if (psi[r] <= 0) return(infeasible)
        Psi_r <- (1 / psi[r]) * Psi_from_H_cpp(psi_mr = psi[-r], H = H)
        B <- Psi_r %*% precomp$ZtZ + Matrix::Diagonal(q)
        # LU of B is cached in B by Matrix, so determinant and solve
        # factorize once. solve() uses a sparsity-exploiting triangular
        # solve, so A is obtained in time proportional to its number of
        # nonzeros for block-structured models.
        d <- tryCatch(Matrix::determinant(B, logarithm = TRUE),
                      error = function(e) NULL)
        if (is.null(d) || !is.finite(d$modulus) || d$sign <= 0 ||
            !.sigma_pd(Psi_r, R = precomp$R, Z = Z, gate = precomp$gate)) {
          return(infeasible)
        }
        A <- Matrix::solve(B, Psi_r)
        ldetB <- as.numeric(d$modulus)
        if (REML) {
          loglik_res(A = A, ldetB = ldetB, psi_r = psi[r], H = H, Y = Y,
                     X = X, Z = Z, XtX = precomp$XtX, XtZ = precomp$XtZ,
                     ZtZ = precomp$ZtZ, XtY = precomp$XtY, ZtY = precomp$ZtY,
                     get_val = get_val, get_score = get_score,
                     get_inf = get_inf, expected = expected)
        } else {
          loglik(A = A, ldetB = ldetB, psi_r = psi[r], H = H, Y = Y, X = X,
                 Z = Z, XtX = precomp$XtX, XtZ = precomp$XtZ,
                 ZtZ = precomp$ZtZ, XtY = precomp$XtY, ZtY = precomp$ZtY,
                 get_val = get_val, get_score = get_score,
                 get_inf = get_inf, expected = expected)
        }
      }
    }
  )
  structure(ll, r = r)
}

#' Profile log-likelihood via the r = 2 spectral formulation
#'
#' Computes the log-likelihood maximized over the fixed effects, and its
#' score vector and information matrix for the covariance parameters, when
#' there is a single random-effect variance (\eqn{r = 2}), so that
#' \eqn{\Sigma = \psi_1 K + \psi_2 I_n} with \eqn{K = Z H_1 Z'}. Because
#' \eqn{\Sigma} shares the eigenvectors of \eqn{K} for every \eqn{\psi}, the
#' one-time eigendecomposition \eqn{K = U D U'} stored by
#' \code{get_precomp(..., method = "spectral")} diagonalizes every
#' evaluation: in the rotated coordinates \eqn{\Sigma} has eigenvalues
#' \eqn{w = \psi_1 d + \psi_2} and all quantities are elementwise operations
#' on \eqn{n}-vectors, costing \eqn{O(n)} up to factors in the fixed
#' dimension \eqn{p}. The q-side and n-side paths factorize a matrix at
#' every evaluation; see \code{?loglik} and \code{?loglik_n}.
#'
#' @param d Vector of length \eqn{n} of eigenvalues of \eqn{K = Z H_1 Z'}.
#' @param Yt Vector of length \eqn{n} of rotated responses \eqn{U'Y}.
#' @param Xt Matrix of size \eqn{n \times p} of rotated predictors \eqn{U'X}.
#' @param psi Vector of length 2 of covariance parameters; the last element
#'        is the error variance.
#' @param get_val If \code{TRUE}, the value of the log-likelihood is computed.
#' @param get_score If \code{TRUE} the score vector is calculated.
#' @param get_inf If \code{TRUE}, an information matrix is calculated.
#' @param expected If \code{TRUE}, the expected information is calculated;
#'        otherwise the observed, or negative Hessian of the profile
#'        log-likelihood.
#'
#' @return A list with components \code{value}, \code{score}, and
#' \code{inf_mat} as in \code{?loglik}.
#'
#' @details The rotation leaves the likelihood unchanged: residuals enter
#' only through \eqn{U'e = Yt - Xt \beta}, and \eqn{\beta}-blocks through
#' \eqn{X'\Sigma^{-1}X = Xt' diag(1/w) Xt}. The fixed effects are profiled
#' out as in \code{?loglik}.
#'
#' Feasibility is decided here as in \code{?loglik_n}: \eqn{\Sigma} is
#' positive definite iff all \eqn{w > 0}, and \eqn{\psi_2 > 0} is checked
#' separately because \eqn{K} alone can be positive definite when
#' \eqn{q \ge n}. At infeasible parameters, or when \eqn{X'\Sigma^{-1}X} is
#' not positive definite, \code{value} is \code{-Inf} (regardless of
#' \code{get_val}) and score and information are zero.
#'
#' @noRd
loglik_spectral <- function(d, Yt, Xt, psi, get_val = TRUE, get_score = TRUE,
                            get_inf = TRUE, expected = TRUE) {
  n <- length(Yt)
  p <- ncol(Xt)

  # Initialize returns; zeros are also the infeasible-parameter returns
  ll <- NA_real_
  S <- numeric(2)
  I <- matrix(0, 2, 2)
  infeasible <- list("value" = -Inf, "score" = S, "inf_mat" = I)

  # Eigenvalues of Sigma; the feasibility gate (see Details)
  w <- psi[1] * d + psi[2]
  if (psi[2] <= 0 || min(w) <= 0) return(infeasible)

  # Profile out beta by generalized least squares in rotated coordinates;
  # all products with Sigma^{-1} are elementwise divisions by w
  Xw <- Xt / w
  et <- Yt
  if (p > 0) {
    ch <- tryCatch(chol(crossprod(Xt, Xw)), error = function(e) NULL)
    if (is.null(ch)) return(infeasible)
    beta <- backsolve(ch, backsolve(ch, crossprod(Xw, Yt), transpose = TRUE))
    et <- as.vector(Yt - Xt %*% beta)
  }
  etilde <- et / w

  if (get_val) {
    ll <- -0.5 * (sum(log(w)) + n * log(2 * pi) + sum(et * etilde))
  }
  if (!get_score && !get_inf) {
    return(list("value" = ll, "score" = S, "inf_mat" = I))
  }

  # Score: with K_1 = K and K_2 = I_n rotating to diag(d) and I_n,
  # S(psi_j) = 0.5 (e'Si K_j Si e - tr(Si K_j))
  S[1] <- 0.5 * (sum(d * etilde^2) - sum(d / w))
  S[2] <- 0.5 * (sum(etilde^2) - sum(1 / w))

  if (get_inf) {
    # Expected information I(psi_j, psi_k) = 0.5 tr(Si K_j Si K_k)
    I <- 0.5 * matrix(c(sum(d^2 / w^2), sum(d / w^2),
                        sum(d / w^2), sum(1 / w^2)), 2, 2)
    if (!expected) {
      # Observed information: flip the sign of the deterministic block and
      # add the stochastic terms u_j'Sigma^{-1}u_k, where u_j = K_j Si e
      Ue <- cbind(d * etilde, etilde, deparse.level = 0)
      I <- -I + crossprod(Ue, Ue / w)
      # Profile out beta: subtract C'U^{-1}C for the cross terms
      # C = X'Sigma^{-1}u_j (zero in expectation), with U = ch'ch
      if (p > 0) I <- I - crossprod(backsolve(ch, crossprod(Xw, Ue),
                                              transpose = TRUE))
    }
  }
  list("value" = ll, "score" = S, "inf_mat" = I)
}

#' Restricted log-likelihood via the r = 2 spectral formulation
#'
#' Computes the restricted log-likelihood, score vector, and information
#' matrix for the covariance parameters using the one-time
#' eigendecomposition of \eqn{K = Z H_1 Z'}; the spectral counterpart of
#' \code{loglik_res}. See \code{?loglik_spectral} for the rotation and when
#' to prefer this path.
#'
#' @param d Vector of length \eqn{n} of eigenvalues of \eqn{K = Z H_1 Z'}.
#' @param Yt Vector of length \eqn{n} of rotated responses \eqn{U'Y}.
#' @param Xt Matrix of size \eqn{n \times p} of rotated predictors \eqn{U'X}.
#' @param psi Vector of length 2 of covariance parameters; the last element
#'        is the error variance.
#' @param get_val If \code{TRUE}, the value of the log-likelihood is computed.
#' @param get_score If \code{TRUE} the score vector is calculated.
#' @param get_inf If \code{TRUE}, an information matrix is calculated.
#' @param expected If \code{TRUE}, the expected information is calculated;
#'        otherwise the observed, or negative Hessian of the restricted
#'        log-likelihood.
#'
#' @return A list with components \code{value}, \code{score}, and
#' \code{inf_mat} as in \code{?loglik_res}.
#'
#' @details In the rotated coordinates the projection
#' \eqn{P = \Sigma^{-1} - \Sigma^{-1}X(X'\Sigma^{-1}X)^{-1}X'\Sigma^{-1}}
#' is \eqn{diag(1/w) - B B'} with \eqn{B = diag(1/w) Xt R^{-1}} and
#' \eqn{R'R = Xt' diag(1/w) Xt}, so every trace is a sum over the
#' \eqn{n}-vector \eqn{m = diag(B B')} and \eqn{p \times p} products.
#'
#' Feasibility is decided here as in \code{?loglik_spectral}: at infeasible
#' parameters, or when \eqn{X'\Sigma^{-1}X} is not positive definite,
#' \code{value} is \code{-Inf} (regardless of \code{get_val}) and score
#' and information are zero.
#'
#' @noRd
loglik_res_spectral <- function(d, Yt, Xt, psi, get_val = TRUE,
                                get_score = TRUE, get_inf = TRUE,
                                expected = TRUE) {
  n <- length(Yt)
  p <- ncol(Xt)

  # Initialize returns; zeros are also the infeasible-parameter returns
  ll <- NA_real_
  s_psi <- numeric(2)
  I_psi <- matrix(0, 2, 2)
  infeasible <- list("value" = -Inf, "score" = s_psi, "inf_mat" = I_psi)

  # Eigenvalues of Sigma; the feasibility gate (see ?loglik_spectral)
  w <- psi[1] * d + psi[2]
  if (psi[2] <= 0 || min(w) <= 0) return(infeasible)

  # GLS quantities in rotated coordinates: Xw = Sigma^{-1}X and
  # U = X'Sigma^{-1}X, whose Cholesky attempt decides positive definiteness
  Xw <- Xt / w
  ch <- tryCatch(chol(crossprod(Xt, Xw)), error = function(e) NULL)
  if (is.null(ch)) return(infeasible)
  beta <- as.vector(backsolve(ch, backsolve(ch, crossprod(Xw, Yt),
                                            transpose = TRUE)))
  et <- as.vector(Yt - Xt %*% beta)
  etilde <- et / w

  if (get_val) {
    ll <- -0.5 * (sum(log(w)) + 2 * sum(log(diag(ch))) + sum(et * etilde) +
                    (n - p) * log(2 * pi))
  }

  if (get_score || get_inf) {
    # P = diag(1/w) - B B' in rotated coordinates (see Details); K_1 = K and
    # K_2 = I_n rotate to diag(d) and I_n
    B <- t(backsolve(ch, t(Xw), transpose = TRUE))
    m <- rowSums(B^2)

    # Score for psi_j: 0.5 (etilde'K_j etilde - tr(P K_j))
    s_psi[1] <- 0.5 * (sum(d * etilde^2) - sum(d / w) + sum(d * m))
    s_psi[2] <- 0.5 * (sum(etilde^2) - sum(1 / w) + sum(m))

    if (get_inf) {
      # Expected information I(psi_j, psi_k) = 0.5 tr(P K_j P K_k)
      #  = 0.5 {sum(d_j d_k / w^2) - 2 sum(d_j d_k m / w) + tr(G_j G_k)}
      # with G_j = B'K_j B of size p x p
      G1 <- crossprod(B, d * B)
      G2 <- crossprod(B)
      I_psi[1, 1] <- sum(d^2 / w^2) - 2 * sum(d^2 * m / w) + sum(G1 * G1)
      I_psi[1, 2] <- sum(d / w^2) - 2 * sum(d * m / w) + sum(G1 * G2)
      I_psi[2, 1] <- I_psi[1, 2]
      I_psi[2, 2] <- sum(1 / w^2) - 2 * sum(m / w) + sum(G2 * G2)
      I_psi <- 0.5 * I_psi
      if (!expected) {
        # Observed information: flip the sign of the deterministic part and
        # add the stochastic terms u_j'P u_k, u_j = K_j etilde (u_2 = etilde)
        Ue <- cbind(d * etilde, etilde, deparse.level = 0)
        BtU <- crossprod(B, Ue)
        I_psi <- crossprod(Ue, Ue / w) - crossprod(BtU) - I_psi
      }
    }
  }
  list("value" = ll, "score" = s_psi, "inf_mat" = I_psi)
}
