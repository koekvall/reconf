#' Log-likelihood
#'
#' Computes log-likelihood, score vector, and information matrix
#' for a linear mixed effects model.
#'
#' @param psi Vector of length \eqn{r} of covariance parameters (see details).
#' @param b Vector of length\eqn{p} of fixed effects parameters
#' @param Y Vector of length \eqn{n} of responses.
#' @param X Dense matrix of size \eqn{n\times p} of regressors.
#' @param Z Sparse design matrix for the random effects of size \eqn{n\times q}.
#' @param Hlist A list of matrices determining how \eqn{\psi} is mapped to \eqn{\Psi} (see details)
#' @param REML If \code{TRUE}, use the restricted likelihood; otherwise the regular likelihood is used
#' @param get_val If \code{TRUE}, the value of the log-likelihood is computed.
#' @param get_score If \code{TRUE} the score vector, or gradient of log-likelihood, is calculated.
#' @param get_inf If \code{TRUE}, an information matrix is calculated.
#' @param get_beta If \code{TRUE} and \code{REML} is \code{FALSE}, return score and
#'  information for \eqn{\theta = [\beta', \psi']'}, otherwise for \eqn{\psi} only.
#' @param expected If \code{TRUE}, the expected information is calculated; otherwise the observed, or negative Hessian of the log-likelihood.
#' @param precomp Optional list of pre-computed quantities. Entries must have the
#' correct names and classes (see details).
#' @param check If \code{TRUE} (default), validate arguments. Internal callers
#' that evaluate the likelihood repeatedly with structurally identical arguments
#' set \code{FALSE} to skip redundant validation.
#' @param method Which computational path to use: \code{"q_side"} works with
#' sparse \eqn{q \times q} matrices via the Woodbury identity;
#' \code{"n_side"} works with dense \eqn{n \times n} matrices and has cost
#' independent of \eqn{q}; \code{"spectral"} applies only when \eqn{r = 2}
#' and evaluates in \eqn{O(n)} time after a one-time eigendecomposition
#' stored in the precomputations (see \code{?loglik_spectral}). The default
#' \code{"auto"} picks a dense path iff \eqn{q \ge n} and \eqn{Z} is dense
#' (more than 10 percent nonzeros), the regime where the sparse path
#' degenerates: \code{"spectral"} if \eqn{r = 2}, otherwise \code{"n_side"}.
#' For sparse \eqn{Z} the q-side is fast even when \eqn{q \gg n}.
#' When \code{precomp} is supplied, its \code{method} tag takes precedence
#' and this argument is ignored.
#'
#'
#' @return A list with components:
#'  \item{value}{The value of the log-likelihood}
#'  \item{score}{By default the score, or gradient of the log-likelihood, for \eqn{\psi}. If \code{get_beta = TRUE}
#'  and \code{REML = FALSE}, the score is for \eqn{\theta = [\beta', \psi']'}}
#'  \item{inf_mat}{By default, an information matrix for \eqn{\psi}. If \code{get_beta = TRUE}
#'  and \code{REML = FALSE}, an information matrix for \eqn{\theta = [\beta', \psi']'}}
#'
#' @details
#' The model is \deqn{Y = X\beta + Z U + E,} where \eqn{U \sim N_q(0, \Psi)}
#' and \eqn{E \sim N_n(0, \psi_r I_n)}. The last element of \eqn{\psi} (or `psi[r]`)
#' is the error variance. The first \eqn{r - 1} elements of \eqn{\psi}
#' are variances and covariances of random effects.
#'
#' \eqn{\Psi = \sum_{j = 1}^{r - 1}\psi_j H_j}, where each \eqn{H_j}
#' is a \eqn{q\times q} matrix of zeros and ones. The argument \code{Hlist} is a list of length
#' \eqn{r - 1} whose \eqn{j}th element is \eqn{H_j}. Each element of
#' \eqn{\Psi} is then one of \eqn{\psi_1, \dots, \psi_{r - 1}}.
#'
#' If \code{precomp} is supplied, it must be a list with elements
#'
#'  - \code{XtX = as.matrix(crossprod(X))}
#'  - \code{XtZ = as.matrix(crossprod(X, Z))}
#'  - \code{ZtZ = methods::as(crossprod(Z), "generalMatrix")}
#'
#' and, if \code{REML} is \code{TRUE}, also
#'
#'  - \code{XtY = as.vector(crossprod(X, Y))}
#'  - \code{ZtY = as.vector(crossprod(Z, Y))}
#'
#' Optional entries, computed internally when absent: \code{H} (the matrix
#' \code{cbind(Hlist)} coerced to \code{generalMatrix}) and \code{R} (sparse
#' Cholesky factor of \code{ZtZ}, used by the feasibility check). A precomp
#' list with \code{method = "n_side"} instead carries the dense matrix
#' \code{K} of concatenated \eqn{Z H_j Z'} (see \code{?loglik_n}), and one
#' with \code{method = "spectral"} carries the eigenvalues \code{d} of
#' \eqn{K = Z H_1 Z'} and the rotated data \code{Yt}, \code{Xt} (see
#' \code{?loglik_spectral}). Use \code{get_precomp} to construct all entries
#' at once for any path.
#'
#' If \code{REML} is \code{FALSE}, the parameter \code{b} must be provided
#' when \code{p > 0} to compute residuals \code{e = Y - X %*% b}. When
#' \code{p = 0}, residuals are computed as \code{e = Y}.
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
#' @useDynLib reconf, .registration=TRUE
#' @import Matrix methods
#' @keywords internal
loglikelihood <-function(psi, b = NULL, Y, X = NULL, Z, Hlist, REML = TRUE, get_val = TRUE,
                         get_score = TRUE, get_inf = TRUE, get_beta = FALSE,
                         expected = TRUE, precomp = NULL, check = TRUE,
                         method = c("auto", "q_side", "n_side", "spectral"))
{
  method <- match.arg(method)
  r <- length(psi)
  n <- length(Y)
  if(is.null(X) || ncol(X) == 0){
    p <- 0
    X <- matrix(0, n, 0)
    REML <- FALSE

    if(!is.null(b)){
      if (check) warning("X has zero columns or is NULL; setting b to NULL")
      b <- NULL
    }
  } else {
    p <- ncol(X)
  }
  if (!is(Z, "generalMatrix")) Z <- as(Z, "generalMatrix")
  q <- ncol(Z)

  if (check) {
    assertthat::assert_that(is.numeric(psi), r > 0,
                            msg = "psi should be a numeric vector of positive length")

    assertthat::assert_that(is.null(b) || is.vector(b, mode = "numeric"),
                            msg = "b should be NULL or a numeric vector")

    assertthat::assert_that(is.vector(Y, mode = "numeric"), n > 0,
                            msg = "Y should be a numeric vector of positive length")

    assertthat::assert_that(is.matrix(X), msg = "X should be NULL or a matrix")

    assertthat::assert_that(is(Z, "sparseMatrix"), q >= 1, nrow(Z) == n,
                            msg = "Z has to be an n x q matrix with q > 0")

    assertthat::assert_that(is.list(Hlist),
                            length(Hlist) == r - 1,
                            all(sapply(Hlist, methods::is, "sparseMatrix")),
                            all(sapply(Hlist, dim) == c(q, q)),
                            msg = "Hlist should be a list of length r - 1 with
                             q x q sparse matrices")

    assertthat::assert_that(is.logical(REML),
                            is.logical(get_val),
                            is.logical(get_score),
                            is.logical(get_inf),
                            is.logical(get_beta),
                            is.logical(expected),
                            msg = "REML, get_val, get_score, get_inf, get_beta,
                            and expected should all be logical")

    if(get_beta && REML){
      warning("Score or information for beta not available for restricted likelihood")
    }
  }

  if(p > 0 && is.null(b) && !REML){
    stop("b cannot be NULL when X has positive number of columns unless using REML")
  }

  if (is.null(precomp)) precomp <- get_precomp(Y = Y, X = X, Z = Z, REML = REML,
                                               Hlist = Hlist, method = method)

  # Dimension of score/information in the returned list
  k <- if (!REML && get_beta && p > 0) p + r else r

  # Feasibility gate. The parameters are feasible iff psi[r] > 0 and
  # Sigma = Z Psi Z' + psi_r I_n is positive definite. On the q side the
  # latter holds iff I_q + F Psi_r F' is positive definite for any F with
  # F'F = Z'Z (the nonzero eigenvalues of F Psi_r F' and Psi_r Z'Z coincide);
  # the determinant sign alone is not sufficient: in balanced designs several
  # eigenvalues cross zero together and the sign may not flip. On the n side
  # the kernels decide it by a dense Cholesky attempt on Sigma itself.
  infeasible <- list("value" = -Inf, "score" = numeric(k),
                     "inf_mat" = matrix(0, k, k))
  if (psi[r] <= 0) return(infeasible)

  # Dense n-by-n path: everything, including the Sigma feasibility gate,
  # happens inside the kernels; see ?loglik_n. A supplied precomp determines
  # the path (untagged lists are q-side ones built by earlier callers).
  if (identical(precomp$method, "n_side")) {
    if (REML) {
      ll_things <- loglik_res_n(K = precomp$K, psi = as.numeric(psi), Y = Y,
                                X = X, get_val = get_val,
                                get_score = get_score, get_inf = get_inf,
                                expected = expected)
    } else {
      e <- if (p == 0) Y else as.vector(Y - X %*% b)
      ll_things <- loglik_n(K = precomp$K, psi = as.numeric(psi), e = e,
                            X = X, get_val = get_val, get_score = get_score,
                            get_inf = get_inf, expected = expected)
      if (!get_beta && p > 0) {
        ll_things$score <- ll_things$score[-(1:p)]
        ll_things$inf_mat <- ll_things$inf_mat[-(1:p), -(1:p), drop = FALSE]
      }
    }
    return(list("value" = ll_things$value,
                "score" = ll_things$score,
                "inf_mat" = ll_things$inf_mat))
  }

  # Spectral r = 2 path: Sigma = psi_1 K + psi_2 I_n shares the eigenvectors
  # of K = Z H_1 Z' for every psi, so the one-time decomposition stored in
  # the precomp makes each evaluation O(n) up to fixed-dimension factors; see
  # ?loglik_spectral. Feasibility (beyond psi_r > 0 above) lives inside the
  # kernels, as on the n side.
  if (identical(precomp$method, "spectral")) {
    if (REML) {
      ll_things <- loglik_res_spectral(d = precomp$d, Yt = precomp$Yt,
                                       Xt = precomp$Xt, psi = as.numeric(psi),
                                       get_val = get_val,
                                       get_score = get_score,
                                       get_inf = get_inf,
                                       expected = expected)
    } else {
      ll_things <- loglik_spectral(d = precomp$d, Yt = precomp$Yt,
                                   Xt = precomp$Xt, psi = as.numeric(psi),
                                   b = b, get_val = get_val,
                                   get_score = get_score, get_inf = get_inf,
                                   expected = expected)
      if (!get_beta && p > 0) {
        ll_things$score <- ll_things$score[-(1:p)]
        ll_things$inf_mat <- ll_things$inf_mat[-(1:p), -(1:p), drop = FALSE]
      }
    }
    return(list("value" = ll_things$value,
                "score" = ll_things$score,
                "inf_mat" = ll_things$inf_mat))
  }

  H <- precomp$H
  if (is.null(H)) H <- methods::as(do.call(cbind, Hlist), "generalMatrix")

  Psi_r <- (1 / psi[r]) * Psi_from_H_cpp(psi_mr = psi[-r], H = H)
  B <- Psi_r %*% precomp$ZtZ + Matrix::Diagonal(q)
  # LU of B is cached in B by Matrix, so determinant and solve factorize once.
  # solve() uses a sparsity-exploiting triangular solve, so A is obtained in
  # time proportional to its number of nonzeros for block-structured models.
  d <- tryCatch(Matrix::determinant(B, logarithm = TRUE),
                error = function(e) NULL)
  if (is.null(d) || !is.finite(d$modulus) || d$sign <= 0 ||
      !.sigma_pd(Psi_r, R = precomp$R, Z = Z, gate = precomp$gate)) {
    return(infeasible)
  }
  A <- Matrix::solve(B, Psi_r)
  ldetB <- as.numeric(d$modulus)

  if(REML){
    ll_things <- loglik_res(A = A,
                            ldetB = ldetB,
                            psi_r = psi[r],
                            H = H,
                            Y = Y,
                            X = X,
                            Z = Z,
                            XtX = precomp$XtX,
                            XtZ = precomp$XtZ,
                            ZtZ = precomp$ZtZ,
                            XtY = precomp$XtY,
                            ZtY = precomp$ZtY,
                            get_val = get_val,
                            get_score = get_score,
                            get_inf = get_inf,
                            expected = expected)
  } else{
    # Always compute e from b to ensure consistency with theta
    e <- if(p == 0) Y else as.vector(Y - X %*% b)
    ll_things <- loglik(A = A,
                        ldetB = ldetB,
                        psi_r = psi[r],
                        H = H,
                        e = e,
                        X = X,
                        Z = Z,
                        XtX = precomp$XtX,
                        XtZ = precomp$XtZ,
                        ZtZ = precomp$ZtZ,
                        get_val = get_val,
                        get_score = get_score,
                        get_inf = get_inf,
                        expected = expected)
    if(!get_beta && p > 0){
      ll_things$score <- ll_things$score[-(1:p)]
      ll_things$inf_mat <- ll_things$inf_mat[-(1:p), -(1:p), drop = FALSE]
    }
  }
  list("value" = ll_things$value,
       "score" = ll_things$score,
       "inf_mat" = ll_things$inf_mat)
}

#' Log-likelihood via the r = 2 spectral formulation
#'
#' Computes the log-likelihood, score vector, and information matrix when
#' there is a single random-effect variance (\eqn{r = 2}), so that
#' \eqn{\Sigma = \psi_1 K + \psi_2 I_n} with \eqn{K = Z H_1 Z'}. Because
#' \eqn{\Sigma} shares the eigenvectors of \eqn{K} for every \eqn{\psi}, the
#' one-time eigendecomposition \eqn{K = U D U'} stored by
#' \code{get_precomp(..., method = "spectral")} diagonalizes every
#' evaluation: in the rotated coordinates \eqn{\Sigma} has eigenvalues
#' \eqn{w = \psi_1 d + \psi_2} and all quantities are elementwise operations
#' on \eqn{n}-vectors, costing \eqn{O(n)} up to factors in the fixed
#' dimension \eqn{p}. The q-side and n-side paths pay a factorization on
#' every evaluation; see \code{?loglik} and \code{?loglik_n}.
#'
#' @param d Vector of length \eqn{n} of eigenvalues of \eqn{K = Z H_1 Z'}.
#' @param Yt Vector of length \eqn{n} of rotated responses \eqn{U'Y}.
#' @param Xt Matrix of size \eqn{n \times p} of rotated predictors \eqn{U'X}.
#' @param psi Vector of length 2 of covariance parameters; the last element
#'        is the error variance.
#' @param b Vector of length \eqn{p} of fixed effects parameters, or
#'        \code{NULL} when \eqn{p = 0}.
#' @param get_val If \code{TRUE}, the value of the log-likelihood is computed.
#' @param get_score If \code{TRUE} the score vector is calculated.
#' @param get_inf If \code{TRUE}, an information matrix is calculated.
#' @param expected If \code{TRUE}, the expected information is calculated;
#'        otherwise the observed, or negative Hessian of the log-likelihood.
#'
#' @return A list with components \code{value}, \code{score}, and
#' \code{inf_mat} as in \code{?loglik}, with score and information of
#' dimension \eqn{p + 2}.
#'
#' @details The rotation leaves the likelihood unchanged: residuals enter
#' only through \eqn{U'e = Yt - Xt b}, and \eqn{\beta}-blocks through
#' \eqn{X'\Sigma^{-1}X = Xt' diag(1/w) Xt}.
#'
#' Feasibility is decided here as in \code{?loglik_n}: \eqn{\Sigma} is
#' positive definite iff all \eqn{w > 0}, and \eqn{\psi_2 > 0} is checked
#' separately because \eqn{K} alone can be positive definite when
#' \eqn{q \ge n}. At infeasible parameters \code{value} is \code{-Inf}
#' (regardless of \code{get_val}) and score and information are zero.
#'
#' @keywords internal
loglik_spectral <- function(d, Yt, Xt, psi, b = NULL, get_val = TRUE,
                            get_score = TRUE, get_inf = TRUE,
                            expected = TRUE) {
  n <- length(Yt)
  p <- ncol(Xt)

  # Initialize returns; zeros are also the infeasible-parameter returns
  ll <- NA_real_
  S <- numeric(p + 2)
  I <- matrix(0, p + 2, p + 2)

  # Eigenvalues of Sigma; the feasibility gate (see Details)
  w <- psi[1] * d + psi[2]
  if (psi[2] <= 0 || min(w) <= 0) {
    return(list("value" = -Inf, "score" = S, "inf_mat" = I))
  }

  # Rotated residuals; all products with Sigma^{-1} become elementwise
  # divisions by w
  et <- if (p > 0) as.vector(Yt - Xt %*% b) else Yt
  etilde <- et / w

  if (get_val) {
    ll <- -0.5 * (sum(log(w)) + n * log(2 * pi) + sum(et * etilde))
  }
  if (!get_score && !get_inf) {
    return(list("value" = ll, "score" = S, "inf_mat" = I))
  }

  # Score: S(beta) = X'Sigma^{-1}e and, with K_1 = K and K_2 = I_n rotating
  # to diag(d) and I_n, S(psi_j) = 0.5 (e'Si K_j Si e - tr(Si K_j))
  if (p > 0) S[1:p] <- as.vector(crossprod(Xt, etilde))
  S[p + 1] <- 0.5 * (sum(d * etilde^2) - sum(d / w))
  S[p + 2] <- 0.5 * (sum(etilde^2) - sum(1 / w))

  if (get_inf) {
    idx <- p + (1:2)
    Xw <- Xt / w
    # Information for beta: X'Sigma^{-1}X (equals the observed block)
    if (p > 0) I[1:p, 1:p] <- crossprod(Xt, Xw)
    # Expected information I(psi_j, psi_k) = 0.5 tr(Si K_j Si K_k)
    I[idx, idx] <- 0.5 * matrix(c(sum(d^2 / w^2), sum(d / w^2),
                                  sum(d / w^2), sum(1 / w^2)), 2, 2)
    if (!expected) {
      # Observed information: flip the sign of the deterministic psi block
      # and add the stochastic terms u_j'Sigma^{-1}u_k, where u_j = K_j Si e
      Ue <- cbind(d * etilde, etilde)
      I[idx, idx] <- -I[idx, idx] + crossprod(Ue, Ue / w)
      # Cross terms with beta: X'Sigma^{-1}u_j (zero in expectation)
      if (p > 0) {
        I[1:p, idx] <- crossprod(Xw, Ue)
        I[idx, 1:p] <- t(I[1:p, idx])
      }
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
#' @return A list with components \code{value}, \code{score}, \code{inf_mat},
#' \code{beta}, and \code{I_b_chol} as in \code{?loglik_res}.
#'
#' @details The projection
#' \eqn{P = \Sigma^{-1} - \Sigma^{-1}X(X'\Sigma^{-1}X)^{-1}X'\Sigma^{-1}}
#' never materializes: in the rotated coordinates
#' \eqn{P = diag(1/w) - B B'} with \eqn{B = diag(1/w) Xt R^{-1}} and
#' \eqn{R'R = Xt' diag(1/w) Xt}, so every trace reduces to sums over the
#' \eqn{n}-vector \eqn{m = diag(B B')} and \eqn{p \times p} products.
#'
#' Feasibility is decided here as in \code{?loglik_spectral}: at infeasible
#' parameters, or when \eqn{X'\Sigma^{-1}X} is not positive definite,
#' \code{value} is \code{-Inf} (regardless of \code{get_val}), score and
#' information are zero, and \code{beta} and \code{I_b_chol} are \code{NA}.
#'
#' @keywords internal
loglik_res_spectral <- function(d, Yt, Xt, psi, get_val = TRUE,
                                get_score = TRUE, get_inf = TRUE,
                                expected = TRUE) {
  n <- length(Yt)
  p <- ncol(Xt)

  # Initialize returns; zeros/NAs are also the infeasible-parameter returns
  ll <- NA_real_
  s_psi <- numeric(2)
  I_psi <- matrix(0, 2, 2)
  infeasible <- list("value" = -Inf, "score" = s_psi, "inf_mat" = I_psi,
                     "beta" = rep(NA_real_, p),
                     "I_b_chol" = matrix(NA_real_, p, p))

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
  list("value" = ll, "score" = s_psi, "inf_mat" = I_psi, "beta" = beta,
       "I_b_chol" = ch)
}
