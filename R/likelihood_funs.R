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
#' independent of \eqn{q}. The default \code{"auto"} picks \code{"n_side"}
#' iff \eqn{q \ge n} and \eqn{Z} is dense (more than 10 percent nonzeros), the
#' regime where the sparse path degenerates; for sparse \eqn{Z} the q-side is
#' fast even when \eqn{q \gg n}. When \code{precomp} is supplied, its
#' \code{method} tag takes precedence and this argument is ignored.
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
#' Specifically, \eqn{\Psi = \sum_{j = 1}^{r - 1}\psi_j H_j} where each \eqn{H_j}
#' is a \eqn{q\times q} matrix of zeros and ones. The argument \code{Hlist} is a list of length
#' \eqn{r - 1} whose \eqn{j}th element is \eqn{H_j}. This specification implies
#' each element of \eqn{\Psi} is one of \eqn{\psi_1, \dots, \psi_{r - 1}}.
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
#' \code{K} of concatenated \eqn{Z H_j Z'} (see \code{?loglik_n}). Use
#' \code{get_precomp} to construct all entries at once for either path.
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
                         method = c("auto", "q_side", "n_side"))
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

    if(!expected && REML){
      warning("Observed information not implemented for restricted likelihood;
              using expected")
    }

    if(get_beta && REML){
      warning("Score or information for beta not available for restricted likelihood")
    }
  }
  if (!expected && REML) expected <- TRUE

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
                                get_score = get_score, get_inf = get_inf)
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
                            get_inf = get_inf)
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
