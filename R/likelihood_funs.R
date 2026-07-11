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
#' Cholesky factor of \code{ZtZ}, used by the feasibility check). Use
#' \code{get_precomp} to construct all entries at once.
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
#' @useDynLib reconf, .registration=TRUE
#' @import Matrix methods
#' @keywords internal
loglikelihood <-function(psi, b = NULL, Y, X = NULL, Z, Hlist, REML = TRUE, get_val = TRUE,
                         get_score = TRUE, get_inf = TRUE, get_beta = FALSE,
                         expected = TRUE, precomp = NULL, check = TRUE)
{
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
                                               Hlist = Hlist)
  H <- precomp$H
  if (is.null(H)) H <- methods::as(do.call(cbind, Hlist), "generalMatrix")

  # Dimension of score/information in the returned list
  k <- if (!REML && get_beta && p > 0) p + r else r

  # Feasibility gate and shared factorization. The parameters are feasible
  # iff psi[r] > 0 and Sigma = psi_r (I_n + Z Psi_r Z') is positive definite,
  # which holds iff I_q + F Psi_r F' is positive definite for any F with
  # F'F = Z'Z (the nonzero eigenvalues of F Psi_r F' and Psi_r Z'Z coincide).
  # The determinant sign alone is not sufficient: in balanced designs several
  # eigenvalues cross zero together and the sign may not flip.
  infeasible <- list("value" = -Inf, "score" = numeric(k),
                     "inf_mat" = matrix(0, k, k))
  if (psi[r] <= 0) return(infeasible)

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

#' Log-likelihood in pure R
#'
#' Calculates the log-likelihood, score vector and Fisher information matrix
#' for the covariance parameter vector \code{psi} in a linear mixed effects model.
#'
#' @param Z A matrix of fixed effects.
#' @param ZtZXe Pre-computed matrix with of \eqn{Z'[Z, X, e]}
#' @param e An error or residual vector (see \code{?loglikelihood}).
#' @param H Matrix of derivatives of Psi with respect to elements of psi.
#'          Assumes \eqn{H = [H_1, \dots , H_{r - 1}]}, where \eqn{H_j} is \eqn{q\times q}.
#' @param Psi_r Covariance matrix of random effects (\eqn{\Psi}) divided by error
#'             variance \eqn{\psi_r}.
#' @param psi_r The error variance, \eqn{\psi_r}.
#' @param get_val If \code{TRUE} (default), the log-likelihood will be calculated.
#' @param get_score If \code{TRUE} (default), the score vector will be calculated.
#' @param get_inf If \code{TRUE} (default), the information matrix will be calculated.
#' @param expected If \code{TRUE} (detault), return expected information;
#' otherwise observed.
#'
#' @return A list with components:
#' \item{value}{The log-likelihood evaluated at supplied parameters}
#' \item{score}{The score at supplied parameters}
#' \item{inf_mat}{The information matrix at supplied parameters}
loglik_psi <- function(Z, ZtZXe, e, H, Psi_r, psi_r, get_val = TRUE,
                       get_score = TRUE, get_inf = TRUE, expected = TRUE)
{
  # Define dimensions
  n <- length(e)
  q <- ncol(Psi_r)
  rm1 <- ncol(H) / q
  r <- rm1 + 1
  p <- ncol(ZtZXe) - q - 1

  # loglikelihood to return
  ll <- NA

  # Score vector to return
  s_psi <- rep(NA, r)

  # Fisher information to return
  I_psi <- matrix(NA, r, r)

  # Pre-compute Psi_rZtZ (columns 1:q), Pzi0ZtX (columns (q + 1):(q + p)),
  # and Psi_rZte (column q + p + 1)
  A <- Matrix::crossprod(Psi_r, ZtZXe)

  # Add loglik term before overwriting
  if(get_val){
    ll <- -0.5 * Matrix::determinant(A[, 1:q] + Matrix::Diagonal(q), logarithm = TRUE)$modulus -
      0.5 * n * log(psi_r)
  }

  # Matrix denoted M in manuscript is A[, 1:q]
  A <- Matrix::solve(A[, 1:q] + Matrix::Diagonal(q), A, sparse = TRUE)

  # Score for error variance psi_r
  # NB: REPLACE e by Sigma^{-1}e
  if(get_val || (get_inf && !expected)){
   e_save <- e
  }
  e <- (1 / psi_r) * (e - Z %*% A[, q + p + 1]) # = Sigma^{-1}e

  if(get_val){
    ll <- ll  - 0.5 * sum(e * e_save)
  }

  trace_M <- sum(Matrix::diag(A[, 1:q]))

  s_psi[r] <- 0.5 * sum(e^2) - (0.5 / psi_r) * (n - trace_M)


  # Use recycling to compute v'H_i v for all Hi
  v <- as.vector(Matrix::crossprod(Z, e)) # sparse matrix does not recycle
  w <- as.vector(Matrix::crossprod(v, H)) # Used later if !expected

  s_psi[-r] <- 0.5 * colSums(matrix(as.vector(w * v), nrow = q))

  # B = Z'Z (M - I_q) in paper notation
  B <- A[, 1:q]
  Matrix::diag(B) <- Matrix::diag(B) - 1
  B <- Matrix::crossprod(ZtZXe[, 1:q], B)

  if(!get_inf){
    # Compute -[ZtZ (I_q - M) * H_1, ..., ZtZ (I_q - M) * H_r] using
    # recycling, where * denotes elementwise multiplication
    # The "if" is because the calculation is a byproduct of a more expensive one
    # (B %*% H) done to get Fisher information
    H <- as.matrix(H) * as.vector(B)

    s_psi[-r] <- s_psi[-r] + (0.5 / psi_r) * colSums(matrix(Matrix::colSums(H), nrow = q))
  } else{
    H <- B %*% H
    I_psi[r, r] <- (0.5 / psi_r^2) * (n - 2 * trace_M +
    sum(Matrix::t(A[, 1:q]) * A[, 1:q]))

    # Subtract identity matrix from M
    Matrix::diag(A[, 1:q]) <- Matrix::diag(A[, 1:q]) - 1

    D <- matrix(Matrix::colSums(as.vector(A[, 1:q])  * as.matrix(H)), nrow = q)

    I_psi[-r, r] <- (0.5 / psi_r^2) * Matrix::colSums(D)

    for(ii in 1:rm1){
      first_idx <- ((ii - 1) * q + 1):(ii * q)
      s_psi[ii] <- s_psi[ii] + (0.5 / psi_r) * sum(Matrix::diag(H[, first_idx]))
      for(jj in ii:rm1){
        second_idx <- ((jj - 1) * q + 1):(jj * q)
        I_psi[ii, jj] <- (0.5 / psi_r^2) * sum(Matrix::t(H[, second_idx]) * H[, first_idx])
      }
    }
    if(!expected){
      I_psi <- -I_psi
      # u = Sigma^{-2}e. Some calculations could be saved from before
      u <- (1 / psi_r^2) * (e_save + Z %*% (-2 * A[, q + p + 1] +
                            (A[, 1:q] + Matrix::Diagonal(q, 1)) %*% A[, q + p + 1]))
      I_psi[r, r] <- I_psi[r, r] + sum(e * u)

      v <- as.vector(Matrix::crossprod(Z, u)) # = Z' Sigma^{-2}e
      I_psi[-r, r] <- I_psi[-r, r] + colSums(matrix(v * w, ncol = rm1))

      for(ii in 1:rm1){
        first_idx <- ((ii - 1) * q + 1):(ii * q)
        for(jj in ii:rm1){
          second_idx <- ((jj - 1) * q + 1):(jj * q)
          I_psi[ii, jj] <- I_psi[ii, jj] - (1 / psi_r) *
            sum(crossprod(w[first_idx], B) * w[second_idx])
        }
      }
    }
  }
  I_psi <- Matrix::forceSymmetric(I_psi, uplo = "U")
  return(list("value" = c(ll),  "score" = s_psi, "inf_mat" = I_psi))
}

chol_solve <- function(U, b)
{
  Matrix::solve(U, Matrix::solve(Matrix::t(U), b))
}


#' Restricted likelihood in pure R
#'
#' Computes the restricted (residual) likelihood, score, and information matrix
#' for the variance parameter \code{psi} a linear mixed effects model
#'
#' @param XtX An n x p matrix of the crossproduct of the design matrix of fixed
#' effects (\eqn{X}) with itself.
#' @param XtY An \eqn{p \times 1} vector of the crossproduct of the design matrix of fixed
#' effects (\eqn{X}) with the response vector (\eqn{Y}).
#' @param XtZ An \eqn{p \times q} matrix of the crossproduct of the design matrix of fixed
#' effects (\eqn{X}) with the design matrix of random effects (\eqn{Z}).
#' @param ZtZ A q x q matrix of the crossproduct of the design matrix of random
#' effects (\eqn{Z}) with itself.
#' @param ZtY A \eqn{q \times 1} vector of the crossproduct of the design matrix of random
#' effects (\eqn{Z}) with the response vector (\eqn{Y}).
#' @param Y An n x 1 vector of the response variable.
#' @param X An n x p matrix of the design matrix of fixed effects.
#' @param Z An n x q matrix of the design matrix of random effects.
#' @param H A q x rq matrix, where \eqn{H = [H_1, \dots , H_r]}, with \eqn{H_j = \partial \Psi / \partial \psi_j}
#' @param Psi_r The covariance matrix of the random effects (\eqn{Psi}) divided by the
#' error variance (\eqn{\psi_r})
#' @param psi_r A scalar value of the error variance.
#' @param get_val If \code{TRUE} (default), the log-likelihood will be computed.
#' @param get_score If \code{TRUE} (default), the score vector will be computed.
#' @param get_inf If \code{TRUE} (default), the Fisher information matrix will be
#' computed.
#'
#' @return A list with elements:
#' \item{value}{The restricted log-likelihood at supplied parameters.}
#' \item{score}{The restricted score at the supplied parameters.}
#' \item{inf_mat}{The expected restricted information matrix at supplied parameters}
#' \item{beta}{Partial maximizer of the likelihood in \eqn{\beta}}
#' \item{I_b_chol}{Cholesky root of \eqn{X'\Sigma^{-1}X}, the expected
#' information matrix for \eqn{\beta}}
res_ll <- function(XtX, XtY, XtZ, ZtZ, ZtY, Y, X, Z, H, Psi_r, psi_r,
                   get_val = TRUE, get_score = FALSE, get_inf = FALSE)
{
  # Define dimensions
  n <- length(Y)
  q <- ncol(Psi_r)
  rm1 <- ncol(H) / q # Assumes H = [H_1, ... , H_{r - 1}], where H_j is q by q
  r <- rm1 + 1
  p <- ncol(XtX)

  # Loglikelihood to return
  ll <- NA

  # Score vector to return
  s_psi <- rep(NA, r)

  # Fisher information to return
  I_psi <- matrix(NA, r, r)
  # Pre-compute A = (I_q + Psi_r Z'Z)^{-1} Psi_r
  A <- Matrix::crossprod(Psi_r, ZtZ) + Matrix::Diagonal(q) # q x q storage

  # Add likelihood term before overwriting
  if(get_val) ll <- Matrix::determinant(A, logarithm = TRUE)$modulus

  A <- Matrix::solve(A, Psi_r)
  B <- XtZ %*% A # q x q

  # Create XtSiX
  U <- Matrix::forceSymmetric((1 / psi_r) * (XtX - Matrix::tcrossprod(B, XtZ))) # p x p, now XtSiX
  U <- try(Matrix::chol(U)) # replace XtSiX by its Cholesky root
  if(inherits(U,"try-error")){
     return(list("value" = -Inf, "score" = s_psi, "inf_mat" = I_psi,
                 "beta" = rep(NA, p), "I_b_chol" = matrix(NA, p, p)))
  }

  # Create XtSiY for use in beta_tilde
  beta_tilde <- (1/ psi_r) * (XtY - XtZ %*% (A %*% ZtY)) # p x 1
  beta_tilde <- chol_solve(U, beta_tilde)

  # Replace Y by residuals
  Y <- Y - X %*% beta_tilde

  # n x 1 vector for storing \Sigma^{-1}e
  a <- (1 / psi_r) * (Y - Z %*% (A %*% Matrix::crossprod(Z, Y))) # n x 1

  if(get_val){
    ll <- ll + 2 * sum(log(Matrix::diag(U)))
    ll <- ll + sum(Y * a) + n * log(2 * pi * psi_r)
    ll <- -0.5 * ll
  }

  if(get_score){
    # Stochastic part of restricted score for psi
    s_psi[r] <- 0.5 * sum(a^2)
    v <- as.vector(Matrix::crossprod(Z, a)) # q x 1 vector storage
    s_psi[-r] <- 0.5 * colSums(matrix(as.vector(Matrix::crossprod(v, H)) * v,
                                      nrow = q))
  }
  #############################################################################
  ## NOTHING BELOW SHOULD DEPEND ON Y.
  #############################################################################

  if(get_inf){
    A <- Matrix::tcrossprod(A, ZtZ) # q x q, called M in manuscript

    s_psi[r] <- s_psi[r] - (0.5 / psi_r) * n + (0.5 / psi_r) * sum(Matrix::diag(A))

    I_psi[r, r] <- (0.5 / psi_r^2) * (n - 2 * sum(Matrix::diag(A)) +
                                       sum(Matrix::t(A) * A))

    E <- Matrix::crossprod(ZtZ, A) # q x q storage
    D <- XtZ %*% A # p x q storage
    XtSiZ <- (1 / psi_r) * (XtZ - D) # p x q
    XtSi2Z <- (1 / psi_r)^2 * (XtZ - 2 * D + D %*% A) # p x q


    C <- Matrix::tcrossprod(B, XtZ) # p x p storage, here XtZA ZtX
    G <- B %*% Matrix::tcrossprod(ZtZ, B) # p x q, here XtZA ZtZ AtZtX
    XtSi2X <- (1 / psi_r)^2 * (XtX - 2 * C + G) # p x p
    XtSi3X <- (1 / psi_r^3) * (XtX - 3 * C + 3 * G -
                                D %*% tcrossprod(A, B)) # p x p
    C <- chol_solve(U, XtSi2X)

    I_psi[r, r] <- I_psi[r, r] + 0.5 * sum(C * Matrix::t(C))
    s_psi[r] <- s_psi[r] + 0.5 * sum(Matrix::diag(C))
    I_psi[r, r] <- I_psi[r, r] - sum(Matrix::diag(chol_solve(U, XtSi3X)))

    # A (q x q), G (p x q) ARE FREE
    A <- (1 / psi_r)^2 * (ZtZ - 2 * E + E %*% A) # ZtSi2Z right now
    E <-  (1/ psi_r) * (ZtZ - E) # Now holds ZtSiZ
    D <- chol_solve(U, XtSiZ)
    A <- A - 2 * Matrix::crossprod(D, XtSi2Z) + Matrix::crossprod(XtSiZ, C %*% D)

    I_psi[-r, r] <- 0.5 * colSums(matrix(Matrix::colSums(as.vector(A) * H), nrow = q))
    s_psi[-r] <- s_psi[-r] - 0.5 * colSums(matrix(Matrix::colSums(
      as.vector(E - Matrix::crossprod(XtSiZ, D)) * H), nrow = q))

    H2 <- Matrix::crossprod(XtSiZ, D %*% H) # Storage can be avoided by
    # muliply in loop
    # Has to come after H2 since H is overwritten
    H <- Matrix::crossprod(E, H) # = ZtSiZ %*% H
    for(jj in 1:rm1){
      idx1 <- ((jj - 1) * q + 1):(jj * q)
      for(ii in 1:jj){
        idx2 <-  ((ii - 1) * q + 1):(ii * q)
        I_psi[ii, jj] <- 0.5 * sum(H[, idx1] * Matrix::t(H[, idx2])) -
          sum(H[, idx1] * Matrix::t(H2[, idx2])) + 0.5 *
          sum(H2[, idx1] * Matrix::t(H2[, idx2]))
      }
    }
  } else if (get_score){
    A <- Matrix::tcrossprod(A, ZtZ) # q x q, called M in manuscript
    s_psi[r] <- s_psi[r] - (0.5 / psi_r) * n + (0.5 / psi_r) * sum(Matrix::diag(A))

    D <- XtZ %*% A # p x q storage

    C <- Matrix::tcrossprod(B, XtZ) # p x p storage
    C <- chol_solve(U, (1 / psi_r)^2 * (XtX - 2 * C + B %*% Matrix::tcrossprod(ZtZ, B)))
    s_psi[r] <- s_psi[r] + 0.5 * sum(Matrix::diag(C))

    A <- Matrix::crossprod(ZtZ, A)
    A <- (1/ psi_r) * (ZtZ - A)

    XtSiZ <- (1 / psi_r) * (XtZ - D) # p x q
    D <- chol_solve(U, XtSiZ) # p x q
    v <- as.vector(A - Matrix::crossprod(XtSiZ, D)) # pq
    s_psi[-r] <- s_psi[-r] - 0.5 * colSums(matrix(Matrix::colSums(v * H), nrow = q))
  }
  I_psi <- Matrix::forceSymmetric(I_psi, uplo = "U")
  return(list("value" = ll[1], "score" = s_psi, "inf_mat" = I_psi, "beta" = beta_tilde,
              "I_b_chol" = U))
}
