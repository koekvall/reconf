#include <RcppEigen.h>

//' Psi_from_H_cpp
//'
//' Constructs the covariance matrix of the random effects
//'
//' @param psi_mr A vector of covariance parameter (see ?loglikelihood)
//' @param H Sparse matrix of derivatives of Psi with respect to elements of psi,
//'        \eqn{H = [H_1, \dots , H_{r - 1}]}, where \eqn{H_j = \partial \Psi / \partial \psi_j}.
//' @return The covariance matrix \eqn{\Psi}
// [[Rcpp::export]]
Eigen::SparseMatrix<double> Psi_from_H_cpp(const Eigen::Map<Eigen::VectorXd> psi_mr,
                                           const Eigen::MappedSparseMatrix<double> H) { //
  int q = H.rows();
  int rm1 = H.cols() / q;
  Eigen::SparseMatrix<double> Psi(q, q);
  for (int ii = 0; ii < rm1; ii++) {
     if(psi_mr(ii) != 0.0){ // avoid initializing elements that are zero anyway
       Psi += psi_mr(ii) * H.middleCols(ii * q, q);
     }
  }
  return Psi;
}

//' Log-likelihood using RcppEigen
//'
//' Computes the log-likelihood, score vector, and information matrix
//' for the covariance parameter vector in a linear mixed effects model.
//'
//' @param A The \eqn{q \times q} sparse matrix
//'        \eqn{A = (I_q + \Psi_r Z'Z)^{-1} \Psi_r}, where \eqn{\Psi_r = \Psi / \psi_r},
//'        precomputed by the caller (see details).
//' @param ldetB Log-determinant of \eqn{B = I_q + \Psi_r Z'Z}, precomputed by
//'        the caller.
//' @param psi_r The error variance \eqn{\psi_r > 0}.
//' @param H Sparse \eqn{q \times (qr - q)} matrix of horizontally concatenated
//'        derivatives of \eqn{\Psi} (see details) of class \code{dgCMatrix}.
//' @param e Vector of length \eqn{n} of errors, or residuals, \eqn{e = Y - X \beta}.
//' @param X Matrix of size \eqn{n \times p} of predictors, of class \code{matrix}.
//' @param Z Sparse \eqn{n \times q} random effect design matrix of class \code{dgCMatrix}.
//' @param XtX Precomputed matrix \code{crossprod(X)} of class \code{matrix}.
//' @param XtZ Precomputed matrix \code{crossprod(X, Z)} of class \code{matrix}.
//' @param ZtZ Precomputed matrix \code{crossprod(Z)} of class \code{dgCMatrix}.
//' @param get_val If \code{TRUE}, the value of the loglikelihood is computed.
//' @param get_score If \code{TRUE} the score vector is calculated.
//' @param get_inf If \code{TRUE}, an information matrix is calculated.
//' @param expected If \code{TRUE}, the expected information is calculated; otherwise
//' the observed, or negative Hessian of the loglikelihood.
//'
//' @return A list with components:
//' \item{value}{The value of the log-likelihood}
//' \item{score}{The score, or gradient of the log-likelihood, for \eqn{\psi}}
//' \item{inf_mat}{The information matrix for \eqn{\psi}}
//'
//' @details The model is \deqn{Y = X\beta + Z U + E,} where \eqn{U \sim N_q(0, \Psi)}
//' and \eqn{E \sim N_n(0, \psi_r I_n)}. The first \eqn{r - 1} elements of \eqn{\psi}
//' parameterize \eqn{\Psi}, while the \eqn{r}th and last element is the error
//' variance. It is assumed that \eqn{H_j = \partial \Psi / \partial \psi_j} is
//' a (usually sparse) matrix of zeros and ones, \eqn{j \in \{1, \dots, r - 1\}},
//' and that \eqn{\Psi = \sum_{j = 1}^{r - 1}\psi_j H_j}. Thus, \eqn{\psi_1, \dots, \psi_{r - 1}}
//' are variances and covariances of random effects.
//' The argument matrix \code{H} is \eqn{H = [H_1, \dots, H_{r - 1}]}.
//'
//' The fixed effects \eqn{\beta} affect the likelihood only through the
//' precomputed \eqn{e = Y - X\beta}.
//'
//' The information matrix includes both \eqn{\beta} and \eqn{\psi} parameters,
//' with dimensions \eqn{(p + r) \times (p + r)}.
//'
//' The caller is responsible for verifying that the parameters are feasible,
//' that is, that \eqn{\psi_r > 0} and \eqn{\Sigma = Z \Psi Z' + \psi_r I_n} is
//' positive definite, and for computing \code{A} and \code{ldetB}. Solving for
//' \code{A} with a solver that exploits sparsity in the right-hand side (for
//' example \code{Matrix::solve}) is much faster than a dense solve when the
//' random effects are block-structured.
//'
//' @useDynLib reconf, .registration=TRUE
//' @import Matrix
// [[Rcpp::export]]

Rcpp::List loglik(const Eigen::MappedSparseMatrix<double> A,
                  const double ldetB,
                  const double psi_r,
                  Eigen::SparseMatrix<double> H,
                  const Eigen::Map<Eigen::VectorXd> e,
                  const Eigen::Map<Eigen::MatrixXd> X,
                  const Eigen::MappedSparseMatrix<double> Z,
                  const Eigen::Map<Eigen::MatrixXd> XtX,
                  const Eigen::Map<Eigen::MatrixXd> XtZ,
                  const Eigen::MappedSparseMatrix<double> ZtZ,
                  const bool get_val = true,
                  const bool get_score = true,
                  const bool get_inf = true,
                  const bool expected = true) {
  // Define dimensions
  int p = X.cols();
  int n = e.size();
  int q = A.cols();
  int r = H.cols() / q + 1;

  // Initialize returns (allocate only when needed)
  double ll = NA_REAL;
  Eigen::VectorXd S = get_score || get_inf ? Eigen::VectorXd::Zero(p + r) : Eigen::VectorXd();
  Eigen::MatrixXd I = get_inf ? Eigen::MatrixXd::Zero(p + r, p + r) : Eigen::MatrixXd();

  // Initialize identity matrix
  Eigen::SparseMatrix<double> Id_q(q, q);
  Id_q.setIdentity();

  // Compute \tilde{e} = \Sigma^{-1}e
  Eigen::VectorXd etilde = (1.0 / psi_r) * (e - Z * (A * (Z.transpose() * e)));

  if (get_val) {
    ll = -0.5 * ldetB - 0.5 * (double)n * log(psi_r) - 0.5 * etilde.dot(e);
  }

  if (!get_score && !get_inf) {
    return Rcpp::List::create(Rcpp::Named("value") = ll,
                              Rcpp::Named("score") = S,
                              Rcpp::Named("inf_mat") = I);
  }

  // Get score for beta
  if(p > 0) {
    S.head(p) = X.transpose() * etilde;
  }

  // Get score for psi
  Eigen::SparseMatrix<double> C = Id_q - A * ZtZ;

  S(p + r - 1) = 0.5 * (etilde.dot(etilde) - (1 / psi_r) * (n - q +
    C.diagonal().sum()));

  // v is Z'Sigma^{-1}e, w is e'Sigma^{-1}Z[H1...Hr]
  Eigen::VectorXd v = (Z.transpose() * etilde);
  Eigen::VectorXd w = H.transpose() * v;
  for (int ii = 0; ii < r - 1; ii++) {
    S(p + ii) = 0.5 * w.middleRows(ii * q, q).dot(v);
  }

  // B holds Z' \Sigma^{-1} Z
  Eigen::SparseMatrix<double> B = (1 / psi_r) * ZtZ * C;

  if (!get_inf) {
    for (int ii = 0; ii < r - 1; ii++) {
      S(p + ii) -= 0.5 * H.middleCols(ii * q, q).cwiseProduct(B).sum();
    }
  } else {
    if(p > 0) {
      I.topLeftCorner(p, p) = (1 / psi_r) * (XtX - XtZ * (A * XtZ.transpose()));
    }
    Eigen::SparseMatrix<double> H1 = C.transpose(); // This is replaced later.
    // Putting instead C.transpose() in next call does not compile
    I(p + r - 1, p + r - 1) = (0.5 / (psi_r * psi_r)) *
      (n - q + C.cwiseProduct(H1).sum());
    // Replace H by Z'\Sigma^{-1}Z H
    H = B * H;
    Eigen::SparseMatrix<double> H2(q, q);
    for (int jj = 0; jj < r - 1; jj++) {
      H1 = H.middleCols(jj * q, q);
      I(p + jj, p + r - 1) = (0.5 / psi_r) * H1.cwiseProduct(C).sum();
      S(p + jj) -= 0.5 * H1.diagonal().sum();
      for (int ii = 0; ii <= jj; ii++) {
        H2 = H.middleCols(ii * q, q).transpose();
        I(p + ii, p + jj) = 0.5 * H1.cwiseProduct(H2).sum();
      }
    }
    if (!expected) {
      I.bottomRightCorner(r, r) = -I.bottomRightCorner(r, r);
      // \check{e} = Sigma^{-2}e = \Sigma^{-1}\tilde{e}
      Eigen::VectorXd e_check = (1 / psi_r) * (etilde - Z * (A * v));
      v = Z.transpose() * e_check;

      if(p > 0) {
        I.topRightCorner(p, 1) = X.transpose() * e_check;
      }

      I(p + r - 1, p + r - 1) += e_check.dot(etilde);

      for (int jj = 0; jj < r - 1; jj++) {
        if(p > 0) {
          I.block(0, p + jj, p, 1) = (1 / psi_r) * (XtZ *(C * w.middleRows(jj * q, q)));
        }
        I(p + jj, p + r - 1) += w.middleRows(jj * q, q).dot(v);
        for (int ii = 0; ii <= jj; ii++) {
          I(p + ii, p + jj) +=  (B * w.middleRows(ii * q, q)).dot(w.middleRows(jj * q, q));
        }
      }
    }
  }
  I = I.selfadjointView<Eigen::Upper>();
  return Rcpp::List::create(Rcpp::Named("value") = ll,
                            Rcpp::Named("score") = S,
                            Rcpp::Named("inf_mat") = I);
}

//' Restricted log-likelihood using RcppEigen
//'
//' Computes the restricted log-likelihood, score vector, and information matrix
//' for the covariance parameter vector in a linear mixed effects model.
//'
//' @param A The \eqn{q \times q} sparse matrix
//'        \eqn{A = (I_q + \Psi_r Z'Z)^{-1} \Psi_r}, where \eqn{\Psi_r = \Psi / \psi_r},
//'        precomputed by the caller (see \code{?loglik}).
//' @param ldetB Log-determinant of \eqn{B = I_q + \Psi_r Z'Z}, precomputed by
//'        the caller.
//' @param psi_r The error variance \eqn{\psi_r > 0}.
//' @param H Sparse \eqn{q \times q(r - 1)} matrix of horizontally concatenated
//'        derivatives of \eqn{\Psi} (see details) of class \code{dgCMatrix}.
//' @param Y Vector of length \eqn{n} of responses, of class \code{numeric}.
//' @param X Matrix of size \eqn{n \times p} of predictors, of class \code{matrix}.
//' @param Z Sparse \eqn{n \times q} random effect design matrix of class \code{dgCMatrix}.
//' @param XtX Precomputed matrix \code{crossprod(X)} of class \code{matrix}.
//' @param XtZ Precomputed matrix \code{crossprod(X, Z)} of class \code{matrix}.
//' @param ZtZ Precomputed matrix \code{crossprod(Z)} of class \code{dgCMatrix}.
//' @param XtY Precomputed vector \code{crossprod(X, Y)} of class \code{numeric}.
//' @param ZtY Precomputed vector \code{crossprod(Z, Y)} of class \code{numeric}.
//' @param get_val If \code{TRUE}, the value of the loglikelihood is computed.
//' @param get_score If \code{TRUE} the score vector is calculated.
//' @param get_inf If \code{TRUE}, an information matrix is calculated.
//'
//' @return A list with components:
//' \item{value}{The value of the restricted log-likelihood}
//' \item{score}{The restricted score, or gradient of the restricted log-likelihood, for \eqn{\psi}}
//' \item{inf_mat}{The restricted information matrix for \eqn{\psi}}
//' \item{beta}{Partial maximizer of the regular likelihood in \eqn{\beta},
//'   \eqn{\tilde{\beta} = (X' \Sigma^{-1} X)^{-1} X' \Sigma^{-1}Y},
//'   where \eqn{\Sigma = Z\Psi Z' + \psi_r I_n}}
//' \item{I_b_chol}{Cholesky root of the expected information matrix
//'   for \eqn{\beta}, \eqn{I(\beta; \psi) = X' \Sigma^{-1} X}}
//'
//' @details The model is \deqn{Y = X\beta + Z U + E,} where \eqn{U \sim N_q(0, \Psi)}
//' and \eqn{E \sim N_n(0, \psi_r I_n)}. The first \eqn{r - 1} elements of \eqn{\psi}
//' parameterize \eqn{\Psi}, while the \eqn{r}th and last element is the error
//' variance. It is assumed that \eqn{H_j = \partial \Psi / \partial \psi_j} is
//' a (usually sparse) matrix of zeros and ones, \eqn{j \in \{1, \dots, r - 1\}},
//' and that \eqn{\Psi = \sum_{j = 1}^{r - 1}\psi_j H_j}. Thus, \eqn{\psi_1, \dots, \psi_{r - 1}}
//' are variances and covariances of random effects.
//' The argument matrix \code{H} is \eqn{H = [H_1, \dots, H_{r - 1}]}.
//'
//' The restricted likelihood integrates out the fixed effects \eqn{\beta}.
//'
//' The caller is responsible for verifying feasibility and computing \code{A}
//' and \code{ldetB}; see \code{?loglik}.
//'
//' \code{A} and \code{H} are received by value (deep copy) rather than as
//' maps: Eigen's products with blocks of mapped sparse matrices fall back to
//' slow generic paths in the information computations, and the copies are
//' cheap relative to the algebra.
//'
//' @useDynLib reconf, .registration=TRUE
//' @import Matrix
// [[Rcpp::export]]
Rcpp::List loglik_res(const Eigen::SparseMatrix<double> A,
                      const double ldetB,
                      const double psi_r,
                      const Eigen::SparseMatrix<double> H,
                      Eigen::VectorXd Y,
                      const Eigen::Map<Eigen::MatrixXd> X,
                      const Eigen::MappedSparseMatrix<double> Z,
                      const Eigen::Map<Eigen::MatrixXd> XtX,
                      const Eigen::Map<Eigen::MatrixXd> XtZ,
                      const Eigen::MappedSparseMatrix<double> ZtZ,
                      const Eigen::Map<Eigen::VectorXd> XtY,
                      const Eigen::Map<Eigen::VectorXd> ZtY,
                      const bool get_val = true,
                      const bool get_score = true,
                      const bool get_inf = true)
{
  // Define dimensions
  int n = Y.size();
  int q = A.cols();
  int r = H.cols() / q + 1;
  int p = X.cols();
  // loglikelihood to return
  double ll = NA_REAL;
  // Score vector to return
  Eigen::VectorXd s_psi = Eigen::VectorXd::Zero(r);
  // Information matrix to return
  Eigen::MatrixXd I_psi =  Eigen::MatrixXd::Zero(r, r);

  Eigen::SparseMatrix<double> Id_q(q, q);
  Id_q.setIdentity();

  //Create XtSiX
  Eigen::MatrixXd U = (1.0 / psi_r) * (XtX - XtZ * A * XtZ.transpose());  //p*p

  // Force symmetric
  U = U.selfadjointView<Eigen::Upper>();
  // llt decomposition
  Eigen::LLT<Eigen::MatrixXd, Eigen::Upper> llt(U);

  if (llt.info() != Eigen::Success) {
    return Rcpp::List::create(
      Rcpp::Named("value") = -R_PosInf,
      Rcpp::Named("score") = s_psi,
      Rcpp::Named("inf_mat") = I_psi,
      Rcpp::Named("beta") = Eigen::VectorXd::Constant(p, NA_REAL),
      Rcpp::Named("I_b_chol") = Eigen::MatrixXd::Constant(p, p, NA_REAL));
  }

  // Create XtSiY and \tilde{\beta}
  // (1/ psi_r) * (XtY - XtZ %*% A %*% ZtY)
  Eigen::VectorXd beta_tilde = (1.0 / psi_r) * (XtY -
    XtZ * (A * ZtY));
  beta_tilde = llt.solve(beta_tilde);

  // Replace response with residuals
  Y = Y - X * beta_tilde;

  // n x 1 vector for storing \tilde{e} = \Sigma^{-1}e
  Eigen::VectorXd etilde = (1.0 / psi_r) * (Y - Z * (A * (Z.transpose() * Y)));


  if (get_val) {
    ll = ldetB;
    ll += 2.0 * llt.matrixLLT().diagonal().array().log().sum();
    ll += Y.dot(etilde) + n * log(2.0 * M_PI * psi_r);
    ll *= -0.5;
  }

  if (get_score) {
    // Stochastic part of the restricted score for psi
    s_psi(r - 1) = 0.5 * etilde.dot(etilde);
    Eigen::VectorXd v = Z.transpose() * etilde;
    for (int ii = 0; ii < r - 1; ii++) {
      s_psi(ii) = 0.5 * v.dot(H.middleCols(ii * q, q) * v);
    }
  }
  /////////////////////////////////////////////////////////////////////////////
  // NOTHING BELOW SHOULD DEPEND ON Y / NONSTOCHASTIC PARTS
  /////////////////////////////////////////////////////////////////////////////
  if (get_inf) {
    // Create matrices used repeatedly
    Eigen::SparseMatrix<double> Id_p(p, p);
    Id_p.setIdentity();
    Eigen::SparseMatrix<double> C = Id_q - A * ZtZ;
    Eigen::MatrixXd XtSiZ = (1.0 / psi_r) * XtZ * C;
    Eigen::MatrixXd E1 = llt.solve(XtSiZ);
    Eigen::MatrixXd D2 = (1.0 / psi_r) * (Id_p - E1 * (A * XtZ.transpose()));
    Eigen::MatrixXd E2 = (1.0 / psi_r) * E1 * C;

    // Score for psi_r is done after this
    s_psi(r - 1) -= (0.5 / psi_r) * (n - q + C.diagonal().sum());
    s_psi(r - 1) += 0.5 * D2.diagonal().sum();

    ////////////////////////////////////////////////////////////////////////////
    // Information for psi_r
    ////////////////////////////////////////////////////////////////////////////
    // The term -tr(D_{(3)})
    I_psi(r - 1, r - 1) = (-1.0 / psi_r) * (D2.diagonal().sum() -
      E2.cwiseProduct(XtZ * A).sum());

    Eigen::SparseMatrix<double> Ct = C.transpose();

    // The term 0.5 tr(\Sigma^{-2})
    I_psi(r - 1, r - 1) += (0.5 / (psi_r * psi_r)) *
                            (n - q + C.cwiseProduct(Ct).sum());

    // The term 0.5 tr(D_{(2)}^2)
    I_psi(r - 1, r - 1) += 0.5 * D2.transpose().cwiseProduct(D2).sum();

    ////////////////////////////////////////////////////////////////////////////
    // Information and score for  psi_{-r}, precompute before loop
    ////////////////////////////////////////////////////////////////////////////

    // Terms for S(psi_j)
    Eigen::SparseMatrix<double> ZtSiZ = (1.0 / psi_r) * ZtZ * C;
    Eigen::MatrixXd ZtSiXE1 = XtSiZ.transpose() * E1;

    // ZtSi2Z
    Eigen::SparseMatrix<double> ZtSi2Z = (1.0 / psi_r) * ZtSiZ * C;

     //Terms for I(psi_j, psi_r)
    Eigen::MatrixXd ZtSiXE2 = XtSiZ.transpose() * E2;
    Eigen::MatrixXd ZtSiXD2E1 = XtSiZ.transpose() * D2 * E1;

    // Terms for I(psi_j, psi_k) not needed; already have ZtSiZ and ZtSiXE1

    // Precompute F = ZtSiZ' - ZtSiXE1 and F * H_j for each j, so the double
    // loop over (jj, kk) only needs dot-product-like operations.
    Eigen::MatrixXd F = ZtSiZ.transpose() - ZtSiXE1;
    std::vector<Eigen::MatrixXd> FH(r - 1);
    for (int jj = 0; jj < r - 1; jj++) {
      FH[jj] = F * H.middleCols(jj * q, q);
    }

    for(int jj = 0; jj < r - 1; jj++) {
      // Score for psi_j
      s_psi(jj) -= 0.5 * ZtSiZ.cwiseProduct(H.middleCols(jj * q, q)).sum() -
        0.5 * ZtSiXE1.cwiseProduct(H.middleCols(jj * q, q)).sum();

      // Information for I(psi_j, psi_r)
      I_psi(jj, r - 1) = 0.5 * ZtSi2Z.cwiseProduct(H.middleCols(jj * q, q)).sum()
       - ZtSiXE2.cwiseProduct(H.middleCols(jj * q, q)).sum()
       + 0.5 * ZtSiXD2E1.cwiseProduct(H.middleCols(jj * q, q)).sum();

      // Information for I(psi_j, psi_k): uses precomputed FH
      // Original: tr(F H_jj F H_kk) = tr(FH[jj] FH[kk])
      //         = sum(FH[jj]' .* FH[kk])
      // (H_jj symmetric, so F H_jj' = F H_jj = FH[jj])
      for(int kk = 0; kk <= jj; kk++) {
        I_psi(kk, jj) = 0.5 * FH[jj].transpose().cwiseProduct(FH[kk]).sum();
      }
    }
  } else if (get_score) {
    // Terms for S(psi_j)
    Eigen::SparseMatrix<double> C = Id_q - A * ZtZ;
    Eigen::MatrixXd XtSiZ = (1.0 / psi_r) * XtZ * C;
    Eigen::MatrixXd E1 = llt.solve(XtSiZ);

    s_psi(r - 1) -= (0.5 / psi_r) * (n - q + C.diagonal().sum());
    s_psi(r - 1) += (0.5 / psi_r) * (p - E1.cwiseProduct(XtZ * A).sum());

    Eigen::SparseMatrix<double> ZtSiZ = (1.0 / psi_r) * ZtZ * C;

    Eigen::MatrixXd ZtSiXE1 = XtSiZ.transpose() * E1;
    for(int jj = 0; jj < r - 1; jj++) {
      // Score for psi_j
      s_psi(jj) -= 0.5 * ZtSiZ.cwiseProduct(H.middleCols(jj * q, q)).sum() -
        0.5 * ZtSiXE1.cwiseProduct(H.middleCols(jj * q, q)).sum();
    }
  }

  Eigen::MatrixXd U_chol = llt.matrixU();
  I_psi = I_psi.selfadjointView<Eigen::Upper>();
  return Rcpp::List::create(
    Rcpp::Named("value") = ll,
    Rcpp::Named("score") = s_psi,
    Rcpp::Named("inf_mat") = I_psi,
    Rcpp::Named("beta") = beta_tilde,
    Rcpp::Named("I_b_chol") = U_chol);
}
