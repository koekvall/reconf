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
//' The caller must verify that the parameters are feasible, that is,
//' \eqn{\psi_r > 0} and \eqn{\Sigma = Z \Psi Z' + \psi_r I_n} positive
//' definite, and must compute \code{A} and \code{ldetB}. Solving for
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
    ll = -0.5 * ldetB - 0.5 * (double)n * log(2.0 * M_PI * psi_r) -
      0.5 * etilde.dot(e);
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

      // B w_j precomputed once per j: the triangular loop below reuses each
      // column r - 1 times, so the stochastic psi block costs O(r q^2)
      // instead of O(r^2 q^2)
      Eigen::MatrixXd Bw(q, r - 1);
      for (int jj = 0; jj < r - 1; jj++) {
        Bw.col(jj) = B * w.middleRows(jj * q, q);
      }
      for (int jj = 0; jj < r - 1; jj++) {
        if(p > 0) {
          I.block(0, p + jj, p, 1) = (1 / psi_r) * (XtZ *(C * w.middleRows(jj * q, q)));
        }
        I(p + jj, p + r - 1) += w.middleRows(jj * q, q).dot(v);
        for (int ii = 0; ii <= jj; ii++) {
          I(p + ii, p + jj) += Bw.col(ii).dot(w.middleRows(jj * q, q));
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
//' @details See \code{?loglik} for the model and the parameterization of
//' \eqn{\Psi} through \code{H}. The restricted likelihood integrates out the
//' fixed effects \eqn{\beta}.
//'
//' The caller must verify feasibility and compute \code{A} and
//' \code{ldetB}; see \code{?loglik}.
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
    ll += Y.dot(etilde) + (n - p) * log(2.0 * M_PI) + n * log(psi_r);
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
    // G = X'Sigma^{-1}Z (p x q, dense)
    Eigen::MatrixXd G = (1.0 / psi_r) * XtZ * C;
    Eigen::MatrixXd E1 = llt.solve(G);
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
    // Information and score for psi_{-r}
    ////////////////////////////////////////////////////////////////////////////
    // All traces use the decomposition F = Z'P Z = S - G'E1, where
    // S = Z'Sigma^{-1}Z is sparse and G'E1 = G'(X'Sigma^{-1}X)^{-1}G has rank
    // p. Expanding tr(F H_j F H_k) then gives one sparse-sparse term, two
    // rank-p cross terms of size p x q, and one p x p term, so no dense
    // q x q matrix is ever formed and the cost is linear in q for
    // block-structured models:
    //   tr(F H_j F H_k) = tr(S H_j S H_k) - tr(S H_j G'E1 H_k)
    //                     - tr(S H_k G'E1 H_j) + tr(E1 H_j G' E1 H_k G')
    // using symmetry of S, H_j, and G'E1, and tr(S H_j G'E1 H_k)
    // = sum((E1 H_k S H_j) .* G).
    Eigen::SparseMatrix<double> S = (1.0 / psi_r) * ZtZ * C;
    Eigen::SparseMatrix<double> ZtSi2Z = (1.0 / psi_r) * S * C;
    Eigen::MatrixXd E3 = D2 * E1;

    std::vector<Eigen::SparseMatrix<double>> SH(r - 1);
    std::vector<Eigen::MatrixXd> Q(r - 1), QS(r - 1), Rp(r - 1);
    for (int jj = 0; jj < r - 1; jj++) {
      SH[jj] = S * H.middleCols(jj * q, q);
      Q[jj]  = E1 * H.middleCols(jj * q, q);   // p x q
      QS[jj] = Q[jj] * S;                      // p x q
      Rp[jj] = Q[jj] * G.transpose();          // p x p
    }

    for(int jj = 0; jj < r - 1; jj++) {
      // Score for psi_j: tr(S H_j) - tr(G'E1 H_j)
      s_psi(jj) -= 0.5 * S.cwiseProduct(H.middleCols(jj * q, q)).sum() -
        0.5 * Q[jj].cwiseProduct(G).sum();

      // Information for I(psi_j, psi_r)
      I_psi(jj, r - 1) = 0.5 * ZtSi2Z.cwiseProduct(H.middleCols(jj * q, q)).sum()
       - (E2 * H.middleCols(jj * q, q)).cwiseProduct(G).sum()
       + 0.5 * (E3 * H.middleCols(jj * q, q)).cwiseProduct(G).sum();

      // Information for I(psi_j, psi_k) = 0.5 tr(F H_j F H_k)
      Eigen::SparseMatrix<double> SHt = SH[jj].transpose();
      for(int kk = 0; kk <= jj; kk++) {
        double t1 = SHt.cwiseProduct(SH[kk]).sum();
        double t2 = (QS[kk] * H.middleCols(jj * q, q)).cwiseProduct(G).sum();
        double t3 = (QS[jj] * H.middleCols(kk * q, q)).cwiseProduct(G).sum();
        double t4 = Rp[jj].transpose().cwiseProduct(Rp[kk]).sum();
        I_psi(kk, jj) = 0.5 * (t1 - t2 - t3 + t4);
      }
    }
  } else if (get_score) {
    // Terms for S(psi_j); same rank-p decomposition as in the get_inf branch
    Eigen::SparseMatrix<double> C = Id_q - A * ZtZ;
    Eigen::MatrixXd G = (1.0 / psi_r) * XtZ * C;
    Eigen::MatrixXd E1 = llt.solve(G);

    s_psi(r - 1) -= (0.5 / psi_r) * (n - q + C.diagonal().sum());
    s_psi(r - 1) += (0.5 / psi_r) * (p - E1.cwiseProduct(XtZ * A).sum());

    Eigen::SparseMatrix<double> S = (1.0 / psi_r) * ZtZ * C;

    for(int jj = 0; jj < r - 1; jj++) {
      // Score for psi_j: tr(S H_j) - tr(G'E1 H_j)
      s_psi(jj) -= 0.5 * S.cwiseProduct(H.middleCols(jj * q, q)).sum() -
        0.5 * (E1 * H.middleCols(jj * q, q)).cwiseProduct(G).sum();
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

//' Log-likelihood via the n-by-n formulation
//'
//' Computes the log-likelihood, score vector, and information matrix for the
//' covariance parameters using dense n-by-n algebra. Intended for models
//' where the number of random effects \eqn{q} exceeds, or is comparable to,
//' the number of observations \eqn{n}; see \code{?loglik} for the model and
//' the q-by-q counterpart.
//'
//' @param K Dense \eqn{n \times n(r - 1)} matrix of horizontally concatenated
//'        \eqn{K_j = Z H_j Z'}. The \eqn{K_j} do not depend on \eqn{\psi} and
//'        are precomputed once by \code{get_precomp}.
//' @param psi Vector of length \eqn{r} of covariance parameters; the last
//'        element is the error variance \eqn{\psi_r}.
//' @param e Vector of length \eqn{n} of errors, or residuals, \eqn{e = Y - X\beta}.
//' @param X Matrix of size \eqn{n \times p} of predictors, of class \code{matrix}.
//' @param get_val If \code{TRUE}, the value of the log-likelihood is computed.
//' @param get_score If \code{TRUE} the score vector is calculated.
//' @param get_inf If \code{TRUE}, an information matrix is calculated.
//' @param expected If \code{TRUE}, the expected information is calculated;
//'        otherwise the observed, or negative Hessian of the log-likelihood.
//'
//' @return A list with components \code{value}, \code{score}, and
//' \code{inf_mat} as in \code{?loglik}, with score and information of
//' dimension \eqn{p + r}.
//'
//' @details Each evaluation forms \eqn{\Sigma = \psi_r I_n + \sum_j \psi_j K_j}
//' and factorizes it densely, so the cost is \eqn{O(r n^3 + r^2 n^2)},
//' independent of \eqn{q}.
//'
//' Unlike \code{loglik}, feasibility is decided here rather than by the
//' caller: \eqn{\Sigma} is positive definite iff its Cholesky factorization
//' succeeds. When \eqn{q \ge n}, \eqn{Z \Psi Z'} alone can be positive
//' definite, so \eqn{\psi_r > 0} is checked separately. At infeasible
//' parameters \code{value} is \code{-Inf} (regardless of \code{get_val}) and
//' score and information are zero.
//'
//' @useDynLib reconf, .registration=TRUE
// [[Rcpp::export]]
Rcpp::List loglik_n(const Eigen::Map<Eigen::MatrixXd> K,
                    const Eigen::Map<Eigen::VectorXd> psi,
                    const Eigen::Map<Eigen::VectorXd> e,
                    const Eigen::Map<Eigen::MatrixXd> X,
                    const bool get_val = true,
                    const bool get_score = true,
                    const bool get_inf = true,
                    const bool expected = true) {
  // Define dimensions
  const int n = e.size();
  const int p = X.cols();
  const int r = psi.size();

  // Initialize returns; zeros are also the infeasible-parameter returns
  double ll = NA_REAL;
  Eigen::VectorXd S = Eigen::VectorXd::Zero(p + r);
  Eigen::MatrixXd I = Eigen::MatrixXd::Zero(p + r, p + r);

  // Form Sigma = psi_r I_n + sum_j psi_j K_j
  Eigen::MatrixXd Sigma = Eigen::MatrixXd::Zero(n, n);
  for (int jj = 0; jj < r - 1; jj++) {
    if (psi(jj) != 0.0) { // avoid touching terms that are zero anyway
      Sigma += psi(jj) * K.middleCols(jj * n, n);
    }
  }
  Sigma.diagonal().array() += psi(r - 1);

  // Feasibility gate: the Cholesky attempt decides positive definiteness of
  // Sigma. psi_r > 0 must be checked separately because Z Psi Z' alone can
  // be positive definite when q >= n.
  Eigen::LLT<Eigen::MatrixXd> llt(Sigma);
  if (psi(r - 1) <= 0.0 || llt.info() != Eigen::Success) {
    return Rcpp::List::create(Rcpp::Named("value") = -R_PosInf,
                              Rcpp::Named("score") = S,
                              Rcpp::Named("inf_mat") = I);
  }

  // log det Sigma from the Cholesky diagonal; etilde = Sigma^{-1} e
  const double ldetS = 2.0 * llt.matrixLLT().diagonal().array().log().sum();
  Eigen::VectorXd etilde = llt.solve(e);

  if (get_val) {
    ll = -0.5 * (ldetS + n * log(2.0 * M_PI) + e.dot(etilde));
  }

  if (!get_score && !get_inf) {
    return Rcpp::List::create(Rcpp::Named("value") = ll,
                              Rcpp::Named("score") = S,
                              Rcpp::Named("inf_mat") = I);
  }

  // Si = Sigma^{-1}; all trace terms reduce to elementwise sums against Si
  // or against the products M_j = Si K_j below
  Eigen::MatrixXd Si = llt.solve(Eigen::MatrixXd::Identity(n, n));

  // Score for beta
  if (p > 0) {
    S.head(p) = X.transpose() * etilde;
  }

  // Score for psi_j: 0.5 (e'Sigma^{-1} K_j Sigma^{-1} e - tr(Sigma^{-1} K_j));
  // K_j and Si are symmetric so the trace is an elementwise product sum
  for (int jj = 0; jj < r - 1; jj++) {
    S(p + jj) = 0.5 * (etilde.dot(K.middleCols(jj * n, n) * etilde) -
      Si.cwiseProduct(K.middleCols(jj * n, n)).sum());
  }
  S(p + r - 1) = 0.5 * (etilde.squaredNorm() - Si.trace());

  if (get_inf) {
    // Information for beta: X'Sigma^{-1}X (equals the observed block)
    if (p > 0) {
      Eigen::MatrixXd W = llt.solve(X);
      I.topLeftCorner(p, p) = X.transpose() * W;
    }

    // M_j = Sigma^{-1} K_j, with K_r = I_n for the error variance, so
    // I(psi_j, psi_k) = 0.5 tr(M_j M_k) uniformly in j, k
    std::vector<Eigen::MatrixXd> M(r);
    for (int jj = 0; jj < r - 1; jj++) {
      M[jj] = Si * K.middleCols(jj * n, n);
    }
    M[r - 1] = Si;
    for (int jj = 0; jj < r; jj++) {
      for (int ii = 0; ii <= jj; ii++) {
        I(p + ii, p + jj) = 0.5 * M[ii].cwiseProduct(M[jj].transpose()).sum();
      }
    }

    if (!expected) {
      // Observed information: flip the sign of the deterministic psi block
      // and add the stochastic terms u_j' Sigma^{-1} u_k, where
      // u_j = K_j Sigma^{-1} e (u_r = Sigma^{-1} e)
      I.bottomRightCorner(r, r) = -I.bottomRightCorner(r, r);
      Eigen::MatrixXd Ue(n, r);
      for (int jj = 0; jj < r - 1; jj++) {
        Ue.col(jj) = K.middleCols(jj * n, n) * etilde;
      }
      Ue.col(r - 1) = etilde;
      Eigen::MatrixXd SUe = llt.solve(Ue);
      I.bottomRightCorner(r, r) += Ue.transpose() * SUe;
      // Cross terms with beta: X'Sigma^{-1} u_j (zero in expectation)
      if (p > 0) {
        I.topRightCorner(p, r) = X.transpose() * SUe;
      }
    }
  }
  I = I.selfadjointView<Eigen::Upper>();
  return Rcpp::List::create(Rcpp::Named("value") = ll,
                            Rcpp::Named("score") = S,
                            Rcpp::Named("inf_mat") = I);
}

//' Restricted log-likelihood via the n-by-n formulation
//'
//' Computes the restricted log-likelihood, score vector, and information
//' matrix for the covariance parameters using dense n-by-n algebra; the
//' n-side counterpart of \code{loglik_res}. See \code{?loglik_n} for when to
//' prefer this path and \code{?loglik} for the model.
//'
//' @param K Dense \eqn{n \times n(r - 1)} matrix of horizontally concatenated
//'        \eqn{K_j = Z H_j Z'}, precomputed by \code{get_precomp}.
//' @param psi Vector of length \eqn{r} of covariance parameters; the last
//'        element is the error variance \eqn{\psi_r}.
//' @param Y Vector of length \eqn{n} of responses.
//' @param X Matrix of size \eqn{n \times p} of predictors, of class \code{matrix}.
//' @param get_val If \code{TRUE}, the value of the log-likelihood is computed.
//' @param get_score If \code{TRUE} the score vector is calculated.
//' @param get_inf If \code{TRUE}, an information matrix is calculated.
//'
//' @return A list with components \code{value}, \code{score}, \code{inf_mat},
//' \code{beta}, and \code{I_b_chol} as in \code{?loglik_res}.
//'
//' @details All quantities are computed from the dense factorization of
//' \eqn{\Sigma = \psi_r I_n + \sum_j \psi_j K_j} and the projection
//' \eqn{P = \Sigma^{-1} - \Sigma^{-1} X (X'\Sigma^{-1}X)^{-1} X'\Sigma^{-1}},
//' formed explicitly as an \eqn{n \times n} matrix: the restricted score is
//' \eqn{0.5\{e'\Sigma^{-1} K_j \Sigma^{-1} e - tr(P K_j)\}} with
//' \eqn{e = Y - X\tilde\beta}, and the expected information is
//' \eqn{0.5 tr(P K_j P K_k)}.
//'
//' Feasibility is decided here as in \code{?loglik_n}: at infeasible
//' parameters, or when \eqn{X'\Sigma^{-1}X} is not positive definite,
//' \code{value} is \code{-Inf} (regardless of \code{get_val}), score and
//' information are zero, and \code{beta} and \code{I_b_chol} are \code{NA}.
//'
//' @useDynLib reconf, .registration=TRUE
// [[Rcpp::export]]
Rcpp::List loglik_res_n(const Eigen::Map<Eigen::MatrixXd> K,
                        const Eigen::Map<Eigen::VectorXd> psi,
                        const Eigen::Map<Eigen::VectorXd> Y,
                        const Eigen::Map<Eigen::MatrixXd> X,
                        const bool get_val = true,
                        const bool get_score = true,
                        const bool get_inf = true) {
  // Define dimensions
  const int n = Y.size();
  const int p = X.cols();
  const int r = psi.size();

  // Initialize returns; zeros/NAs are also the infeasible-parameter returns
  double ll = NA_REAL;
  Eigen::VectorXd s_psi = Eigen::VectorXd::Zero(r);
  Eigen::MatrixXd I_psi = Eigen::MatrixXd::Zero(r, r);
  Rcpp::List infeasible = Rcpp::List::create(
    Rcpp::Named("value") = -R_PosInf,
    Rcpp::Named("score") = s_psi,
    Rcpp::Named("inf_mat") = I_psi,
    Rcpp::Named("beta") = Eigen::VectorXd::Constant(p, NA_REAL),
    Rcpp::Named("I_b_chol") = Eigen::MatrixXd::Constant(p, p, NA_REAL));

  // Form Sigma = psi_r I_n + sum_j psi_j K_j
  Eigen::MatrixXd Sigma = Eigen::MatrixXd::Zero(n, n);
  for (int jj = 0; jj < r - 1; jj++) {
    if (psi(jj) != 0.0) {
      Sigma += psi(jj) * K.middleCols(jj * n, n);
    }
  }
  Sigma.diagonal().array() += psi(r - 1);

  // Feasibility gate; see ?loglik_n
  Eigen::LLT<Eigen::MatrixXd> llt(Sigma);
  if (psi(r - 1) <= 0.0 || llt.info() != Eigen::Success) {
    return infeasible;
  }

  const double ldetS = 2.0 * llt.matrixLLT().diagonal().array().log().sum();

  // W = Sigma^{-1} X and U = X'Sigma^{-1}X
  Eigen::MatrixXd W = llt.solve(X);
  Eigen::MatrixXd U = X.transpose() * W;
  U = U.selfadjointView<Eigen::Upper>();
  Eigen::LLT<Eigen::MatrixXd, Eigen::Upper> llt_U(U);
  if (llt_U.info() != Eigen::Success) {
    return infeasible;
  }

  // GLS coefficient, residuals, and etilde = Sigma^{-1} e = P Y
  Eigen::VectorXd beta_tilde = llt_U.solve(W.transpose() * Y);
  Eigen::VectorXd e = Y - X * beta_tilde;
  Eigen::VectorXd etilde = llt.solve(e);

  if (get_val) {
    ll = ldetS;
    ll += 2.0 * llt_U.matrixLLT().diagonal().array().log().sum();
    ll += e.dot(etilde) + (n - p) * log(2.0 * M_PI);
    ll *= -0.5;
  }

  if (get_score || get_inf) {
    // P = Sigma^{-1} - W U^{-1} W', formed explicitly (n x n, rank-p update)
    Eigen::MatrixXd P = llt.solve(Eigen::MatrixXd::Identity(n, n));
    P.noalias() -= W * llt_U.solve(W.transpose());

    // Score for psi_j: 0.5 (etilde' K_j etilde - tr(P K_j)); P and K_j are
    // symmetric so the trace is an elementwise product sum
    for (int jj = 0; jj < r - 1; jj++) {
      s_psi(jj) = 0.5 * (etilde.dot(K.middleCols(jj * n, n) * etilde) -
        P.cwiseProduct(K.middleCols(jj * n, n)).sum());
    }
    s_psi(r - 1) = 0.5 * (etilde.squaredNorm() - P.trace());

    if (get_inf) {
      // Expected information I(psi_j, psi_k) = 0.5 tr(P K_j P K_k), with
      // K_r = I_n so that the error variance is handled uniformly
      std::vector<Eigen::MatrixXd> PK(r);
      for (int jj = 0; jj < r - 1; jj++) {
        PK[jj] = P * K.middleCols(jj * n, n);
      }
      PK[r - 1] = P;
      for (int jj = 0; jj < r; jj++) {
        for (int ii = 0; ii <= jj; ii++) {
          I_psi(ii, jj) = 0.5 * PK[ii].cwiseProduct(PK[jj].transpose()).sum();
        }
      }
      I_psi = I_psi.selfadjointView<Eigen::Upper>();
    }
  }

  Eigen::MatrixXd U_chol = llt_U.matrixU();
  return Rcpp::List::create(
    Rcpp::Named("value") = ll,
    Rcpp::Named("score") = s_psi,
    Rcpp::Named("inf_mat") = I_psi,
    Rcpp::Named("beta") = beta_tilde,
    Rcpp::Named("I_b_chol") = U_chol);
}
