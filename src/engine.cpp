// Exact kriging cross-validation from one Cholesky factorisation of C.
//
// With C = R'R (R upper triangular) and W = R^{-1}, C^{-1} = W W'. For
// ordinary/universal kriging, B = W'F, A = F'C^{-1}F = B'B = S'S and
// G = S^{-T} (C^{-1}F)', so that
//   Q = C^{-1} - C^{-1}F A^{-1} F'C^{-1} = W W' - G'G.
// Only the diagonal blocks Q[I, I] = W[I, ] W[I, ]' - G[, I]' G[, I] are
// formed, never the full n x n matrix Q. For simple kriging G is empty.
//
// Errors are returned as a status code (and the position of the offending
// fold) so that the R wrapper raises exactly the same conditions as the R
// engine:
//   1 C not positive definite        4 fold leaves the trend unidentifiable
//   2 trend matrix rank deficient    5 Q[I, I] numerically singular
//   3 fold contains every observation 6 out of memory

#include <RcppArmadillo.h>

#include "lapack.h"

namespace {

// Folds whose remaining trend Gram matrix (in the orthonormal basis of F)
// has an eigenvalue below this are declared unidentifiable.
const double ident_tol = 1e-12;

Rcpp::List fail(int status, int fold) {
  return Rcpp::List::create(Rcpp::Named("status") = status,
                            Rcpp::Named("fold") = fold);
}

// Inverse of a small symmetric positive definite matrix through its
// Cholesky factor; false if the factorisation fails.
bool fold_inverse(const arma::mat& Qii, arma::mat& V) {
  arma::mat Rq;
  if (!arma::chol(Rq, arma::symmatu(Qii))) return false;
  arma::mat Ri = arma::inv(arma::trimatu(Rq));
  V = Ri * Ri.t();
  return true;
}

// Is the trend identifiable from the rows of F outside `idx`? U is an
// orthonormal basis of the column space of F, so this holds exactly when
// U[-idx, ]' U[-idx, ] is nonsingular. Computed from the smaller of the fold
// and its complement.
bool identifiable(const arma::mat& U, const arma::uvec& idx) {
  const arma::uword n = U.n_rows, p = U.n_cols;
  arma::mat Gm;
  if (2 * idx.n_elem <= n) {
    arma::mat Ui = U.rows(idx);
    Gm = arma::eye(p, p) - Ui.t() * Ui;
  } else {
    std::vector<char> in(n, 0);
    for (arma::uword i : idx) in[i] = 1;
    arma::uvec rest(n - idx.n_elem);
    arma::uword r = 0;
    for (arma::uword i = 0; i < n; ++i) {
      if (!in[i]) rest[r++] = i;
    }
    arma::mat Ur = U.rows(rest);
    Gm = Ur.t() * Ur;
  }
  arma::vec ev = arma::eig_sym(arma::symmatu(Gm));
  return ev.min() >= ident_tol;
}

Rcpp::List engine(const Rcpp::NumericMatrix& Cr, const arma::mat& F,
                  const arma::mat& Y, const arma::mat& U,
                  const Rcpp::List& fold_idx, bool simple,
                  const arma::vec& beta) {
  const int n = Cr.nrow();
  const arma::uword T = Y.n_cols;
  const int nf = fold_idx.size();

  // W <- R^{-1}, in place, from a single copy of C.
  arma::mat W(Cr.begin(), n, n);
  if (chol_upper_inplace(W.memptr(), n) != 0) return fail(1, 0);
  if (tri_inverse_upper_inplace(W.memptr(), n) != 0) return fail(1, 0);
  for (int j = 0; j < n; ++j) {
    double* col = W.colptr(j);
    for (int i = j + 1; i < n; ++i) col[i] = 0.0;
  }

  arma::mat Qy, G;
  if (simple) {
    arma::mat Yt = Y;
    Yt.each_col() -= F * beta;
    Qy = W * (W.t() * Yt);
  } else {
    arma::mat B = W.t() * F;
    arma::mat S;
    if (!arma::chol(S, arma::symmatu(B.t() * B))) return fail(2, 0);
    arma::mat H = W * B;  // C^{-1} F
    G = arma::solve(arma::trimatl(S.t()), H.t());
    arma::mat Z = W.t() * Y;
    arma::mat FtCiY = B.t() * Z;
    Qy = W * Z - G.t() * arma::solve(arma::trimatl(S.t()), FtCiY);
  }

  arma::mat pred(n, T);
  arma::vec var(n);

  bool loo = true;
  for (int f = 0; f < nf && loo; ++f) {
    loo = Rcpp::as<Rcpp::IntegerVector>(fold_idx[f]).size() == 1;
  }

  if (loo) {
    // Q_ii = ||W[i, ]||^2 - ||G[, i]||^2, accumulated column by column.
    arma::vec qd(n, arma::fill::zeros);
    for (int j = 0; j < n; ++j) {
      const double* col = W.colptr(j);
      for (int i = 0; i <= j; ++i) qd[i] += col[i] * col[i];
    }
    if (!simple) qd -= arma::sum(arma::square(G), 0).t();
    arma::vec lev;
    if (!simple) lev = arma::sum(arma::square(U), 1);
    for (int f = 0; f < nf; ++f) {
      int i = Rcpp::as<Rcpp::IntegerVector>(fold_idx[f])[0];
      if (n == 1) return fail(3, f + 1);
      if (!simple && 1.0 - lev[i] < ident_tol) return fail(4, f + 1);
      if (!(qd[i] > 0.0) || !std::isfinite(qd[i])) return fail(5, f + 1);
      var[i] = 1.0 / qd[i];
      pred.row(i) = Y.row(i) - Qy.row(i) / qd[i];
    }
  } else {
    for (int f = 0; f < nf; ++f) {
      Rcpp::IntegerVector iv = fold_idx[f];
      arma::uvec idx = Rcpp::as<arma::uvec>(iv);
      if (idx.n_elem == static_cast<arma::uword>(n)) return fail(3, f + 1);
      if (!simple && !identifiable(U, idx)) return fail(4, f + 1);
      arma::mat Wi = W.rows(idx);
      arma::mat Qii = Wi * Wi.t();
      if (!simple) {
        arma::mat Gi = G.cols(idx);
        Qii -= Gi.t() * Gi;
      }
      arma::mat V;
      if (!fold_inverse(Qii, V)) return fail(5, f + 1);
      pred.rows(idx) = Y.rows(idx) - V * Qy.rows(idx);
      var.elem(idx) = V.diag();
      Rcpp::checkUserInterrupt();
    }
  }
  return Rcpp::List::create(Rcpp::Named("status") = 0,
                            Rcpp::Named("pred") = pred,
                            Rcpp::Named("var") = Rcpp::NumericVector(
                                var.begin(), var.end()));
}

}  // namespace

// [[Rcpp::export]]
Rcpp::List kcv_engine_cpp(const Rcpp::NumericMatrix& C, const arma::mat& F,
                          const arma::mat& Y, const arma::mat& U,
                          const Rcpp::List& fold_idx, bool simple,
                          const arma::vec& beta) {
  try {
    return engine(C, F, Y, U, fold_idx, simple, beta);
  } catch (std::bad_alloc&) {
    return fail(6, 0);
  }
}

// Test hook for the numerically singular Q[I, I] branch.
// [[Rcpp::export]]
SEXP fold_inverse_cpp(const arma::mat& Qii) {
  arma::mat V;
  if (!fold_inverse(Qii, V)) return R_NilValue;
  return Rcpp::wrap(V);
}
