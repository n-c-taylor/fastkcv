// Exact kriging cross-validation from one Cholesky factorisation of C.
//
// With C = R'R (R upper triangular) and W = R^{-1}, C^{-1} = W W'. For
// ordinary/universal kriging, B = W'F, A = F'C^{-1}F = B'B = S'S and
// G = S^{-T} (C^{-1}F)', so that
//   Q = C^{-1} - C^{-1}F A^{-1} F'C^{-1} = W W' - G'G.
// Only the diagonal blocks Q[I, I] = W[I, ] W[I, ]' - G[, I]' G[, I] are
// formed, never the full n x n matrix Q. For simple kriging B, S and G are
// empty.
//
// The work is split in two stages so that one factorisation can serve many
// fold schemes and responses:
//   prepare: W (overwriting a single n x n buffer), B, S, G and the diagonal
//            of Q used by leave-one-out;
//   run:     Qy and the per-fold solves (run_cpp() for a partition into
//            folds, run_buffer_cpp() for buffered leave-one-out, where the
//            hold-out sets overlap and only one point per set is scored).
//
// Errors are returned as a status code (and the position of the offending
// fold) so that the R wrapper raises exactly the same conditions as the R
// engine:
//   1 C not positive definite        4 fold leaves the trend unidentifiable
//   2 trend matrix rank deficient    5 Q[I, I] numerically singular
//   3 fold contains every observation 6 out of memory

#include <RcppArmadillo.h>

#include <algorithm>
#include <cmath>
#include <vector>

#include "covariance.h"
#include "lapack.h"

namespace {

// Folds whose remaining trend Gram matrix (in the orthonormal basis of F)
// has an eigenvalue below this are declared unidentifiable.
const double ident_tol = 1e-12;

Rcpp::List fail(int status, int fold) {
  return Rcpp::List::create(Rcpp::Named("status") = status,
                            Rcpp::Named("fold") = fold);
}

Rcpp::List fail_buffer(int status, arma::uword i, arma::uword m) {
  return Rcpp::List::create(Rcpp::Named("status") = status,
                            Rcpp::Named("fold") = int(i + 1),
                            Rcpp::Named("size") = int(m));
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

// Q y for ordinary/universal kriging, Q (y - F beta) for simple kriging.
arma::mat compute_Qy(const arma::mat& W, const arma::mat& B,
                     const arma::mat& S, const arma::mat& G,
                     const arma::mat& F, const arma::mat& Y, bool simple,
                     const arma::vec& beta) {
  if (simple) {
    arma::mat Yt = Y;
    Yt.each_col() -= F * beta;
    return W * (W.t() * Yt);
  }
  arma::mat Z = W.t() * Y;
  return W * Z - G.t() * arma::solve(arma::trimatl(S.t()), B.t() * Z);
}

// Solves Q[I, I] x = e_pos through the Cholesky factor of Q[I, I], giving
// the row of Q[I, I]^{-1} for the scored point; false if the factorisation
// fails.
bool solve_unit(const arma::mat& Qii, arma::uword pos, arma::vec& x) {
  arma::mat Rq;
  if (!arma::chol(Rq, arma::symmatu(Qii))) return false;
  arma::vec e(Qii.n_rows, arma::fill::zeros);
  e[pos] = 1.0;
  arma::vec z = arma::solve(arma::trimatl(Rq.t()), e);
  x = arma::solve(arma::trimatu(Rq), z);
  return true;
}

// Factorises the covariance matrix held in `Wr` in place (Wr <- R^{-1}) and
// returns the prepared object, or a failure status.
Rcpp::List factor_inplace(Rcpp::NumericMatrix& Wr, const arma::mat& F,
                          bool simple) {
  const int n = Wr.nrow();
  if (chol_upper_inplace(Wr.begin(), n) != 0) return fail(1, 0);
  if (tri_inverse_upper_inplace(Wr.begin(), n) != 0) return fail(1, 0);
  arma::mat W(Wr.begin(), n, n, false, true);
  for (int j = 0; j < n; ++j) {
    double* col = W.colptr(j);
    for (int i = j + 1; i < n; ++i) col[i] = 0.0;
  }

  // Diagonal of C^{-1}: squared row norms of W, accumulated column by column.
  arma::vec qdiag(n, arma::fill::zeros);
  for (int j = 0; j < n; ++j) {
    const double* col = W.colptr(j);
    for (int i = 0; i <= j; ++i) qdiag[i] += col[i] * col[i];
  }

  arma::mat B, S, G;
  if (!simple) {
    B = W.t() * F;
    if (!arma::chol(S, arma::symmatu(B.t() * B))) return fail(2, 0);
    G = arma::solve(arma::trimatl(S.t()), (W * B).t());
    qdiag -= arma::sum(arma::square(G), 0).t();
  }
  return Rcpp::List::create(
      Rcpp::Named("status") = 0, Rcpp::Named("W") = Wr,
      Rcpp::Named("B") = B, Rcpp::Named("S") = S, Rcpp::Named("G") = G,
      Rcpp::Named("qdiag") =
          Rcpp::NumericVector(qdiag.begin(), qdiag.end()));
}

}  // namespace

// Builds the covariance matrix of `coords` directly in the buffer that is
// then factorised, so only one n x n matrix is ever allocated.
// [[Rcpp::export]]
Rcpp::List prepare_model_cpp(const Rcpp::NumericMatrix& coords,
                             const Rcpp::IntegerVector& type,
                             const Rcpp::NumericVector& psill,
                             const Rcpp::NumericVector& range,
                             const Rcpp::NumericVector& kappa,
                             const Rcpp::NumericVector& ang1,
                             const Rcpp::NumericVector& anis1,
                             const arma::mat& F, bool simple) {
  try {
    std::vector<fastkcv::Component> comp =
        fastkcv::make_components(type, psill, range, kappa, ang1, anis1);
    Rcpp::NumericMatrix W(Rcpp::no_init(coords.nrow(), coords.nrow()));
    fastkcv::fill_cov(W.begin(), coords, coords, comp, true);
    return factor_inplace(W, F, simple);
  } catch (std::bad_alloc&) {
    return fail(6, 0);
  }
}

// Factorises a user-supplied covariance matrix; C itself is not modified.
// [[Rcpp::export]]
Rcpp::List prepare_matrix_cpp(const Rcpp::NumericMatrix& C,
                              const arma::mat& F, bool simple) {
  try {
    Rcpp::NumericMatrix W = Rcpp::clone(C);
    return factor_inplace(W, F, simple);
  } catch (std::bad_alloc&) {
    return fail(6, 0);
  }
}

// Cross-validation residuals and variances for one fold scheme, from a
// prepared factorisation. `fold_idx` is a list of 0-based row indices, one
// element per fold, in the order in which folds are checked.
// [[Rcpp::export]]
Rcpp::List run_cpp(const Rcpp::NumericMatrix& Wr, const arma::mat& B,
                   const arma::mat& S, const arma::mat& G,
                   const arma::mat& U, const arma::vec& qdiag,
                   const arma::mat& F, const arma::mat& Y,
                   const Rcpp::List& fold_idx, bool simple,
                   const arma::vec& beta) {
  try {
    const int n = Wr.nrow();
    const arma::uword T = Y.n_cols;
    const int nf = fold_idx.size();
    const arma::mat W(const_cast<double*>(Wr.begin()), n, n, false, true);

    const arma::mat Qy = compute_Qy(W, B, S, G, F, Y, simple, beta);

    arma::mat pred(n, T);
    arma::vec var(n);

    bool loo = true;
    for (int f = 0; f < nf && loo; ++f) {
      loo = Rcpp::as<Rcpp::IntegerVector>(fold_idx[f]).size() == 1;
    }

    if (loo) {
      arma::vec lev;
      if (!simple) lev = arma::sum(arma::square(U), 1);
      for (int f = 0; f < nf; ++f) {
        int i = Rcpp::as<Rcpp::IntegerVector>(fold_idx[f])[0];
        if (n == 1) return fail(3, f + 1);
        if (!simple && 1.0 - lev[i] < ident_tol) return fail(4, f + 1);
        if (!(qdiag[i] > 0.0) || !std::isfinite(qdiag[i])) {
          return fail(5, f + 1);
        }
        var[i] = 1.0 / qdiag[i];
        pred.row(i) = Y.row(i) - Qy.row(i) / qdiag[i];
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
    return Rcpp::List::create(
        Rcpp::Named("status") = 0, Rcpp::Named("pred") = pred,
        Rcpp::Named("var") = Rcpp::NumericVector(var.begin(), var.end()));
  } catch (std::bad_alloc&) {
    return fail(6, 0);
  }
}

// Buffered leave-one-out: for every point i, the hold-out set is
// I = {j : ||s_j - s_i|| <= radius} (Euclidean distance on the coordinates)
// and only point i is scored. The sets overlap, so this is not a partition
// and cannot go through run_cpp(). Neighbours are found on the fly from the
// points sorted by x, so no neighbour lists are stored. A set of size one
// uses the leave-one-out formulas, so radius 0 reproduces run_cpp() exactly.
// On failure, `fold` is the 1-based point and `size` the size of its set.
// [[Rcpp::export]]
Rcpp::List run_buffer_cpp(const Rcpp::NumericMatrix& Wr, const arma::mat& B,
                          const arma::mat& S, const arma::mat& G,
                          const arma::mat& U, const arma::vec& qdiag,
                          const arma::mat& F, const arma::mat& Y,
                          const arma::mat& coords, double radius,
                          bool simple, const arma::vec& beta) {
  try {
    const arma::uword n = Wr.nrow();
    const arma::uword T = Y.n_cols;
    const arma::mat W(const_cast<double*>(Wr.begin()), n, n, false, true);
    const arma::mat Qy = compute_Qy(W, B, S, G, F, Y, simple, beta);

    arma::vec lev;
    if (!simple) lev = arma::sum(arma::square(U), 1);
    const arma::vec x = coords.col(0), y = coords.col(1);
    const arma::uvec ord = arma::sort_index(x);
    const arma::vec xs = x.elem(ord);

    arma::mat pred(n, T);
    arma::vec var(n), size(n);
    std::vector<arma::uword> members;
    for (arma::uword i = 0; i < n; ++i) {
      // Points with x in [x_i - radius, x_i + radius], then the distance test.
      const double* lo =
          std::lower_bound(xs.begin(), xs.end(), x[i] - radius);
      const double* hi =
          std::upper_bound(xs.begin(), xs.end(), x[i] + radius);
      members.clear();
      for (arma::uword k = lo - xs.begin(); k < arma::uword(hi - xs.begin());
           ++k) {
        const arma::uword j = ord[k];
        const double dx = x[j] - x[i], dy = y[j] - y[i];
        if (j == i || std::sqrt(dx * dx + dy * dy) <= radius) {
          members.push_back(j);
        }
      }
      std::sort(members.begin(), members.end());
      const arma::uword m = members.size();
      size[i] = m;

      if (m == n) return fail_buffer(3, i, m);
      if (m == 1) {
        if (!simple && 1.0 - lev[i] < ident_tol) return fail_buffer(4, i, m);
        if (!(qdiag[i] > 0.0) || !std::isfinite(qdiag[i])) {
          return fail_buffer(5, i, m);
        }
        var[i] = 1.0 / qdiag[i];
        pred.row(i) = Y.row(i) - Qy.row(i) / qdiag[i];
        continue;
      }

      arma::uvec idx(m);
      arma::uword pos = 0;
      for (arma::uword a = 0; a < m; ++a) {
        idx[a] = members[a];
        if (members[a] == i) pos = a;
      }
      if (!simple && !identifiable(U, idx)) return fail_buffer(4, i, m);
      // Row idx[a] of the upper triangular W is zero left of column idx[a],
      // so columns before min(idx) = idx[0] contribute nothing.
      arma::mat Wi = W.submat(idx, arma::regspace<arma::uvec>(idx[0], n - 1));
      arma::mat Qii = Wi * Wi.t();
      if (!simple) {
        arma::mat Gi = G.cols(idx);
        Qii -= Gi.t() * Gi;
      }
      arma::vec v;
      if (!solve_unit(Qii, pos, v)) return fail_buffer(5, i, m);
      var[i] = v[pos];
      pred.row(i) = Y.row(i) - v.t() * Qy.rows(idx);
      if (i % 64 == 0) Rcpp::checkUserInterrupt();
    }
    return Rcpp::List::create(
        Rcpp::Named("status") = 0, Rcpp::Named("pred") = pred,
        Rcpp::Named("var") = Rcpp::NumericVector(var.begin(), var.end()),
        Rcpp::Named("size") = Rcpp::IntegerVector(size.begin(), size.end()));
  } catch (std::bad_alloc&) {
    return fail(6, 0);
  }
}

// Test hook for the numerically singular Q[I, I] branch of
// run_buffer_cpp(): the row of Q[I, I]^{-1} for position `pos` (0-based),
// or NULL.
// [[Rcpp::export]]
SEXP buffer_solve_cpp(const arma::mat& Qii, int pos) {
  arma::vec x;
  if (!solve_unit(Qii, pos, x)) return R_NilValue;
  return Rcpp::wrap(x);
}

// Test hook for the numerically singular Q[I, I] branch.
// [[Rcpp::export]]
SEXP fold_inverse_cpp(const arma::mat& Qii) {
  arma::mat V;
  if (!fold_inverse(Qii, V)) return R_NilValue;
  return Rcpp::wrap(V);
}
