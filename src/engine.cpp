// Exact kriging cross-validation from one Cholesky factorization of C.
//
// With C = R'R (R upper triangular), W = R^{-1} and Wt = W' = R^{-T}
// (lower triangular, the matrix that is stored), C^{-1} = W W' = Wt' Wt. For
// ordinary/universal kriging, B = W'F = Wt F, A = F'C^{-1}F = B'B = S'S and
// G = S^{-T} (C^{-1}F)', so that
//   Q = C^{-1} - C^{-1}F A^{-1} F'C^{-1} = Wt' Wt - G'G.
// Q depends on F only through its column space, and the R side passes an
// orthonormal basis of it (the raw F only for simple kriging, as F beta).
// Only the diagonal blocks Q[I, I] = Wt[, I]' Wt[, I] - G[, I]' G[, I] are
// formed, never the full n x n matrix Q. Storing Wt rather than W makes the
// columns Wt[, I] contiguous in memory. For simple kriging B, S and G are
// empty.
//
// The work is split in two stages so that one factorization can serve many
// fold schemes and responses:
//   prepare: Wt (overwriting a single n x n buffer), B, S, G and the diagonal
//            of Q used by leave-one-out;
//   run:     Qy and the per-fold solves (run_cpp() for a partition into
//            folds, run_buffer_cpp() for buffered leave-one-out, where the
//            hold-out sets overlap and only one point per set is scored).
//
// Errors are returned as a status code (and the position of the offending
// fold) so that the R wrapper raises exactly the same conditions as the R
// engine:
//   1 C not positive definite        4 fold leaves the trend unidentifiable
//   2 F'C^{-1}F numerically singular 5 Q[I, I] numerically singular
//   3 fold contains every observation 6 out of memory
//   7 C not symmetric                 8 C has a non-finite entry

#include <RcppArmadillo.h>

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <cstddef>
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
// Cholesky factor; false if the factorization fails.
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

// Q y for ordinary/universal kriging, Q (y - F beta) for simple kriging,
// with C^{-1} = Wt' Wt. `Beta` is p x 1 (shared by every column of Y) or
// p x T (column t for response t); F * Beta.col(t) is a gemv per column, as
// for a single response.
arma::mat compute_Qy(const arma::mat& Wt, const arma::mat& B,
                     const arma::mat& S, const arma::mat& G,
                     const arma::mat& F, const arma::mat& Y, bool simple,
                     const arma::mat& Beta) {
  if (simple) {
    arma::mat Yt = Y;
    if (Beta.n_cols == 1) {
      Yt.each_col() -= F * arma::vec(Beta.unsafe_col(0));
    } else {
      for (arma::uword t = 0; t < Yt.n_cols; ++t) {
        Yt.col(t) -= F * arma::vec(Beta.unsafe_col(t));
      }
    }
    return Wt.t() * (Wt * Yt);
  }
  arma::mat Z = Wt * Y;
  return Wt.t() * Z - G.t() * arma::solve(arma::trimatl(S.t()), B.t() * Z);
}

// Columns `idx` (sorted) of the lower triangular Wt, rows c..n-1 only: with
// c = min(idx), the rows above c are zero in every gathered column. Each
// column segment is contiguous.
arma::mat gather_cols(const arma::mat& Wt, const arma::uvec& idx) {
  const arma::uword n = Wt.n_rows, c = idx[0], m = idx.n_elem;
  arma::mat Wi(n - c, m);
  for (arma::uword a = 0; a < m; ++a) {
    std::copy(Wt.colptr(idx[a]) + c, Wt.colptr(idx[a]) + n, Wi.colptr(a));
  }
  return Wi;
}

// Q[I, I] = Wt[, I]' Wt[, I] - G[, I]' G[, I].
arma::mat q_block(const arma::mat& Wt, const arma::mat& G,
                  const arma::uvec& idx, bool simple) {
  arma::mat Wi = gather_cols(Wt, idx);
  arma::mat Qii = Wi.t() * Wi;
  if (!simple) {
    arma::mat Gi = G.cols(idx);
    Qii -= Gi.t() * Gi;
  }
  return Qii;
}

// Single entry Q[j, k], j <= k: the columns of Wt are zero above their own
// index, so only rows k..n-1 contribute.
double q_entry(const arma::mat& Wt, const arma::mat& G, arma::uword j,
               arma::uword k, bool simple) {
  const arma::uword n = Wt.n_rows;
  const double* a = Wt.colptr(j);
  const double* b = Wt.colptr(k);
  // Four independent partial sums: a single running sum is limited by the
  // latency of each add (the compiler may not reorder a floating-point
  // reduction), which made this loop the bottleneck of the pair cache.
  double s0 = 0.0, s1 = 0.0, s2 = 0.0, s3 = 0.0;
  arma::uword r = std::max(j, k);
  for (; r + 4 <= n; r += 4) {
    s0 += a[r] * b[r];
    s1 += a[r + 1] * b[r + 1];
    s2 += a[r + 2] * b[r + 2];
    s3 += a[r + 3] * b[r + 3];
  }
  for (; r < n; ++r) s0 += a[r] * b[r];
  double s = (s0 + s1) + (s2 + s3);
  if (!simple) s -= arma::dot(G.col(j), G.col(k));
  return s;
}

// Solves Q[I, I] x = e_pos through the Cholesky factor of Q[I, I], giving
// the row of Q[I, I]^{-1} for the scored point; false if the factorization
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

// In-place transpose of an n x n column-major matrix, in cache blocks.
// std::size_t, not arma::uword (32 bits here), so that i + j n cannot wrap.
void transpose_inplace(double* a, std::size_t n) {
  const std::size_t bs = 64;
  for (std::size_t jb = 0; jb < n; jb += bs) {
    for (std::size_t ib = 0; ib <= jb; ib += bs) {
      const std::size_t jmax = std::min(jb + bs, n);
      for (std::size_t j = jb; j < jmax; ++j) {
        const std::size_t imax = std::min(ib + bs, j);
        for (std::size_t i = ib; i < imax; ++i) {
          std::swap(a[i + j * n], a[j + i * n]);
        }
      }
    }
  }
}

// Factorizes the covariance matrix held in `Wr` in place and returns the
// prepared object, or a failure status. With C = R'R, the buffer ends up
// holding Wt = R^{-T} (lower triangular), so that C^{-1} = Wt' Wt and the
// rows of W = R^{-1} needed for Q[I, I] are contiguous columns of Wt.
Rcpp::List factor_inplace(Rcpp::NumericMatrix& Wr, const arma::mat& F,
                          bool simple) {
  const int n = Wr.nrow();
  if (chol_upper_inplace(Wr.begin(), n) != 0) return fail(1, 0);
  if (tri_inverse_upper_inplace(Wr.begin(), n) != 0) return fail(1, 0);
  // Zero the (unreferenced) lower triangle of R^{-1}, then transpose.
  double* w = Wr.begin();
  for (int j = 0; j < n; ++j) {
    for (int i = j + 1; i < n; ++i) w[i + std::size_t(j) * n] = 0.0;
  }
  transpose_inplace(w, n);
  arma::mat Wt(w, n, n, false, true);

  // Diagonal of C^{-1}: squared norms of the columns of Wt.
  arma::vec qdiag(n);
  for (int j = 0; j < n; ++j) {
    const double* col = Wt.colptr(j);
    double s = 0.0;
    for (int i = j; i < n; ++i) s += col[i] * col[i];
    qdiag[j] = s;
  }

  arma::mat B, S, G;
  if (!simple) {
    // F is an orthonormal basis of the trend here (trend_basis() in R), so
    // cond(B'B) <= cond(C). Below rcond eps, Armadillo's triangular solve()
    // with S would print a warning and switch to an approximate (SVD)
    // solution; a raw quadratic trend in UTM coordinates had rcond(S) =
    // 3e-18 and predictions off by 0.2. Refuse instead.
    B = Wt * F;
    if (!arma::chol(S, arma::symmatu(B.t() * B))) return fail(2, 0);
    // rcond of S' (1-norm), the matrix that solve() below works with.
    const double rc = arma::rcond(arma::mat(S.t()));
    if (!(rc >= DBL_EPSILON)) return fail(2, 0);
    G = arma::solve(arma::trimatl(S.t()), (Wt.t() * B).t());
    qdiag -= arma::sum(arma::square(G), 0).t();
  }
  return Rcpp::List::create(
      Rcpp::Named("status") = 0, Rcpp::Named("Wt") = Wr,
      Rcpp::Named("B") = B, Rcpp::Named("S") = S, Rcpp::Named("G") = G,
      Rcpp::Named("qdiag") =
          Rcpp::NumericVector(qdiag.begin(), qdiag.end()));
}

// Points within `radius` of point i (Euclidean), from candidates sorted by
// x: `ord` holds point indices and `xs` their x coordinates. Appends to
// `out` in arbitrary order. The x window is widened by a few ulps: x[i] +
// radius is rounded independently of the distance test, and without the
// slack a point at distance about `radius` could be dropped here although
// the test accepts it (and accepts i from that point's side).
void neighbors(const arma::vec& x, const arma::vec& y, const arma::uvec& ord,
                const arma::vec& xs, arma::uword i, double radius,
                std::vector<arma::uword>& out) {
  const double w = radius + 4.0 * DBL_EPSILON * (std::fabs(x[i]) + radius);
  const double* lo = std::lower_bound(xs.begin(), xs.end(), x[i] - w);
  const double* hi = std::upper_bound(xs.begin(), xs.end(), x[i] + w);
  for (arma::uword k = lo - xs.begin(); k < arma::uword(hi - xs.begin());
       ++k) {
    const arma::uword j = ord[k];
    const double dx = x[j] - x[i], dy = y[j] - y[i];
    if (j == i || std::sqrt(dx * dx + dy * dy) <= radius) out.push_back(j);
  }
}

// Checks the n x n column-major matrix `a`: 0 if every entry is finite and
// it is symmetric up to 1e-8 times its largest diagonal entry, else the
// failure status (8 a non-finite entry, 7 asymmetric). dpotrf reads only
// the upper triangle, so an asymmetric matrix would otherwise be used
// silently, and some LAPACK builds (OpenBLAS) do not flag a NaN on the
// diagonal. Visits the pairs in cache blocks, as transpose_inplace() does,
// without allocating.
int check_matrix(const double* a, std::size_t n) {
  double dmax = 0.0;
  for (std::size_t j = 0; j < n; ++j) {
    const double d = a[j + j * n];
    if (!std::isfinite(d)) return 8;
    dmax = std::max(dmax, std::fabs(d));
  }
  const double tol = 1e-8 * dmax;
  bool symmetric = true;
  const std::size_t bs = 64;
  for (std::size_t jb = 0; jb < n; jb += bs) {
    for (std::size_t ib = 0; ib <= jb; ib += bs) {
      const std::size_t jmax = std::min(jb + bs, n);
      for (std::size_t j = jb; j < jmax; ++j) {
        const std::size_t imax = std::min(ib + bs, j);
        for (std::size_t i = ib; i < imax; ++i) {
          const double u = a[i + j * n], l = a[j + i * n];
          if (!std::isfinite(u) || !std::isfinite(l)) return 8;
          if (!(std::fabs(u - l) <= tol)) symmetric = false;
        }
      }
    }
  }
  return symmetric ? 0 : 7;
}

}  // namespace

// Builds the covariance matrix of `coords` directly in the buffer that is
// then factorized, so only one n x n matrix is ever allocated.
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

// Factorizes a user-supplied covariance matrix; C itself is not modified.
// [[Rcpp::export]]
Rcpp::List prepare_matrix_cpp(const Rcpp::NumericMatrix& C,
                              const arma::mat& F, bool simple) {
  try {
    // One explicit copy: `W = Rcpp::clone(C)` duplicated the matrix twice
    // (measured with tracemem), so the peak was three n x n matrices.
    const int n = C.nrow();
    const int status = check_matrix(C.begin(), n);
    if (status != 0) return fail(status, 0);
    Rcpp::NumericMatrix W(Rcpp::no_init(n, n));
    std::copy(C.begin(), C.end(), W.begin());
    return factor_inplace(W, F, simple);
  } catch (std::bad_alloc&) {
    return fail(6, 0);
  }
}

// Cross-validation residuals and variances for one fold scheme, from a
// prepared factorization. `fold_idx` is a list of 0-based row indices, one
// element per fold, in the order in which folds are checked.
// [[Rcpp::export]]
Rcpp::List run_cpp(const Rcpp::NumericMatrix& Wr, const arma::mat& B,
                   const arma::mat& S, const arma::mat& G,
                   const arma::mat& U, const arma::vec& qdiag,
                   const arma::mat& F, const arma::mat& Y,
                   const Rcpp::List& fold_idx, bool simple,
                   const arma::mat& beta) {
  try {
    const int n = Wr.nrow();
    const arma::uword T = Y.n_cols;
    const int nf = fold_idx.size();
    const arma::mat Wt(const_cast<double*>(Wr.begin()), n, n, false, true);

    const arma::mat Qy = compute_Qy(Wt, B, S, G, F, Y, simple, beta);

    arma::mat pred(n, T);
    arma::vec var(n);

    bool loo = true;
    for (int f = 0; f < nf && loo; ++f) {
      loo = Rcpp::as<Rcpp::IntegerVector>(fold_idx[f]).size() == 1;
    }

    if (loo) {
      if (n == 1) return fail(3, 1);
      arma::vec lev;
      if (!simple) lev = arma::sum(arma::square(U), 1);
      for (int f = 0; f < nf; ++f) {
        const int i = Rcpp::as<Rcpp::IntegerVector>(fold_idx[f])[0];
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
        arma::uvec idx = arma::sort(Rcpp::as<arma::uvec>(iv));
        if (idx.n_elem == static_cast<arma::uword>(n)) return fail(3, f + 1);
        if (!simple && !identifiable(U, idx)) return fail(4, f + 1);
        arma::mat Qii = q_block(Wt, G, idx, simple);
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

// Buffered leave-one-out: for every target point i, the hold-out set is
// I = {j : ||s_j - s_i|| <= radius} (Euclidean distance on the coordinates,
// over all points) and only i is scored. The sets overlap, so this is not a
// partition and cannot go through run_cpp(). A set of size one uses the
// leave-one-out formulas, so radius 0 reproduces run_cpp() exactly.
//
// Q[I, I] comes from one of two strategies (`method`: 0 auto, 1 pairs,
// 2 points):
//   points: per target, gather the m columns of Wt and form the block,
//           about m^2 (n - min I) flops each;
//   pairs:  every pair of points that share a hold-out set lies within
//           2 radius, so Q[j, k] is computed once for all such pairs (about
//           n - k flops each) and cached (12 bytes per pair); each block is
//           then assembled from the cache by a merge walk per row.
// Auto picks the cheaper by flop count, using pairs only if the cache fits
// in `max_cache_bytes`.
//
// Errors are reported for the first failing target in `targets` order. On
// failure, `fold` is the 1-based point and `size` the size of its set.
// [[Rcpp::export]]
Rcpp::List run_buffer_cpp(const Rcpp::NumericMatrix& Wr, const arma::mat& B,
                          const arma::mat& S, const arma::mat& G,
                          const arma::mat& U, const arma::vec& qdiag,
                          const arma::mat& F, const arma::mat& Y,
                          const arma::mat& coords, double radius,
                          const arma::uvec& targets, bool simple,
                          const arma::mat& beta, int method,
                          double max_cache_bytes) {
  try {
    const arma::uword n = Wr.nrow();
    const arma::uword T = Y.n_cols;
    const arma::uword nt = targets.n_elem;
    const arma::mat Wt(const_cast<double*>(Wr.begin()), n, n, false, true);

    arma::vec lev;
    if (!simple) lev = arma::sum(arma::square(U), 1);
    const arma::vec x = coords.col(0), y = coords.col(1);
    const arma::uvec ord = arma::sort_index(x);
    const arma::vec xs = x.elem(ord);

    // Hold-out sets of all targets (sorted), as offsets into `members`.
    std::vector<arma::uword> members, buf;
    std::vector<std::size_t> off(nt + 1, 0);
    for (arma::uword t = 0; t < nt; ++t) {
      buf.clear();
      neighbors(x, y, ord, xs, targets[t], radius, buf);
      std::sort(buf.begin(), buf.end());
      members.insert(members.end(), buf.begin(), buf.end());
      off[t + 1] = members.size();
      if (t % 256 == 0) Rcpp::checkUserInterrupt();
    }
    auto set_of = [&](arma::uword t) {
      arma::uvec idx(off[t + 1] - off[t]);
      for (std::size_t a = off[t]; a < off[t + 1]; ++a) {
        idx[a - off[t]] = members[a];
      }
      return idx;
    };

    // Statuses 3 and 4 (and 5 for sets of size one) for every target, so
    // that the first failing target in order is reported.
    arma::uword first_fail = nt;
    int fail_status = 0;
    double cost_points = 0.0;
    for (arma::uword t = 0; t < nt && first_fail == nt; ++t) {
      const arma::uword i = targets[t], m = off[t + 1] - off[t];
      if (m == n) {
        first_fail = t;
        fail_status = 3;
      } else if (m == 1) {
        if (!simple && 1.0 - lev[i] < ident_tol) {
          first_fail = t;
          fail_status = 4;
        } else if (!(qdiag[i] > 0.0) || !std::isfinite(qdiag[i])) {
          first_fail = t;
          fail_status = 5;
        }
      } else {
        if (!simple && !identifiable(U, set_of(t))) {
          first_fail = t;
          fail_status = 4;
        }
        cost_points += double(m) * m * (n - members[off[t]]);
      }
    }
    const arma::uword n_score = first_fail;

    // Pair cache over the points that occur in some hold-out set of a
    // target that is scored with a block (|I| > 1).
    bool use_pairs = false;
    std::vector<std::size_t> poff;
    std::vector<arma::uword> pidx;   // sorted partner k >= j, per j
    std::vector<double> pval;
    std::vector<std::size_t> pstart(n, 0), pend(n, 0);
    if (method != 2 && radius > 0.0) {
      std::vector<char> in_union(n, 0);
      for (arma::uword t = 0; t < n_score; ++t) {
        if (off[t + 1] - off[t] > 1) {
          for (std::size_t a = off[t]; a < off[t + 1]; ++a) {
            in_union[members[a]] = 1;
          }
        }
      }
      arma::uvec uord_all(n);
      arma::uword nu = 0;
      for (arma::uword k = 0; k < n; ++k) {
        if (in_union[ord[k]]) uord_all[nu++] = ord[k];
      }
      if (nu > 0) {
        const arma::uvec uord = uord_all.head(nu);
        const arma::vec uxs = x.elem(uord);
        // Slightly above 2 radius, so rounding in the distances cannot drop
        // a pair that shares a hold-out set; a missing pair would still be
        // computed on the fly below.
        const double r2 = 2.0 * radius * (1.0 + 1e-12);
        const double max_pairs = max_cache_bytes / 12.0;
        // Pass 1 counts the pairs and their cost, so that the cache is
        // allocated once at its exact size (growing it by push_back could
        // briefly need twice the memory that the budget allows for).
        double cost_pairs = 0.0;
        std::size_t n_pairs = 0;
        bool fits = true;
        for (arma::uword k = 0; k < nu && fits; ++k) {
          const arma::uword j = uord[k];
          buf.clear();
          neighbors(x, y, uord, uxs, j, r2, buf);
          pstart[j] = n_pairs;
          for (arma::uword kk : buf) {
            if (kk >= j) {
              ++n_pairs;
              cost_pairs += double(n - kk);
            }
          }
          pend[j] = n_pairs;
          if (method == 0 && double(n_pairs) > max_pairs) fits = false;
          if (k % 256 == 0) Rcpp::checkUserInterrupt();
        }
        use_pairs = method == 1 || (fits && cost_pairs < cost_points);
        if (use_pairs) {
          // Pass 2 stores the sorted partners and the entries Q[j, k].
          pidx.resize(n_pairs);
          pval.resize(n_pairs);
          for (arma::uword k = 0; k < nu; ++k) {
            const arma::uword j = uord[k];
            buf.clear();
            neighbors(x, y, uord, uxs, j, r2, buf);
            std::sort(buf.begin(), buf.end());
            std::size_t a = pstart[j];
            for (arma::uword kk : buf) {
              if (kk >= j) {
                pidx[a] = kk;
                pval[a] = q_entry(Wt, G, j, kk, simple);
                ++a;
              }
            }
            if (k % 64 == 0) Rcpp::checkUserInterrupt();
          }
        }
      }
    }
    const arma::mat Qy = compute_Qy(Wt, B, S, G, F, Y, simple, beta);
    arma::mat pred(nt, T);
    arma::vec var(nt);
    arma::ivec size(nt);
    for (arma::uword t = 0; t < n_score; ++t) {
      const arma::uword i = targets[t], m = off[t + 1] - off[t];
      size[t] = m;
      if (m == 1) {
        var[t] = 1.0 / qdiag[i];
        pred.row(t) = Y.row(i) - Qy.row(i) / qdiag[i];
        continue;
      }
      const arma::uvec idx = set_of(t);
      arma::uword pos = 0;
      while (idx[pos] != i) ++pos;
      arma::mat Qii;
      if (use_pairs) {
        // Row a: the partners of idx[a] (sorted) are walked in step with
        // idx[a..m) (also sorted), a linear merge instead of one binary
        // search per entry (which was the bulk of the time for large m).
        Qii.set_size(m, m);
        for (arma::uword a = 0; a < m; ++a) {
          const arma::uword j = idx[a];
          std::size_t c = pstart[j];
          const std::size_t e = pend[j];
          for (arma::uword b = a; b < m; ++b) {
            const arma::uword k = idx[b];
            while (c < e && pidx[c] < k) ++c;
            const double q = (c < e && pidx[c] == k)
                                 ? pval[c]
                                 : q_entry(Wt, G, j, k, simple);
            Qii(a, b) = Qii(b, a) = q;
          }
        }
      } else {
        Qii = q_block(Wt, G, idx, simple);
      }
      arma::vec v;
      if (!solve_unit(Qii, pos, v)) return fail_buffer(5, i, m);
      var[t] = v[pos];
      pred.row(t) = Y.row(i) - v.t() * Qy.rows(idx);
      if (t % 64 == 0) Rcpp::checkUserInterrupt();
    }
    if (first_fail < nt) {
      return fail_buffer(fail_status, targets[first_fail],
                         off[first_fail + 1] - off[first_fail]);
    }
    return Rcpp::List::create(
        Rcpp::Named("status") = 0, Rcpp::Named("pred") = pred,
        Rcpp::Named("var") = Rcpp::NumericVector(var.begin(), var.end()),
        Rcpp::Named("size") = Rcpp::IntegerVector(size.begin(), size.end()),
        Rcpp::Named("method") = use_pairs ? "pairs" : "points");
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
