// High-precision oracle for bench-conditioning.R: kriging cross-validation
// in quadruple precision (__float128, a 113-bit significand, machine
// epsilon about 1.9e-34), for small n. It takes the double-precision
// covariance matrix C as exact input and computes, with the algebra of
// kcv_engine(),
//   Q = C^-1 - C^-1 F (F'C^-1 F)^-1 F'C^-1   (ordinary/universal kriging)
// and, for each fold I, residual Q[I, I]^-1 (Q y)[I] and variance
// diag(Q[I, I]^-1). The rounding error is about cond(C) * 1e-34, far
// below double precision for any matrix whose double Cholesky succeeds
// (cond(C) up to about 1e16 to 1e19).
//
// Plain arithmetic on __float128 needs only GCC's libgcc; the square root
// is a few Newton steps from the double square root, so libquadmath is not
// needed. Compiled with Rcpp::sourceCpp(); not part of the package.

#include <Rcpp.h>

#include <cmath>
#include <vector>

typedef __float128 qd;

namespace {

qd qsqrt(qd a) {
  if (a <= 0) return 0;
  qd x = std::sqrt(static_cast<double>(a));
  for (int k = 0; k < 4; ++k) x = (x + a / x) / 2;
  return x;
}

// Dense column-major matrix.
struct Mat {
  int nr, nc;
  std::vector<qd> v;
  Mat(int r, int c) : nr(r), nc(c), v(std::size_t(r) * c, 0) {}
  qd& operator()(int i, int j) { return v[i + std::size_t(j) * nr]; }
  qd operator()(int i, int j) const { return v[i + std::size_t(j) * nr]; }
};

// Lower Cholesky factor of a symmetric positive definite matrix; false if
// a pivot is not positive.
bool chol(const Mat& A, Mat& L) {
  const int n = A.nr;
  L = Mat(n, n);
  for (int j = 0; j < n; ++j) {
    qd s = A(j, j);
    for (int k = 0; k < j; ++k) s -= L(j, k) * L(j, k);
    if (!(s > 0)) return false;
    L(j, j) = qsqrt(s);
    for (int i = j + 1; i < n; ++i) {
      qd t = A(i, j);
      for (int k = 0; k < j; ++k) t -= L(i, k) * L(j, k);
      L(i, j) = t / L(j, j);
    }
  }
  return true;
}

// Solves A X = B given the lower Cholesky factor L of A (B overwritten).
void chol_solve(const Mat& L, Mat& B) {
  const int n = L.nr;
  for (int c = 0; c < B.nc; ++c) {
    for (int i = 0; i < n; ++i) {
      qd t = B(i, c);
      for (int k = 0; k < i; ++k) t -= L(i, k) * B(k, c);
      B(i, c) = t / L(i, i);
    }
    for (int i = n - 1; i >= 0; --i) {
      qd t = B(i, c);
      for (int k = i + 1; k < n; ++k) t -= L(k, i) * B(k, c);
      B(i, c) = t / L(i, i);
    }
  }
}

Mat from_r(const Rcpp::NumericMatrix& M) {
  Mat A(M.nrow(), M.ncol());
  for (int j = 0; j < M.ncol(); ++j) {
    for (int i = 0; i < M.nrow(); ++i) A(i, j) = M(i, j);
  }
  return A;
}

}  // namespace

// [[Rcpp::export]]
Rcpp::List quad_kcv(const Rcpp::NumericMatrix& C, const Rcpp::NumericMatrix& F,
                    const Rcpp::NumericVector& y, const Rcpp::List& folds,
                    bool simple, const Rcpp::NumericVector& beta) {
  const int n = C.nrow(), p = F.ncol();
  Mat L(0, 0);
  if (!chol(from_r(C), L)) Rcpp::stop("quad Cholesky of C failed");
  // Q = C^-1, then the GLS correction for ordinary/universal kriging.
  Mat Q(n, n);
  for (int i = 0; i < n; ++i) Q(i, i) = 1;
  chol_solve(L, Q);
  std::vector<qd> r(n);
  for (int i = 0; i < n; ++i) {
    r[i] = y[i];
    if (simple) {
      for (int k = 0; k < p; ++k) r[i] -= qd(F(i, k)) * qd(beta[k]);
    }
  }
  if (!simple) {
    Mat CiF = from_r(F);
    chol_solve(L, CiF);
    Mat A(p, p);
    for (int a = 0; a < p; ++a) {
      for (int b = 0; b < p; ++b) {
        qd s = 0;
        for (int i = 0; i < n; ++i) s += qd(F(i, a)) * CiF(i, b);
        A(a, b) = s;
      }
    }
    Mat LA(0, 0);
    if (!chol(A, LA)) Rcpp::stop("quad Cholesky of F'C^-1 F failed");
    // X = A^-1 CiF', so Q -= CiF X.
    Mat X(p, n);
    for (int i = 0; i < n; ++i) {
      for (int a = 0; a < p; ++a) X(a, i) = CiF(i, a);
    }
    chol_solve(LA, X);
    for (int j = 0; j < n; ++j) {
      for (int i = 0; i < n; ++i) {
        qd s = 0;
        for (int a = 0; a < p; ++a) s += CiF(i, a) * X(a, j);
        Q(i, j) -= s;
      }
    }
  }
  std::vector<qd> Qy(n, 0);
  for (int j = 0; j < n; ++j) {
    for (int i = 0; i < n; ++i) Qy[i] += Q(i, j) * r[j];
  }
  Rcpp::NumericVector pred(n), var(n);
  for (int f = 0; f < folds.size(); ++f) {
    Rcpp::IntegerVector I = folds[f];  // 1-based
    const int m = I.size();
    Mat Qii(m, m);
    for (int a = 0; a < m; ++a) {
      for (int b = 0; b < m; ++b) Qii(a, b) = Q(I[a] - 1, I[b] - 1);
    }
    Mat Lq(0, 0);
    if (!chol(Qii, Lq)) Rcpp::stop("quad Cholesky of Q[I, I] failed");
    Mat V(m, m);
    for (int a = 0; a < m; ++a) V(a, a) = 1;
    chol_solve(Lq, V);
    for (int a = 0; a < m; ++a) {
      qd e = 0;
      for (int b = 0; b < m; ++b) e += V(a, b) * Qy[I[b] - 1];
      const int i = I[a] - 1;
      pred[i] = static_cast<double>(qd(y[i]) - e);
      var[i] = static_cast<double>(V(a, a));
    }
  }
  return Rcpp::List::create(Rcpp::Named("pred") = pred,
                            Rcpp::Named("var") = var);
}
