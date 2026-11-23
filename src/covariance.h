#ifndef FASTKCV_COVARIANCE_H
#define FASTKCV_COVARIANCE_H

// Covariance matrices C(h) = total sill - gamma(h) for gstat variogram models.
//
// Parameterizations follow gstat (checked against
// gstat::variogramLine(covariance = TRUE)); with a = range, h the
// (anisotropy-reduced) lag distance and k = kappa, the correlation of each
// component is
//   Nug  1 if the lag vector is exactly zero, else 0
//   Exp  exp(-h / a)
//   Sph  1 - 1.5 h/a + 0.5 (h/a)^3 for h < a, else 0
//   Gau  exp(-(h / a)^2)
//   Mat  x^k K_k(x) / (2^(k-1) Gamma(k)),  x = h / a
//   Ste  as Mat with x = 2 sqrt(k) h / a  (Stein's parameterization)
//   Cir  1 - (2/pi) (r sqrt(1 - r^2) + asin r) for r = h/a < 1, else 0
//   Pen  1 - (15/8 r - 5/4 r^3 + 3/8 r^5) for r = h/a < 1, else 0
//   Bes  r K_1(r), r = h / a (1 at r = 0)
//   Wav  sin(pi r) / (pi r), r = h / a (1 at r = 0)
//   Hol  sin(r) / r, r = h / a (1 at r = 0)
// All are positive (semi)definite in two dimensions. Bounded Lin is not
// (its triangular covariance is valid in one dimension only), so it is
// rejected in R.
// Mat and Ste are evaluated on the log scale with the exponentially scaled
// Bessel function, which stays finite where the direct product
// x^k * K_k(x) underflows times overflows (large kappa, small lags); where
// even the scaled Bessel function overflows, the small-x series is used.

#include <Rcpp.h>
#include <Rmath.h>
#include <cmath>
#include <vector>

namespace fastkcv {

// Codes are positions in R's `supported_models` minus 1; new models are
// appended so that the existing codes never change.
enum ModelType {
  NUG = 0, EXP = 1, SPH = 2, GAU = 3, MAT = 4, STE = 5,
  CIR = 6, PEN = 7, BES = 8, WAV = 9, HOL = 10
};

// Regular part of the small-x series of the Matern correlation,
//   sum_j (-1)^j (x/2)^(2j) / (j! (k-1)(k-2)...(k-j)),  j < k.
// Used only where K_k(x) overflows a double, i.e. tiny x relative to a large
// kappa; the neglected singular part is then of order (x/2)^(2k), far below
// double precision.
inline double matern_series(double x, double k) {
  double q = 0.25 * x * x, term = 1.0, sum = 1.0;
  for (int j = 1; j < k; ++j) {
    term *= -q / (j * (k - j));
    sum += term;
    if (std::fabs(term) < 1e-17 * std::fabs(sum)) break;
  }
  return sum;
}

inline double matern_corr(double x, double k) {
  if (x == 0.0) return 1.0;
  double kb = R::bessel_k(x, k, 2.0);  // exp(x) * K_k(x)
  if (!std::isfinite(kb)) return matern_series(x, k);
  if (kb <= 0.0) return 0.0;
  double lv = k * std::log(x) - (k - 1.0) * M_LN2 - R::lgammafn(k) +
              std::log(kb) - x;
  return std::exp(lv);
}

inline double corr(int type, double h, double a, double k) {
  switch (type) {
  case EXP:
    return std::exp(-h / a);
  case SPH: {
    if (h >= a) return 0.0;
    double r = h / a;
    return 1.0 - 1.5 * r + 0.5 * r * r * r;
  }
  case GAU: {
    double r = h / a;
    return std::exp(-r * r);
  }
  case MAT:
    return matern_corr(h / a, k);
  case STE:
    return matern_corr(2.0 * std::sqrt(k) * h / a, k);
  case CIR: {
    if (h >= a) return 0.0;
    double r = h / a;
    return 1.0 - M_2_PI * (r * std::sqrt(1.0 - r * r) + std::asin(r));
  }
  case PEN: {
    if (h >= a) return 0.0;
    double r = h / a, r2 = r * r;
    return 1.0 - r * (1.875 - r2 * (1.25 - 0.375 * r2));
  }
  case BES: {
    double r = h / a;
    if (r == 0.0) return 1.0;
    // exp(r) K_1(r) (scaled) times r exp(-r); r K_1(r) -> 1 as r -> 0.
    return r * R::bessel_k(r, 1.0, 2.0) * std::exp(-r);
  }
  case WAV: {
    double r = M_PI * h / a;
    return r == 0.0 ? 1.0 : std::sin(r) / r;
  }
  case HOL: {
    double r = h / a;
    return r == 0.0 ? 1.0 : std::sin(r) / r;
  }
  default:
    return 0.0;
  }
}

struct Component {
  int type;
  double psill, range, kappa, sin_a, cos_a, ratio;
  bool iso;
};

inline double pair_cov(const std::vector<Component>& comp, double dx, double dy) {
  bool zero = (dx == 0.0 && dy == 0.0);
  double d_iso = std::sqrt(dx * dx + dy * dy);
  double out = 0.0;
  for (const Component& c : comp) {
    if (c.type == NUG) {
      if (zero) out += c.psill;
      continue;
    }
    double h;
    if (c.iso) {
      h = d_iso;
    } else {
      double h_major = dx * c.sin_a + dy * c.cos_a;
      double h_minor = (dx * c.cos_a - dy * c.sin_a) / c.ratio;
      h = std::sqrt(h_major * h_major + h_minor * h_minor);
    }
    out += c.psill * corr(c.type, h, c.range, c.kappa);
  }
  return out;
}

// Builds the covariance components from the vgm columns passed from R.
inline std::vector<Component> make_components(const Rcpp::IntegerVector& type,
                                              const Rcpp::NumericVector& psill,
                                              const Rcpp::NumericVector& range,
                                              const Rcpp::NumericVector& kappa,
                                              const Rcpp::NumericVector& ang1,
                                              const Rcpp::NumericVector& anis1) {
  const int K = type.size();
  std::vector<Component> comp(K);
  for (int j = 0; j < K; ++j) {
    double a = ang1[j] * M_PI / 180.0;
    comp[j].type = type[j];
    comp[j].psill = psill[j];
    comp[j].range = range[j];
    comp[j].kappa = kappa[j];
    comp[j].sin_a = std::sin(a);
    comp[j].cos_a = std::cos(a);
    comp[j].ratio = anis1[j];
    comp[j].iso = (ang1[j] == 0.0 && anis1[j] == 1.0);
  }
  return comp;
}

// Writes the n x n covariance matrix of `coords` (column-major) into `out`:
// the diagonal is the total sill, and only i < j is evaluated and mirrored,
// so the matrix is exactly symmetric.
inline void fill_cov(double* out, const Rcpp::NumericMatrix& coords,
                     const std::vector<Component>& comp) {
  const R_xlen_t n = coords.nrow();
  double total = 0.0;
  for (const Component& c : comp) total += c.psill;
  for (R_xlen_t j = 0; j < n; ++j) {
    out[j + j * n] = total;
    for (R_xlen_t i = j + 1; i < n; ++i) {
      double v = pair_cov(comp, coords(i, 0) - coords(j, 0),
                          coords(i, 1) - coords(j, 1));
      out[i + j * n] = v;
      out[j + i * n] = v;
    }
    if (j % 256 == 0) Rcpp::checkUserInterrupt();
  }
}

}  // namespace fastkcv

#endif
