// Covariance matrices C(h) = total sill - gamma(h) for gstat variogram models.
//
// Parameterisations follow gstat (checked against
// gstat::variogramLine(covariance = TRUE)); with a = range, h the
// (anisotropy-reduced) lag distance and k = kappa, the correlation of each
// component is
//   Nug  1 if the lag vector is exactly zero, else 0
//   Exp  exp(-h / a)
//   Sph  1 - 1.5 h/a + 0.5 (h/a)^3 for h < a, else 0
//   Gau  exp(-(h / a)^2)
//   Mat  x^k K_k(x) / (2^(k-1) Gamma(k)),  x = h / a
//   Ste  as Mat with x = 2 sqrt(k) h / a  (Stein's parameterisation)
// Mat and Ste are evaluated on the log scale with the exponentially scaled
// Bessel function, which stays finite where the direct product
// x^k * K_k(x) underflows times overflows (large kappa, small lags); where
// even the scaled Bessel function overflows, the small-x series is used.

#include <Rcpp.h>
#include <Rmath.h>
#include <cmath>
#include <vector>

namespace {

enum ModelType { NUG = 0, EXP = 1, SPH = 2, GAU = 3, MAT = 4, STE = 5 };

// Regular part of the small-x series of the Matern correlation,
//   sum_j (-1)^j (x/2)^(2j) / (j! (k-1)(k-2)...(k-j)),  j < k.
// Used only where K_k(x) overflows a double, i.e. tiny x relative to a large
// kappa; the neglected singular part is then of order (x/2)^(2k), far below
// double precision.
double matern_series(double x, double k) {
  double q = 0.25 * x * x, term = 1.0, sum = 1.0;
  for (int j = 1; j < k; ++j) {
    term *= -q / (j * (k - j));
    sum += term;
    if (std::fabs(term) < 1e-17 * std::fabs(sum)) break;
  }
  return sum;
}

double matern_corr(double x, double k) {
  if (x == 0.0) return 1.0;
  double kb = R::bessel_k(x, k, 2.0);  // exp(x) * K_k(x)
  if (!std::isfinite(kb)) return matern_series(x, k);
  if (kb <= 0.0) return 0.0;
  double lv = k * std::log(x) - (k - 1.0) * M_LN2 - R::lgammafn(k) +
              std::log(kb) - x;
  return std::exp(lv);
}

double corr(int type, double h, double a, double k) {
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
  default:
    return 0.0;
  }
}

struct Component {
  int type;
  double psill, range, kappa, sin_a, cos_a, ratio;
  bool iso;
};

double pair_cov(const std::vector<Component>& comp, double dx, double dy) {
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

}  // namespace

// [[Rcpp::export]]
Rcpp::NumericMatrix cov_matrix_cpp(Rcpp::NumericMatrix coords1,
                                   Rcpp::NumericMatrix coords2,
                                   Rcpp::IntegerVector type,
                                   Rcpp::NumericVector psill,
                                   Rcpp::NumericVector range,
                                   Rcpp::NumericVector kappa,
                                   Rcpp::NumericVector ang1,
                                   Rcpp::NumericVector anis1,
                                   bool symmetric) {
  const int K = type.size();
  std::vector<Component> comp(K);
  double total = 0.0;
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
    total += psill[j];
  }

  const int n1 = coords1.nrow(), n2 = coords2.nrow();
  Rcpp::NumericMatrix out(n1, n2);
  if (symmetric) {
    for (int j = 0; j < n1; ++j) {
      out(j, j) = total;
      for (int i = j + 1; i < n1; ++i) {
        double v = pair_cov(comp, coords1(i, 0) - coords1(j, 0),
                            coords1(i, 1) - coords1(j, 1));
        out(i, j) = v;
        out(j, i) = v;
      }
      if (j % 256 == 0) Rcpp::checkUserInterrupt();
    }
  } else {
    for (int j = 0; j < n2; ++j) {
      for (int i = 0; i < n1; ++i) {
        out(i, j) = pair_cov(comp, coords1(i, 0) - coords2(j, 0),
                             coords1(i, 1) - coords2(j, 1));
      }
      if (j % 256 == 0) Rcpp::checkUserInterrupt();
    }
  }
  return out;
}
