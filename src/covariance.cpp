// Covariance matrix of a gstat variogram model at a set of locations; see
// covariance.h. Used by the tests (against gstat::variogramLine()) and the
// benchmarks; the package itself builds the matrix in place
// (prepare_model_cpp() in engine.cpp).

#include "covariance.h"

// [[Rcpp::export]]
Rcpp::NumericMatrix cov_matrix_cpp(Rcpp::NumericMatrix coords,
                                   Rcpp::IntegerVector type,
                                   Rcpp::NumericVector psill,
                                   Rcpp::NumericVector range,
                                   Rcpp::NumericVector kappa,
                                   Rcpp::NumericVector ang1,
                                   Rcpp::NumericVector anis1) {
  std::vector<fastkcv::Component> comp =
      fastkcv::make_components(type, psill, range, kappa, ang1, anis1);
  Rcpp::NumericMatrix out(coords.nrow(), coords.nrow());
  fastkcv::fill_cov(out.begin(), coords, comp);
  return out;
}
