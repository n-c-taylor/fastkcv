// Covariance matrices for gstat variogram models; see covariance.h.

#include "covariance.h"

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
  std::vector<fastkcv::Component> comp =
      fastkcv::make_components(type, psill, range, kappa, ang1, anis1);
  Rcpp::NumericMatrix out(coords1.nrow(),
                          symmetric ? coords1.nrow() : coords2.nrow());
  fastkcv::fill_cov(out.begin(), coords1, coords2, comp, symmetric);
  return out;
}
