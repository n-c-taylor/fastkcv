// In-place LAPACK calls on a dense n x n matrix, kept apart from Armadillo
// so that R's and Armadillo's LAPACK prototypes are never seen together.

#define USE_FC_LEN_T
#include <R_ext/Lapack.h>
#ifndef FCONE
#define FCONE
#endif

#include <cstddef>
#include <vector>

#include "lapack.h"

// Overwrites the upper triangle of `a` with R, where A = R'R. Returns the
// LAPACK info code (0 on success).
int chol_upper_inplace(double* a, int n) {
  int info = 0;
  F77_CALL(dpotrf)("U", &n, a, &n, &info FCONE);
  return info;
}

// Overwrites the upper triangular matrix in `a` with its inverse.
int tri_inverse_upper_inplace(double* a, int n) {
  int info = 0;
  F77_CALL(dtrtri)("U", "N", &n, a, &n, &info FCONE FCONE);
  return info;
}

// Reciprocal condition number (1-norm) of A = R'R from its Cholesky factor
// R (upper triangle of `r`) and `anorm`, the 1-norm of A. An estimate
// (dpocon): within a small factor of the true value, O(n^2) operations.
// Returns -1 if LAPACK reports an error.
double chol_rcond_upper(const double* r, int n, double anorm) {
  int info = 0;
  double rcond = 0.0;
  std::vector<double> work(3 * static_cast<std::size_t>(n));
  std::vector<int> iwork(n);
  F77_CALL(dpocon)("U", &n, r, &n, &anorm, &rcond, work.data(), iwork.data(),
                   &info FCONE);
  return info == 0 ? rcond : -1.0;
}
