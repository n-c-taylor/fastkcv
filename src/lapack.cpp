// In-place LAPACK calls on a dense n x n matrix, kept apart from Armadillo
// so that R's and Armadillo's LAPACK prototypes are never seen together.

#define USE_FC_LEN_T
#include <R_ext/Lapack.h>
#ifndef FCONE
#define FCONE
#endif

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
