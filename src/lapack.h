#ifndef FASTKCV_LAPACK_H
#define FASTKCV_LAPACK_H

int chol_upper_inplace(double* a, int n);
int tri_inverse_upper_inplace(double* a, int n);
double chol_rcond_upper(const double* r, int n, double anorm);

#endif
