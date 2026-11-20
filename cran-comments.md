## Submission

This is a new submission.

## Test environments

* Local: Linux Mint 22.3 (Ubuntu 24.04 base), R 4.6.2, with the reference
  BLAS/LAPACK and with OpenBLAS 0.3.26.
* GitHub Actions: ubuntu-latest, R release.

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new submission.

The local check also reports the non-portable compilation flag
`-mno-omit-leaf-frame-pointer`. It comes from the system R's `Makeconf`
(Ubuntu `r-base` build flags), not from the package's `Makevars`.
