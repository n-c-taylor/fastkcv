# fastkcv <img src="man/figures/logo.png" align="right" height="139" alt="fastkcv logo" />

<!-- badges: start -->
[![R-CMD-check](https://github.com/n-c-taylor/fastkcv/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/n-c-taylor/fastkcv/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

Documentation: <https://n-c-taylor.github.io/fastkcv/>

Fast, exact cross-validation for global kriging, as a companion to
[gstat](https://CRAN.R-project.org/package=gstat). `kcv()` takes the same
inputs as `gstat::krige.cv()` and returns the same output. It factorizes the
covariance matrix of the observations once, instead of refitting the
kriging system for every fold, and the results agree with `krige.cv()` to
numerical precision.

* Leave-one-out, k-fold and buffered (spatial) leave-one-out
  cross-validation for simple, ordinary and universal kriging.
* gstat variogram models `Nug`, `Exp`, `Sph`, `Gau`, `Mat`, `Ste`, `Cir`,
  `Pen`, `Bes`, `Wav` and `Hol`, nested sums of them and 2D geometric
  anisotropy. data.frame, `sf` and `sp` data.
* Random, spatial block and k-means fold builders, and `kcv_metrics()` to
  compare fold designs (RMSE, bias, z-scores, interval coverage).
* Several responses at the same locations (`cbind(z1, z2) ~ 1`, or
  `kcv_run(prepared, y = )`) cross-validated from one factorization.
* `kcv_prepare()` / `kcv_run()` to factorize once and run many fold schemes;
  `kcv_engine()` for a precomputed covariance matrix.

## Installation

Using `remotes`:
```r
# install.packages("remotes")
remotes::install_github("n-c-taylor/fastkcv", build_vignettes = TRUE)
```

Using `pak` (faster):
```r
# install.packages("pak")
pak::pkg_install("n-c-taylor/fastkcv")
```

The package contains C++ code, so installing from GitHub needs a compiler
toolchain: Rtools on Windows, the Xcode command line tools on macOS.

## Example

```r
library(fastkcv)
library(gstat)
data(meuse, package = "sp")
v <- vgm(0.59, "Sph", 897, 0.05)

# Leave-one-out, as krige.cv(log(zinc) ~ 1, ~ x + y, meuse, model = v)
cv <- kcv(log(zinc) ~ 1, ~ x + y, meuse, model = v)
kcv_metrics(cv)

# One factorization, several fold designs
p <- kcv_prepare(log(zinc) ~ sqrt(dist), ~ x + y, meuse, model = v)
set.seed(1)
kcv_metrics(kcv_run(p, folds = fold_random(~ x + y, meuse, k = 10)))
kcv_metrics(kcv_run(p, folds = fold_blocks(~ x + y, meuse, blocks = 3)))
kcv_metrics(kcv_run(p, buffer = 300))
```

## Performance

Ordinary kriging on simulated data, 4-core laptop with OpenBLAS (median
seconds; from the scripts in `inst/bench/`):

| n     | LOO: fastkcv / `krige.cv()` | 5-fold: fastkcv / `krige.cv()` |
|------:|:----------------------------|:-------------------------------|
| 1,000 | 0.06 / 67 (1,100x)          | 0.07 / 0.8 (12x)               |
| 2,000 | 0.25 / 368 (1,500x)         | 0.29 / 5.7 (19x)               |
| 5,000 | 2.1 / not run               | 2.6 / 85 (32x)                 |

The gain is largest for leave-one-out and smallest for a few large folds.
`kcv()` handles n = 20,000 in about two minutes.

The vignette (`vignette("fastkcv")`) covers the fold designs, buffered
leave-one-out, precomputed covariance matrices, and the performance and
memory notes.

## Scope

Global kriging only: every observation is in the kriging neighborhood.
Local neighborhoods (`nmax`, `nmin`, `maxdist`), block kriging, unbounded
variograms, geographic coordinates and re-fitting the variogram in each
fold are not supported, and give an error. Memory is one n x n matrix
(8 n^2 bytes, 3.2 GB at n = 20,000), checked before anything is allocated.

## Related packages

fastkcv is a companion to gstat, not a fork or a replacement: `kcv()` takes
the inputs of `krige.cv()` and returns its output, and fitting variograms
and kriging new locations stay with gstat. It is a separate package for
three reasons. The single factorization works only for global kriging,
while gstat's engine sets up a kriging system for each prediction location,
which is what makes local neighborhoods (`nmax`, `maxdist`) possible. It
needs compiled C++ code (Rcpp, RcppArmadillo) and a dense n x n matrix
(8 n^2 bytes, checked up front). And most of the package is a
cross-validation toolkit rather than a faster `krige.cv()`: buffered
leave-one-out, spatial fold builders, `kcv_metrics()`, `kcv_prepare()` /
`kcv_run()` to reuse one factorization, and `kcv_engine()` for any
covariance matrix.

`leaveOneOut.km()` in DiceKriging uses the same closed-form leave-one-out
algebra, without a nugget term, with tensor-product kernels on a design
matrix rather than gstat variogram models on spatial data, and without
k-fold or buffered designs. From the same covariance matrix, fastkcv is
1.6 to 4.1 times faster end to end (OpenBLAS) and about level on the
leave-one-out step alone. blockCV and spatialsample build spatial folds
for any model but do not fit one; a partition from them, such as
`folds_ids` from blockCV's `cv_spatial()`, can be passed to
`kcv_run(folds = )`.

|                                | LOO              | k-fold | Buffered LOO | Block / cluster folds | Nugget, gstat models | Input                  | Exact, no refit per fold |
|:-------------------------------|:-----------------|:-------|:-------------|:----------------------|:---------------------|:-----------------------|:-------------------------|
| fastkcv                        | yes              | yes    | yes          | yes                   | yes                  | data.frame, `sf`, `sp` | yes                      |
| gstat `krige.cv()`             | yes              | yes    | no           | any fold vector       | yes                  | data.frame, `sf`, `sp` | no (refits)              |
| DiceKriging `leaveOneOut.km()` | yes              | no     | no           | no                    | no                   | design matrix          | yes                      |
| blockCV                        | NNDM variant     | yes    | yes          | yes                   | (folds only)         | `sf`, `sp` points      | (folds only)             |
| spatialsample                  | yes (`v = Inf`)  | yes    | yes          | yes                   | (folds only)         | `sf`                   | (folds only)             |

## Use of AI tools

This package was developed with the help of an AI coding assistant
(Anthropic's Claude Opus 5.5), used for writing and reviewing code, tests and
documentation. The author directed the work and is responsible for the
package; its results are checked against `gstat::krige.cv()` and
independent R implementations by the test suite.
