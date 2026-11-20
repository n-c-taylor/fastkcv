# fastkcv

<!-- badges: start -->
[![R-CMD-check](https://github.com/n-c-taylor/fastkcv/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/n-c-taylor/fastkcv/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

Fast, exact cross-validation for global kriging, as a companion to
[gstat](https://CRAN.R-project.org/package=gstat). `kcv()` takes the same
inputs as `gstat::krige.cv()` and returns the same output. It factorizes the
covariance matrix of the observations once, instead of refitting the
kriging system for every fold, and the results agree with `krige.cv()` to
numerical precision.

* Leave-one-out, k-fold and buffered (spatial) leave-one-out
  cross-validation for simple, ordinary and universal kriging.
* gstat variogram models `Nug`, `Exp`, `Sph`, `Gau`, `Mat` and `Ste`, nested
  sums of them and 2D geometric anisotropy. data.frame, `sf` and `sp` data.
* Random, spatial block and k-means fold builders, and `kcv_metrics()` to
  compare fold designs (RMSE, bias, z-scores, interval coverage).
* `kcv_prepare()` / `kcv_run()` to factorize once and run many fold schemes;
  `kcv_engine()` for a precomputed covariance matrix.

## Installation

```r
# install.packages("remotes")
remotes::install_github("n-c-taylor/fastkcv")
```

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
