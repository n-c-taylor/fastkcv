# fastkcv <img src="man/figures/logo.png" align="right" height="139" alt="fastkcv logo" />

<!-- badges: start -->
[![R-CMD-check](https://github.com/n-c-taylor/fastkcv/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/n-c-taylor/fastkcv/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

Fast, exact cross-validation for global kriging, as a companion to
[gstat](https://CRAN.R-project.org/package=gstat). `kcv()` takes the same
inputs as `gstat::krige.cv()` and returns the same output. It factorizes the
covariance matrix of the observations once, instead of refitting the
kriging system for every fold, and the results agree with `krige.cv()` to
numerical precision. A covariance matrix too ill-conditioned for that
(typically a Gaussian model without a nugget) gives a warning.

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

DiceKriging computes leave-one-out (`leaveOneOut.km()`) and, since version
1.6.0, multiple-fold cross-validation (`cv()`) with the same closed-form
formulas (Dubrule 1983; Ginsbourger and Schärer 2024), for simple kriging
and for kriging with the trend re-estimated in each fold. Both handle a
homogeneous nugget; observations with individual noise variances
(`noise.var`) are not supported. `cv()` accepts folds that share points
and also returns the covariance matrices of the cross-validation errors,
which fastkcv does not. DiceKriging works with tensor-product kernels on a
design matrix rather than gstat variogram models on spatial data. From the
same covariance matrix, with or without a nugget, the two give the same
variances to about 1e-15 and predictions to about 1e-13. The exception is the smoother
`matern5_2` kernel without a nugget and with the trend re-estimated: its
covariance matrix is badly conditioned, and the predictions agree only to
between 1e-9 and 1e-7. On leave-one-out without a nugget, fastkcv is 1.6
to 4.1 times faster end to end (model setup plus leave-one-out) with
OpenBLAS and 1.4 to 2.1 times with the reference BLAS; on the
leave-one-out step alone the two are about level (0.9 to 1.8 times with
OpenBLAS, 0.75 to 0.9 times with the reference BLAS).

spmodel fits spatial linear models with its own covariance families
(including a nugget and anisotropy) and computes exact leave-one-out
(`loocv()`) and k-fold cross-validation (`kcv()`, with user folds through
`folds_index`) from one inverse of the covariance matrix, updated for each
fold, with the trend re-estimated by generalized least squares. Above
n = 5,000 it switches to local approximations by default. It reports
bias, MSPE, RMSPE, the squared correlation and interval coverage. spmodel
also exports a function named `kcv()`: when both packages are attached,
the one attached last masks the other, so call `fastkcv::kcv()` or
`spmodel::kcv()` explicitly. Given a fitted spmodel model, fastkcv's
`kcv()` stops with a pointer to `spmodel::kcv()`.

For computer experiments, hetGP (`LOO_preds()`, for homoscedastic and
heteroscedastic Gaussian process models) and RobustGaSP
(`leave_one_out_rgasp()`) also give closed-form leave-one-out predictions.

What fastkcv adds to these is gstat's variogram models and input formats,
spatial fold builders and buffered leave-one-out computed from the same
factorization, several responses from one factorization, and the split
into `kcv_prepare()` and `kcv_run()` that reuses one factorization across
fold designs.

blockCV and spatialsample build spatial folds (blocks, clusters, buffers,
nearest-neighbor distance matching) for any model, but do not fit or
predict with one. A partition from them, such as `folds_ids` from
blockCV's `cv_spatial()`, `cv_cluster()` or `cv_group()` (leave-group-out),
can be passed to `kcv_run(folds = )`. blockCV's `cv_buffer()` and
`cv_nndm()` return overlapping hold-out sets (`folds_list`) rather than a
partition; fastkcv's own `buffer` argument covers buffered leave-one-out.

|                                          | LOO             | k-fold | Buffered LOO                   | Block / cluster folds       | Models                                     | Input                    | One factorization, exact     | CV error covariance |
|:-----------------------------------------|:----------------|:-------|:-------------------------------|:----------------------------|:-------------------------------------------|:-------------------------|:-----------------------------|:--------------------|
| fastkcv                                  | yes             | yes    | yes                            | builders, any fold vector   | gstat models, nugget, anisotropy           | data.frame, `sf`, `sp`   | yes                          | no                  |
| gstat `krige.cv()`                       | yes             | yes    | no                             | any fold vector             | gstat models, nugget, anisotropy           | data.frame, `sf`, `sp`   | no (refits each fold)        | no                  |
| DiceKriging `leaveOneOut.km()`, `cv()`   | yes             | yes    | overlapping folds, no builder  | any index list              | tensor-product kernels, homogeneous nugget | design matrix            | yes                          | yes                 |
| spmodel `loocv()`, `kcv()`               | yes             | yes    | no                             | any fold vector             | own families, nugget, anisotropy           | data.frame, `sf`         | yes (local above n = 5,000)  | no                  |
| blockCV                                  | NNDM, kNNDM     | yes    | yes                            | yes, and leave-group-out    | (folds only)                               | `sf`, `sp` points        | (folds only)                 | (folds only)        |
| spatialsample                            | yes (`v = Inf`) | yes    | yes                            | yes                         | (folds only)                               | `sf`                     | (folds only)                 | (folds only)        |

"(folds only)": the package returns fold assignments and fits no model.
NNDM is nearest-neighbor distance matching (`cv_nndm()`), a variant of
leave-one-out; blockCV's buffered leave-one-out is `cv_buffer()`.
spatialsample's buffered leave-one-out is `spatial_buffer_vfold_cv()` with
`v = Inf` and a `buffer`.

## Use of AI tools

This package was developed with the help of an AI coding assistant
(Anthropic's Claude Opus 5.5), used for writing and reviewing code, tests and
documentation. The author directed the work and is responsible for the
package; its results are checked against `gstat::krige.cv()` and
independent R implementations by the test suite.
