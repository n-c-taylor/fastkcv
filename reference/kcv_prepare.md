# Factorize once, cross-validate many times

`kcv_prepare()` does the expensive part of
[`kcv()`](https://n-c-taylor.github.io/fastkcv/reference/kcv.md)
(building and factorizing the covariance matrix, which depends only on
the locations, the variogram model and the trend) and returns it as an
object. `kcv_run()` then computes the cross-validation for one fold
scheme from it, at a fraction of the cost. `kcv(...)` is equivalent to
`kcv_run(kcv_prepare(...), ...)`.

## Usage

``` r
kcv_prepare(formula, locations, data, model, beta = NULL, ...)

kcv_run(
  prepared,
  nfold,
  folds = NULL,
  buffer = NULL,
  targets = NULL,
  y = NULL,
  beta = NULL
)
```

## Arguments

- formula:

  Formula defining the response and the trend, e.g. `log(zinc) ~ 1` or
  `log(zinc) ~ sqrt(dist)`. A matrix response, such as
  `cbind(log(zinc), log(copper)) ~ 1` or a matrix column of `data`,
  cross-validates each column with the same model, trend and folds; see
  "Several responses" below.

- locations:

  A formula giving the coordinates (e.g. `~x + y`) when `data` is a
  data.frame, or an `sf` object or `sp` `Spatial*DataFrame` holding both
  the locations and the variables (then `data` is not used).

- data:

  A data.frame, when `locations` is a formula.

- model:

  A `gstat` variogram model
  ([`gstat::vgm()`](https://r-spatial.github.io/gstat/reference/vgm.html))
  with fixed parameters. Supported: `Nug`, `Exp`, `Sph`, `Gau`, `Mat`,
  `Ste`, `Cir`, `Pen`, `Bes`, `Wav`, `Hol` and nested sums of them, with
  optional 2D geometric anisotropy.

- beta:

  Optional, for a `prepared` object made for simple kriging (with
  `beta`): new trend coefficients, a vector of length p or a p x T
  matrix as in
  [`kcv()`](https://n-c-taylor.github.io/fastkcv/reference/kcv.md). The
  factorization of simple kriging does not depend on them.

- ...:

  As in
  [`kcv()`](https://n-c-taylor.github.io/fastkcv/reference/kcv.md).

- prepared:

  A `kcv_prepared` object from `kcv_prepare()`.

- nfold:

  As in
  [`gstat::krige.cv()`](https://r-spatial.github.io/gstat/reference/krige.cv.html):
  missing or `nrow(data)` for leave-one-out; a single integer `k < n`
  for random k-fold assignment (drawn exactly as
  [`krige.cv()`](https://r-spatial.github.io/gstat/reference/krige.cv.html)
  does, so the same seed gives the same folds); or a vector of length
  `n` of fold identifiers.

- folds:

  Optional vector of length `n` of fold identifiers, as an alternative
  to a vector `nfold`. Giving both `nfold` and `folds` is an error.

- buffer:

  Optional radius for buffered leave-one-out, in coordinate units (\\\ge
  0\\). Each observation is predicted from the observations farther than
  `buffer` from it; see "Buffered leave-one-out" below. Cannot be
  combined with `nfold` or `folds`.

- targets:

  Optional, with `buffer`: the observations to predict, as row numbers
  or a logical vector of length `n`. The hold-out sets are still formed
  from all observations. Useful to estimate buffered leave-one-out
  statistics from a sample of points when the radius is large. By
  default every observation is predicted.

- y:

  Optional new responses, replacing those of the data given to
  `kcv_prepare()`: a numeric vector of length n (output as from a single
  response) or an n x T matrix (T responses, output with a `response`
  column), with rows in the order of that data. The locations, trend
  covariates and model are those of `prepared`. Missing or non-finite
  values are an error.

## Value

`kcv_prepare()`: a `kcv_prepared` object. `kcv_run()`: the same output
as [`kcv()`](https://n-c-taylor.github.io/fastkcv/reference/kcv.md).

## Details

The prepared object holds an \\n \times n\\ matrix (\\8 n^2\\ bytes) and
can be saved with [`saveRDS()`](https://rdrr.io/r/base/readRDS.html).

`kcv_run(prepared, y = )` cross-validates new responses at the same
locations, with the trend matrix and model stored in `prepared`, without
another factorization: for example, many fields simulated at the
observation locations. See "Several responses" in
[`kcv()`](https://n-c-taylor.github.io/fastkcv/reference/kcv.md).

## See also

[`kcv_engine_prepare()`](https://n-c-taylor.github.io/fastkcv/reference/kcv_engine.md)
for a precomputed covariance matrix;
[`vignette("fastkcv")`](https://n-c-taylor.github.io/fastkcv/articles/fastkcv.md)
for timings.

## Examples

``` r
data(meuse, package = "sp")
v <- gstat::vgm(0.59, "Sph", 897, 0.05)
p <- kcv_prepare(log(zinc) ~ sqrt(dist), ~ x + y, meuse, model = v)
p
#> <kcv_prepared> universal kriging, n = 155, p = 2, 187.7 Kb
#> Estimated reciprocal condition number of C: 0.0051
loo <- kcv_run(p)
set.seed(1)
rmse_10fold <- replicate(5, sqrt(mean(kcv_run(p, nfold = 10)$residual^2)))
rmse_10fold
#> [1] 0.3876193 0.3918699 0.3886007 0.3739821 0.3781501
# Buffered leave-one-out at several radii from the same factorization.
sapply(c(0, 200, 400), function(r) kcv_metrics(kcv_run(p, buffer = r))$rmse)
#> [1] 0.3771621 0.4349331 0.4640929

# New responses at the same locations: here 20 random permutations of
# the observations, cross-validated from the same factorization.
set.seed(2)
perm <- replicate(20, sample(log(meuse$zinc)))
cv_perm <- kcv_run(p, nfold = 10, y = perm)
summary(kcv_metrics(cv_perm)$rmse)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>  0.7316  0.8015  0.8146  0.8090  0.8259  0.8700 
```
