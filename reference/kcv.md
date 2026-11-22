# Fast exact cross-validation for global kriging

A drop-in companion to
[`gstat::krige.cv()`](https://r-spatial.github.io/gstat/reference/krige.cv.html)
for global kriging (every observation in the kriging neighborhood).
Instead of refitting the kriging system once per fold, `kcv()`
factorizes the covariance matrix once and obtains every fold's
predictions and variances from it. Results match
[`krige.cv()`](https://r-spatial.github.io/gstat/reference/krige.cv.html)
to numerical precision.

## Usage

``` r
kcv(
  formula,
  locations,
  data,
  model,
  nfold,
  folds = NULL,
  beta = NULL,
  buffer = NULL,
  targets = NULL,
  ...
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

- nfold:

  As in
  [`gstat::krige.cv()`](https://r-spatial.github.io/gstat/reference/krige.cv.html):
  missing or `nrow(data)` for leave-one-out; a single integer `k < n`
  for random k-fold assignment (drawn exactly as
  [`krige.cv()`](https://r-spatial.github.io/gstat/reference/krige.cv.html)
  does, so the same seed gives the same folds); or a vector of length
  `n` of fold identifiers.

- folds:

  Optional vector of length `n` of fold identifiers. Takes precedence
  over `nfold`.

- beta:

  Optional trend coefficients; if given, simple kriging with this known
  trend is used. A vector of length p (the number of trend coefficients)
  is shared by all responses; with several responses, a p x T matrix
  gives one column of coefficients per response.

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

- ...:

  Accepted for compatibility with
  [`krige.cv()`](https://r-spatial.github.io/gstat/reference/krige.cv.html):
  `verbose` and `debug.level` are ignored. Local-neighborhood or block
  arguments (`nmax`, `nmin`, `maxdist`, `block`) are not supported and
  raise an error.

## Value

An object of the same kind as the input (a data.frame with the
coordinate columns appended, an `sf` object, or an `sp`
`SpatialPointsDataFrame`) with columns `var1.pred`, `var1.var`,
`observed`, `residual` (observed minus predicted), `zscore` (residual
divided by the kriging standard error) and `fold`, as in
[`krige.cv()`](https://r-spatial.github.io/gstat/reference/krige.cv.html).
With `buffer`, `fold` is the row number (each observation is its own
hold-out set) and a column `n_heldout` gives the number of observations
held out for it, itself included; with `targets`, only the target rows
are returned, in increasing row order. With a matrix response, the
results of the responses are stacked (all rows of the first response,
then the second, ...) and a factor column `response`, after `fold` (and
`n_heldout`), names the response of each row.

## Details

Ordinary kriging (`formula = z ~ 1`) and universal kriging (`z ~ x + y`,
external drift, ...) re-estimate the trend by generalized least squares
in each fold, as gstat does. Simple kriging is used when `beta` is
supplied.

## Several responses

The kriging weights depend only on the locations, the variogram model,
the trend and the folds, not on the observed values. A matrix response
with T columns is therefore cross-validated from one factorization: each
additional response costs a few matrix products (no further
factorization), and in buffered leave-one-out the hold-out sets and
their factorizations are shared too. Typical uses are many simulated
fields, several variables or time slices observed at the same locations,
or bootstrap replicates.
[`kcv_run()`](https://n-c-taylor.github.io/fastkcv/reference/kcv_prepare.md)
with `y` cross-validates new responses from an existing prepared object.

A response with one column (such as `scale(z) ~ 1`, or a one-column
matrix column) is a single response, as in
[`stats::model.response()`](https://rdrr.io/r/stats/model.extract.html),
unless the left-hand side is an explicit
[`cbind()`](https://rdrr.io/r/base/cbind.html).

Every response is predicted at every location, so a missing or
non-finite value in any column is an error. Response names come from the
column names, else from the arguments of
[`cbind()`](https://rdrr.io/r/base/cbind.html), else they are `"1"`,
`"2"`, .... Each response's rows are the output of the single-response
call, up to rounding in the last bits (with several columns the BLAS
uses matrix-matrix products).
[`kcv_metrics()`](https://n-c-taylor.github.io/fastkcv/reference/kcv_metrics.md)
summarizes each response separately. The output holds about \\10 n\\
numbers per response; several responses are checked against the
`fastkcv.max_gb` budget (see
[`kcv_engine()`](https://n-c-taylor.github.io/fastkcv/reference/kcv_engine.md))
left after the factorization, and the error message gives a batch size
that fits.

## Buffered leave-one-out

With `buffer = r`, observation \\i\\ is held out together with every
observation within Euclidean distance \\r\\ of it (distance \\\le r\\,
computed on the coordinates as given, without any variogram anisotropy),
and only \\i\\ is predicted and reported. This removes the spatially
correlated neighbors that make ordinary leave-one-out optimistic.
`buffer = 0` gives ordinary leave-one-out exactly.

Hold-out sets of neighboring points overlap, so each point is a separate
kriging problem, solved from the same factorization with a hold-out set
\\I\\ of size \\m\\, which needs the block \\Q\_{II}\\ of the matrix
\\Q\\ described in
[`kcv_engine()`](https://n-c-taylor.github.io/fastkcv/reference/kcv_engine.md).
Two strategies are used, whichever needs fewer operations:

- per point: form \\Q\_{II}\\ directly, about \\m^2 n\\ operations per
  point and \\8 m n\\ bytes of temporary memory (160 MB for \\m =
  1000\\, \\n = 20000\\);

- shared pairs: two points in a common hold-out set are at most \\2r\\
  apart, so the entries of \\Q\\ for all pairs of points within \\2r\\
  are computed once (about \\n\\ operations each, roughly \\2 \bar m
  n^2\\ in total) and each \\Q\_{II}\\ is assembled from them. The cache
  takes 12 bytes per pair (about \\24 \bar m n\\ bytes) and is used only
  if it fits, together with the factorization, in the `fastkcv.max_gb`
  budget (see
  [`kcv_engine()`](https://n-c-taylor.github.io/fastkcv/reference/kcv_engine.md)).

Each point also needs a Cholesky factorization of \\Q\_{II}\\ (\\m^3 /
3\\ operations). Since \\m\\ grows with the area of the buffer, large
radii are expensive at large \\n\\; `targets` restricts the predictions
to a sample of points, at a proportional saving. These computations
consist of many small matrix operations, which a multithreaded BLAS does
not speed up and can slow down: with OpenBLAS, setting the environment
variable `OPENBLAS_NUM_THREADS=1` before starting R is usually as fast
or faster. Use
[`kcv_prepare()`](https://n-c-taylor.github.io/fastkcv/reference/kcv_prepare.md)
and
[`kcv_run()`](https://n-c-taylor.github.io/fastkcv/reference/kcv_prepare.md)
to try several radii from one factorization.

An error names the first point that cannot be predicted: a buffer that
holds out every observation, or (universal kriging) one that removes all
observations that identify the trend, for example every observation of
one level of a factor.

## Performance

Almost all of the time goes into one dense Cholesky factorization and
triangular inverse, which R delegates to its BLAS/LAPACK libraries. The
reference BLAS shipped with many R builds is single-threaded and slow;
an optimized BLAS is typically many times faster on these operations.
Check which one R uses with `extSoftVersion()["BLAS"]` or
[`sessionInfo()`](https://rdrr.io/r/utils/sessionInfo.html).

- Linux (Debian/Ubuntu): install an optimized BLAS such as OpenBLAS
  (`sudo apt install libopenblas0-pthread`); the system's alternatives
  mechanism switches R to it without reinstalling anything.

- macOS: the CRAN build of R can use Apple's Accelerate framework; see
  the R for macOS FAQ for how to switch.

- Windows: see "R Installation and Administration" for replacing R's
  BLAS.

Multithreaded BLAS libraries use all cores by default; when running
several R processes in parallel, limit the threads (for OpenBLAS, the
environment variable `OPENBLAS_NUM_THREADS`). OpenBLAS starts one thread
per logical CPU, and on CPUs with simultaneous multithreading
(hyper-threading) one thread per physical core can be much faster: on a
4-core, 8-thread laptop, `OPENBLAS_NUM_THREADS=4` made 10-fold `kcv()`
1.6 to 2 times faster than the default at n = 1,500 to 5,000.

The factorization depends only on the locations, the model and the
trend, not on the folds. To cross-validate several fold schemes for the
same data, call
[`kcv_prepare()`](https://n-c-taylor.github.io/fastkcv/reference/kcv_prepare.md)
once and
[`kcv_run()`](https://n-c-taylor.github.io/fastkcv/reference/kcv_prepare.md)
for each scheme. Each fold of size \\m\\ adds about \\m^2 n\\
operations, so a few large folds cost the most: with 2 folds, the folds
take about as long as the factorization.

See
[`kcv_engine()`](https://n-c-taylor.github.io/fastkcv/reference/kcv_engine.md)
for memory requirements.

## Examples

``` r
library(sp)
data(meuse, package = "sp")
v <- gstat::vgm(0.59, "Sph", 897, 0.05)
cv <- kcv(log(zinc) ~ 1, ~ x + y, meuse, model = v)
head(cv)
#>   var1.pred  var1.var observed    residual      zscore fold      x      y
#> 1  6.769182 0.1800190 6.929517  0.16033461  0.37789233    1 181072 333611
#> 2  6.767296 0.1747339 7.039660  0.27236448  0.65157117    2 181025 333558
#> 3  6.296517 0.1818894 6.461468  0.16495146  0.38676967    3 181165 333537
#> 4  6.039217 0.2281430 5.549076 -0.49014098 -1.02616546    4 181298 333484
#> 5  5.581987 0.1743017 5.594711  0.01272471  0.03047874    5 181307 333330
#> 6  5.425952 0.2311933 5.638355  0.21240263  0.44174562    6 181390 333260
sqrt(mean(cv$residual^2))
#> [1] 0.3917495

# Buffered leave-one-out: neighbors within 300 m are held out too.
cv_buf <- kcv(log(zinc) ~ 1, ~ x + y, meuse, model = v, buffer = 300)
summary(cv_buf$n_heldout)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>   1.000   7.000  10.000   9.981  13.000  19.000 
kcv_metrics(cv_buf)
#>     n      rmse       mae       bias     mean_z     var_z coverage_50
#> 1 155 0.5313131 0.4242791 0.04478428 0.06573575 0.7752399    0.516129
#>   coverage_80 coverage_90 coverage_95
#> 1   0.8903226   0.9290323   0.9612903

# The same statistics from a sample of 50 points.
set.seed(1)
kcv_metrics(kcv(log(zinc) ~ 1, ~ x + y, meuse, model = v, buffer = 300,
  targets = sample(nrow(meuse), 50)))
#>    n     rmse       mae         bias      mean_z     var_z coverage_50
#> 1 50 0.482048 0.3723308 -0.005066683 -0.01554473 0.6595127        0.58
#>   coverage_80 coverage_90 coverage_95
#> 1        0.92        0.96        0.98

# Several responses with the same model, trend and folds.
cv3 <- kcv(cbind(log(zinc), log(copper), log(lead)) ~ 1, ~ x + y, meuse,
  model = v, nfold = 5)
kcv_metrics(cv3)
#>      response   n      rmse       mae        bias      mean_z     var_z
#> 1   log(zinc) 155 0.4116874 0.3111180 0.004152948 0.008842182 0.8627269
#> 2 log(copper) 155 0.3143176 0.2505060 0.002310917 0.004644471 0.5196108
#> 3   log(lead) 155 0.4108739 0.3031579 0.005709850 0.011264441 0.8696293
#>   coverage_50 coverage_80 coverage_90 coverage_95
#> 1   0.5741935   0.8193548   0.9290323   0.9677419
#> 2   0.6258065   0.9354839   0.9677419   1.0000000
#> 3   0.6258065   0.8451613   0.9032258   0.9612903
```
