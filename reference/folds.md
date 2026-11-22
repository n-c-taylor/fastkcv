# Fold builders for cross-validation

Assign observations to cross-validation folds. Each function returns a
plain integer vector with one fold identifier per observation, to be
passed as `folds` to
[`kcv()`](https://n-c-taylor.github.io/fastkcv/reference/kcv.md) or
[`kcv_run()`](https://n-c-taylor.github.io/fastkcv/reference/kcv_prepare.md)
(or as `nfold` to
[`gstat::krige.cv()`](https://r-spatial.github.io/gstat/reference/krige.cv.html)).

## Usage

``` r
fold_random(locations, data = NULL, k)

fold_blocks(locations, data = NULL, blocks = NULL, size = NULL, k = NULL)

fold_kmeans(locations, data = NULL, k, nstart = 10, iter.max = 100, ...)
```

## Arguments

- locations:

  A numeric two-column coordinate matrix, a formula giving the
  coordinates (e.g. `~x + y`) with a data.frame in `data`, an `sf`
  object with POINT geometries, or an `sp` `Spatial*` object.

- data:

  A data.frame, when `locations` is a formula.

- k:

  Number of folds. For `fold_blocks()`, optional: the number of folds
  into which the blocks are grouped at random.

- blocks:

  Number of grid cells along x and y: a single integer (used for both)
  or a vector `c(nx, ny)`. The cells divide the bounding box of the
  locations into equal parts.

- size:

  Cell side length along x and y, in coordinate units: a single number
  (square cells) or `c(dx, dy)`. The grid starts at the minimum corner
  of the bounding box.

- nstart, iter.max:

  Passed to [`stats::kmeans()`](https://rdrr.io/r/stats/kmeans.html).

- ...:

  Further arguments to
  [`stats::kmeans()`](https://rdrr.io/r/stats/kmeans.html), such as
  `algorithm`.

## Value

An integer vector of fold identifiers, one per observation.

## Details

- `fold_random()` draws folds exactly as
  [`gstat::krige.cv()`](https://r-spatial.github.io/gstat/reference/krige.cv.html)
  does for a scalar `nfold`: `sample(k, n, replace = TRUE)` when
  `k < n`, so the same seed gives the same folds in both packages. Fold
  sizes are random and a fold may be empty (then it is simply absent).
  With `k >= n` it returns `1:n` (leave-one-out).

- `fold_blocks()` lays a regular grid over the bounding box of the
  locations and makes each non-empty grid cell a fold. Specify the grid
  either by the number of cells (`blocks`) or by the cell size (`size`).
  Empty cells are dropped and the remaining cells are numbered
  consecutively (x fastest, then y, starting at the minimum corner), so
  the number of folds can be smaller than the number of cells. With `k`,
  the blocks are grouped at random into `k` folds, each receiving the
  same number of blocks (up to one).

- `fold_kmeans()` clusters the coordinates with
  [`stats::kmeans()`](https://rdrr.io/r/stats/kmeans.html) and makes
  each cluster a fold. The result depends on the random starts, so set a
  seed for reproducibility. On large, strongly clustered data the
  default Hartigan-Wong algorithm can reach its internal step limit
  ("Quick-TRANSfer stage steps exceeded"); the clustering is then still
  a valid set of folds, possibly not a locally optimal one, so this
  warning is not shown unless `algorithm` is given explicitly.

Coordinates are used as given (projected units); all three functions
reject geographic (longitude/latitude) coordinates, as
[`kcv()`](https://n-c-taylor.github.io/fastkcv/reference/kcv.md) does.

## See also

[`kcv()`](https://n-c-taylor.github.io/fastkcv/reference/kcv.md),
[`kcv_run()`](https://n-c-taylor.github.io/fastkcv/reference/kcv_prepare.md),
[`kcv_metrics()`](https://n-c-taylor.github.io/fastkcv/reference/kcv_metrics.md);
for buffered leave-one-out, the `buffer` argument of
[`kcv()`](https://n-c-taylor.github.io/fastkcv/reference/kcv.md).

## Examples

``` r
data(meuse, package = "sp")
v <- gstat::vgm(0.59, "Sph", 897, 0.05)
set.seed(1)
f_rand <- fold_random(~ x + y, meuse, k = 10)
f_block <- fold_blocks(~ x + y, meuse, blocks = 4)
f_km <- fold_kmeans(~ x + y, meuse, k = 10)
table(f_block)
#> f_block
#>  1  2  3  4  5  6  7  8  9 10 11 
#> 23  9  3 14 20  9 12 25  6  3 31 

p <- kcv_prepare(log(zinc) ~ 1, ~ x + y, meuse, model = v)
rbind(
  random = kcv_metrics(kcv_run(p, folds = f_rand)),
  blocks = kcv_metrics(kcv_run(p, folds = f_block)),
  kmeans = kcv_metrics(kcv_run(p, folds = f_km))
)
#>          n      rmse       mae         bias      mean_z     var_z coverage_50
#> random 155 0.4095691 0.3047347  0.013372260  0.02630708 0.8672858   0.6000000
#> blocks 155 0.5346710 0.4262636 -0.097951373 -0.15705485 0.7134334   0.5548387
#> kmeans 155 0.4845936 0.3866041 -0.008659043 -0.02156850 0.5862294   0.5935484
#>        coverage_80 coverage_90 coverage_95
#> random   0.8129032   0.9096774   0.9483871
#> blocks   0.8774194   0.9419355   0.9806452
#> kmeans   0.9225806   0.9612903   0.9870968
```
