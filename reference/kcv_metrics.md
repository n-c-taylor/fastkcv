# Summary statistics of a cross-validation

Summarizes the output of
[`kcv()`](https://n-c-taylor.github.io/fastkcv/reference/kcv.md) or
[`kcv_run()`](https://n-c-taylor.github.io/fastkcv/reference/kcv_prepare.md)
(or of
[`gstat::krige.cv()`](https://r-spatial.github.io/gstat/reference/krige.cv.html),
which has the same columns), overall or per fold.

## Usage

``` r
kcv_metrics(
  x,
  by_fold = FALSE,
  level = c(0.5, 0.8, 0.9, 0.95),
  by_response = TRUE
)
```

## Arguments

- x:

  Output of
  [`kcv()`](https://n-c-taylor.github.io/fastkcv/reference/kcv.md) or
  [`kcv_run()`](https://n-c-taylor.github.io/fastkcv/reference/kcv_prepare.md):
  a data.frame, an `sf` object or an `sp` `SpatialPointsDataFrame`, with
  columns `residual` and `var1.var` (and `fold` when `by_fold = TRUE`).

- by_fold:

  If `TRUE`, one row per fold (in `sort(unique(fold))` order) instead of
  one overall row.

- level:

  Nominal coverage levels of the prediction intervals, in (0, 1).

- by_response:

  If `TRUE` (the default) and `x` has a `response` column (the output of
  a cross-validation of several responses), one row per response (in the
  order of its levels), or per response and fold with `by_fold = TRUE`.
  `FALSE` pools the responses. Ignored without a `response` column.

## Value

A data.frame with columns `n`, `rmse`, `mae`, `bias`, `mean_z`, `var_z`
and one `coverage_*` column per level, preceded by `fold` when
`by_fold = TRUE`, and by `response` for several responses.

## Details

With residuals \\e_i\\ (observed minus predicted), kriging variances
\\\sigma_i^2\\ and z-scores \\z_i = e_i / \sigma_i\\:

- `rmse`: \\\sqrt{\mathrm{mean}(e_i^2)}\\;

- `mae`: \\\mathrm{mean}(\|e_i\|)\\;

- `bias`: \\\mathrm{mean}(e_i)\\ (positive when the predictions are too
  low on average);

- `mean_z`, `var_z`: mean and variance of the z-scores, ideally close to
  0 and 1. `var_z` uses the \\n - 1\\ denominator of
  [`stats::var()`](https://rdrr.io/r/stats/cor.html), so it is `NA` for
  a fold with a single observation;

- `coverage_<100 L>` for each `level` \\L\\: the proportion of
  observations inside the Gaussian prediction interval \\\hat y_i \pm q
  \sigma_i\\, \\q\\ = `qnorm((1 + L) / 2)`, i.e. with \\\|z_i\| \le q\\.
  Well-calibrated variances give coverage close to \\L\\. The default
  levels 50, 80, 90 and 95% probe both the center and the tails of the
  predictive distribution.

The z-scores are recomputed from the `residual` and `var1.var` columns.

## See also

[`kcv()`](https://n-c-taylor.github.io/fastkcv/reference/kcv.md),
[`fold_blocks()`](https://n-c-taylor.github.io/fastkcv/reference/folds.md)
and the other fold builders.

## Examples

``` r
data(meuse, package = "sp")
v <- gstat::vgm(0.59, "Sph", 897, 0.05)
cv <- kcv(log(zinc) ~ 1, ~ x + y, meuse, model = v)
kcv_metrics(cv)
#>     n      rmse       mae          bias       mean_z     var_z coverage_50
#> 1 155 0.3917495 0.2921011 -1.256051e-05 0.0001815253 0.8281059   0.6064516
#>   coverage_80 coverage_90 coverage_95
#> 1   0.8258065   0.9225806   0.9677419

cv5 <- kcv(log(zinc) ~ 1, ~ x + y, meuse, model = v,
  folds = fold_blocks(~ x + y, meuse, blocks = 2))
kcv_metrics(cv5, by_fold = TRUE)
#>   fold  n      rmse       mae          bias     mean_z     var_z coverage_50
#> 1    1 66 0.6304688 0.5340387 -0.1237597885 -0.1662439 0.6162833   0.5606061
#> 2    2 12 0.6730908 0.5421966 -0.0009212365 -0.0179022 1.0958447   0.4166667
#> 3    3 12 0.6107563 0.5330891  0.1781119209  0.2709644 0.9429501   0.2500000
#> 4    4 65 0.6439722 0.5533069 -0.2784011327 -0.3697079 0.5745268   0.4769231
#>   coverage_80 coverage_90 coverage_95
#> 1   0.9090909   0.9848485   1.0000000
#> 2   0.8333333   0.9166667   0.9166667
#> 3   0.8333333   1.0000000   1.0000000
#> 4   0.8923077   0.9692308   1.0000000
```
