# Package index

## Cross-validation

[`kcv()`](https://n-c-taylor.github.io/fastkcv/reference/kcv.md) is
[`gstat::krige.cv()`](https://r-spatial.github.io/gstat/reference/krige.cv.html)
from one factorization;
[`kcv_prepare()`](https://n-c-taylor.github.io/fastkcv/reference/kcv_prepare.md)
and
[`kcv_run()`](https://n-c-taylor.github.io/fastkcv/reference/kcv_prepare.md)
reuse that factorization for any number of fold designs and responses.

- [`kcv()`](https://n-c-taylor.github.io/fastkcv/reference/kcv.md) :
  Fast exact cross-validation for global kriging
- [`kcv_prepare()`](https://n-c-taylor.github.io/fastkcv/reference/kcv_prepare.md)
  [`kcv_run()`](https://n-c-taylor.github.io/fastkcv/reference/kcv_prepare.md)
  : Factorize once, cross-validate many times

## Fold designs

Random, spatial block and k-means folds for `nfold` or `folds`.

- [`fold_random()`](https://n-c-taylor.github.io/fastkcv/reference/folds.md)
  [`fold_blocks()`](https://n-c-taylor.github.io/fastkcv/reference/folds.md)
  [`fold_kmeans()`](https://n-c-taylor.github.io/fastkcv/reference/folds.md)
  : Fold builders for cross-validation

## Summaries

- [`kcv_metrics()`](https://n-c-taylor.github.io/fastkcv/reference/kcv_metrics.md)
  : Summary statistics of a cross-validation

## Precomputed covariance matrix

The same computation from any covariance matrix.

- [`kcv_engine()`](https://n-c-taylor.github.io/fastkcv/reference/kcv_engine.md)
  [`kcv_engine_prepare()`](https://n-c-taylor.github.io/fastkcv/reference/kcv_engine.md)
  [`kcv_engine_run()`](https://n-c-taylor.github.io/fastkcv/reference/kcv_engine.md)
  : Exact kriging cross-validation from a covariance matrix
