# Changelog

## fastkcv 0.1.0

First release.

- [`kcv()`](https://n-c-taylor.github.io/fastkcv/reference/kcv.md):
  exact leave-one-out and k-fold cross-validation for global simple,
  ordinary and universal kriging from one factorization of the
  covariance matrix. Same inputs (`formula`, `locations`, `data`,
  `model`, `nfold`, `beta`) and output as
  [`gstat::krige.cv()`](https://r-spatial.github.io/gstat/reference/krige.cv.html),
  for data.frame, `sf` and `sp` data; variogram models `Nug`, `Exp`,
  `Sph`, `Gau`, `Mat`, `Ste`, `Cir`, `Pen`, `Bes`, `Wav`, `Hol`, nested
  sums and 2D geometric anisotropy.
- [`kcv_prepare()`](https://n-c-taylor.github.io/fastkcv/reference/kcv_prepare.md)
  and
  [`kcv_run()`](https://n-c-taylor.github.io/fastkcv/reference/kcv_prepare.md):
  factorize once, cross-validate any number of fold schemes.
- Buffered leave-one-out (`buffer`), optionally for a subset of
  `targets`.
- Several responses from one factorization: a matrix left-hand side
  (`cbind(z1, z2) ~ 1`, or a matrix column) in
  [`kcv()`](https://n-c-taylor.github.io/fastkcv/reference/kcv.md) and
  [`kcv_prepare()`](https://n-c-taylor.github.io/fastkcv/reference/kcv_prepare.md),
  and new responses for a prepared object with
  `kcv_run(prepared, y = )`, for every fold scheme and buffered
  leave-one-out. The output stacks one block per response with a
  `response` column;
  [`kcv_metrics()`](https://n-c-taylor.github.io/fastkcv/reference/kcv_metrics.md)
  summarizes each response. Simple kriging takes one `beta` for all
  responses or a p x T matrix (also in `kcv_run(beta = )` and
  `kcv_engine_run(beta = )`).
- [`kcv_engine()`](https://n-c-taylor.github.io/fastkcv/reference/kcv_engine.md),
  [`kcv_engine_prepare()`](https://n-c-taylor.github.io/fastkcv/reference/kcv_engine.md)
  and
  [`kcv_engine_run()`](https://n-c-taylor.github.io/fastkcv/reference/kcv_engine.md):
  the same computation from a precomputed covariance matrix, also for a
  matrix-valued response.
- Fold builders
  [`fold_random()`](https://n-c-taylor.github.io/fastkcv/reference/folds.md)
  (the same draw as `krige.cv(nfold = k)`),
  [`fold_blocks()`](https://n-c-taylor.github.io/fastkcv/reference/folds.md)
  and
  [`fold_kmeans()`](https://n-c-taylor.github.io/fastkcv/reference/folds.md).
- [`kcv_metrics()`](https://n-c-taylor.github.io/fastkcv/reference/kcv_metrics.md):
  RMSE, MAE, bias, z-score mean and variance, and prediction-interval
  coverage, overall or per fold.
- Memory requirements are checked up front against
  `getOption("fastkcv.max_gb")`.
- Ordinary and universal kriging estimate the trend in an orthonormal
  basis of the trend matrix, so trends in large raw coordinates (such as
  a quadratic trend in UTM meters) stay accurate.
