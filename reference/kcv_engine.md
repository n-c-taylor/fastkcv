# Exact kriging cross-validation from a covariance matrix

Low-level engine behind
[`kcv()`](https://n-c-taylor.github.io/fastkcv/reference/kcv.md),
implemented in C++. Given the covariance matrix `C` of the observations,
the trend matrix `F` and fold assignments, it computes cross-validation
predictions and variances for every fold from a single factorization of
`C`:

## Usage

``` r
kcv_engine(C, F, y, folds, beta = NULL)

kcv_engine_prepare(C, F, beta = NULL)

kcv_engine_run(prepared, y, folds, beta = NULL)
```

## Arguments

- C:

  n x n covariance matrix of the observations (symmetric positive
  definite; an asymmetric matrix or a non-finite entry is an error).

- F:

  n x p trend matrix (a column of ones for ordinary kriging).

- y:

  Response: a numeric vector or an n x T matrix. Each column is
  cross-validated with the same kriging weights.

- folds:

  Vector of length n of fold identifiers.

- beta:

  Optional trend coefficients for simple kriging: a length-p vector
  shared by every column of `y`, or a p x T matrix with one column per
  column of `y`. In `kcv_engine_run()`, replaces the coefficients given
  to `kcv_engine_prepare()`, which must have been made with `beta` (the
  factorization of simple kriging does not depend on `beta`).

- prepared:

  A `kcv_prepared` object from `kcv_engine_prepare()` or
  [`kcv_prepare()`](https://n-c-taylor.github.io/fastkcv/reference/kcv_prepare.md).

## Value

`kcv_engine()` and `kcv_engine_run()` return a list with `pred` (n x T
matrix), `var` (length-n vector; the same for every column of `y`) and
`folds`. `kcv_engine_prepare()` returns a `kcv_prepared` object holding
the factorization, to be passed to `kcv_engine_run()` any number of
times with different `y` or `folds`.

## Details

- ordinary/universal kriging (`beta = NULL`, trend re-estimated by GLS
  in each fold): \\Q = C^{-1} - C^{-1}F(F'C^{-1}F)^{-1}F'C^{-1}\\;

- simple kriging (known `beta`): \\Q = C^{-1}\\, applied to \\y -
  F\beta\\.

For a fold \\I\\, the residual is \\e_I = Q\_{II}^{-1}(Qy)\_I\\, the
prediction \\y_I - e_I\\, and the prediction-error variance
\\\mathrm{diag}(Q\_{II}^{-1})\\.

Only the diagonal blocks \\Q\_{II}\\ are formed, never the full \\n
\times n\\ matrix \\Q\\. Leave-one-out (every fold of size one) uses
\\e_i = (Qy)\_i / Q\_{ii}\\ and \\1 / Q\_{ii}\\ directly.

## Memory

Factorizing `C` needs one \\n \times n\\ workspace in addition to `C`,
so the peak memory of `kcv_engine()` is about \\16 n^2\\ bytes. (The
model-based
[`kcv()`](https://n-c-taylor.github.io/fastkcv/reference/kcv.md) builds
the covariance directly in the workspace and needs half of that.)
Because operating systems often over-commit memory, an allocation that
is too large may not fail cleanly, so the estimate is checked up front
against a budget of `getOption("fastkcv.max_gb", 6.4)` GB, which allows
\\n\\ up to about 20,000 here and 28,000 in
[`kcv()`](https://n-c-taylor.github.io/fastkcv/reference/kcv.md). Raise
it with `options(fastkcv.max_gb = )` on machines with enough memory. A
prepared object keeps its \\n \times n\\ factor (\\8 n^2\\ bytes) for as
long as it exists.

## See also

[`kcv_prepare()`](https://n-c-taylor.github.io/fastkcv/reference/kcv_prepare.md)
and
[`kcv_run()`](https://n-c-taylor.github.io/fastkcv/reference/kcv_prepare.md)
for the model-based equivalent.

## Examples

``` r
set.seed(1)
xy <- matrix(runif(100), ncol = 2)
C <- exp(-as.matrix(dist(xy)) / 0.3)
y <- rnorm(50)
res <- kcv_engine(C, matrix(1, 50, 1), y, folds = seq_len(50))
head(res$pred)
#>             [,1]
#> [1,] -1.25265287
#> [2,]  0.29321589
#> [3,]  0.43726382
#> [4,]  0.01225811
#> [5,] -0.04605711
#> [6,] -0.27491546

# Factorize once, then cross-validate several fold schemes.
p <- kcv_engine_prepare(C, matrix(1, 50, 1))
loo <- kcv_engine_run(p, y, seq_len(50))
k5 <- kcv_engine_run(p, y, rep(1:5, 10))
```
