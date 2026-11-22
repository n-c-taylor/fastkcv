<!-- Site home page (pkgdown). The sections from "Example" on are copied
     from README.md; keep the two in step. -->

```{=html}
<div class="hero">
  <div class="hero-text">
    <span class="eyebrow">R package &middot; companion to gstat</span>
    <h1>fastkcv</h1>
    <p class="hero-tagline">Fast, exact cross-validation for <span class="accent">global kriging</span></p>
    <p class="hero-summary"><code>kcv()</code> takes the same inputs as <code>gstat::krige.cv()</code> and returns the same output. It factorizes the covariance matrix of the observations once, instead of refitting the kriging system for every fold, and the results agree with <code>krige.cv()</code> to numerical precision.</p>
    <div class="hero-actions">
      <a class="btn-accent" href="articles/fastkcv.html"><span class="fa fa-rocket"></span> Get started</a>
      <a class="btn-ghost" href="reference/index.html"><span class="fa fa-book"></span> Reference</a>
      <a class="btn-ghost" href="https://github.com/n-c-taylor/fastkcv"><span class="fab fa-github"></span> GitHub</a>
    </div>
    <div class="hero-links">
      <a href="https://github.com/n-c-taylor/fastkcv/actions/workflows/R-CMD-check.yaml"><img src="https://github.com/n-c-taylor/fastkcv/actions/workflows/R-CMD-check.yaml/badge.svg" alt="R-CMD-check"></a>
      <span class="pill pill-neutral">MIT license</span>
      <span class="pill">R / C++</span>
    </div>
  </div>
  <div class="hero-photo"><img src="reference/figures/logo.png" alt="fastkcv logo"></div>
</div>

<div class="stat-strip">
  <div class="stat"><span class="stat-num">1,500x</span><span class="stat-label">faster leave-one-out than <code>krige.cv()</code> at n = 2,000</span></div>
  <div class="stat"><span class="stat-num">1</span><span class="stat-label">factorization for every fold, no refit per fold</span></div>
  <div class="stat"><span class="stat-num">11</span><span class="stat-label">gstat variogram models, nested and anisotropic</span></div>
  <div class="stat"><span class="stat-num">20,000</span><span class="stat-label">observations in about two minutes</span></div>
</div>

<div class="card-grid">
  <div class="feature-card"><span class="icon-box"><span class="fa fa-right-left"></span></span><span class="card-title">A drop-in for krige.cv()</span><span class="card-text">Same formula, locations, data, model, nfold and beta; the same output columns for data.frame, sf and sp data.</span></div>
  <div class="feature-card"><span class="icon-box"><span class="fa fa-bolt"></span></span><span class="card-title">Factorize once</span><span class="card-text"><code>kcv_prepare()</code> and <code>kcv_run()</code> reuse one factorization for any number of fold designs.</span></div>
  <div class="feature-card"><span class="icon-box"><span class="fa fa-table-cells"></span></span><span class="card-title">Spatial fold designs</span><span class="card-text">Random, block and k-means folds, and <code>kcv_metrics()</code> to compare them: RMSE, bias, z-scores, coverage.</span></div>
  <div class="feature-card"><span class="icon-box"><span class="fa fa-bullseye"></span></span><span class="card-title">Buffered leave-one-out</span><span class="card-text">Hold out every observation within a radius of each point, for all points or a subset of targets.</span></div>
  <div class="feature-card"><span class="icon-box"><span class="fa fa-layer-group"></span></span><span class="card-title">Several responses</span><span class="card-text"><code>cbind(z1, z2) ~ 1</code> or <code>kcv_run(prepared, y = )</code>: many responses at the same locations from one factorization.</span></div>
  <div class="feature-card"><span class="icon-box"><span class="fa fa-table"></span></span><span class="card-title">Any covariance matrix</span><span class="card-text"><code>kcv_engine()</code> runs the same exact computation from a precomputed covariance matrix.</span></div>
</div>
```

## Installation

Using `pak`:

```{=html}
<div class="install-box"><code>pak::pkg_install("n-c-taylor/fastkcv")</code><button type="button" class="copy-btn" data-copy='pak::pkg_install("n-c-taylor/fastkcv")' aria-label="Copy install command"><span class="fa fa-copy"></span><span>Copy</span></button></div>
```

Using `remotes`:
```r
# install.packages("remotes")
remotes::install_github("n-c-taylor/fastkcv", build_vignettes = TRUE)
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

The vignette ([Get started](articles/fastkcv.html)) covers the fold designs, buffered
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
