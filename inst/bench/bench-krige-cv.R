# Benchmark: kcv() (C++ engine) vs the pure-R engine vs gstat::krige.cv(),
# global ordinary kriging, leave-one-out and 5/10/25-fold.
#
# Install the package first so the compiled code is optimised
# (devtools::load_all() compiles without optimisation):
#   R CMD INSTALL .
#   Rscript inst/bench/bench-krige-cv.R
# On Debian/Ubuntu with several BLAS libraries installed, a particular one
# can be selected for a single run through R_LD_LIBRARY_PATH (R prepends it
# to LD_LIBRARY_PATH), e.g.
#   R_LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu/openblas-pthread:/usr/lib/R/lib:/usr/lib/x86_64-linux-gnu \
#     Rscript inst/bench/bench-krige-cv.R
#   R_LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu/blas:/usr/lib/x86_64-linux-gnu/lapack:/usr/lib/R/lib:/usr/lib/x86_64-linux-gnu \
#     Rscript inst/bench/bench-krige-cv.R
# Not run on CRAN (inst/bench is in .Rbuildignore). krige.cv() is run only
# up to n = 1000; at larger n only the two fast engines are compared.
# bench's mem_alloc counts R allocations only, not memory allocated in C++.

library(fastkcv)
library(gstat)

# The same pipeline as kcv(), with the variogramLine() covariance and the
# pure-R engine.
kcv_r <- function(formula, locations, data, model, folds) {
  inp <- fastkcv:::parse_inputs(formula, locations, data)
  C <- fastkcv:::covariance_matrix_r(model, inp$coords)
  res <- fastkcv:::kcv_engine_r(C, inp$F, inp$y, folds)
  fastkcv:::build_output(inp, res)
}

set.seed(2026)
data(meuse, package = "sp")

sim_data <- function(n, model) {
  d <- data.frame(x = runif(n, 0, 1000), y = runif(n, 0, 1000))
  C <- fastkcv:::covariance_matrix(model, as.matrix(d))
  d$z <- drop(crossprod(chol(C), rnorm(n))) + 5
  d
}

v_sim <- vgm(1, "Exp", 200, 0.1)
datasets <- list(
  list(name = "meuse", formula = log(zinc) ~ 1, data = meuse,
    model = vgm(0.5, "Sph", 900, 0.1)),
  list(name = "sim", formula = z ~ 1, data = sim_data(1000, v_sim), model = v_sim),
  list(name = "sim", formula = z ~ 1, data = sim_data(2000, v_sim), model = v_sim),
  list(name = "sim", formula = z ~ 1, data = sim_data(5000, v_sim), model = v_sim)
)

# `expr` is captured unevaluated and re-evaluated on every iteration.
median_s <- function(expr, iterations) {
  e <- substitute(expr)
  env <- parent.frame()
  b <- bench::mark(eval(e, env), iterations = iterations, check = FALSE,
    filter_gc = FALSE)
  list(time = as.numeric(b$median), mem = format(b$mem_alloc))
}

rows <- list()
for (ds in datasets) {
  n <- nrow(ds$data)
  run_gstat <- n <= 1000
  it <- if (n <= 155) 10 else if (n <= 1000) 5 else if (n <= 2000) 3 else 1
  it_gstat <- if (n <= 155) 5 else 1
  for (k in c(n, 25, 10, 5)) {
    folds <- if (k == n) seq_len(n) else sample(rep(seq_len(k), length.out = n))
    label <- if (k == n) "LOO" else paste0(k, "-fold")
    cat(sprintf("n = %d, %s ...\n", n, label))

    a <- median_s(kcv(ds$formula, ~ x + y, ds$data, model = ds$model,
      folds = folds), it)
    b <- median_s(kcv_r(ds$formula, ~ x + y, ds$data, ds$model, folds),
      max(1, it %/% 2))
    g <- list(time = NA_real_)
    max_diff <- NA_real_
    if (run_gstat) {
      g <- median_s(krige.cv(ds$formula, ~ x + y, ds$data, model = ds$model,
        nfold = folds, verbose = FALSE), it_gstat)
      ref <- krige.cv(ds$formula, ~ x + y, ds$data, model = ds$model,
        nfold = folds, verbose = FALSE)
      out <- kcv(ds$formula, ~ x + y, ds$data, model = ds$model, folds = folds)
      max_diff <- max(abs(out$var1.pred - ref$var1.pred),
        abs(out$var1.var - ref$var1.var))
    }
    rows[[length(rows) + 1]] <- data.frame(
      data = ds$name, n = n, folds = label,
      kcv_s = signif(a$time, 3),
      kcv_r_s = signif(b$time, 3),
      krige.cv_s = signif(g$time, 3),
      `R/C++` = round(b$time / a$time, 1),
      `krige.cv/C++` = round(g$time / a$time, 1),
      max_abs_diff = signif(max_diff, 2),
      kcv_R_mem = a$mem,
      check.names = FALSE
    )
  }
}
res <- do.call(rbind, rows)
print(res, row.names = FALSE)

# Several fold schemes on the same data: kcv() for each scheme vs one
# kcv_prepare() followed by kcv_run() for each scheme.
rows <- list()
for (ds in datasets[-(1:2)]) {
  n <- nrow(ds$data)
  schemes <- c(list(seq_len(n)), lapply(c(25, 10, 5), function(k) {
    sample(rep(seq_len(k), length.out = n))
  }))
  cat(sprintf("n = %d, %d fold schemes ...\n", n, length(schemes)))
  t_kcv <- system.time(for (f in schemes) {
    kcv(ds$formula, ~ x + y, ds$data, model = ds$model, folds = f)
  })[["elapsed"]]
  t_prep <- system.time(
    p <- kcv_prepare(ds$formula, ~ x + y, ds$data, model = ds$model)
  )[["elapsed"]]
  t_run <- system.time(for (f in schemes) kcv_run(p, folds = f))[["elapsed"]]
  rows[[length(rows) + 1]] <- data.frame(
    n = n, schemes = "LOO, 25, 10, 5-fold",
    kcv_each_s = signif(t_kcv, 3),
    prepare_s = signif(t_prep, 3),
    runs_s = signif(t_run, 3),
    speedup = round(t_kcv / (t_prep + t_run), 1)
  )
}
cat("\n")
print(do.call(rbind, rows), row.names = FALSE)
cat(
  "\nR", as.character(getRversion()),
  "| gstat", as.character(packageVersion("gstat")),
  "| RcppArmadillo", as.character(packageVersion("RcppArmadillo")),
  "\nBLAS:  ", extSoftVersion()[["BLAS"]],
  "\nLAPACK:", La_library(), "\n"
)
