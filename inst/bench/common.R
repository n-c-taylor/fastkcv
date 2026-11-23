# Shared helpers for the benchmark scripts in inst/bench/.
#
# Every script is run from the package root, against an optimized install:
#   rm -f src/*.o src/*.so && R CMD INSTALL --preclean .
#   Rscript inst/bench/bench-folds.R
# (devtools::load_all() compiles without optimization, and an in-place
# install silently reuses those object files.)
#
# BLAS: on Debian/Ubuntu with several BLAS libraries installed, one is
# selected for a single run through R_LD_LIBRARY_PATH (R puts it before
# LD_LIBRARY_PATH; plain LD_LIBRARY_PATH has no effect):
#   OpenBLAS:
#   R_LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu/openblas-pthread:/usr/lib/R/lib:/usr/lib/x86_64-linux-gnu
#   reference BLAS:
#   R_LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu/blas:/usr/lib/x86_64-linux-gnu/lapack:/usr/lib/R/lib:/usr/lib/x86_64-linux-gnu
# inst/bench/run-all.sh runs every script under both.
#
# Results go to inst/bench/results/<script>-<blas>.csv (and .rds);
# inst/bench/summarize.R turns them into the summary tables. Simulated data
# sets are cached in inst/bench/cache/ (not tracked), so both BLAS runs use
# identical data.

suppressPackageStartupMessages({
  library(fastkcv)
  library(gstat)
})

bench_dir <- file.path("inst", "bench")
results_dir <- file.path(bench_dir, "results")
cache_dir <- file.path(bench_dir, "cache")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)

# Sample sizes of a script: `default`, or FASTKCV_BENCH_NS (comma separated)
# for a smoke run, e.g. FASTKCV_BENCH_NS="250,500" FASTKCV_BENCH_CAP=5. A
# script whose cost grows too fast for those sizes passes its own `smoke`.
bench_ns <- function(default, smoke = NULL) {
  env <- Sys.getenv("FASTKCV_BENCH_NS")
  if (!nzchar(env)) return(default)
  if (!is.null(smoke)) smoke else as.numeric(strsplit(env, ",")[[1]])
}

# Competitors predicted to take longer than this many seconds are skipped.
bench_cap <- function() as.numeric(Sys.getenv("FASTKCV_BENCH_CAP", "600"))

# "openblas" or "refblas", from the BLAS library actually loaded.
blas_tag <- function() {
  lib <- paste(extSoftVersion()[["BLAS"]], La_library())
  if (grepl("openblas", lib, ignore.case = TRUE)) "openblas" else "refblas"
}

session_stamp <- function() {
  data.frame(
    r = as.character(getRversion()),
    fastkcv = as.character(packageVersion("fastkcv")),
    gstat = as.character(packageVersion("gstat")),
    blas = blas_tag(),
    blas_lib = basename(extSoftVersion()[["BLAS"]]),
    lapack_lib = basename(La_library()),
    cpus = parallel::detectCores(),
    date = format(Sys.Date())
  )
}

# The simulation model: exponential covariance, range 200 on [0, 1000]^2,
# sill 1, nugget 0.1. The field is drawn exactly from this covariance
# (Cholesky), so the data are a realization of the model used for kriging.
sim_model <- function() vgm(1, "Exp", 200, 0.1)

# Simulated data set with n points, cached on disk (generating it needs a
# Cholesky factor of the n x n covariance, slow with the reference BLAS).
sim_data <- function(n, seed = 1000 + n, model = sim_model()) {
  path <- file.path(cache_dir, sprintf("sim-%d-%d.rds", n, seed))
  if (file.exists(path)) return(readRDS(path))
  d <- withr::with_seed(seed, {
    d <- data.frame(x = stats::runif(n, 0, 1000), y = stats::runif(n, 0, 1000))
    C <- fastkcv:::covariance_matrix(model, as.matrix(d))
    d$z <- drop(crossprod(chol(C), stats::rnorm(n))) + 5
    d
  })
  saveRDS(d, path)
  d
}

meuse_data <- function() {
  e <- new.env()
  utils::data("meuse", package = "sp", envir = e)
  e$meuse
}

# The pure-R pipelines that the competitors stand for: the covariance from
# gstat::variogramLine() and either the explicit-Q engine (kcv_engine_r) or
# the per-fold refit with solve() (kcv_reference).
kcv_r <- function(formula, locations, data, model, folds,
                  engine = c("engine_r", "reference")) {
  engine <- match.arg(engine)
  inp <- fastkcv:::parse_inputs(formula, locations, data)
  C <- fastkcv:::covariance_matrix_r(model, inp$coords)
  f <- switch(engine,
    engine_r = fastkcv:::kcv_engine_r,
    reference = fastkcv:::kcv_reference
  )
  fastkcv:::build_output(inp, f(C, inp$F, inp$y, folds))
}

# Times an expression (captured unevaluated, re-evaluated per iteration).
# One bench::mark() iteration first; if that is quick, more iterations
# (up to `max_iter`, at least `min_time` seconds in total) and the median.
# Returns seconds, R-side allocation in bytes (bench's mem_alloc: memory
# allocated in C++ is not counted) and the number of iterations.
time_expr <- function(expr, env = parent.frame(), quick = 2, min_time = 2,
                      max_iter = 20) {
  force(env) # before bench::mark() changes what parent.frame() means
  e <- substitute(expr)
  b <- bench::mark(eval(e, env), iterations = 1, check = FALSE,
    filter_gc = FALSE, memory = TRUE)
  t1 <- as.numeric(b$median)
  if (t1 < quick) {
    b <- bench::mark(eval(e, env), min_time = min_time,
      min_iterations = 3, max_iterations = max_iter, check = FALSE,
      filter_gc = FALSE, memory = TRUE)
  }
  list(
    time = as.numeric(b$median),
    mem = as.numeric(b$mem_alloc),
    iter = as.integer(b$n_itr)
  )
}

# Predicted run time at n from earlier measurements of the same series
# (`hist`, a list of list(n, time)), t ~ n^power. With two or more points
# the exponent is estimated from the last two, clamped to [2, max_power]
# (fixed overheads make small-n growth look slower than it is); with one
# point max_power is used. NA when there is no measurement.
predict_time <- function(hist, n, max_power) {
  k <- length(hist)
  if (k == 0) return(NA_real_)
  b <- hist[[k]]
  power <- max_power
  if (k >= 2) {
    a <- hist[[k - 1]]
    power <- min(max_power, max(2, log(b$time / a$time) / log(b$n / a$n)))
  }
  b$time * (n / b$n)^power
}

write_results <- function(df, name) {
  tag <- blas_tag()
  df$blas <- tag
  base <- file.path(results_dir, sprintf("%s-%s", name, tag))
  utils::write.csv(df, paste0(base, ".csv"), row.names = FALSE)
  saveRDS(list(results = df, session = session_stamp()), paste0(base, ".rds"))
  message("wrote ", base, ".csv")
  invisible(df)
}

max_abs_diff <- function(a, b) {
  max(abs(a$var1.pred - b$var1.pred), abs(a$var1.var - b$var1.var))
}

# Runs rss-worker.R in a fresh R process under /usr/bin/time -v (the child
# inherits this process's environment, so the same BLAS) and returns the
# worker's RESULT fields plus the peak resident set size in MB.
run_rss <- function(config, n, ...) {
  extra <- list(...)
  extra <- if (length(extra)) paste0(names(extra), "=", unlist(extra)) else NULL
  err <- tempfile()
  out <- system2("/usr/bin/time",
    c("-v", file.path(R.home("bin"), "Rscript"),
      file.path(bench_dir, "rss-worker.R"), config, n, extra),
    stdout = TRUE, stderr = err)
  e <- readLines(err)
  unlink(err)
  rss <- grep("Maximum resident set size", e, value = TRUE)
  res <- grep("^RESULT", out, value = TRUE)
  if (length(rss) != 1 || length(res) != 1) {
    stop("worker failed (", config, ", n = ", n, "):\n",
      paste(utils::tail(c(out, e), 15), collapse = "\n"))
  }
  f <- strsplit(sub("^RESULT ", "", trimws(res)), " ")[[1]]
  vals <- stats::setNames(sub("^[^=]*=", "", f), sub("=.*$", "", f))
  c(as.list(vals), list(
    peak_mb = as.numeric(sub(".*: *", "", rss)) / 1024
  ))
}
