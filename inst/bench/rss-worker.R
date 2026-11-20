# One configuration for the peak-memory measurements, run in a fresh R
# process under /usr/bin/time -v by bench-rss.R and bench-buffer-large.R
# (bench::mark() only sees R allocations, not the C++ ones).
#
#   Rscript inst/bench/rss-worker.R <config> <n> [key=value ...]
#
# Configs (all start by loading the packages and the cached data set):
#   base            nothing more (baseline for the data-based configs)
#   kcv             kcv() leave-one-out
#   prepare         kcv_prepare() only
#   engine_base     also builds the covariance matrix C (baseline for engine)
#   engine          C, then kcv_engine() leave-one-out
#   engine_r        the pure-R pipeline (variogramLine covariance and
#                   kcv_engine_r), leave-one-out
#   krige.cv5       gstat::krige.cv(nfold = 5)
#   buffer_base     readRDS() of a prepared object (path=...), nothing more
#   buffer          readRDS(), then run_buffer_prepared() with radius=,
#                   method= (auto, pairs, points), targets= (number of
#                   targets drawn with seed 1; all if absent), max_gb=,
#                   save= (an .rds path for the predictions and variances)
# Prints one line "RESULT key=value ..." with the elapsed time of the
# measured step and, for buffer, the method used.

source(file.path("inst", "bench", "common.R"))

args <- commandArgs(trailingOnly = TRUE)
config <- args[1]
n <- as.integer(args[2])
kv <- list()
for (a in args[-(1:2)]) {
  p <- strsplit(a, "=", fixed = TRUE)[[1]]
  kv[[p[1]]] <- p[2]
}
if (!is.null(kv$max_gb)) options(fastkcv.max_gb = as.numeric(kv$max_gb))

out <- list(config = config, n = n)
if (config %in% c("buffer_base", "buffer")) {
  p <- readRDS(kv$path)
  if (config == "buffer") {
    tg <- if (is.null(kv$targets)) NULL else
      withr::with_seed(1, sort(sample(n, as.integer(kv$targets))))
    t <- system.time(res <- fastkcv:::run_buffer_prepared(p, p$inp$y,
      as.numeric(kv$radius), tg, method = kv$method))[["elapsed"]]
    out$time <- t
    out$method_used <- res$method
    out$mean_heldout <- mean(res$n_heldout)
    if (!is.null(kv$save)) saveRDS(res[c("pred", "var", "n_heldout")], kv$save)
  }
} else {
  d <- sim_data(n)
  model <- sim_model()
  t <- switch(config,
    base = 0,
    kcv = system.time(kcv(z ~ 1, ~ x + y, d, model = model))[["elapsed"]],
    prepare = system.time(kcv_prepare(z ~ 1, ~ x + y, d,
      model = model))[["elapsed"]],
    engine_base = ,
    engine = {
      C <- fastkcv:::covariance_matrix(model, as.matrix(d[, c("x", "y")]))
      if (config == "engine") {
        system.time(kcv_engine(C, matrix(1, n, 1), d$z,
          seq_len(n)))[["elapsed"]]
      } else {
        0
      }
    },
    engine_r = system.time(kcv_r(z ~ 1, ~ x + y, d, model,
      seq_len(n)))[["elapsed"]],
    krige.cv5 = system.time(krige.cv(z ~ 1, ~ x + y, d, model = model,
      nfold = withr::with_seed(5, fold_random(~ x + y, d, k = 5)),
      verbose = FALSE))[["elapsed"]],
    stop("unknown config ", config)
  )
  out$time <- t
}
cat("RESULT", paste0(names(out), "=", unlist(out)), "\n")
