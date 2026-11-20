# Peak resident memory (measured with /usr/bin/time -v in a fresh R process
# per configuration; see rss-worker.R) of kcv(), kcv_prepare(),
# kcv_engine() and two competitors, against the n x n matrices each holds.
#
# extra_mb = peak - baseline, the baseline being a process that has loaded
# the packages and the data set (and, for kcv_engine, built the user's
# covariance matrix C). matrices = extra_mb / (8 n^2 bytes): kcv() and
# kcv_prepare() build the covariance in place, so about 1; kcv_engine()
# copies C before factorizing it, so about 1 on top of the user's C.
#
# kcv_engine() is measured up to n = 10,000: at n = 20,000 the user's C and
# the engine's copy (2 x 3.2 GB) reach the default 6.4 GB budget, and on a
# 15 GB machine the run pushed the system into swap. Under the reference
# BLAS the n = 20,000 runs are skipped (each would take about 45 minutes);
# memory does not depend on the BLAS.
# Writes results/rss-<blas>.csv. See common.R for how to run it.

source(file.path("inst", "bench", "common.R"))

ns <- c(1000, 2000, 5000, 10000, 20000)
if (nzchar(Sys.getenv("FASTKCV_BENCH_NS"))) {
  ns <- as.numeric(strsplit(Sys.getenv("FASTKCV_BENCH_NS"), ",")[[1]])
}
max_n_ref <- 10000

# Competitors only where they finish in reasonable time and memory.
configs <- list(
  kcv = c(base = "base"),
  prepare = c(base = "base"),
  engine = c(base = "engine_base", max_n = 10000),
  engine_r = c(base = "base", max_n = 5000),
  krige.cv5 = c(base = "base", max_n = 2000)
)

rows <- list()
for (n in ns) {
  if (blas_tag() == "refblas" && n > max_n_ref) next
  sim_data(n) # make sure the cache exists before measuring
  base <- list()
  for (b in c("base", "engine_base")) {
    message(sprintf("n = %d, %s", n, b))
    base[[b]] <- run_rss(b, n)$peak_mb
  }
  for (cfg in names(configs)) {
    mx <- configs[[cfg]]["max_n"]
    if (!is.na(mx) && n > as.numeric(mx)) next
    message(sprintf("n = %d, %s", n, cfg))
    r <- run_rss(cfg, n)
    b <- base[[configs[[cfg]][["base"]]]]
    rows[[length(rows) + 1]] <- data.frame(
      n = n, config = cfg, time_s = as.numeric(r$time),
      peak_mb = r$peak_mb, baseline_mb = b, extra_mb = r$peak_mb - b,
      matrix_mb = 8 * n^2 / 2^20,
      matrices = (r$peak_mb - b) / (8 * n^2 / 2^20)
    )
  }
}
write_results(do.call(rbind, rows), "rss")
