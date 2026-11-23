# Buffered leave-one-out near the size limit: n = 20,000 (the factor Wt
# alone is 8 n^2 = 3.2 GB of the default fastkcv.max_gb = 6.4 GB budget).
#
# For several radii and target samples, each run is a fresh R process under
# /usr/bin/time -v (rss-worker.R) that reads one saved kcv_prepared object
# and runs the buffered LOO with the "auto", "pairs" or "points" strategy
# (the internal run_buffer_prepared(), which reports the method used).
# extra_mb = peak - (peak of a process that only reads the object), the
# memory of the buffered run itself. pairs is the number of pairs within
# 2 r among the points of the scored hold-out sets (counted here in R), and
# cache_est_mb = 12 pairs bytes, the budget estimate of the pair cache.
#
# The fallback check reruns r = 25 with fastkcv.max_gb set just above
# 8 n^2, so the cache does not fit: "auto" must then pick "points" and give
# bitwise identical results to a forced "points" run.
#
# The prepared object (and the data) are created once, under whichever BLAS
# runs first, and cached in inst/bench/cache/ (uncompressed, 3.2 GB); the
# reference-BLAS run reuses them (factorizing n = 20,000 with the reference
# BLAS takes about 45 minutes) and runs a reduced set of configurations.
# Writes results/buffer-large-<blas>.csv. See common.R for how to run it.

source(file.path("inst", "bench", "common.R"))

# A smoke run (FASTKCV_BENCH_NS) uses its largest n instead.
n <- max(bench_ns(20000))
path <- file.path(cache_dir, sprintf("prepared-%d.rds", n))
d <- sim_data(n)
if (!file.exists(path)) {
  message("preparing n = ", n)
  p <- kcv_prepare(z ~ 1, ~ x + y, d, model = sim_model())
  saveRDS(p, path, compress = FALSE)
  rm(p)
  gc()
}
coords <- as.matrix(d[, c("x", "y")])

# Pairs j <= k within 2 r (1 + 1e-12) among the points of the hold-out sets
# of `targets` that have more than one point (as run_buffer_cpp() caches).
count_pairs <- function(coords, r, targets = NULL) {
  n <- nrow(coords)
  if (is.null(targets)) targets <- seq_len(n)
  ord <- order(coords[, 1])
  xs <- coords[ord, 1]
  near <- function(i, rad) {
    lo <- findInterval(coords[i, 1] - rad, xs, left.open = TRUE) + 1
    hi <- findInterval(coords[i, 1] + rad, xs)
    if (hi < lo) return(integer(0))
    cand <- ord[lo:hi]
    cand[(coords[cand, 1] - coords[i, 1])^2 +
      (coords[cand, 2] - coords[i, 2])^2 <= rad^2]
  }
  in_union <- logical(n)
  for (i in targets) {
    s <- near(i, r)
    if (length(s) > 1) in_union[s] <- TRUE
  }
  u <- which(in_union)
  sub <- coords[u, , drop = FALSE]
  total <- 0
  ord_u <- order(sub[, 1])
  xs_u <- sub[ord_u, 1]
  r2 <- 2 * r * (1 + 1e-12)
  for (a in seq_along(u)) {
    lo <- findInterval(sub[a, 1] - r2, xs_u, left.open = TRUE) + 1
    hi <- findInterval(sub[a, 1] + r2, xs_u)
    cand <- ord_u[lo:hi]
    cand <- cand[u[cand] >= u[a]]
    total <- total + sum((sub[cand, 1] - sub[a, 1])^2 +
      (sub[cand, 2] - sub[a, 2])^2 <= r2^2)
  }
  total
}

targets_of <- function(k) {
  if (is.na(k)) NULL else withr::with_seed(1, sort(sample(n, k)))
}

ref <- blas_tag() == "refblas"
runs <- if (!ref) {
  data.frame(
    radius = c(10, 25, 25, 25, 50, 50, 100, 100, 100, 200),
    method = c("auto", "auto", "pairs", "points", "auto", "points",
      "auto", "pairs", "points", "auto"),
    targets = c(NA, NA, NA, NA, NA, 500, 500, 500, 500, 100)
  )
} else {
  data.frame(
    radius = c(10, 25, 25, 25, 100, 100),
    method = c("auto", "auto", "pairs", "points", "pairs", "points"),
    targets = c(NA, NA, NA, NA, 500, 500)
  )
}

base_mb <- run_rss("buffer_base", n, path = path)$peak_mb
rows <- list()
for (i in seq_len(nrow(runs))) {
  ru <- runs[i, ]
  message(sprintf("r = %g, %s, targets = %s", ru$radius, ru$method, ru$targets))
  extra <- list(path = path, radius = ru$radius, method = ru$method)
  if (!is.na(ru$targets)) extra$targets <- ru$targets
  r <- do.call(run_rss, c(list("buffer", n), extra))
  P <- count_pairs(coords, ru$radius, targets_of(ru$targets))
  rows[[length(rows) + 1]] <- data.frame(
    n = n, radius = ru$radius, targets = ru$targets, method = ru$method,
    method_used = r$method_used, max_gb = 6.4,
    mean_heldout = as.numeric(r$mean_heldout), time_s = as.numeric(r$time),
    peak_mb = r$peak_mb, baseline_mb = base_mb,
    extra_mb = r$peak_mb - base_mb, pairs = P, cache_est_mb = 12 * P / 2^20
  )
}

# Fallback: a budget that leaves beside Wt only half of what the pair cache
# for r = 25 needs (about 10 MB at n = 20,000), so the cache cannot fit at
# any n.
P <- count_pairs(coords, 25)
max_gb <- (8 * n^2 + 6 * P) / 1e9
fb <- tempfile(fileext = ".rds")
pt <- tempfile(fileext = ".rds")
r_fb <- run_rss("buffer", n, path = path, radius = 25, method = "auto",
  max_gb = max_gb, save = fb)
r_pt <- run_rss("buffer", n, path = path, radius = 25, method = "points",
  save = pt)
same <- identical(readRDS(fb), readRDS(pt))
unlink(c(fb, pt))
rows[[length(rows) + 1]] <- data.frame(
  n = n, radius = 25, targets = NA, method = "auto (small budget)",
  method_used = r_fb$method_used, max_gb = max_gb,
  mean_heldout = as.numeric(r_fb$mean_heldout),
  time_s = as.numeric(r_fb$time), peak_mb = r_fb$peak_mb,
  baseline_mb = base_mb, extra_mb = r_fb$peak_mb - base_mb, pairs = P,
  cache_est_mb = 12 * P / 2^20
)
res <- do.call(rbind, rows)
res$fallback_identical_to_points <- ifelse(
  res$method == "auto (small budget)", same, NA
)
write_results(res, "buffer-large")
if (!same || r_fb$method_used != "points") {
  stop("fallback check failed: method ", r_fb$method_used,
    ", identical = ", same)
}
