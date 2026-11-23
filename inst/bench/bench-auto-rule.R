# The automatic choice between the two ways of forming the buffered
# leave-one-out blocks Q[I, I] (run_buffer_cpp(), method "auto"):
#   points: q_block() per target, cost_points = sum_t m_t^2 (n - min I_t);
#   pairs:  cache Q[j, k] for all pairs within 2 r, cost_pairs = sum (n - k);
# auto picks pairs when cost_pairs < cost_points (and the cache fits in
# the memory budget). The two costs count flops, not BLAS efficiency: the
# points path runs as matrix products, the pairs path as dot products.
#
# Grid: n in {1000, 2000, 5000, 10000} and radii giving a mean hold-out
# size m of about 5, 20, 50, 150 and 500 on [0, 1000]^2 for uniform points
# (sim_data()). A clustered design (20 clusters of normal spread 30 around
# uniform centers, plus 20% uniform background points; seed 2000 + n) uses
# the same radii: there m differs by orders of magnitude between points in
# a cluster and background points. Runs with n m > 1e6 use 500 random
# targets (seed 1) instead of all points. Both strategies are
# forced and timed; the costs are recomputed here in R exactly as in C++,
# and auto's choice is cost_pairs < cost_points (checked against the C++
# choice on the cheaper runs).
#
# Writes results/auto-rule-<blas>.csv; summarize.R reports where auto picks
# the slower strategy and by how much, per BLAS.
# See common.R for how to run it.

source(file.path("inst", "bench", "common.R"))

ns <- bench_ns(c(1000, 2000, 5000, 10000))
ms <- c(5, 20, 50, 150, 500)

# Clustered locations; the response does not affect the timings.
clustered_data <- function(n, seed = 2000 + n) {
  withr::with_seed(seed, {
    n_bg <- round(0.2 * n)
    centers <- matrix(stats::runif(40, 100, 900), ncol = 2)
    id <- sample(20, n - n_bg, replace = TRUE)
    xy <- rbind(centers[id, ] + matrix(stats::rnorm(2 * (n - n_bg), 0, 30),
      ncol = 2), matrix(stats::runif(2 * n_bg, 0, 1000), ncol = 2))
    data.frame(x = xy[, 1], y = xy[, 2], z = stats::rnorm(n))
  })
}

# Hold-out sets (sorted), as in run_buffer_cpp().
sets_of <- function(coords, r, targets) {
  ord <- order(coords[, 1])
  xs <- coords[ord, 1]
  lapply(targets, function(i) {
    lo <- findInterval(coords[i, 1] - r, xs, left.open = TRUE) + 1
    hi <- findInterval(coords[i, 1] + r, xs)
    cand <- ord[seq_len(max(0, hi - lo + 1)) + lo - 1]
    sort(cand[(coords[cand, 1] - coords[i, 1])^2 +
      (coords[cand, 2] - coords[i, 2])^2 <= r^2])
  })
}

# cost_points, cost_pairs and the number of cached pairs, as in C++
# (0-based indices there: n - k0 = n - k + 1 here).
rule_costs <- function(coords, r, targets) {
  n <- nrow(coords)
  sets <- sets_of(coords, r, targets)
  m <- lengths(sets)
  big <- m > 1
  cost_points <- sum(as.double(m[big])^2 * (n - vapply(sets[big], min, 1) + 1))
  u <- sort(unique(unlist(sets[big])))
  sub <- coords[u, , drop = FALSE]
  ord <- order(sub[, 1])
  xs <- sub[ord, 1]
  r2 <- 2 * r * (1 + 1e-12)
  cost_pairs <- 0
  pairs <- 0
  for (a in seq_along(u)) {
    lo <- findInterval(sub[a, 1] - r2, xs, left.open = TRUE) + 1
    hi <- findInterval(sub[a, 1] + r2, xs)
    cand <- ord[lo:hi]
    cand <- cand[(sub[cand, 1] - sub[a, 1])^2 +
      (sub[cand, 2] - sub[a, 2])^2 <= r2^2]
    k <- u[cand]
    k <- k[k >= u[a]]
    pairs <- pairs + length(k)
    cost_pairs <- cost_pairs + sum(n - k + 1)
  }
  list(cost_points = cost_points, cost_pairs = cost_pairs, pairs = pairs,
    mean_m = mean(m), max_m = max(m))
}

rows <- list()
for (design in c("uniform", "clustered")) for (n in ns) {
  d <- if (design == "uniform") sim_data(n) else clustered_data(n)
  coords <- as.matrix(d[, c("x", "y")])
  p <- kcv_prepare(z ~ 1, ~ x + y, d, model = sim_model())
  y <- p$inp$y
  # A mean hold-out set of n / 4 or more is not buffered LOO any more (and
  # at smoke sizes would hold out every point).
  for (m in ms[ms < n / 4]) {
    r <- round(sqrt(m * 1e6 / (pi * n)), 1)
    tg <- if (n * m > 1e6) withr::with_seed(1, sort(sample(n, 500))) else NULL
    targets <- if (is.null(tg)) seq_len(n) else tg
    message(sprintf("%s, n = %d, r = %g, targets = %d", design, n, r,
      length(targets)))
    cs <- rule_costs(coords, r, targets)
    tpa <- time_expr(a <- fastkcv:::run_buffer_prepared(p, y, r, tg,
      method = "pairs"), quick = 1)
    tpo <- time_expr(b <- fastkcv:::run_buffer_prepared(p, y, r, tg,
      method = "points"), quick = 1)
    auto_r <- if (cs$cost_pairs < cs$cost_points) "pairs" else "points"
    auto_cpp <- if (min(tpa$time, tpo$time) < 5) {
      fastkcv:::run_buffer_prepared(p, y, r, tg)$method
    } else {
      NA_character_
    }
    rows[[length(rows) + 1]] <- data.frame(
      design = design, n = n, target_m = m, radius = r, n_targets = length(targets),
      mean_m = cs$mean_m, max_m = cs$max_m, pairs = cs$pairs, cost_pairs = cs$cost_pairs,
      cost_points = cs$cost_points, pairs_s = tpa$time, points_s = tpo$time,
      auto = auto_r, auto_cpp = auto_cpp,
      faster = if (tpa$time < tpo$time) "pairs" else "points",
      max_abs_diff = max(abs(a$pred - b$pred), abs(a$var - b$var))
    )
  }
  rm(p)
  gc()
}
res <- do.call(rbind, rows)
if (any(!is.na(res$auto_cpp) & res$auto_cpp != res$auto)) {
  warning("the R replica of the auto rule disagrees with C++")
}
write_results(res, "auto-rule")
