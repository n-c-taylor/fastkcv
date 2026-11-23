# Several responses from one factorization: T = 1, 10, 100 responses at
# n = 1000 and 5000 on sim_data(), ordinary kriging, for LOO, 10-fold
# (random, fixed seed) and buffered LOO (r = 50).
#
# Methods (seconds, time_expr() median):
#   kcv_multi  kcv(Y ~ 1) with an n x T matrix column Y (end to end);
#   kcv_loop   T single-response kcv() calls (each one factorizes again);
#   run_multi  kcv_run(p, y = Y) on one prepared object (run stage only);
#   run_loop   T kcv_run(p, y = Y[, t]) calls (run stage only).
# prepare_s is the time of kcv_prepare() alone. kcv_loop is skipped (NA)
# when T times its T = 1 time exceeds FASTKCV_BENCH_CAP (default 600 s).
# max_abs_diff compares every response block of kcv_multi with the
# corresponding single run (var1.pred and var1.var).
# Writes results/multi-<blas>.csv. See common.R for how to run it.

source(file.path("inst", "bench", "common.R"))

cap <- bench_cap()
ns <- bench_ns(c(1000, 5000))
Ts <- c(1, 10, 100)
model <- sim_model()

rows <- list()
for (n in ns) {
  d <- sim_data(n)
  # Response t is the simulated field plus independent noise; the values do
  # not affect the time.
  Y_all <- withr::with_seed(n, d$z + matrix(stats::rnorm(n * max(Ts), 0,
    0.3), n, max(Ts)))
  folds10 <- withr::with_seed(1, sample(10, n, replace = TRUE))
  schemes <- list(
    LOO = list(),
    `10-fold` = list(folds = folds10),
    `buffer 50` = list(buffer = 50)
  )
  t_prep <- time_expr(kcv_prepare(z ~ 1, ~ x + y, d, model = model))$time
  p <- kcv_prepare(z ~ 1, ~ x + y, d, model = model)
  for (sn in names(schemes)) {
    sch <- schemes[[sn]]
    t1_loop <- NA_real_
    for (nT in Ts) {
      message(sprintf("n = %d, %s, T = %d", n, sn, nT))
      d$Y <- Y_all[, seq_len(nT), drop = FALSE]
      Y <- d$Y
      multi <- do.call(kcv, c(list(Y ~ 1, ~ x + y, d, model = model), sch))
      kcv_multi <- time_expr(do.call(kcv, c(list(Y ~ 1, ~ x + y, d,
        model = model), sch)))$time
      run_multi <- time_expr(do.call(kcv_run, c(list(p, y = Y), sch)))$time
      run_loop <- time_expr(for (t in seq_len(nT)) {
        do.call(kcv_run, c(list(p, y = Y[, t]), sch))
      })$time
      loop_ok <- is.na(t1_loop) || nT * t1_loop <= cap
      kcv_loop <- NA_real_
      if (loop_ok) {
        kcv_loop <- time_expr(for (t in seq_len(nT)) {
          d$yt <- Y[, t]
          do.call(kcv, c(list(yt ~ 1, ~ x + y, d, model = model), sch))
        })$time
        if (nT == 1) t1_loop <- kcv_loop
      }
      # Agreement of every block with the single-response run.
      diff <- 0
      for (t in seq_len(nT)) {
        one <- do.call(kcv_run, c(list(p, y = Y[, t]), sch))
        # T = 1: a one-column matrix column is a single response.
        b <- if (nT == 1) multi else
          multi[multi$response == levels(multi$response)[t], ]
        diff <- max(diff, max_abs_diff(b, one))
      }
      rows[[length(rows) + 1]] <- data.frame(
        n = n, scheme = sn, T = nT, prepare_s = t_prep,
        kcv_multi = kcv_multi, kcv_loop = kcv_loop,
        run_multi = run_multi, run_loop = run_loop,
        max_abs_diff = diff,
        status = if (loop_ok) "ok" else "kcv_loop skipped"
      )
      print(rows[[length(rows)]])
    }
  }
}
write_results(do.call(rbind, rows), "multi")
