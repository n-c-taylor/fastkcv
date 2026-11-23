# Global kriging with fastkcv vs local kriging with gstat
# (krige.cv(nmax = 50)), ordinary kriging on sim_data().
#
# Local kriging solves a different, approximate kriging system (the 50
# nearest observations of each prediction location), so this is a
# comparison of cost, not of results. Its time grows about linearly in n
# (one small system per prediction) while fastkcv's single factorization
# grows as n^3, so the two cross over at some n that depends on the BLAS.
#
# Schemes: random 10-fold (fold_random(k = 10)) and leave-one-out. fastkcv
# is timed end to end (one kcv() call). Under OpenBLAS each case runs with
# 4 BLAS threads (one per physical core of the 4-core benchmark laptop) and
# with OpenBLAS's default (one per logical CPU); the reference BLAS is
# single-threaded. Local gstat does not use the BLAS for its small systems
# in any way that threads help, so its time is taken once per n.
#
# Writes results/local-<blas>.csv. See common.R for how to run it.

source(file.path("inst", "bench", "common.R"))
if (!requireNamespace("RhpcBLASctl", quietly = TRUE)) {
  stop("bench-local.R needs the RhpcBLASctl package.")
}

cap <- bench_cap()
ns <- bench_ns(c(1000, 2000, 3000, 4000))
tag <- blas_tag()
default_threads <- RhpcBLASctl::blas_get_num_procs()
threads <- if (tag == "openblas") unique(c(4, default_threads)) else 1

rmse <- function(cv) sqrt(mean(cv$residual^2))
hist <- list()
rows <- list()
for (n in ns) {
  d <- sim_data(n)
  sch <- list(
    `random 10` = withr::with_seed(10, fold_random(~ x + y, d, k = 10)),
    LOO = seq_len(n)
  )
  for (s in names(sch)) {
    f <- sch[[s]]
    message(sprintf("n = %d, %s", n, s))
    pred <- predict_time(hist[[s]], n, 1.5)
    if (!is.na(pred) && pred > cap) {
      tl <- list(time = NA_real_)
      status <- "skipped"
      rmse_l <- NA_real_
    } else {
      tl <- time_expr(l <- krige.cv(z ~ 1, ~ x + y, d, model = sim_model(),
        nmax = 50, nfold = f, verbose = FALSE), max_iter = 3)
      hist[[s]] <- c(hist[[s]], list(list(n = n, time = tl$time)))
      status <- "ok"
      rmse_l <- rmse(l)
    }
    for (th in threads) {
      RhpcBLASctl::blas_set_num_threads(th)
      tf <- time_expr(cv <- kcv(z ~ 1, ~ x + y, d, model = sim_model(),
        folds = f))
      rows[[length(rows) + 1]] <- data.frame(
        n = n, scheme = s, threads = th, status = status,
        fastkcv_s = tf$time, local_s = tl$time,
        local_over_fastkcv = tl$time / tf$time, local_predicted_s = pred,
        rmse_fastkcv = rmse(cv), rmse_local = rmse_l
      )
    }
    RhpcBLASctl::blas_set_num_threads(default_threads)
  }
}
write_results(do.call(rbind, rows), "local")
