# Leave-one-out without a nugget: fastkcv's kcv_engine() vs
# DiceKriging::leaveOneOut.km(), which implements the same closed-form
# formulas (Dubrule 1983) but has no nugget / measurement-error term.
#
# Both start from the same covariance: a DiceKriging km() object with all
# parameters fixed (covtype "exp" and "matern5_2"; DiceKriging's kernels are
# tensor products over the coordinates, so the matrix is taken from
# DiceKriging itself and passed to kcv_engine()). Two settings:
#   SK: known mean, leaveOneOut.km(type = "SK", trend.reestim = FALSE)
#       vs kcv_engine(beta = mean);
#   OK: mean re-estimated by GLS in every fold,
#       leaveOneOut.km(type = "UK", trend.reestim = TRUE) vs kcv_engine().
# (The two mixed combinations of `type` and `trend.reestim` pair the mean of
# one setting with the variance of the other, so they match neither.)
#
# Timed: dice_s = leaveOneOut.km() on a fitted km object; fastkcv_s =
# kcv_engine() from the covariance matrix; plus km_s (building the km object,
# which factorizes the covariance) and cov_s (DiceKriging's covariance
# matrix alone), so the end-to-end totals are km_s + dice_s and
# cov_s + fastkcv_s.
#
# DiceKriging is not a dependency of fastkcv; the script needs it installed.
# Writes results/dice-<blas>.csv. See common.R for how to run it.

source(file.path("inst", "bench", "common.R"))
if (!requireNamespace("DiceKriging", quietly = TRUE)) {
  stop("bench-dice.R needs the DiceKriging package.")
}

cap <- as.numeric(Sys.getenv("FASTKCV_BENCH_CAP", "600"))
ns <- c(250, 500, 1000, 2000, 5000)
if (nzchar(Sys.getenv("FASTKCV_BENCH_NS"))) {
  ns <- as.numeric(strsplit(Sys.getenv("FASTKCV_BENCH_NS"), ",")[[1]])
}
mu <- 5
sigma2 <- 1.1
theta <- 200

hist <- list()
rows <- list()
for (n in ns) {
  d <- sim_data(n)
  X <- d[, c("x", "y")]
  for (covtype in c("exp", "matern5_2")) {
    tk <- time_expr(m <- DiceKriging::km(~1, design = X, response = d$z,
      covtype = covtype, coef.trend = mu, coef.cov = c(theta, theta),
      coef.var = sigma2))
    tc <- time_expr(C <- DiceKriging::covMatrix(m@covariance,
      as.matrix(X))$C)
    F <- matrix(1, n, 1)
    for (setting in c("SK", "OK")) {
      message(sprintf("n = %d, %s, %s", n, covtype, setting))
      key <- paste(covtype, setting)
      pred <- predict_time(hist[[key]], n, 3)
      tf <- time_expr(e <- kcv_engine(C, F, d$z, seq_len(n),
        beta = if (setting == "SK") mu))
      if (!is.na(pred) && pred > cap) {
        status <- "skipped"
        td <- list(time = NA_real_, mem = NA_real_)
        diff_pred <- diff_var <- NA_real_
      } else {
        status <- "ok"
        td <- time_expr(l <- DiceKriging::leaveOneOut.km(m,
          type = if (setting == "SK") "SK" else "UK",
          trend.reestim = setting == "OK"))
        hist[[key]] <- c(hist[[key]], list(list(n = n, time = td$time)))
        diff_pred <- max(abs(l$mean - e$pred))
        diff_var <- max(abs(l$sd^2 - e$var))
      }
      rows[[length(rows) + 1]] <- data.frame(
        n = n, covtype = covtype, setting = setting, status = status,
        fastkcv_s = tf$time, dice_s = td$time, km_s = tk$time,
        cov_s = tc$time, speedup = td$time / tf$time,
        speedup_total = (tk$time + td$time) / (tc$time + tf$time),
        fastkcv_mem = tf$mem, dice_mem = td$mem, dice_predicted_s = pred,
        max_abs_diff_pred = diff_pred, max_abs_diff_var = diff_var
      )
    }
    rm(m, C)
    gc()
  }
}
write_results(do.call(rbind, rows), "dice")
