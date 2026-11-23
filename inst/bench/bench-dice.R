# fastkcv's kcv_engine() vs DiceKriging, which implements the same
# closed-form formulas: leaveOneOut.km() for leave-one-out (Dubrule 1983)
# and cv() for multiple folds (Ginsbourger and Schärer, DiceKriging >=
# 1.6.0). Both handle a homogeneous nugget (km(nugget = )); DiceKriging
# excludes only observations with individual noise variances (noise.var).
# Arms: nugget 0 and 0.1 (the variances of both packages then include the
# nugget); schemes leave-one-out and random 5- and 10-fold (fold_random();
# cv(light = TRUE), which skips the covariances between folds that fastkcv
# does not compute either).
#
# Both start from the same covariance: a DiceKriging km() object with all
# parameters fixed (covtype "exp" and "matern5_2"; DiceKriging's kernels are
# tensor products over the coordinates, so the matrix is taken from
# DiceKriging itself and passed to kcv_engine()). Two settings:
#   SK: known mean, type = "SK", trend.reestim = FALSE
#       vs kcv_engine(beta = mean);
#   OK: mean re-estimated by GLS in every fold,
#       type = "UK", trend.reestim = TRUE vs kcv_engine().
# (The two mixed combinations of `type` and `trend.reestim` pair the mean of
# one setting with the variance of the other, so they match neither.)
#
# Timed: dice_s = leaveOneOut.km() or cv() on a fitted km object;
# fastkcv_s = kcv_engine() from the covariance matrix; plus km_s (building
# the km object, which factorizes the covariance) and cov_s (DiceKriging's
# covariance matrix alone), so the end-to-end totals are km_s + dice_s and
# cov_s + fastkcv_s.
#
# DiceKriging is not a dependency of fastkcv; the script needs it installed.
# Writes results/dice-<blas>.csv. See common.R for how to run it.

source(file.path("inst", "bench", "common.R"))
if (!requireNamespace("DiceKriging", quietly = TRUE)) {
  stop("bench-dice.R needs the DiceKriging package.")
}

cap <- bench_cap()
ns <- bench_ns(c(250, 500, 1000, 2000, 5000))
mu <- 5
sigma2 <- 1.1
theta <- 200

hist <- list()
rows <- list()
for (n in ns) {
  d <- sim_data(n)
  X <- d[, c("x", "y")]
  F <- matrix(1, n, 1)
  schemes <- list(
    LOO = seq_len(n),
    `random 10` = withr::with_seed(10, fold_random(~ x + y, d, k = 10)),
    `random 5` = withr::with_seed(5, fold_random(~ x + y, d, k = 5))
  )
  for (nugget in c(0, 0.1)) {
    for (covtype in c("exp", "matern5_2")) {
      tk <- time_expr(m <- DiceKriging::km(~1, design = X, response = d$z,
        covtype = covtype, coef.trend = mu, coef.cov = c(theta, theta),
        coef.var = sigma2, nugget = if (nugget > 0) nugget,
        control = list(trace = FALSE)))
      tc <- time_expr(C <- DiceKriging::covMatrix(m@covariance,
        as.matrix(X))$C)
      for (sname in names(schemes)) {
        folds <- schemes[[sname]]
        fold_list <- unname(split(seq_len(n), folds))
        for (setting in c("SK", "OK")) {
          message(sprintf("n = %d, nugget %g, %s, %s, %s", n, nugget,
            covtype, sname, setting))
          key <- paste(nugget, covtype, sname, setting)
          pred <- predict_time(hist[[key]], n, 3)
          tf <- time_expr(e <- kcv_engine(C, F, d$z, folds,
            beta = if (setting == "SK") mu))
          type <- if (setting == "SK") "SK" else "UK"
          if (!is.na(pred) && pred > cap) {
            status <- "skipped"
            td <- list(time = NA_real_, mem = NA_real_)
            diff_pred <- diff_var <- NA_real_
          } else {
            status <- "ok"
            if (sname == "LOO") {
              td <- time_expr(l <- DiceKriging::leaveOneOut.km(m,
                type = type, trend.reestim = setting == "OK"))
              dp <- l$mean
              dv <- l$sd^2
            } else {
              td <- time_expr(l <- DiceKriging::cv(m, fold_list, type = type,
                trend.reestim = setting == "OK", fast = TRUE, light = TRUE))
              idx <- unlist(fold_list)
              dp <- dv <- numeric(n)
              dp[idx] <- unlist(l$mean)
              dv[idx] <- unlist(lapply(l$cvcov.list, diag))
            }
            hist[[key]] <- c(hist[[key]], list(list(n = n, time = td$time)))
            diff_pred <- max(abs(dp - e$pred))
            diff_var <- max(abs(dv - e$var))
          }
          rows[[length(rows) + 1]] <- data.frame(
            n = n, nugget = nugget, covtype = covtype, scheme = sname,
            setting = setting, status = status,
            fastkcv_s = tf$time, dice_s = td$time, km_s = tk$time,
            cov_s = tc$time, speedup = td$time / tf$time,
            speedup_total = (tk$time + td$time) / (tc$time + tf$time),
            fastkcv_mem = tf$mem, dice_mem = td$mem, dice_predicted_s = pred,
            max_abs_diff_pred = diff_pred, max_abs_diff_var = diff_var
          )
        }
      }
      rm(m, C)
      gc()
    }
  }
}
write_results(do.call(rbind, rows), "dice")
