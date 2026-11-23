# fastkcv's kcv() vs spmodel's loocv() and kcv(), ordinary kriging with
# fixed covariance parameters on sim_data().
#
# spmodel's "exponential" covariance with de = 1, ie = 0.1, range = 200 is
# exactly gstat's vgm(1, "Exp", 200, 0.1) (checked: identical matrices), and
# spcov_initial(known = "given") fixes every parameter, so splm() only
# estimates the trend by generalized least squares. spmodel's
# cross-validation is exact too: one inverse of the covariance matrix, then
# for each fold an update of it for the remaining data, with the trend
# re-estimated (about n^2 operations per fold for the predictions and more
# for the standard errors). Above n = 5,000 spmodel switches to local
# approximations by default; the sizes here stay at or below that.
#
# Timed end to end: spmodel_s = splm() + loocv()/kcv() with predictions and
# standard errors (the output kcv() always gives); spmodel_pred_s = the same
# with predictions only; fastkcv_s = one kcv() call. Schemes:
# leave-one-out and random 10- and 5-fold (fold_random(), passed to
# spmodel as `folds_index`).
#
# spmodel is not a dependency of fastkcv; the script needs it installed.
# Writes results/spmodel-<blas>.csv. See common.R for how to run it.

source(file.path("inst", "bench", "common.R"))
if (!requireNamespace("spmodel", quietly = TRUE)) {
  stop("bench-spmodel.R needs the spmodel package.")
}

cap <- bench_cap()
ns <- bench_ns(c(250, 500, 1000, 2000, 5000))

spmodel_cv <- function(d, folds, se) {
  fit <- spmodel::splm(z ~ 1, data = d, xcoord = "x", ycoord = "y",
    spcov_initial = spmodel::spcov_initial("exponential", de = 1, ie = 0.1,
      range = 200, known = "given"))
  if (length(unique(folds)) == nrow(d)) {
    spmodel::loocv(fit, cv_predict = TRUE, se.fit = se, local = FALSE)
  } else {
    spmodel::kcv(fit, folds_index = folds, cv_predict = TRUE, se.fit = se,
      local = FALSE)
  }
}

hist <- list()
rows <- list()
for (n in ns) {
  d <- sim_data(n)
  schemes <- list(
    LOO = seq_len(n),
    `random 10` = withr::with_seed(10, fold_random(~ x + y, d, k = 10)),
    `random 5` = withr::with_seed(5, fold_random(~ x + y, d, k = 5))
  )
  for (s in names(schemes)) {
    f <- schemes[[s]]
    message(sprintf("n = %d, %s", n, s))
    tf <- time_expr(cv <- kcv(z ~ 1, ~ x + y, d, model = sim_model(),
      folds = f))
    pred <- predict_time(hist[[s]], n, 3)
    if (!is.na(pred) && pred > cap) {
      status <- "skipped"
      ts <- tp <- list(time = NA_real_)
      diff_pred <- diff_var <- NA_real_
    } else {
      status <- "ok"
      ts <- time_expr(sp <- spmodel_cv(d, f, TRUE), max_iter = 5)
      tp <- time_expr(spmodel_cv(d, f, FALSE), max_iter = 5)
      hist[[s]] <- c(hist[[s]], list(list(n = n, time = ts$time)))
      diff_pred <- max(abs(sp$cv_predict - cv$var1.pred))
      diff_var <- max(abs(sp$se.fit^2 - cv$var1.var))
    }
    rows[[length(rows) + 1]] <- data.frame(
      n = n, scheme = s, status = status, fastkcv_s = tf$time,
      spmodel_s = ts$time, spmodel_pred_s = tp$time,
      speedup = ts$time / tf$time, speedup_pred = tp$time / tf$time,
      spmodel_predicted_s = pred,
      max_abs_diff_pred = diff_pred, max_abs_diff_var = diff_var
    )
  }
}
write_results(do.call(rbind, rows), "spmodel")
