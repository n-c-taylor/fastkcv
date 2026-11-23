# Accuracy against a high-precision oracle as the covariance matrix becomes
# ill-conditioned: kcv(), gstat::krige.cv() and kcv_reference() (the
# per-fold refit with solve()) against the same cross-validation computed
# in quadruple precision (quad-oracle.cpp, compiled with Rcpp::sourceCpp()).
#
# The oracle takes the double-precision covariance matrix C as exact input,
# so the comparison isolates the rounding error of each method. Data: n
# points uniform on [0, 1000]^2 (seed 2), z = sin(x / 150) + cos(y / 200),
# n = 100 and 200. Models, sill 1:
#   Gau range: Gaussian without nugget, range 100 to 300;
#   Gau nugget, Ste nugget: Gaussian (range 200) and Stein (kappa 5,
#     range 100) with nugget 1e-1 down to 0.
# Ordinary kriging (z ~ 1) and universal kriging (z ~ x + y);
# leave-one-out and random 5-fold.
#
# Recorded: the reciprocal condition number estimated by fastkcv (dpocon),
# whether kcv() warned, the largest absolute errors in var1.pred and
# var1.var, and the error in var1.pred relative to the leave-one-out RMSE
# of the exact predictions (the error a user would compare with). The rcond
# warning threshold (getOption("fastkcv.rcond_warn")) is the largest rcond
# at which kcv()'s relative error exceeded 1e-3, rounded up to a power of
# ten, the larger over the two BLAS runs (1e-14 with OpenBLAS, 1e-13 with
# the reference BLAS, so 1e-13); the script prints its own.
#
# Writes results/conditioning-<blas>.csv. See common.R for how to run it.

source(file.path("inst", "bench", "common.R"))
Rcpp::sourceCpp(file.path(bench_dir, "quad-oracle.cpp"))

# The quadruple-precision oracle is slow (software floating point).
ns <- bench_ns(c(100, 200), smoke = 40)
nuggets <- c(1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 0)
cases <- c(
  lapply(c(100, 130, 160, 200, 230, 250, 270, 300), function(r) {
    list(family = "Gau range", par = r, model = vgm(1, "Gau", r))
  }),
  lapply(nuggets, function(g) {
    list(family = "Gau nugget", par = g,
      model = if (g > 0) vgm(1, "Gau", 200, g) else vgm(1, "Gau", 200))
  }),
  lapply(nuggets, function(g) {
    list(family = "Ste nugget", par = g,
      model = if (g > 0) vgm(1, "Ste", 100, g, kappa = 5) else
        vgm(1, "Ste", 100, kappa = 5))
  })
)
formulas <- list(OK = z ~ 1, UK = z ~ x + y)

attempt <- function(expr) {
  warned <- FALSE
  val <- tryCatch(
    withCallingHandlers(expr, warning = function(w) {
      warned <<- TRUE
      invokeRestart("muffleWarning")
    }),
    error = function(e) e
  )
  list(value = val, warned = warned, failed = inherits(val, "error"))
}
errs <- function(res, truth) {
  if (res$failed) return(c(NA_real_, NA_real_))
  v <- res$value
  pred <- if (is.data.frame(v)) v$var1.pred else drop(v$pred)
  var <- if (is.data.frame(v)) v$var1.var else v$var
  c(max(abs(pred - truth$pred)), max(abs(var - truth$var)))
}

rows <- list()
for (n in ns) {
  d <- withr::with_seed(2, data.frame(x = runif(n, 0, 1000),
    y = runif(n, 0, 1000)))
  d$z <- sin(d$x / 150) + cos(d$y / 200)
  xy <- as.matrix(d[, c("x", "y")])
  schemes <- list(
    LOO = seq_len(n),
    `random 5` = withr::with_seed(5, fold_random(xy, k = 5))
  )
  for (cs in cases) {
    C <- fastkcv:::covariance_matrix(cs$model, xy)
    for (fn in names(formulas)) {
      F <- model.matrix(formulas[[fn]], d)
      p <- attempt(kcv_prepare(formulas[[fn]], ~ x + y, d, model = cs$model))
      rc <- if (p$failed) NA_real_ else p$value$rcond
      for (s in names(schemes)) {
        f <- schemes[[s]]
        message(sprintf("n = %d, %s %g, %s, %s", n, cs$family, cs$par, fn, s))
        fold_list <- unname(split(seq_len(n), f))
        truth <- tryCatch(
          quad_kcv(C, F, d$z, fold_list, FALSE, numeric(0)),
          error = function(e) NULL
        )
        if (is.null(truth)) next
        loo_truth <- if (s == "LOO") truth else
          quad_kcv(C, F, d$z, as.list(seq_len(n)), FALSE, numeric(0))
        rmse <- sqrt(mean((d$z - loo_truth$pred)^2))
        fk <- attempt(kcv(formulas[[fn]], ~ x + y, d, model = cs$model,
          folds = f))
        gs <- attempt(krige.cv(formulas[[fn]], ~ x + y, d, model = cs$model,
          nfold = f, verbose = FALSE))
        rf <- attempt(fastkcv:::kcv_reference(C, F, d$z, f))
        e_fk <- errs(fk, truth)
        e_gs <- errs(gs, truth)
        e_rf <- errs(rf, truth)
        rows[[length(rows) + 1]] <- data.frame(
          n = n, family = cs$family, par = cs$par, trend = fn, scheme = s,
          rcond = rc, kcv_warned = fk$warned, kcv_failed = fk$failed,
          krige_cv_failed = gs$failed, reference_failed = rf$failed,
          loo_rmse = rmse,
          kcv_pred_err = e_fk[1], kcv_var_err = e_fk[2],
          krige_cv_pred_err = e_gs[1], krige_cv_var_err = e_gs[2],
          reference_pred_err = e_rf[1], reference_var_err = e_rf[2],
          kcv_rel_err = e_fk[1] / rmse, krige_cv_rel_err = e_gs[1] / rmse
        )
      }
    }
  }
}
res <- do.call(rbind, rows)
write_results(res, "conditioning")

bad <- res[!res$kcv_failed & res$kcv_rel_err > 1e-3, ]
if (nrow(bad)) {
  message(sprintf(
    "Largest rcond with kcv() error above 1e-3 of the LOO RMSE: %.3g; ",
    max(bad$rcond)
  ), sprintf("threshold 10^%d", ceiling(log10(max(bad$rcond)))))
}
