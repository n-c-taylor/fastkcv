# Leave-one-out and k-fold cross-validation (global ordinary kriging):
# fastkcv vs the pure-R explicit-Q engine (kcv_engine_r), the naive per-fold
# refit (kcv_reference) and gstat::krige.cv().
#
# n in {250, 500, 1000, 2000, 5000, 10000} (simulated, fixed seeds) plus
# meuse (n = 155, OK and UK); fold schemes: LOO, random k-fold (k = 25, 10,
# 5; fold_random(), the same draw as krige.cv(nfold = k)) and spatial block
# folds (a 10 x 10 grid of blocks grouped at random into k folds,
# fold_blocks(blocks = 10, k = k)).
#
# fastkcv is timed in its two stages: kcv_prepare() once per data set
# (covariance and factorization, independent of the folds) and kcv_run()
# per fold scheme; kcv_s = prepare_s + run_s is the cost of one kcv() call.
# The competitors are timed end to end. A competitor is skipped when the
# time predicted from its previous (smaller) n exceeds `cap` seconds; the
# row then records status "skipped" and the predicted time.
#
# Writes results/folds-<blas>.csv and results/designs-<blas>.csv.
# See common.R for how to run it.

source(file.path("inst", "bench", "common.R"))

cap <- as.numeric(Sys.getenv("FASTKCV_BENCH_CAP", "600"))
ns <- c(250, 500, 1000, 2000, 5000, 10000)
# Smoke test: FASTKCV_BENCH_NS="250,500" FASTKCV_BENCH_CAP=5
if (nzchar(Sys.getenv("FASTKCV_BENCH_NS"))) {
  ns <- as.numeric(strsplit(Sys.getenv("FASTKCV_BENCH_NS"), ",")[[1]])
}
ks <- c(25, 10, 5)

datasets <- c(
  list(
    list(name = "meuse OK", formula = log(zinc) ~ 1, data = meuse_data(),
      model = vgm(0.5, "Sph", 900, 0.1)),
    list(name = "meuse UK", formula = log(zinc) ~ sqrt(dist),
      data = meuse_data(), model = vgm(0.5, "Sph", 900, 0.1))
  ),
  lapply(ns, function(n) {
    list(name = "sim", formula = z ~ 1, data = sim_data(n), model = sim_model())
  })
)

schemes <- function(d) {
  n <- nrow(d)
  out <- list(LOO = list(type = "LOO", k = n, folds = seq_len(n)))
  for (k in ks) {
    out[[paste0("random ", k)]] <- list(type = "random", k = k,
      folds = withr::with_seed(k, fold_random(~ x + y, d, k = k)))
    out[[paste0("blocks ", k)]] <- list(type = "blocks", k = k,
      folds = withr::with_seed(k, fold_blocks(~ x + y, d, blocks = 10, k = k)))
  }
  out
}

# Time growth per competitor: LOO refits n systems of size n - 1 (n^4);
# k-fold refits k systems of size about n (n^3); kcv_engine_r is n^3.
powers <- list(
  engine_r = c(LOO = 3, random = 3, blocks = 3),
  reference = c(LOO = 4, random = 3, blocks = 3),
  krige.cv = c(LOO = 4, random = 3, blocks = 3)
)

last <- list() # measurements per competitor and scheme, for predictions
rows <- list()
for (ds in datasets) {
  d <- ds$data
  n <- nrow(d)
  message(sprintf("%s, n = %d: prepare", ds$name, n))
  env <- environment()
  tp <- time_expr(p <- kcv_prepare(ds$formula, ~ x + y, d, model = ds$model),
    env)
  for (sname in names(sc <- schemes(d))) {
    s <- sc[[sname]]
    folds <- s$folds
    message(sprintf("%s, n = %d, %s", ds$name, n, sname))
    tr <- time_expr(out <- kcv_run(p, folds = folds), env)
    base <- data.frame(
      data = ds$name, n = n, scheme = sname, type = s$type, k = s$k,
      n_folds = length(unique(folds)), largest_fold = max(table(folds))
    )
    rows[[length(rows) + 1]] <- cbind(base, engine = "fastkcv",
      status = "ok", time_s = tp$time + tr$time, prepare_s = tp$time,
      run_s = tr$time, mem_bytes = tp$mem + tr$mem, iterations = tr$iter,
      predicted_s = NA_real_, max_abs_diff = NA_real_)

    for (eng in c("engine_r", "reference", "krige.cv")) {
      key <- paste(eng, s$type, if (s$type == "LOO") 0 else s$k)
      pred <- if (ds$name == "sim") {
        predict_time(last[[key]], n, powers[[eng]][[s$type]])
      } else {
        NA_real_
      }
      if (!is.na(pred) && pred > cap) {
        rows[[length(rows) + 1]] <- cbind(base, engine = eng,
          status = "skipped", time_s = NA_real_, prepare_s = NA_real_,
          run_s = NA_real_, mem_bytes = NA_real_, iterations = 0L,
          predicted_s = pred, max_abs_diff = NA_real_)
        next # later predictions still start from the last measured n
      }
      t <- switch(eng,
        engine_r = time_expr(ref <- kcv_r(ds$formula, ~ x + y, d, ds$model,
          folds, "engine_r"), env),
        reference = time_expr(ref <- kcv_r(ds$formula, ~ x + y, d, ds$model,
          folds, "reference"), env),
        krige.cv = time_expr(ref <- krige.cv(ds$formula, ~ x + y, d,
          model = ds$model, nfold = folds, verbose = FALSE), env)
      )
      if (ds$name == "sim") {
        last[[key]] <- c(last[[key]], list(list(n = n, time = t$time)))
      }
      rows[[length(rows) + 1]] <- cbind(base, engine = eng, status = "ok",
        time_s = t$time, prepare_s = NA_real_, run_s = NA_real_,
        mem_bytes = t$mem, iterations = t$iter, predicted_s = pred,
        max_abs_diff = max_abs_diff(out, ref))
    }
  }
  rm(p)
  gc()
}
res <- do.call(rbind, rows)
# max_abs_diff: largest absolute difference in var1.pred or var1.var
# between fastkcv and that competitor.
write_results(res, "folds")

# Fold designs from one factorization: what each design costs and what it
# says about the model (RMSE, variance of z-scores, 90% coverage).
rows <- list()
for (n in intersect(c(1000, 5000), ns)) {
  d <- sim_data(n)
  p <- kcv_prepare(z ~ 1, ~ x + y, d, model = sim_model())
  designs <- withr::with_seed(1, list(
    `LOO` = seq_len(n),
    `random k = 10` = fold_random(~ x + y, d, k = 10),
    `random k = 25` = fold_random(~ x + y, d, k = 25),
    `blocks 5 x 5` = fold_blocks(~ x + y, d, blocks = 5),
    `blocks 10 x 10, k = 10` = fold_blocks(~ x + y, d, blocks = 10, k = 10),
    `kmeans k = 10` = fold_kmeans(~ x + y, d, k = 10),
    `kmeans k = 25` = fold_kmeans(~ x + y, d, k = 25)
  ))
  for (nm in names(designs)) {
    f <- designs[[nm]]
    t <- time_expr(out <- kcv_run(p, folds = f))
    m <- kcv_metrics(out)
    rows[[length(rows) + 1]] <- data.frame(
      n = n, design = nm, n_folds = length(unique(f)),
      largest_fold = max(table(f)), run_s = t$time,
      rmse = m$rmse, mean_z = m$mean_z, var_z = m$var_z,
      coverage_90 = m$coverage_90
    )
  }
  for (r in c(25, 50, 100)) {
    t <- time_expr(out <- kcv_run(p, buffer = r))
    m <- kcv_metrics(out)
    rows[[length(rows) + 1]] <- data.frame(
      n = n, design = sprintf("buffered LOO r = %g", r), n_folds = n,
      largest_fold = max(out$n_heldout), run_s = t$time,
      rmse = m$rmse, mean_z = m$mean_z, var_z = m$var_z,
      coverage_90 = m$coverage_90
    )
  }
}
write_results(do.call(rbind, rows), "designs")
