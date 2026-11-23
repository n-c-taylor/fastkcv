# Summary tables from the raw results in inst/bench/results/ (one CSV per
# benchmark script and BLAS). From the package root:
#   Rscript inst/bench/summarize.R            # tables to the console
#   Rscript inst/bench/summarize.R summary.md # also written as Markdown

out_file <- commandArgs(trailingOnly = TRUE)[1]
results_dir <- file.path("inst", "bench", "results")

read_all <- function(name) {
  files <- list.files(results_dir, sprintf("^%s-(openblas|refblas)\\.csv$", name),
    full.names = TRUE)
  if (!length(files)) return(NULL)
  do.call(rbind, lapply(files, utils::read.csv, stringsAsFactors = FALSE))
}

lines <- character()
emit <- function(title, df, note = NULL) {
  if (is.null(df) || !nrow(df)) return(invisible())
  tab <- knitr::kable(df, format = "pipe", row.names = FALSE)
  lines <<- c(lines, "", paste("###", title), if (!is.null(note)) c("", note),
    "", tab)
  cat("\n", title, "\n", sep = "")
  if (!is.null(note)) cat(note, "\n")
  print(df, row.names = FALSE)
}
sig <- function(x, d = 3) signif(x, d)

# 1. Fold schemes: seconds per engine and speed-ups.
f <- read_all("folds")
if (!is.null(f)) {
  key <- c("blas", "data", "n", "scheme")
  wide <- function(eng, col = "time_s") {
    x <- f[f$engine == eng, c(key, col)]
    names(x)[ncol(x)] <- eng
    x
  }
  fk <- f[f$engine == "fastkcv", c(key, "time_s", "prepare_s", "run_s")]
  names(fk)[5:7] <- c("fastkcv", "prepare", "run")
  t1 <- Reduce(function(a, b) merge(a, b, by = key, sort = FALSE),
    list(fk, wide("engine_r"), wide("reference"), wide("refit"),
      wide("krige.cv")))
  t1$`engine_r/fastkcv` <- round(t1$engine_r / t1$fastkcv, 1)
  t1$`reference/fastkcv` <- round(t1$reference / t1$fastkcv, 1)
  t1$`refit/fastkcv` <- round(t1$refit / t1$fastkcv, 2)
  t1$`krige.cv/fastkcv` <- round(t1$krige.cv / t1$fastkcv, 1)
  t1$`prepare share` <- round(t1$prepare / t1$fastkcv, 2)
  t1 <- t1[order(t1$blas, t1$data != "sim", t1$n,
    match(sub(" .*", "", t1$scheme), c("LOO", "random", "blocks")),
    -as.numeric(sub("^\\D*", "", sub("LOO", "1e9", t1$scheme)))), ]
  for (col in c("fastkcv", "prepare", "run", "engine_r", "reference",
    "refit", "krige.cv")) t1[[col]] <- sig(t1[[col]])
  emit("Fold schemes: seconds (bench::mark median) and speed-ups", t1,
    "fastkcv = kcv_prepare() + kcv_run(); NA = skipped (predicted above the time cap).")

  d <- max(f$max_abs_diff, na.rm = TRUE)
  cat(sprintf("\nLargest |difference| in var1.pred/var1.var vs any competitor: %.2g\n", d))
  lines <- c(lines, "", sprintf(
    "Largest absolute difference in var1.pred or var1.var against any competitor: %.2g.", d))

  # Advantage vs number of folds: the speed-up over krige.cv and over the
  # per-fold refit, random folds, as a function of k.
  adv <- t1[t1$data == "sim" & !grepl("blocks", t1$scheme),
    c("blas", "n", "scheme", "krige.cv/fastkcv", "reference/fastkcv",
      "refit/fastkcv", "prepare share")]
  emit("Where the advantage shrinks (simulated data, random folds)", adv,
    paste("Speed-up of fastkcv over the per-fold refits by n and number of",
      "folds; refit/fastkcv < 1 where the Cholesky refit is faster."))
  sk <- f[f$status == "skipped", c("blas", "n", "scheme", "engine",
    "predicted_s")]
  sk$predicted_s <- sig(sk$predicted_s, 2)
  emit("Skipped competitor runs (predicted seconds)", sk)
}

# 2. Fold designs.
g <- read_all("designs")
if (!is.null(g)) {
  for (col in c("run_s", "rmse", "mean_z", "var_z", "coverage_90")) {
    g[[col]] <- sig(g[[col]])
  }
  emit("Fold designs from one factorization", g)
}

# 3. DiceKriging.
dk <- read_all("dice")
if (!is.null(dk)) {
  # Results from before the nugget arm and the k-fold schemes.
  if (is.null(dk$nugget)) dk$nugget <- 0
  if (is.null(dk$scheme)) dk$scheme <- "LOO"
  dk <- dk[, c("blas", "n", "nugget", "covtype", "scheme", "setting",
    "status", "fastkcv_s",
    "dice_s", "km_s", "cov_s", "speedup", "speedup_total",
    "max_abs_diff_pred", "max_abs_diff_var")]
  for (col in c("fastkcv_s", "dice_s", "km_s", "cov_s")) dk[[col]] <- sig(dk[[col]])
  for (col in c("speedup", "speedup_total")) dk[[col]] <- round(dk[[col]], 2)
  for (col in c("max_abs_diff_pred", "max_abs_diff_var")) {
    dk[[col]] <- formatC(dk[[col]], format = "e", digits = 1)
  }
  emit("DiceKriging leaveOneOut.km() (LOO) and cv() (k-fold)", dk,
    paste("speedup = DiceKriging / kcv_engine (both from a ready",
      "covariance; km() has already factorized it for DiceKriging);",
        "speedup_total = (km + DiceKriging CV) / (covariance + kcv_engine)."))
}

# 4. Peak memory.
r <- read_all("rss")
if (!is.null(r)) {
  r <- r[, c("blas", "n", "config", "time_s", "peak_mb", "extra_mb",
    "matrix_mb", "matrices")]
  for (col in c("peak_mb", "extra_mb", "matrix_mb")) r[[col]] <- round(r[[col]])
  r$matrices <- round(r$matrices, 2)
  emit("Peak resident memory (/usr/bin/time -v)", r,
    "extra_mb = peak minus a process with the packages and data loaded; matrices = extra / (8 n^2).")
}

# 5. Buffered LOO at n = 20,000.
b <- read_all("buffer-large")
if (!is.null(b)) {
  b <- b[, c("blas", "radius", "targets", "method", "method_used",
    "mean_heldout", "time_s", "peak_mb", "extra_mb", "pairs",
    "cache_est_mb", "fallback_identical_to_points")]
  for (col in c("time_s", "mean_heldout")) b[[col]] <- sig(b[[col]])
  for (col in c("peak_mb", "extra_mb", "cache_est_mb")) b[[col]] <- round(b[[col]])
  emit("Buffered LOO at n = 20,000", b,
    "extra_mb = peak minus a process that has only read the prepared object (3.2 GB).")
}

# 6. The automatic pairs/points rule.
a <- read_all("auto-rule")
if (!is.null(a)) {
  a$auto_s <- ifelse(a$auto == "pairs", a$pairs_s, a$points_s)
  a$best_s <- pmin(a$pairs_s, a$points_s)
  a$loss <- round(a$auto_s / a$best_s - 1, 2)
  a$cost_ratio <- sig(a$cost_points / a$cost_pairs, 2)
  if (is.null(a$design)) a$design <- "uniform"
  a <- a[, c("blas", "design", "n", "radius", "n_targets", "mean_m",
    "cost_ratio",
    "pairs_s", "points_s", "auto", "faster", "loss")]
  for (col in c("pairs_s", "points_s", "mean_m")) a[[col]] <- sig(a[[col]])
  emit("Automatic pairs/points choice", a,
    "cost_ratio = cost_points / cost_pairs (auto picks pairs when > 1); loss = time of auto's pick / faster - 1.")
}

# 7. Several responses from one factorization.
m <- read_all("multi")
if (!is.null(m)) {
  # Time per extra response, (t(T) - t(1)) / (T - 1), for each method.
  key <- paste(m$blas, m$n, m$scheme)
  for (col in c("kcv_multi", "kcv_loop", "run_multi", "run_loop")) {
    t1 <- m[[col]][m$T == 1][match(key, key[m$T == 1])]
    m[[paste0(col, "_extra")]] <- ifelse(m$T > 1,
      (m[[col]] - t1) / (m$T - 1), NA)
  }
  m$`loop/multi` <- round(m$kcv_loop / m$kcv_multi, 1)
  m <- m[order(m$blas, m$n, match(m$scheme, c("LOO", "10-fold",
    "buffer 50")), m$T), c("blas", "n", "scheme", "T", "prepare_s",
    "kcv_multi", "kcv_loop", "loop/multi", "kcv_multi_extra",
    "kcv_loop_extra", "run_multi", "run_loop", "run_multi_extra",
    "run_loop_extra", "max_abs_diff")]
  for (col in grep("_s$|kcv_|run_", names(m), value = TRUE)) {
    m[[col]] <- sig(m[[col]])
  }
  m$max_abs_diff <- sig(m$max_abs_diff, 2)
  emit("Several responses from one factorization (seconds)", m,
    "*_extra = (t(T) - t(1)) / (T - 1), the time per extra response; kcv_* end to end, run_* run stage only; NA = skipped.")
}

# 8. spmodel.
sm <- read_all("spmodel")
if (!is.null(sm)) {
  for (col in c("fastkcv_s", "spmodel_s", "spmodel_pred_s")) {
    sm[[col]] <- sig(sm[[col]])
  }
  for (col in c("speedup", "speedup_pred")) sm[[col]] <- round(sm[[col]], 1)
  for (col in c("max_abs_diff_pred", "max_abs_diff_var")) {
    sm[[col]] <- formatC(sm[[col]], format = "e", digits = 1)
  }
  emit("spmodel loocv() / kcv() with fixed covariance parameters",
    sm[, c("blas", "n", "scheme", "status", "fastkcv_s", "spmodel_s",
      "spmodel_pred_s", "speedup", "speedup_pred", "max_abs_diff_pred",
      "max_abs_diff_var")],
    paste("speedup = spmodel (splm + CV, with standard errors) / kcv();",
      "speedup_pred = the same with predictions only."))
}

# 9. Local gstat (nmax = 50).
lo <- read_all("local")
if (!is.null(lo)) {
  for (col in c("fastkcv_s", "local_s", "rmse_fastkcv", "rmse_local")) {
    lo[[col]] <- sig(lo[[col]])
  }
  lo$local_over_fastkcv <- round(lo$local_over_fastkcv, 2)
  emit("Global fastkcv vs local krige.cv(nmax = 50)",
    lo[, c("blas", "threads", "n", "scheme", "status", "fastkcv_s",
      "local_s", "local_over_fastkcv", "rmse_fastkcv", "rmse_local")],
    "local_over_fastkcv > 1 where fastkcv is faster; threads = BLAS threads for fastkcv.")
}

# 10. Accuracy against the quadruple-precision oracle.
co <- read_all("conditioning")
if (!is.null(co)) {
  co <- co[order(co$blas, co$n, co$family, -co$par, co$trend, co$scheme), ]
  for (col in c("rcond", "kcv_pred_err", "kcv_var_err", "krige_cv_pred_err",
    "reference_pred_err", "kcv_rel_err", "krige_cv_rel_err")) {
    co[[col]] <- formatC(co[[col]], format = "e", digits = 1)
  }
  emit("Accuracy against a quadruple-precision oracle",
    co[, c("blas", "n", "family", "par", "trend", "scheme", "rcond",
      "kcv_warned", "kcv_pred_err", "kcv_var_err", "kcv_rel_err",
      "krige_cv_pred_err", "krige_cv_rel_err", "reference_pred_err")],
    paste("*_err = largest absolute error in var1.pred (var: var1.var);",
      "*_rel_err = that error / the exact LOO RMSE; NA = the method failed."))
}

if (!is.na(out_file)) {
  writeLines(c("# Benchmark summary", lines), out_file)
  message("wrote ", out_file)
}
