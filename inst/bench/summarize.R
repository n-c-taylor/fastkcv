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
    list(fk, wide("engine_r"), wide("reference"), wide("krige.cv")))
  t1$`engine_r/fastkcv` <- round(t1$engine_r / t1$fastkcv, 1)
  t1$`reference/fastkcv` <- round(t1$reference / t1$fastkcv, 1)
  t1$`krige.cv/fastkcv` <- round(t1$krige.cv / t1$fastkcv, 1)
  t1$`prepare share` <- round(t1$prepare / t1$fastkcv, 2)
  t1 <- t1[order(t1$blas, t1$data != "sim", t1$n,
    match(sub(" .*", "", t1$scheme), c("LOO", "random", "blocks")),
    -as.numeric(sub("^\\D*", "", sub("LOO", "1e9", t1$scheme)))), ]
  for (col in c("fastkcv", "prepare", "run", "engine_r", "reference",
    "krige.cv")) t1[[col]] <- sig(t1[[col]])
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
      "prepare share")]
  emit("Where the advantage shrinks (simulated data, random folds)", adv,
    "Speed-up of fastkcv over the per-fold refits by n and number of folds.")
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
  dk <- dk[, c("blas", "n", "covtype", "setting", "status", "fastkcv_s",
    "dice_s", "km_s", "cov_s", "speedup", "speedup_total",
    "max_abs_diff_pred", "max_abs_diff_var")]
  for (col in c("fastkcv_s", "dice_s", "km_s", "cov_s")) dk[[col]] <- sig(dk[[col]])
  for (col in c("speedup", "speedup_total")) dk[[col]] <- round(dk[[col]], 2)
  for (col in c("max_abs_diff_pred", "max_abs_diff_var")) {
    dk[[col]] <- formatC(dk[[col]], format = "e", digits = 1)
  }
  emit("No-nugget LOO vs DiceKriging::leaveOneOut.km", dk,
    paste("speedup = leaveOneOut.km / kcv_engine (both from a ready",
      "covariance; km() has already factorized it for DiceKriging);",
      "speedup_total = (km + leaveOneOut.km) / (covariance + kcv_engine)."))
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
  a <- a[, c("blas", "n", "radius", "n_targets", "mean_m", "cost_ratio",
    "pairs_s", "points_s", "auto", "faster", "loss")]
  for (col in c("pairs_s", "points_s", "mean_m")) a[[col]] <- sig(a[[col]])
  emit("Automatic pairs/points choice", a,
    "cost_ratio = cost_points / cost_pairs (auto picks pairs when > 1); loss = time of auto's pick / faster - 1.")
}

if (!is.na(out_file)) {
  writeLines(c("# Benchmark summary", lines), out_file)
  message("wrote ", out_file)
}
