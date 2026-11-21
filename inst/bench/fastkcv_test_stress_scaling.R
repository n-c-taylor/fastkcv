# ==============================================================================
# Stress & Scaling Comparison: fastkcv::kcv vs. gstat::krige.cv
# Covers:
#   0. Environment & BLAS/LAPACK metadata
#   1. Variogram families, ill-conditioning, UTM-scale UK (raw & scaled),
#      duplicates, SK, sf, unequal/unsorted/character folds, rank-deficient folds
#   2. Fold-count sweep at fixed n (K = 2 .. n)
#   3. Sample-size scaling (10-fold UK CV) vs. global & local (nmax = 50) gstat
#   4. Isolated subprocess peak RSS measurement (C-level + R-heap memory)
# ==============================================================================

library(fastkcv)
library(gstat)
library(sp)
library(bench)

# ==============================================================================
# PART 0: Environment & Linear Algebra Metadata
# ==============================================================================

fmt_lib <- function(x) if (length(x) == 1 && nzchar(x)) x else "<bundled/default>"

cat("=== Environment & BLAS/LAPACK Metadata ===\n")
cat("R version:   ", R.version.string, "\n")
cat("Platform:    ", R.version$platform, "\n")
cat("BLAS:        ", fmt_lib(sessionInfo()$BLAS), "\n")
cat("LAPACK:      ", fmt_lib(La_library()), "\n")
cat("fastkcv:     ", as.character(packageVersion("fastkcv")), "\n")
cat("gstat:       ", as.character(packageVersion("gstat")), "\n")
thread_vars <- Sys.getenv(
  c("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"),
  unset = "<unset>"
)
cat("Threads:     ", paste(names(thread_vars), thread_vars, sep = "=", collapse = ", "), "\n\n")

# ==============================================================================
# Helpers
# ==============================================================================

make_spatial_data <- function(n,
                              seed = 123,
                              origin = c(0, 0),
                              scale = 1,
                              add_close_pair = FALSE,
                              add_duplicate = FALSE) {
  set.seed(seed)
  u <- cbind(x = runif(n, 0, 1), y = runif(n, 0, 1))
  if (add_close_pair && n >= 2) {
    u[2, ] <- u[1, ] + c(1e-4, 1e-4)
  }
  if (add_duplicate && n >= 2) {
    u[2, ] <- u[1, ]
  }

  d_unit <- as.matrix(dist(u))
  cov_mat <- 1.0 * exp(-d_unit / 0.3) + diag(0.1, n)
  eps <- as.numeric(t(chol(cov_mat)) %*% rnorm(n))

  df <- data.frame(
    x = origin[1] + u[, "x"] * scale,
    y = origin[2] + u[, "y"] * scale,
    z = 2 + 1.5 * u[, "x"] - 0.8 * u[, "y"] + 0.5 * u[, "x"]^2 + eps
  )
  coordinates(df) <- ~ x + y
  df
}

extract_cv <- function(x) {
  df <- as.data.frame(x)
  pred_col <- intersect(c("var1.pred", "pred", "predicted"), names(df))[1]
  var_col <- intersect(c("var1.var", "var", "variance"), names(df))[1]
  if (is.na(pred_col) || is.na(var_col)) {
    stop("Could not locate prediction/variance columns in output.")
  }
  list(pred = as.numeric(df[[pred_col]]), var = as.numeric(df[[var_col]]))
}

safe_max <- function(x) {
  if (length(x) == 0 || all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)
}

expr_labels <- function(x) {
  nms <- names(x)
  if (!is.null(nms) && all(nzchar(nms))) nms else as.character(x)
}

run_comparison <- function(label, expr_g, expr_f, tol_pred = 1e-10, tol_var = 1e-10, core = TRUE) {
  res_g <- tryCatch(eval.parent(substitute(expr_g)), error = function(e) e)
  res_f <- tryCatch(eval.parent(substitute(expr_f)), error = function(e) e)

  err_g <- inherits(res_g, "error")
  err_f <- inherits(res_f, "error")

  if (err_g || err_f) {
    status <- if (err_g && err_f) "BOTH_ERR" else if (err_g) "GSTAT_ERR" else "FASTKCV_ERR"
    msg <- paste(
      if (err_g) paste0("gstat: ", conditionMessage(res_g)),
      if (err_f) paste0("fastkcv: ", conditionMessage(res_f)),
      sep = " | "
    )
    return(data.frame(
      case = label,
      core = core,
      max_abs_pred = NA_real_,
      max_rel_pred = NA_real_,
      max_abs_var = NA_real_,
      na_parity = NA,
      neg_var_f = NA_integer_,
      status = status,
      note = trimws(gsub("^\\|\\s*|\\s*\\|$", "", msg)),
      stringsAsFactors = FALSE
    ))
  }

  g <- extract_cv(res_g)
  f <- extract_cv(res_f)

  if (length(g$pred) != length(f$pred)) {
    return(data.frame(
      case = label, core = core,
      max_abs_pred = NA_real_, max_rel_pred = NA_real_, max_abs_var = NA_real_,
      na_parity = FALSE, neg_var_f = NA_integer_,
      status = "FAIL", note = "Output length mismatch",
      stringsAsFactors = FALSE
    ))
  }

  na_parity <- identical(is.na(g$pred), is.na(f$pred)) &&
    identical(is.na(g$var), is.na(f$var))

  diff_pred <- abs(g$pred - f$pred)
  diff_var <- abs(g$var - f$var)

  max_abs_pred <- safe_max(diff_pred)
  max_rel_pred <- safe_max(diff_pred / pmax(abs(g$pred), 1e-12))
  max_abs_var <- safe_max(diff_var)
  neg_var_f <- sum(f$var < -1e-14, na.rm = TRUE)

  within_tol <- (is.na(max_abs_pred) || max_abs_pred <= tol_pred) &&
    (is.na(max_abs_var) || max_abs_var <= tol_var)

  status <- if (na_parity && within_tol) "PASS" else "FAIL"

  data.frame(
    case = label,
    core = core,
    max_abs_pred = max_abs_pred,
    max_rel_pred = max_rel_pred,
    max_abs_var = max_abs_var,
    na_parity = na_parity,
    neg_var_f = neg_var_f,
    status = status,
    note = "",
    stringsAsFactors = FALSE
  )
}

# ==============================================================================
# PART 1: Numerical Stability, Geostatistical Features & Edge Cases
# ==============================================================================

dat_base <- make_spatial_data(n = 150, seed = 42)
dat_close <- make_spatial_data(n = 150, seed = 43, add_close_pair = TRUE)
dat_dup <- make_spatial_data(n = 150, seed = 44, add_duplicate = TRUE)
dat_utm <- make_spatial_data(n = 150, seed = 45, origin = c(500000, 4500000), scale = 1e5)
dat_utm$xs <- as.numeric(scale(coordinates(dat_utm)[, 1]))
dat_utm$ys <- as.numeric(scale(coordinates(dat_utm)[, 2]))

v_exp <- vgm(psill = 1.0, model = "Exp", range = 0.3, nugget = 0.1)
v_mat <- vgm(psill = 1.0, model = "Mat", range = 0.3, nugget = 0.1, kappa = 1.5)
v_zero <- vgm(psill = 1.0, model = "Exp", range = 0.3, nugget = 0.0)
v_gau_ill <- vgm(psill = 1.0, model = "Gau", range = 0.3, nugget = 1e-6)
v_nested <- vgm(
  psill = 0.6, model = "Exp", range = 0.2, nugget = 0.1,
  add.to = vgm(psill = 0.4, model = "Sph", range = 0.6)
)
v_anis <- vgm(psill = 1.0, model = "Exp", range = 0.3, nugget = 0.1, anis = c(45, 0.5))
v_utm <- vgm(psill = 1.0, model = "Exp", range = 30000, nugget = 0.1)

uk_quad <- z ~ x + y + I(x^2) + I(y^2) + I(x * y)
uk_quad_std <- z ~ xs + ys + I(xs^2) + I(ys^2) + I(xs * ys)

# Unequal spatial block folds + non-contiguous & character fold labels
set.seed(42)
blk_int <- kmeans(coordinates(dat_base), centers = 5)$cluster
blk_nonseq <- c(10L, 2L, 99L, 4L, 50L)[blk_int]
blk_char <- paste0("zone_", LETTERS[blk_int])

# Rank-deficient fold: 2 training points left for 3-parameter linear UK
dat_small <- make_spatial_data(n = 20, seed = 46)
folds_degen <- c(rep(1L, 18), rep(2L, 2))

rows <- list(
  # Core variogram families & zero-nugget
  run_comparison(
    "OK LOO: Exp (nugget)",
    krige.cv(z ~ 1, dat_base, model = v_exp, nfold = nrow(dat_base), verbose = FALSE),
    fastkcv::kcv(z ~ 1, dat_base, model = v_exp, nfold = nrow(dat_base))
  ),
  run_comparison(
    "OK LOO: Mat (kappa=1.5)",
    krige.cv(z ~ 1, dat_base, model = v_mat, nfold = nrow(dat_base), verbose = FALSE),
    fastkcv::kcv(z ~ 1, dat_base, model = v_mat, nfold = nrow(dat_base))
  ),
  run_comparison(
    "OK LOO: Exp (zero nugget)",
    krige.cv(z ~ 1, dat_base, model = v_zero, nfold = nrow(dat_base), verbose = FALSE),
    fastkcv::kcv(z ~ 1, dat_base, model = v_zero, nfold = nrow(dat_base))
  ),
  # Ill-conditioned Gaussian variogram + close point pair
  run_comparison(
    "OK LOO: Gau (nug=1e-6, close pair)",
    krige.cv(z ~ 1, dat_close, model = v_gau_ill, nfold = nrow(dat_close), verbose = FALSE),
    fastkcv::kcv(z ~ 1, dat_close, model = v_gau_ill, nfold = nrow(dat_close)),
    tol_pred = 1e-6, tol_var = 1e-6
  ),
  # Quadratic UK & realistic UTM coordinate magnitudes (raw diagnostic vs. scaled core)
  run_comparison(
    "UK Quadratic LOO (unit square)",
    krige.cv(uk_quad, dat_base, model = v_exp, nfold = nrow(dat_base), verbose = FALSE),
    fastkcv::kcv(uk_quad, dat_base, model = v_exp, nfold = nrow(dat_base))
  ),
  run_comparison(
    "UK Quadratic LOO (UTM raw, unscaled)",
    krige.cv(uk_quad, dat_utm, model = v_utm, nfold = nrow(dat_utm), verbose = FALSE),
    fastkcv::kcv(uk_quad, dat_utm, model = v_utm, nfold = nrow(dat_utm)),
    tol_pred = 1e-6, tol_var = 1e-6, core = FALSE
  ),
  run_comparison(
    "UK Quadratic LOO (UTM centered/scaled)",
    krige.cv(uk_quad_std, dat_utm, model = v_utm, nfold = nrow(dat_utm), verbose = FALSE),
    fastkcv::kcv(uk_quad_std, dat_utm, model = v_utm, nfold = nrow(dat_utm))
  ),
  # Spatial block folds & fold-label handling
  run_comparison(
    "UK Quad Block-5 (contiguous int)",
    krige.cv(uk_quad, dat_base, model = v_exp, nfold = blk_int, verbose = FALSE),
    fastkcv::kcv(uk_quad, dat_base, model = v_exp, folds = blk_int)
  ),
  run_comparison(
    "UK Quad Block-5 (non-seq int IDs)",
    krige.cv(uk_quad, dat_base, model = v_exp, nfold = blk_nonseq, verbose = FALSE),
    fastkcv::kcv(uk_quad, dat_base, model = v_exp, folds = blk_nonseq),
    core = FALSE
  ),
  run_comparison(
    "UK Quad Block-5 (character IDs)",
    krige.cv(uk_quad, dat_base, model = v_exp, nfold = blk_char, verbose = FALSE),
    fastkcv::kcv(uk_quad, dat_base, model = v_exp, folds = blk_char),
    core = FALSE
  ),
  # Extended geostatistical features & edge cases (non-core: test support cleanly)
  run_comparison(
    "OK LOO: Duplicate coordinates",
    krige.cv(z ~ 1, dat_dup, model = v_exp, nfold = nrow(dat_dup), verbose = FALSE),
    fastkcv::kcv(z ~ 1, dat_dup, model = v_exp, nfold = nrow(dat_dup)),
    core = FALSE
  ),
  run_comparison(
    "OK LOO: Nested vgm (Exp + Sph)",
    krige.cv(z ~ 1, dat_base, model = v_nested, nfold = nrow(dat_base), verbose = FALSE),
    fastkcv::kcv(z ~ 1, dat_base, model = v_nested, nfold = nrow(dat_base)),
    core = FALSE
  ),
  run_comparison(
    "OK LOO: Anisotropic vgm (45, 0.5)",
    krige.cv(z ~ 1, dat_base, model = v_anis, nfold = nrow(dat_base), verbose = FALSE),
    fastkcv::kcv(z ~ 1, dat_base, model = v_anis, nfold = nrow(dat_base)),
    core = FALSE
  ),
  run_comparison(
    "Simple Kriging LOO (known beta)",
    krige.cv(z ~ 1, dat_base, model = v_exp, beta = 2.0, nfold = nrow(dat_base), verbose = FALSE),
    fastkcv::kcv(z ~ 1, dat_base, model = v_exp, beta = 2.0, nfold = nrow(dat_base)),
    core = FALSE
  ),
  run_comparison(
    "UK Linear: Rank-deficient fold",
    krige.cv(z ~ x + y, dat_small, model = v_exp, nfold = folds_degen, verbose = FALSE),
    fastkcv::kcv(z ~ x + y, dat_small, model = v_exp, folds = folds_degen),
    core = FALSE
  )
)

if (requireNamespace("sf", quietly = TRUE)) {
  dat_sf <- sf::st_as_sf(dat_base)
  rows[[length(rows) + 1]] <- run_comparison(
    "OK LOO: sf input object",
    krige.cv(z ~ 1, dat_sf, model = v_exp, nfold = nrow(dat_sf), verbose = FALSE),
    fastkcv::kcv(z ~ 1, dat_sf, model = v_exp, nfold = nrow(dat_sf)),
    core = FALSE
  )
}

res_edge <- do.call(rbind, rows)
cat("=== PART 1: Numerical Stability & Feature Compatibility ===\n")
print(res_edge[, c("case", "max_abs_pred", "max_rel_pred", "max_abs_var", "na_parity", "neg_var_f", "status")],
  row.names = FALSE
)

err_notes <- res_edge[nzchar(res_edge$note), c("case", "status", "note")]
if (nrow(err_notes) > 0) {
  cat("\nEdge-case notes / unsupported features:\n")
  print(err_notes, row.names = FALSE)
}

core_failed <- res_edge[res_edge$core & res_edge$status != "PASS", ]
if (nrow(core_failed) > 0) {
  stop("Core numerical agreement check(s) failed: ", paste(core_failed$case, collapse = ", "))
}

# ==============================================================================
# PART 2: Fold-Count Sweep at Fixed n (n = 300, K = 2 .. 300)
# ==============================================================================

cat("\n=== PART 2: Fold-Count Sweep (n = 300, Linear UK) ===\n")
n_sweep <- 300
dat_sweep <- make_spatial_data(n_sweep, seed = 300)

# Verify agreement at K = 2 (largest Schur-complement blocks) before timing
set.seed(42)
f_k2 <- sample(rep(1:2, length.out = n_sweep))
chk_k2 <- run_comparison(
  "Sweep check: n=300, K=2",
  krige.cv(z ~ x + y, dat_sweep, model = v_exp, nfold = f_k2, verbose = FALSE),
  fastkcv::kcv(z ~ x + y, dat_sweep, model = v_exp, folds = f_k2)
)
stopifnot(chk_k2$status == "PASS")

bm_fold_sweep <- bench::press(
  k = c(2, 5, 10, 25, 50, 100, 300),
  {
    set.seed(42)
    f_k <- if (k == n_sweep) seq_len(n_sweep) else sample(rep(seq_len(k), length.out = n_sweep))
    gc(FALSE)

    bench::mark(
      gstat = krige.cv(z ~ x + y, dat_sweep, model = v_exp, nfold = f_k, verbose = FALSE),
      fastkcv = if (k == n_sweep) {
        fastkcv::kcv(z ~ x + y, dat_sweep, model = v_exp, nfold = n_sweep)
      } else {
        fastkcv::kcv(z ~ x + y, dat_sweep, model = v_exp, folds = f_k)
      },
      min_time = 0.5,
      min_iterations = 2,
      check = FALSE
    )
  }
)

sweep_tbl <- bm_fold_sweep[, c("expression", "k", "min", "median", "itr/sec", "mem_alloc", "n_itr")]
print(sweep_tbl)

# Speedup summary table
lbl_sweep <- expr_labels(bm_fold_sweep$expression)
g_sw <- bm_fold_sweep[lbl_sweep == "gstat", ]
f_sw <- bm_fold_sweep[lbl_sweep == "fastkcv", ]
idx_sw <- match(g_sw$k, f_sw$k)

cat("\nFold-Sweep Speedup (gstat / fastkcv):\n")
print(data.frame(
  k = g_sw$k,
  gstat_ms = round(as.numeric(g_sw$median) * 1e3, 2),
  fastkcv_ms = round(as.numeric(f_sw$median[idx_sw]) * 1e3, 2),
  speedup = round(as.numeric(g_sw$median) / as.numeric(f_sw$median[idx_sw]), 2)
), row.names = FALSE)

# ==============================================================================
# PART 3: Sample-Size Scaling — 10-Fold UK CV (Global vs. Local nmax = 50)
# ==============================================================================

cat("\n=== PART 3: Sample-Size Scaling (10-Fold UK CV, Global vs. Local nmax=50) ===\n")
n_grid <- c(100, 250, 500, 1000, 1500)

bm_scale <- bench::press(
  n = n_grid,
  {
    dat <- make_spatial_data(n, seed = 500 + n)
    set.seed(42)
    f_k <- sample(rep(seq_len(10), length.out = n))
    gc(FALSE)

    bench::mark(
      gstat_global = krige.cv(z ~ x + y, dat, model = v_exp, nfold = f_k, verbose = FALSE),
      gstat_local = krige.cv(z ~ x + y, dat, model = v_exp, nmax = 50, nfold = f_k, verbose = FALSE),
      fastkcv = fastkcv::kcv(z ~ x + y, dat, model = v_exp, folds = f_k),
      min_time = 0.5,
      min_iterations = 2,
      check = FALSE
    )
  }
)
print(bm_scale[, c("expression", "n", "min", "median", "itr/sec", "mem_alloc", "n_itr")])

lbl_scale <- expr_labels(bm_scale$expression)
sg <- bm_scale[lbl_scale == "gstat_global", ]
sl <- bm_scale[lbl_scale == "gstat_local", ]
sf <- bm_scale[lbl_scale == "fastkcv", ]
idx_l <- match(sg$n, sl$n)
idx_f <- match(sg$n, sf$n)

cat("\nSample-Size Scaling Speedup (10-Fold UK CV):\n")
print(data.frame(
  n = sg$n,
  gstat_global_ms = round(as.numeric(sg$median) * 1e3, 2),
  gstat_local_ms = round(as.numeric(sl$median[idx_l]) * 1e3, 2),
  fastkcv_ms = round(as.numeric(sf$median[idx_f]) * 1e3, 2),
  vs_global_x = round(as.numeric(sg$median) / as.numeric(sf$median[idx_f]), 2),
  vs_local_x = round(as.numeric(sl$median[idx_l]) / as.numeric(sf$median[idx_f]), 2)
), row.names = FALSE)

# Numerical verification + global-vs-local accuracy check at largest n
n_max <- max(n_grid)
dat_max <- make_spatial_data(n_max, seed = 500 + n_max)
set.seed(42)
f_max <- sample(rep(seq_len(10), length.out = n_max))

res_g_max <- krige.cv(z ~ x + y, dat_max, model = v_exp, nfold = f_max, verbose = FALSE)
res_l_max <- krige.cv(z ~ x + y, dat_max, model = v_exp, nmax = 50, nfold = f_max, verbose = FALSE)
res_f_max <- fastkcv::kcv(z ~ x + y, dat_max, model = v_exp, folds = f_max)

chk_max <- run_comparison(
  sprintf("Max-n check: n=%d, 10-fold UK", n_max),
  res_g_max, res_f_max
)
rmse <- function(obs, pred) sqrt(mean((obs - pred)^2, na.rm = TRUE))

cat(sprintf("\nVerification at n = %d (10-fold UK):\n", n_max))
print(chk_max[, c("case", "max_abs_pred", "max_rel_pred", "max_abs_var", "na_parity", "status")],
  row.names = FALSE
)
cat(sprintf(
  "RMSE — gstat_global: %.6f | fastkcv: %.6f | gstat_local (nmax=50): %.6f\n",
  rmse(dat_max$z, extract_cv(res_g_max)$pred),
  rmse(dat_max$z, extract_cv(res_f_max)$pred),
  rmse(dat_max$z, extract_cv(res_l_max)$pred)
))
stopifnot(chk_max$status == "PASS")

# ==============================================================================
# PART 4: Isolated Subprocess Peak Process Memory (RSS)
# Captures C-level allocations invisible to bench::mark()$mem_alloc
# Pre-generates inputs in parent so data setup does not inflate subprocess HWM
# ==============================================================================

cat("\n=== PART 4: Isolated Subprocess Memory (Peak / Post-Call RSS) ===\n")

if (!requireNamespace("callr", quietly = TRUE)) {
  cat("Skipped: install 'callr' (and optionally 'ps' on macOS/Windows) for subprocess RSS measurement.\n")
} else {
  rss_ns <- c(750, 1500)
  rss_inputs <- setNames(lapply(rss_ns, function(n) {
    dat <- make_spatial_data(n, seed = 500 + n)
    set.seed(42)
    list(
      dat = dat,
      v_mod = vgm(psill = 1.0, model = "Exp", range = 0.3, nugget = 0.1),
      f_k = sample(rep(seq_len(10), length.out = n))
    )
  }), as.character(rss_ns))

  measure_subprocess_rss <- function(method, n) {
    inp <- rss_inputs[[as.character(n)]]
    callr::r(
      function(method, n, dat, v_mod, f_k) {
        suppressPackageStartupMessages({
          library(fastkcv)
          library(gstat)
          library(sp)
        })

        get_rss_mb <- function() {
          # Linux: exact kernel High-Water Mark (VmHWM)
          if (file.exists("/proc/self/status")) {
            st <- readLines("/proc/self/status", warn = FALSE)
            hwm <- grep("^VmHWM:", st, value = TRUE)
            if (length(hwm) == 1) {
              return(list(metric = "VmHWM_MB", value = as.numeric(gsub("[^0-9]", "", hwm)) / 1024))
            }
          }
          if (requireNamespace("ps", quietly = TRUE)) {
            mi <- ps::ps_memory_info()
            if ("peak_wset" %in% names(mi)) {
              return(list(metric = "PeakWSet_MB", value = unname(mi["peak_wset"]) / 1024^2))
            }
            return(list(metric = "RSS_MB", value = unname(mi["rss"]) / 1024^2))
          }
          list(metric = "NA", value = NA_real_)
        }

        gc(FALSE)
        before <- get_rss_mb()
        t0 <- proc.time()[["elapsed"]]
        if (method == "gstat_global") {
          out <- krige.cv(z ~ x + y, dat, model = v_mod, nfold = f_k, verbose = FALSE)
        } else if (method == "gstat_local") {
          out <- krige.cv(z ~ x + y, dat, model = v_mod, nmax = 50, nfold = f_k, verbose = FALSE)
        } else {
          out <- fastkcv::kcv(z ~ x + y, dat, model = v_mod, folds = f_k)
        }
        elapsed <- proc.time()[["elapsed"]] - t0
        after <- get_rss_mb()

        data.frame(
          method = method,
          n = n,
          elapsed_s = round(elapsed, 4),
          metric = after$metric,
          before_mb = round(before$value, 2),
          after_mb = round(after$value, 2),
          delta_mb = round(after$value - before$value, 2),
          stringsAsFactors = FALSE
        )
      },
      args = list(method = method, n = n, dat = inp$dat, v_mod = inp$v_mod, f_k = inp$f_k)
    )
  }

  rss_grid <- expand.grid(
    method = c("gstat_global", "gstat_local", "fastkcv"),
    n = rss_ns,
    stringsAsFactors = FALSE
  )
  rss_results <- do.call(rbind, Map(measure_subprocess_rss, rss_grid$method, rss_grid$n))
  print(rss_results, row.names = FALSE)
}
