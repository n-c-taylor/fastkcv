# ==============================================================================
# Extended Comparison: Universal Kriging & Buffered CV
# fastkcv::kcv vs. gstat::krige.cv / gstat::krige
# ==============================================================================

library(fastkcv)
library(gstat)
library(sp)
library(bench)

# 1. Data & Residual Variogram for Universal Kriging ---------------------------
data(meuse)
coordinates(meuse) <- ~ x + y

uk_formula <- log(zinc) ~ sqrt(dist)

# Fit residual variogram (sill is smaller once sqrt(dist) trend is removed)
v_emp_uk <- variogram(uk_formula, meuse)
v_mod_uk <- fit.variogram(
  v_emp_uk,
  vgm(psill = 0.2, model = "Sph", range = 900, nugget = 0.05)
)

# Shared folds for 10-fold CV
set.seed(42)
folds_k <- sample(rep(seq_len(10), length.out = nrow(meuse)))

# Helper to compare numerical agreement
compare_results <- function(res_gstat, res_fast, obs, label) {
  pred_g <- res_gstat$var1.pred
  var_g <- res_gstat$var1.var

  df_f <- as.data.frame(res_fast)
  pred_f <- df_f[[intersect(c("var1.pred", "pred", "predicted"), names(df_f))[1]]]
  var_f <- df_f[[intersect(c("var1.var", "var", "variance"), names(df_f))[1]]]

  rmse <- function(o, p) sqrt(mean((o - p)^2, na.rm = TRUE))

  cat(sprintf("=== Numerical Agreement: %s ===\n", label))
  cat("Max abs diff (pred):", max(abs(pred_g - pred_f), na.rm = TRUE), "\n")
  cat("Max abs diff (var): ", max(abs(var_g - var_f), na.rm = TRUE), "\n")
  cat("Correlation (pred): ", cor(pred_g, pred_f, use = "complete.obs"), "\n")
  cat("RMSE (gstat):       ", rmse(obs, pred_g), "\n")
  cat("RMSE (fastkcv):     ", rmse(obs, pred_f), "\n\n")
}

# ==============================================================================
# PART 1: Universal Kriging (LOO-CV & 10-Fold CV)
# ==============================================================================

run_gstat_uk_loo <- function() {
  krige.cv(uk_formula,
    locations = meuse, model = v_mod_uk,
    nfold = nrow(meuse), verbose = FALSE
  )
}

run_fastkcv_uk_loo <- function() {
  fastkcv::kcv(uk_formula,
    locations = meuse, model = v_mod_uk,
    nfold = nrow(meuse)
  )
}

run_gstat_uk_kfold <- function() {
  krige.cv(uk_formula,
    locations = meuse, model = v_mod_uk,
    nfold = folds_k, verbose = FALSE
  )
}

run_fastkcv_uk_kfold <- function() {
  fastkcv::kcv(uk_formula,
    locations = meuse, model = v_mod_uk,
    folds = folds_k
  )
}

obs <- log(meuse$zinc)
compare_results(run_gstat_uk_loo(), run_fastkcv_uk_loo(), obs, "UK LOO-CV")
compare_results(run_gstat_uk_kfold(), run_fastkcv_uk_kfold(), obs, "UK 10-Fold CV")

cat("=== Performance Benchmark: Universal Kriging (meuse, n = 155) ===\n")
bm_uk <- bench::mark(
  gstat_uk_loo      = run_gstat_uk_loo(),
  fastkcv_uk_loo    = run_fastkcv_uk_loo(),
  gstat_uk_10fold   = run_gstat_uk_kfold(),
  fastkcv_uk_10fold = run_fastkcv_uk_kfold(),
  iterations        = 10,
  check             = FALSE
)
print(bm_uk[, c("expression", "min", "median", "itr/sec", "mem_alloc")])

# ==============================================================================
# PART 2: Buffered LOO-CV (buffer = 300 m)
# ==============================================================================
# gstat::krige.cv() has no buffer parameter, so we exclude points within
# `buffer_dist` manually using gstat::krige() in a loop.

buffer_dist <- 300 # meters

run_gstat_buffered_loo <- function() {
  d_mat <- sp::spDists(meuse)
  n <- nrow(meuse)
  pred <- numeric(n)
  var <- numeric(n)

  for (i in seq_len(n)) {
    train_idx <- which(d_mat[i, ] > buffer_dist)
    res <- gstat::krige(
      formula     = uk_formula,
      locations   = meuse[train_idx, ],
      newdata     = meuse[i, ],
      model       = v_mod_uk,
      debug.level = 0
    )
    pred[i] <- res$var1.pred
    var[i] <- res$var1.var
  }

  data.frame(var1.pred = pred, var1.var = var)
}

run_fastkcv_buffered_loo <- function() {
  fastkcv::kcv(
    formula   = uk_formula,
    locations = meuse,
    model     = v_mod_uk,
    buffer    = buffer_dist
  )
}

cat("\n")
compare_results(
  run_gstat_buffered_loo(),
  run_fastkcv_buffered_loo(),
  obs,
  sprintf("UK Buffered LOO-CV (buffer = %d m)", buffer_dist)
)

cat(sprintf("=== Performance Benchmark: UK Buffered LOO-CV (buffer = %d m) ===\n", buffer_dist))
bm_buf <- bench::mark(
  gstat_buffered_loo   = run_gstat_buffered_loo(),
  fastkcv_buffered_loo = run_fastkcv_buffered_loo(),
  iterations           = 5,
  check                = FALSE
)
print(bm_buf[, c("expression", "min", "median", "itr/sec", "mem_alloc")])
