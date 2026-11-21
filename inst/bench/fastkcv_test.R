# ==============================================================================
# Benchmark & Accuracy Comparison: fastkcv::kcv vs. gstat::krige.cv
# ==============================================================================

library(fastkcv)
library(gstat)
library(sp)
library(bench)

# 1. Data & Variogram Model ----------------------------------------------------
data(meuse)
coordinates(meuse) <- ~ x + y

v_emp <- variogram(log(zinc) ~ 1, meuse)
v_mod <- fit.variogram(v_emp, vgm(psill = 0.5, model = "Sph", range = 900, nugget = 0.05))

# 2. Shared Folds for 10-Fold CV -----------------------------------------------
set.seed(42)
k <- 10
folds_k <- sample(rep(seq_len(k), length.out = nrow(meuse)))

# 3. Wrapper Functions ---------------------------------------------------------
run_gstat_loo <- function() {
  krige.cv(
    formula   = log(zinc) ~ 1,
    locations = meuse,
    model     = v_mod,
    nfold     = nrow(meuse),
    verbose   = FALSE
  )
}

run_fastkcv_loo <- function() {
  fastkcv::kcv(
    formula   = log(zinc) ~ 1,
    locations = meuse,
    model     = v_mod,
    nfold     = nrow(meuse)
  )
}

run_gstat_kfold <- function() {
  krige.cv(
    formula   = log(zinc) ~ 1,
    locations = meuse,
    model     = v_mod,
    nfold     = folds_k,
    verbose   = FALSE
  )
}

run_fastkcv_kfold <- function() {
  fastkcv::kcv(
    formula   = log(zinc) ~ 1,
    locations = meuse,
    model     = v_mod,
    folds     = folds_k
  )
}

# 4. Numerical Agreement Checks (LOO-CV & 10-Fold CV) --------------------------
compare_results <- function(res_gstat, res_fast, label) {
  # Extract predictions & variances (handles common column naming conventions)
  pred_g <- res_gstat$var1.pred
  var_g <- res_gstat$var1.var

  df_f <- as.data.frame(res_fast)
  pred_f <- df_f[[intersect(c("var1.pred", "pred", "predicted"), names(df_f))[1]]]
  var_f <- df_f[[intersect(c("var1.var", "var", "variance"), names(df_f))[1]]]

  obs <- log(meuse$zinc)
  rmse <- function(o, p) sqrt(mean((o - p)^2, na.rm = TRUE))

  cat(sprintf("=== Numerical Agreement (%s) ===\n", label))
  cat("Max abs diff (pred):", max(abs(pred_g - pred_f), na.rm = TRUE), "\n")
  cat("Max abs diff (var): ", max(abs(var_g - var_f), na.rm = TRUE), "\n")
  cat("Correlation (pred): ", cor(pred_g, pred_f, use = "complete.obs"), "\n")
  cat("RMSE (gstat):       ", rmse(obs, pred_g), "\n")
  cat("RMSE (fastkcv):     ", rmse(obs, pred_f), "\n\n")
}

compare_results(run_gstat_loo(), run_fastkcv_loo(), "LOO-CV")
compare_results(run_gstat_kfold(), run_fastkcv_kfold(), "10-Fold CV")

# 5. Speed & Memory Benchmark (meuse, n = 155) ---------------------------------
cat("=== Performance Benchmark (meuse, n = 155) ===\n")
bm_meuse <- bench::mark(
  gstat_loo      = run_gstat_loo(),
  fastkcv_loo    = run_fastkcv_loo(),
  gstat_10fold   = run_gstat_kfold(),
  fastkcv_10fold = run_fastkcv_kfold(),
  iterations     = 10,
  check          = FALSE
)
print(bm_meuse[, c("expression", "min", "median", "itr/sec", "mem_alloc")])

# 6. Scaling Benchmark (Synthetic, n = 1,000) ----------------------------------
set.seed(42)
n_large <- 1000
sim_sp <- data.frame(
  x = runif(n_large, 0, 1000),
  y = runif(n_large, 0, 1000),
  z = rnorm(n_large, mean = 10, sd = 2)
)
coordinates(sim_sp) <- ~ x + y
v_sim <- vgm(psill = 3, model = "Exp", range = 250, nugget = 0.5)

cat("\n=== Scaling Benchmark (Synthetic LOO-CV, n =", n_large, ") ===\n")
bm_large <- bench::mark(
  gstat_loo   = krige.cv(z ~ 1, locations = sim_sp, model = v_sim, nfold = n_large, verbose = FALSE),
  fastkcv_loo = fastkcv::kcv(z ~ 1, locations = sim_sp, model = v_sim, nfold = n_large),
  iterations  = 3,
  check       = FALSE
)
print(bm_large[, c("expression", "min", "median", "itr/sec", "mem_alloc")])
