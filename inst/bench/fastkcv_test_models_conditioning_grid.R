# ==============================================================================
# Extended Diagnostics & Workflow Comparison: fastkcv::kcv vs. gstat::krige.cv
# File: fastkcv_test_models_conditioning_grid.R
# Covers:
#   1. Variogram model zoo (bounded, Matern-Stein, hole-effect, pure nugget,
#      intrinsic/unbounded models)
#   2. 3D kriging, categorical/interaction UK covariates, and NA handling
#   3. Coordinate-offset conditioning sweep (raw vs. centered/scaled UK)
#   4. End-to-end variogram hyperparameter grid search (LOO-CV & 10-fold CV)
# ==============================================================================

library(fastkcv)
library(gstat)
library(sp)
library(bench)

# ==============================================================================
# Helpers
# ==============================================================================

extract_cv <- function(x) {
  df <- as.data.frame(x)
  pred_col <- intersect(c("var1.pred", "pred", "predicted"), names(df))[1]
  var_col <- intersect(c("var1.var", "var", "variance"), names(df))[1]
  if (is.na(pred_col) || is.na(var_col)) {
    stop("Could not locate prediction/variance columns.")
  }
  list(pred = as.numeric(df[[pred_col]]), var = as.numeric(df[[var_col]]))
}

safe_max <- function(x) {
  if (length(x) == 0 || all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)
}

eval_capture <- function(expr, env) {
  warns <- character(0)
  val <- tryCatch(
    withCallingHandlers(
      eval(expr, envir = env),
      warning = function(w) {
        warns <<- c(warns, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) e
  )
  list(value = val, warnings = warns, is_err = inherits(val, "error"))
}

compare_pair <- function(label, expr_g, expr_f, tol_pred = 1e-10, tol_var = 1e-10) {
  env <- parent.frame()
  rg <- eval_capture(substitute(expr_g), env)
  rf <- eval_capture(substitute(expr_f), env)

  if (rg$is_err || rf$is_err) {
    status <- if (rg$is_err && rf$is_err) "BOTH_ERR" else if (rg$is_err) "GSTAT_ERR" else "FASTKCV_ERR"
    msg <- paste(
      if (rg$is_err) paste0("gstat: ", conditionMessage(rg$value)),
      if (rf$is_err) paste0("fastkcv: ", conditionMessage(rf$value)),
      sep = " | "
    )
    return(data.frame(
      case = label, max_abs_pred = NA_real_, max_abs_var = NA_real_,
      na_parity = NA, warn_g = length(rg$warnings), warn_f = length(rf$warnings),
      status = status, note = trimws(gsub("^\\|\\s*|\\s*\\|$", "", msg)),
      stringsAsFactors = FALSE
    ))
  }

  g <- extract_cv(rg$value)
  f <- extract_cv(rf$value)
  na_parity <- identical(is.na(g$pred), is.na(f$pred)) &&
    identical(is.na(g$var), is.na(f$var))

  max_abs_pred <- safe_max(abs(g$pred - f$pred))
  max_abs_var <- safe_max(abs(g$var - f$var))
  ok <- na_parity &&
    (is.na(max_abs_pred) || max_abs_pred <= tol_pred) &&
    (is.na(max_abs_var) || max_abs_var <= tol_var)

  data.frame(
    case = label,
    max_abs_pred = max_abs_pred,
    max_abs_var = max_abs_var,
    na_parity = na_parity,
    warn_g = length(rg$warnings),
    warn_f = length(rf$warnings),
    status = if (ok) "PASS" else "FAIL",
    note = "",
    stringsAsFactors = FALSE
  )
}

make_2d_data <- function(n = 120, seed = 101, offset = 0) {
  set.seed(seed)
  u <- cbind(x = runif(n, 0, 1), y = runif(n, 0, 1))
  d_mat <- as.matrix(dist(u))
  cov_mat <- 1.0 * exp(-d_mat / 0.3) + diag(0.1, n)
  eps <- as.numeric(t(chol(cov_mat)) %*% rnorm(n))

  df <- data.frame(
    x = u[, "x"] + offset,
    y = u[, "y"] + offset,
    xs = as.numeric(scale(u[, "x"] + offset)),
    ys = as.numeric(scale(u[, "y"] + offset)),
    z = 2 + 1.5 * u[, "x"] - 0.8 * u[, "y"] + 0.5 * u[, "x"]^2 + eps
  )
  coordinates(df) <- ~ x + y
  df
}

make_3d_data <- function(n = 120, seed = 202) {
  set.seed(seed)
  u <- cbind(x = runif(n, 0, 1), y = runif(n, 0, 1), z = runif(n, 0, 1))
  d_mat <- as.matrix(dist(u))
  cov_mat <- 1.0 * exp(-d_mat / 0.35) + diag(0.1, n)
  eps <- as.numeric(t(chol(cov_mat)) %*% rnorm(n))

  df <- data.frame(
    x = u[, "x"],
    y = u[, "y"],
    z = u[, "z"],
    val = 1.0 + 2.0 * u[, "x"] - 1.2 * u[, "y"] + 0.9 * u[, "z"] + eps
  )
  coordinates(df) <- ~ x + y + z
  df
}

# ==============================================================================
# PART 1: Variogram Model Zoo (Bounded, Pure Nugget & Unbounded/Intrinsic)
# ==============================================================================

dat_zoo <- make_2d_data(n = 100, seed = 101)

vgm_zoo <- list(
  Nug_pure      = vgm(psill = 0.5, model = "Nug", range = 0),
  Sph           = vgm(psill = 1.0, model = "Sph", range = 0.4, nugget = 0.1),
  Cir           = vgm(psill = 1.0, model = "Cir", range = 0.4, nugget = 0.1),
  Pen           = vgm(psill = 1.0, model = "Pen", range = 0.4, nugget = 0.1),
  Ste_k2.5      = vgm(psill = 1.0, model = "Ste", range = 0.3, nugget = 0.1, kappa = 2.5),
  Bes           = vgm(psill = 1.0, model = "Bes", range = 0.3, nugget = 0.1),
  Wav           = vgm(psill = 1.0, model = "Wav", range = 0.25, nugget = 0.1),
  Hol           = vgm(psill = 1.0, model = "Hol", range = 0.25, nugget = 0.1),
  Lin_bounded   = vgm(psill = 1.0, model = "Lin", range = 0.5, nugget = 0.1),
  Lin_unbounded = vgm(psill = 1.0, model = "Lin", range = 0.0, nugget = 0.1),
  Pow_1.5       = vgm(psill = 1.0, model = "Pow", range = 1.5, nugget = 0.1),
  Log_intrinsic = vgm(psill = 1.0, model = "Log", range = 0.3, nugget = 0.1)
)

zoo_rows <- lapply(names(vgm_zoo), function(nm) {
  m <- vgm_zoo[[nm]]
  compare_pair(
    paste("OK LOO:", nm),
    krige.cv(z ~ 1, dat_zoo, model = m, nfold = nrow(dat_zoo), verbose = FALSE),
    fastkcv::kcv(z ~ 1, dat_zoo, model = m, nfold = nrow(dat_zoo))
  )
})
res_zoo <- do.call(rbind, zoo_rows)

cat("=== PART 1: Variogram Model Zoo (n = 100, OK LOO-CV) ===\n")
print(res_zoo[, c("case", "max_abs_pred", "max_abs_var", "na_parity", "warn_g", "warn_f", "status")],
  row.names = FALSE
)
zoo_notes <- res_zoo[nzchar(res_zoo$note), c("case", "status", "note")]
if (nrow(zoo_notes) > 0) {
  cat("\nVariogram zoo error notes:\n")
  print(zoo_notes, row.names = FALSE)
}

# ==============================================================================
# PART 2: 3D Kriging, Categorical/Interaction UK & Missing Values (NA)
# ==============================================================================

dat_3d <- make_3d_data(n = 120, seed = 202)
v_3d <- vgm(psill = 1.0, model = "Exp", range = 0.35, nugget = 0.1)
set.seed(42)
folds_3d <- sample(rep(seq_len(10), length.out = nrow(dat_3d)))

data(meuse)
coordinates(meuse) <- ~ x + y
v_meuse <- fit.variogram(
  variogram(log(zinc) ~ sqrt(dist), meuse),
  vgm(psill = 0.2, model = "Sph", range = 900, nugget = 0.05)
)
set.seed(42)
folds_meuse <- sample(rep(seq_len(10), length.out = nrow(meuse)))

# Leave-one-factor-level-out is rank-deficient for intercept + factor;
# test balanced 10-fold CV and LOO-CV with factors & interactions
part2_rows <- list(
  compare_pair(
    "3D OK LOO-CV (x, y, z)",
    krige.cv(val ~ 1, dat_3d, model = v_3d, nfold = nrow(dat_3d), verbose = FALSE),
    fastkcv::kcv(val ~ 1, dat_3d, model = v_3d, nfold = nrow(dat_3d))
  ),
  compare_pair(
    "3D UK 10-Fold CV (val ~ x + y + z)",
    krige.cv(val ~ x + y + z, dat_3d, model = v_3d, nfold = folds_3d, verbose = FALSE),
    fastkcv::kcv(val ~ x + y + z, dat_3d, model = v_3d, folds = folds_3d)
  ),
  compare_pair(
    "meuse UK LOO: factors (ffreq + soil)",
    krige.cv(log(zinc) ~ sqrt(dist) + ffreq + soil, meuse,
      model = v_meuse,
      nfold = nrow(meuse), verbose = FALSE
    ),
    fastkcv::kcv(log(zinc) ~ sqrt(dist) + ffreq + soil, meuse,
      model = v_meuse,
      nfold = nrow(meuse)
    )
  ),
  compare_pair(
    "meuse UK 10-Fold: interaction (sqrt(dist) * ffreq)",
    krige.cv(log(zinc) ~ sqrt(dist) * ffreq, meuse,
      model = v_meuse,
      nfold = folds_meuse, verbose = FALSE
    ),
    fastkcv::kcv(log(zinc) ~ sqrt(dist) * ffreq, meuse,
      model = v_meuse,
      folds = folds_meuse
    )
  ),
  compare_pair(
    "meuse OK LOO: NA in response (om ~ 1)",
    krige.cv(om ~ 1, meuse, model = v_meuse, nfold = nrow(meuse), verbose = FALSE),
    fastkcv::kcv(om ~ 1, meuse, model = v_meuse, nfold = nrow(meuse))
  ),
  compare_pair(
    "meuse UK LOO: NA in covariate (log(zinc) ~ om)",
    krige.cv(log(zinc) ~ om, meuse, model = v_meuse, nfold = nrow(meuse), verbose = FALSE),
    fastkcv::kcv(log(zinc) ~ om, meuse, model = v_meuse, nfold = nrow(meuse))
  )
)
res_part2 <- do.call(rbind, part2_rows)

cat("\n=== PART 2: 3D Kriging, Factors/Interactions & NA Handling ===\n")
print(res_part2[, c("case", "max_abs_pred", "max_abs_var", "na_parity", "warn_g", "warn_f", "status")],
  row.names = FALSE
)
p2_notes <- res_part2[nzchar(res_part2$note), c("case", "status", "note")]
if (nrow(p2_notes) > 0) {
  cat("\nPart 2 error notes:\n")
  print(p2_notes, row.names = FALSE)
}

# ==============================================================================
# PART 3: Coordinate-Offset Conditioning Sweep (Quadratic UK LOO-CV)
# Compares both solvers against the exact offset = 0 reference solution
# ==============================================================================

v_exp <- vgm(psill = 1.0, model = "Exp", range = 0.3, nugget = 0.1)
f_raw <- z ~ x + y + I(x^2) + I(y^2) + I(x * y)
f_std <- z ~ xs + ys + I(xs^2) + I(ys^2) + I(xs * ys)

dat_ref <- make_2d_data(n = 120, seed = 45, offset = 0)
ref_pred <- extract_cv(
  krige.cv(f_raw, dat_ref, model = v_exp, nfold = nrow(dat_ref), verbose = FALSE)
)$pred

offsets <- c(0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6)

cond_rows <- lapply(offsets, function(off) {
  dat_off <- make_2d_data(n = 120, seed = 45, offset = off)
  df_off <- as.data.frame(dat_off)
  X_raw <- model.matrix(f_raw, data = df_off)
  kappa_raw <- kappa(X_raw, exact = TRUE)

  rg_raw <- eval_capture(
    quote(krige.cv(f_raw, dat_off, model = v_exp, nfold = nrow(dat_off), verbose = FALSE)),
    environment()
  )
  rf_raw <- eval_capture(
    quote(fastkcv::kcv(f_raw, dat_off, model = v_exp, nfold = nrow(dat_off))),
    environment()
  )
  rf_std <- eval_capture(
    quote(fastkcv::kcv(f_std, dat_off, model = v_exp, nfold = nrow(dat_off))),
    environment()
  )

  err_vs_ref <- function(res) {
    if (res$is_err) NA_real_ else safe_max(abs(extract_cv(res$value)$pred - ref_pred))
  }

  data.frame(
    offset = off,
    kappa_X = kappa_raw,
    gstat_raw_err = err_vs_ref(rg_raw),
    fastkcv_raw_err = err_vs_ref(rf_raw),
    fastkcv_std_err = err_vs_ref(rf_std),
    fastkcv_raw_warns = length(rf_raw$warnings),
    stringsAsFactors = FALSE
  )
})
res_cond <- do.call(rbind, cond_rows)

cat("\n=== PART 3: Quadratic UK Conditioning Sweep vs. Exact Reference (n = 120) ===\n")
print(res_cond, row.names = FALSE)

# ==============================================================================
# PART 4: Real-World Workflow — Variogram Hyperparameter Grid Search (meuse)
# ==============================================================================

uk_meuse <- log(zinc) ~ sqrt(dist)
obs_meuse <- log(meuse$zinc)
rmse <- function(o, p) sqrt(mean((o - p)^2, na.rm = TRUE))

grid_10f <- expand.grid(
  model = c("Exp", "Sph", "Mat"),
  range = c(500, 900, 1500),
  nugget = c(0.01, 0.05, 0.15),
  stringsAsFactors = FALSE
)
grid_loo <- expand.grid(
  model = c("Exp", "Sph", "Mat"),
  range = c(600, 1200),
  nugget = c(0.02, 0.10),
  stringsAsFactors = FALSE
)

eval_grid <- function(grid_df, engine, mode = c("10fold", "loo")) {
  mode <- match.arg(mode)
  rmses <- numeric(nrow(grid_df))
  for (i in seq_len(nrow(grid_df))) {
    m <- vgm(
      psill = 0.25 - grid_df$nugget[i],
      model = grid_df$model[i],
      range = grid_df$range[i],
      nugget = grid_df$nugget[i],
      kappa = if (grid_df$model[i] == "Mat") 1.5 else 0.5
    )
    res <- if (engine == "gstat") {
      if (mode == "loo") {
        krige.cv(uk_meuse, meuse, model = m, nfold = nrow(meuse), verbose = FALSE)
      } else {
        krige.cv(uk_meuse, meuse, model = m, nfold = folds_meuse, verbose = FALSE)
      }
    } else {
      if (mode == "loo") {
        fastkcv::kcv(uk_meuse, meuse, model = m, nfold = nrow(meuse))
      } else {
        fastkcv::kcv(uk_meuse, meuse, model = m, folds = folds_meuse)
      }
    }
    rmses[i] <- rmse(obs_meuse, extract_cv(res)$pred)
  }
  rmses
}

# Verify identical model selection & grid-wide RMSE parity
rmse_10f_g <- eval_grid(grid_10f, "gstat", "10fold")
rmse_10f_f <- eval_grid(grid_10f, "fastkcv", "10fold")

rmse_loo_g <- eval_grid(grid_loo, "gstat", "loo")
rmse_loo_f <- eval_grid(grid_loo, "fastkcv", "loo")

cat("\n=== PART 4: Variogram Grid Search Agreement & Workflow Speedup (meuse) ===\n")
cat(sprintf(
  "10-Fold Grid (%d models): max |RMSE diff| = %.3e | same best model = %s (idx %d vs %d)\n",
  nrow(grid_10f), max(abs(rmse_10f_g - rmse_10f_f)),
  identical(which.min(rmse_10f_g), which.min(rmse_10f_f)),
  which.min(rmse_10f_g), which.min(rmse_10f_f)
))
cat(sprintf(
  "LOO-CV Grid  (%d models): max |RMSE diff| = %.3e | same best model = %s (idx %d vs %d)\n",
  nrow(grid_loo), max(abs(rmse_loo_g - rmse_loo_f)),
  identical(which.min(rmse_loo_g), which.min(rmse_loo_f)),
  which.min(rmse_loo_g), which.min(rmse_loo_f)
))

cat("\nBest model (LOO-CV grid):\n")
best_idx <- which.min(rmse_loo_f)
print(cbind(grid_loo[best_idx, , drop = FALSE], rmse = rmse_loo_f[best_idx]), row.names = FALSE)

cat("\nWorkflow Benchmarks:\n")
bm_grid_10f <- bench::mark(
  gstat_10f_27models = eval_grid(grid_10f, "gstat", "10fold"),
  fastkcv_10f_27models = eval_grid(grid_10f, "fastkcv", "10fold"),
  iterations = 3,
  check = FALSE
)
bm_grid_loo <- bench::mark(
  gstat_loo_12models = eval_grid(grid_loo, "gstat", "loo"),
  fastkcv_loo_12models = eval_grid(grid_loo, "fastkcv", "loo"),
  iterations = 1,
  check = FALSE
)

print(bm_grid_10f[, c("expression", "min", "median", "itr/sec", "mem_alloc")])
print(bm_grid_loo[, c("expression", "min", "median", "itr/sec", "mem_alloc")])
