# ==============================================================================
# Adversarial Stress & High-K Scaling Suite: fastkcv::kcv vs. gstat::krige.cv
# File: fastkcv_test_adversarial_stress.R
# Covers:
#   0. Environment & BLAS/LAPACK metadata
#   1. Pathological geometries, extreme anisotropy/scales, 3-level nested vgm,
#      high-dimensional/collinear UK (p = 20), leave-one-cluster-out, mixed
#      extreme fold imbalance, row-order invariance & unused factor levels
#   2. Covariance conditioning sweep (kappa(C)): Gau & Ste(k=5) LOO vs. 5-fold
#   3. High-K crossover (50-fold & LOO) vs. local gstat (nmax = 50) on a
#      long-range spatial field (n up to 2000)
#   4. 20x repeated 10-fold CV endurance & isolated subprocess peak RSS
#      ceiling (n = 1500, 2000, 2500; parent-generated inputs)
# ==============================================================================

library(fastkcv)
library(gstat)
library(sp)
library(bench)

# ==============================================================================
# PART 0: Environment & Linear Algebra Metadata
# ==============================================================================

fmt_lib <- function(x) if (length(x) == 1 && nzchar(x)) x else "<bundled/default>"

cat("=== PART 0: Environment & BLAS/LAPACK Metadata ===\n")
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

expr_labels <- function(expr_col) {
  nms <- names(expr_col)
  if (!is.null(nms) && all(nzchar(nms))) nms else as.character(expr_col)
}

rmse <- function(obs, pred) sqrt(mean((obs - pred)^2, na.rm = TRUE))

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

# Note: Armadillo prints C++ rcond warnings directly to stderr rather than
# Rf_warning(), so warn_f only counts R-level warnings.
run_comparison <- function(label, expr_g, expr_f,
                           tol_pred = 1e-10, tol_var = 1e-10,
                           tol_rel = NULL, core = TRUE) {
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
      case = label, core = core,
      max_abs_pred = NA_real_, max_rel_pred = NA_real_,
      max_abs_var = NA_real_, max_rel_var = NA_real_,
      na_parity = NA, neg_var_f = NA_integer_,
      warn_g = length(rg$warnings), warn_f = length(rf$warnings),
      status = status, note = trimws(gsub("^\\|\\s*|\\s*\\|$", "", msg)),
      stringsAsFactors = FALSE
    ))
  }

  g <- extract_cv(rg$value)
  f <- extract_cv(rf$value)
  if (length(g$pred) != length(f$pred)) {
    return(data.frame(
      case = label, core = core,
      max_abs_pred = NA_real_, max_rel_pred = NA_real_,
      max_abs_var = NA_real_, max_rel_var = NA_real_,
      na_parity = FALSE, neg_var_f = NA_integer_,
      warn_g = length(rg$warnings), warn_f = length(rf$warnings),
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
  max_rel_var <- safe_max(diff_var / pmax(abs(g$var), 1e-12))
  neg_var_f <- sum(f$var < -1e-14, na.rm = TRUE)

  within_tol <- if (!is.null(tol_rel)) {
    (is.na(max_rel_pred) || max_rel_pred <= tol_rel) &&
      (is.na(max_rel_var) || max_rel_var <= tol_rel)
  } else {
    (is.na(max_abs_pred) || max_abs_pred <= tol_pred) &&
      (is.na(max_abs_var) || max_abs_var <= tol_var)
  }

  data.frame(
    case = label, core = core,
    max_abs_pred = max_abs_pred, max_rel_pred = max_rel_pred,
    max_abs_var = max_abs_var, max_rel_var = max_rel_var,
    na_parity = na_parity, neg_var_f = neg_var_f,
    warn_g = length(rg$warnings), warn_f = length(rf$warnings),
    status = if (na_parity && within_tol) "PASS" else "FAIL",
    note = "",
    stringsAsFactors = FALSE
  )
}

simulate_sp_data <- function(coords, range = 0.3, psill = 1.0, nugget = 0.1,
                             seed = 123, z_mult = 1.0, extra_df = NULL) {
  set.seed(seed)
  n <- nrow(coords)
  d_mat <- as.matrix(dist(coords))
  cov_mat <- psill * exp(-d_mat / range) + diag(nugget, n)
  eps <- as.numeric(t(chol(cov_mat)) %*% rnorm(n))

  xs <- as.numeric(scale(coords[, 1]))
  ys <- as.numeric(scale(coords[, 2]))
  df <- data.frame(
    x = coords[, 1],
    y = coords[, 2],
    xs = xs,
    ys = ys,
    z = (2.0 + 1.2 * xs - 0.7 * ys + eps) * z_mult
  )
  if (!is.null(extra_df)) df <- cbind(df, extra_df)
  coordinates(df) <- ~ x + y
  df
}

make_uniform_sp <- function(n, range = 0.3, psill = 1.0, nugget = 0.1, seed = 123) {
  set.seed(seed)
  coords <- cbind(x = runif(n, 0, 1), y = runif(n, 0, 1))
  simulate_sp_data(coords, range = range, psill = psill, nugget = nugget, seed = seed + 1000)
}

# ==============================================================================
# PART 1: Adversarial Geometries, High-p Trends & Extreme Fold Structures
# ==============================================================================

# 1a. Regular 12x12 lattice (exact distance symmetries)
grid_xy <- as.matrix(expand.grid(
  x = seq(0, 1, length.out = 12),
  y = seq(0, 1, length.out = 12)
))
dat_grid <- simulate_sp_data(grid_xy, seed = 101)

# 1b. Multi-scale tight clusters (6 clusters x 25 pts, intra-cluster sd = 0.002)
set.seed(102)
centers <- cbind(c(0.1, 0.5, 0.9, 0.2, 0.7, 0.85), c(0.15, 0.2, 0.1, 0.8, 0.75, 0.9))
clust_id <- rep(seq_len(nrow(centers)), each = 25L)
clust_xy <- cbind(
  x = centers[clust_id, 1] + rnorm(length(clust_id), 0, 0.002),
  y = centers[clust_id, 2] + rnorm(length(clust_id), 0, 0.002)
)
dat_clust <- simulate_sp_data(clust_xy, seed = 103)

# 1c. Extreme aspect-ratio strip (10,000:1 domain elongation)
set.seed(104)
strip_xy <- cbind(x = runif(150, 0, 10), y = runif(150, 0, 1e-3))
dat_strip <- simulate_sp_data(strip_xy, range = 2.0, seed = 105)

# 1d. Base dataset + extreme response/variance magnitudes (1e-6 and 1e+6)
set.seed(106)
base_xy <- cbind(x = runif(150, 0, 1), y = runif(150, 0, 1))
dat_base <- simulate_sp_data(base_xy, seed = 107)
dat_tiny_z <- simulate_sp_data(base_xy, seed = 107, z_mult = 1e-6)
dat_huge_z <- simulate_sp_data(base_xy, seed = 107, z_mult = 1e6)

# 1e. High-dimensional UK (p = 20 orthogonal vs. p = 20 collinear rho = 0.99)
set.seed(108)
n_hp <- 150
V_orth <- matrix(rnorm(n_hp * 20), nrow = n_hp, ncol = 20)
colnames(V_orth) <- paste0("v", seq_len(20))
u_shared <- rnorm(n_hp)
V_coll <- sqrt(0.99) * u_shared + sqrt(0.01) * matrix(rnorm(n_hp * 20), nrow = n_hp, ncol = 20)
colnames(V_coll) <- paste0("c", seq_len(20))
dat_highp <- simulate_sp_data(
  base_xy,
  seed = 109,
  extra_df = as.data.frame(cbind(V_orth, V_coll))
)
f_uk_orth <- reformulate(colnames(V_orth), response = "z")
f_uk_coll <- reformulate(colnames(V_coll), response = "z")

# Variograms for Part 1
v_exp <- vgm(psill = 1.0, model = "Exp", range = 0.3, nugget = 0.1)
v_strip <- vgm(psill = 1.0, model = "Exp", range = 2.0, nugget = 0.1)
v_anis_ext <- vgm(psill = 1.0, model = "Exp", range = 0.5, nugget = 0.1, anis = c(37, 0.02))
v_short <- vgm(psill = 1.0, model = "Exp", range = 0.01, nugget = 0.05)
v_long <- vgm(psill = 1.0, model = "Exp", range = 10.0, nugget = 0.05)
v_tiny <- vgm(psill = 1e-12, model = "Exp", range = 0.3, nugget = 1e-13)
v_huge <- vgm(psill = 1e12, model = "Exp", range = 0.3, nugget = 1e11)
v_nested3 <- vgm(
  psill = 0.40, model = "Exp", range = 0.15, nugget = 0.05, anis = c(30, 0.5),
  add.to = vgm(
    psill = 0.35, model = "Sph", range = 0.50, anis = c(120, 0.4),
    add.to = vgm(psill = 0.25, model = "Gau", range = 0.80)
  )
)

# Extreme fold partitions on dat_base (n = 150)
set.seed(42)
# 80% holdout block (120 pts) + 25-pt block + five 1-pt LOO folds in one call
folds_imbal <- sample(c(rep(1L, 120), rep(2L, 25), 3L, 4L, 5L, 6L, 7L))
folds_blk5 <- kmeans(coordinates(dat_base), centers = 5)$cluster
folds_unused <- factor(folds_blk5, levels = c(1:5, 99))

# Permuted row order check
set.seed(99)
perm <- sample(seq_len(nrow(dat_base)))
dat_perm <- dat_base[perm, ]
folds_perm <- folds_blk5[perm]

p1_rows <- list(
  run_comparison(
    "Geometry: Regular 12x12 lattice (OK LOO)",
    krige.cv(z ~ 1, dat_grid, model = v_exp, nfold = nrow(dat_grid), verbose = FALSE),
    fastkcv::kcv(z ~ 1, dat_grid, model = v_exp, nfold = nrow(dat_grid))
  ),
  run_comparison(
    "Geometry: Tight 6x25 clusters (OK LOO)",
    krige.cv(z ~ 1, dat_clust, model = v_exp, nfold = nrow(dat_clust), verbose = FALSE),
    fastkcv::kcv(z ~ 1, dat_clust, model = v_exp, nfold = nrow(dat_clust))
  ),
  run_comparison(
    "Geometry: Leave-one-cluster-out (6 blocks)",
    krige.cv(z ~ xs + ys, dat_clust, model = v_exp, nfold = clust_id, verbose = FALSE),
    fastkcv::kcv(z ~ xs + ys, dat_clust, model = v_exp, folds = clust_id)
  ),
  run_comparison(
    "Geometry: 10000:1 strip (UK LOO)",
    krige.cv(z ~ xs + ys, dat_strip, model = v_strip, nfold = nrow(dat_strip), verbose = FALSE),
    fastkcv::kcv(z ~ xs + ys, dat_strip, model = v_strip, nfold = nrow(dat_strip))
  ),
  run_comparison(
    "Variogram: Extreme anisotropy (37 deg, 0.02)",
    krige.cv(z ~ 1, dat_base, model = v_anis_ext, nfold = nrow(dat_base), verbose = FALSE),
    fastkcv::kcv(z ~ 1, dat_base, model = v_anis_ext, nfold = nrow(dat_base))
  ),
  run_comparison(
    "Variogram: 3-level nested + mixed anis",
    krige.cv(z ~ xs + ys, dat_base, model = v_nested3, nfold = nrow(dat_base), verbose = FALSE),
    fastkcv::kcv(z ~ xs + ys, dat_base, model = v_nested3, nfold = nrow(dat_base))
  ),
  run_comparison(
    "Variogram: Ultra-short range (r = 0.01)",
    krige.cv(z ~ 1, dat_base, model = v_short, nfold = nrow(dat_base), verbose = FALSE),
    fastkcv::kcv(z ~ 1, dat_base, model = v_short, nfold = nrow(dat_base))
  ),
  run_comparison(
    "Variogram: Ultra-long range (r = 10.0)",
    krige.cv(z ~ 1, dat_base, model = v_long, nfold = nrow(dat_base), verbose = FALSE),
    fastkcv::kcv(z ~ 1, dat_base, model = v_long, nfold = nrow(dat_base))
  ),
  run_comparison("Scale: Tiny response/sill (z * 1e-6)",
    krige.cv(z ~ xs + ys, dat_tiny_z, model = v_tiny, nfold = nrow(dat_tiny_z), verbose = FALSE),
    fastkcv::kcv(z ~ xs + ys, dat_tiny_z, model = v_tiny, nfold = nrow(dat_tiny_z)),
    tol_rel = 1e-9
  ),
  run_comparison("Scale: Huge response/sill (z * 1e+6)",
    krige.cv(z ~ xs + ys, dat_huge_z, model = v_huge, nfold = nrow(dat_huge_z), verbose = FALSE),
    fastkcv::kcv(z ~ xs + ys, dat_huge_z, model = v_huge, nfold = nrow(dat_huge_z)),
    tol_rel = 1e-9
  ),
  run_comparison(
    "High-p UK: p = 20 orthogonal (10-fold)",
    krige.cv(f_uk_orth, dat_highp, model = v_exp, nfold = folds_blk5, verbose = FALSE),
    fastkcv::kcv(f_uk_orth, dat_highp, model = v_exp, folds = folds_blk5)
  ),
  run_comparison("High-p UK: p = 20 collinear rho=0.99 (LOO)",
    krige.cv(f_uk_coll, dat_highp, model = v_exp, nfold = nrow(dat_highp), verbose = FALSE),
    fastkcv::kcv(f_uk_coll, dat_highp, model = v_exp, nfold = nrow(dat_highp)),
    tol_pred = 1e-8, tol_var = 1e-8
  ),
  run_comparison(
    "Folds: Extreme imbalance (80% blk + LOO)",
    krige.cv(z ~ xs + ys, dat_base, model = v_exp, nfold = folds_imbal, verbose = FALSE),
    fastkcv::kcv(z ~ xs + ys, dat_base, model = v_exp, folds = folds_imbal)
  ),
  run_comparison(
    "Folds: Shuffled row order (Block-5)",
    krige.cv(z ~ xs + ys, dat_perm, model = v_exp, nfold = folds_perm, verbose = FALSE),
    fastkcv::kcv(z ~ xs + ys, dat_perm, model = v_exp, folds = folds_perm)
  ),
  run_comparison("Folds: Factor with unused level",
    krige.cv(z ~ xs + ys, dat_base, model = v_exp, nfold = folds_unused, verbose = FALSE),
    fastkcv::kcv(z ~ xs + ys, dat_base, model = v_exp, folds = folds_unused),
    core = FALSE
  )
)

res_p1 <- do.call(rbind, p1_rows)
cat("=== PART 1: Adversarial Geometries, High-p UK & Extreme Folds ===\n")
print(
  res_p1[, c("case", "max_abs_pred", "max_rel_pred", "max_abs_var", "na_parity", "neg_var_f", "status")],
  row.names = FALSE
)

# Explicit permutation-invariance check against unpermuted baseline
ref_unperm <- extract_cv(fastkcv::kcv(z ~ xs + ys, dat_base, model = v_exp, folds = folds_blk5))
res_perm_f <- extract_cv(fastkcv::kcv(z ~ xs + ys, dat_perm, model = v_exp, folds = folds_perm))
inv_perm <- order(perm)
perm_diff <- max(abs(res_perm_f$pred[inv_perm] - ref_unperm$pred))
cat(sprintf("\nRow-permutation invariance (fastkcv unpermuted vs. baseline): %.3e\n", perm_diff))
stopifnot(perm_diff < 1e-12)

p1_notes <- res_p1[nzchar(res_p1$note), c("case", "status", "note")]
if (nrow(p1_notes) > 0) {
  cat("\nPart 1 non-core notes / errors:\n")
  print(p1_notes, row.names = FALSE)
}

p1_core_fail <- res_p1[res_p1$core & res_p1$status != "PASS", ]
if (nrow(p1_core_fail) > 0) {
  warning("Part 1 core check(s) failed: ", paste(p1_core_fail$case, collapse = ", "))
}

# ==============================================================================
# PART 2: Covariance Conditioning Sweep (kappa(C) Breakdown Curve)
# Tests Schur-complement block updates (5-fold) vs. LOO vs. direct gstat solves
# ==============================================================================

cat("\n=== PART 2: Covariance Conditioning Sweep (n = 120, Gau & Ste k=5) ===\n")
dat_cond <- make_uniform_sp(n = 120, seed = 201)
d_cond <- as.matrix(dist(coordinates(dat_cond)))
set.seed(42)
folds_cond5 <- kmeans(coordinates(dat_cond), centers = 5)$cluster

nug_grid <- c(1e-1, 1e-3, 1e-5, 1e-7, 1e-9, 1e-11, 0)
cond_models <- list(
  Gau     = function(nug) vgm(psill = 1.0, model = "Gau", range = 0.4, nugget = nug),
  Ste_k5  = function(nug) vgm(psill = 1.0, model = "Ste", range = 0.4, nugget = nug, kappa = 5.0)
)

cond_rows <- list()
for (mnm in names(cond_models)) {
  for (nug in nug_grid) {
    vm <- cond_models[[mnm]](nug)

    # Exact covariance matrix condition number kappa(C)
    C_mat <- if (mnm == "Gau") {
      M <- 1.0 * exp(-(d_cond / 0.4)^2)
      diag(M) <- 1.0 + nug
      M
    } else {
      vl <- gstat::variogramLine(vm, dist_vector = as.numeric(d_cond), covariance = TRUE)$gamma
      M <- matrix(vl, nrow = nrow(d_cond), ncol = ncol(d_cond))
      diag(M) <- 1.0 + nug
      M
    }
    kap_C <- tryCatch(kappa(C_mat, exact = TRUE), error = function(e) Inf)

    r_loo <- run_comparison(
      paste0(mnm, "_LOO"),
      krige.cv(z ~ 1, dat_cond, model = vm, nfold = nrow(dat_cond), verbose = FALSE),
      fastkcv::kcv(z ~ 1, dat_cond, model = vm, nfold = nrow(dat_cond)),
      tol_pred = 1e-6, tol_var = 1e-6, core = FALSE
    )
    r_k5 <- run_comparison(
      paste0(mnm, "_K5"),
      krige.cv(z ~ 1, dat_cond, model = vm, nfold = folds_cond5, verbose = FALSE),
      fastkcv::kcv(z ~ 1, dat_cond, model = vm, folds = folds_cond5),
      tol_pred = 1e-6, tol_var = 1e-6, core = FALSE
    )

    cond_rows[[length(cond_rows) + 1]] <- data.frame(
      model = mnm,
      nugget = nug,
      kappa_C = kap_C,
      loo_pred_err = r_loo$max_abs_pred,
      loo_var_err = r_loo$max_abs_var,
      loo_status = r_loo$status,
      k5_pred_err = r_k5$max_abs_pred,
      k5_var_err = r_k5$max_abs_var,
      k5_neg_var_f = r_k5$neg_var_f,
      k5_status = r_k5$status,
      stringsAsFactors = FALSE
    )
  }
}

res_cond <- do.call(rbind, cond_rows)
print(res_cond, row.names = FALSE)

# ==============================================================================
# PART 3: High-K Scaling Crossover — 50-Fold & LOO vs. Local gstat (nmax = 50)
# Long-range field (range = 1.5 on [0,1]^2) stresses local-neighborhood accuracy
# ==============================================================================

v_long_scale <- vgm(psill = 1.0, model = "Exp", range = 1.5, nugget = 0.05)

# 3a. 50-Fold UK CV across n = 250 .. 2000
cat("\n=== PART 3a: 50-Fold UK CV Crossover (Long-Range Field, n = 250..2000) ===\n")
n_grid_50 <- c(250, 500, 1000, 1500, 2000)

bm_50f <- bench::press(
  n = n_grid_50,
  {
    dat <- make_uniform_sp(n, range = 1.5, psill = 1.0, nugget = 0.05, seed = 300 + n)
    set.seed(42)
    f_50 <- sample(rep(seq_len(50), length.out = n))
    gc(FALSE)

    bench::mark(
      gstat_local = krige.cv(z ~ xs + ys, dat,
        model = v_long_scale, nmax = 50,
        nfold = f_50, verbose = FALSE
      ),
      fastkcv = fastkcv::kcv(z ~ xs + ys, dat, model = v_long_scale, folds = f_50),
      min_time = 0.5,
      min_iterations = 2,
      check = FALSE
    )
  }
)

lbl_50 <- expr_labels(bm_50f$expression)
s50_l <- bm_50f[lbl_50 == "gstat_local", ]
s50_f <- bm_50f[lbl_50 == "fastkcv", ]
idx_50 <- match(s50_l$n, s50_f$n)
print(data.frame(
  n = s50_l$n,
  gstat_local_ms = round(as.numeric(s50_l$median) * 1e3, 2),
  fastkcv_ms = round(as.numeric(s50_f$median[idx_50]) * 1e3, 2),
  fastkcv_vs_local_x = round(as.numeric(s50_l$median) / as.numeric(s50_f$median[idx_50]), 2)
), row.names = FALSE)

# 3b. LOO-CV across n = 250 .. 1500 (tests whether local gstat ever catches up)
cat("\n=== PART 3b: LOO-CV Crossover (Long-Range Field, n = 250..1500) ===\n")
n_grid_loo <- c(250, 500, 1000, 1500)

bm_loo <- bench::press(
  n = n_grid_loo,
  {
    dat <- make_uniform_sp(n, range = 1.5, psill = 1.0, nugget = 0.05, seed = 400 + n)
    gc(FALSE)

    bench::mark(
      gstat_local = krige.cv(z ~ xs + ys, dat,
        model = v_long_scale, nmax = 50,
        nfold = n, verbose = FALSE
      ),
      fastkcv = fastkcv::kcv(z ~ xs + ys, dat, model = v_long_scale, nfold = n),
      min_time = 0.5,
      min_iterations = 1,
      check = FALSE
    )
  }
)

lbl_loo <- expr_labels(bm_loo$expression)
sloo_l <- bm_loo[lbl_loo == "gstat_local", ]
sloo_f <- bm_loo[lbl_loo == "fastkcv", ]
idx_loo <- match(sloo_l$n, sloo_f$n)
print(data.frame(
  n = sloo_l$n,
  gstat_local_ms = round(as.numeric(sloo_l$median) * 1e3, 2),
  fastkcv_ms = round(as.numeric(sloo_f$median[idx_loo]) * 1e3, 2),
  fastkcv_vs_local_x = round(as.numeric(sloo_l$median) / as.numeric(sloo_f$median[idx_loo]), 2)
), row.names = FALSE)

# Explicit accuracy verification + global vs. local RMSE check outside timing
cat("\n=== PART 3c: Accuracy Verification (50-Fold UK, n = 1000 & n = 1500) ===\n")
dat_v1000 <- make_uniform_sp(1000, range = 1.5, psill = 1.0, nugget = 0.05, seed = 1300)
set.seed(42)
f_v1000 <- sample(rep(seq_len(50), length.out = 1000))

res_g1000 <- krige.cv(z ~ xs + ys, dat_v1000, model = v_long_scale, nfold = f_v1000, verbose = FALSE)
res_l1000 <- krige.cv(z ~ xs + ys, dat_v1000, model = v_long_scale, nmax = 50, nfold = f_v1000, verbose = FALSE)
res_f1000 <- fastkcv::kcv(z ~ xs + ys, dat_v1000, model = v_long_scale, folds = f_v1000)

chk_1000 <- run_comparison("Exact check: n=1000, 50-fold UK", res_g1000, res_f1000)
print(chk_1000[, c("case", "max_abs_pred", "max_rel_pred", "max_abs_var", "na_parity", "status")],
  row.names = FALSE
)
cat(sprintf(
  "Long-range RMSE (n=1000) — gstat_global: %.6f | fastkcv: %.6f | gstat_local(50): %.6f\n",
  rmse(dat_v1000$z, extract_cv(res_g1000)$pred),
  rmse(dat_v1000$z, extract_cv(res_f1000)$pred),
  rmse(dat_v1000$z, extract_cv(res_l1000)$pred)
))
stopifnot(chk_1000$status == "PASS")

# ==============================================================================
# PART 4a: Repeated 10-Fold CV Endurance (20 Repeats, n = 300)
# ==============================================================================

cat("\n=== PART 4a: Repeated 10-Fold CV Endurance (20 Repeats, n = 300) ===\n")
n_rep <- 300
R_reps <- 20L
dat_rep <- make_uniform_sp(n_rep, seed = 601)
rep_folds <- lapply(seq_len(R_reps), function(r) {
  set.seed(700 + r)
  sample(rep(seq_len(10), length.out = n_rep))
})

run_repeated <- function(engine) {
  vapply(rep_folds, function(fk) {
    res <- if (engine == "gstat") {
      krige.cv(z ~ xs + ys, dat_rep, model = v_exp, nfold = fk, verbose = FALSE)
    } else {
      fastkcv::kcv(z ~ xs + ys, dat_rep, model = v_exp, folds = fk)
    }
    rmse(dat_rep$z, extract_cv(res)$pred)
  }, numeric(1))
}

rmses_rep_g <- run_repeated("gstat")
rmses_rep_f <- run_repeated("fastkcv")
max_rep_diff <- max(abs(rmses_rep_g - rmses_rep_f))
cat(sprintf(
  "20-repeat 10-fold RMSE parity: max |diff| = %.3e | mean RMSE = %.6f (sd = %.6f)\n",
  max_rep_diff, mean(rmses_rep_f), sd(rmses_rep_f)
))
stopifnot(max_rep_diff < 1e-12)

gc(FALSE)
bm_rep <- bench::mark(
  gstat_20x10f = run_repeated("gstat"),
  fastkcv_20x10f = run_repeated("fastkcv"),
  iterations = 2,
  check = FALSE
)
print(bm_rep[, c("expression", "min", "median", "itr/sec", "mem_alloc")])

# ==============================================================================
# PART 4b: High-n Isolated Subprocess Peak RSS (n = 1500, 2000, 2500)
# Inputs pre-generated in parent process so setup matrices do not taint peak RSS
# ==============================================================================

cat("\n=== PART 4b: Subprocess Peak Process Memory Ceiling (10-Fold UK) ===\n")

if (!requireNamespace("callr", quietly = TRUE)) {
  cat("Skipped: install 'callr' (and 'ps' on Windows/macOS) for subprocess RSS.\n")
} else {
  measure_subprocess_rss <- function(method, dat, v_mod, f_k) {
    callr::r(
      function(method, dat, v_mod, f_k) {
        suppressPackageStartupMessages({
          library(fastkcv)
          library(gstat)
          library(sp)
        })

        get_rss_mb <- function() {
          if (file.exists("/proc/self/status")) {
            hwm <- grep("^VmHWM:", readLines("/proc/self/status", warn = FALSE), value = TRUE)
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
        out <- switch(method,
          gstat_global = krige.cv(z ~ xs + ys, dat, model = v_mod, nfold = f_k, verbose = FALSE),
          gstat_local  = krige.cv(z ~ xs + ys, dat, model = v_mod, nmax = 50, nfold = f_k, verbose = FALSE),
          fastkcv      = fastkcv::kcv(z ~ xs + ys, dat, model = v_mod, folds = f_k)
        )
        elapsed <- proc.time()[["elapsed"]] - t0
        after <- get_rss_mb()

        data.frame(
          method = method,
          n = nrow(dat),
          elapsed_s = round(elapsed, 4),
          metric = after$metric,
          before_mb = round(before$value, 2),
          after_mb = round(after$value, 2),
          delta_mb = round(after$value - before$value, 2),
          stringsAsFactors = FALSE
        )
      },
      args = list(method = method, dat = dat, v_mod = v_mod, f_k = f_k)
    )
  }

  rss_ns <- c(1500L, 2000L, 2500L)
  parent_inputs <- lapply(rss_ns, function(n) {
    set.seed(42)
    list(
      dat = make_uniform_sp(n, seed = 800 + n),
      f_k = sample(rep(seq_len(10), length.out = n))
    )
  })
  names(parent_inputs) <- as.character(rss_ns)

  rss_plan <- data.frame(
    method = c(
      "gstat_global", "gstat_local", "fastkcv",
      "gstat_global", "gstat_local", "fastkcv",
      "gstat_local",  "fastkcv"
    ),
    n = c(1500L, 1500L, 1500L, 2000L, 2000L, 2000L, 2500L, 2500L),
    stringsAsFactors = FALSE
  )

  rss_rows <- lapply(seq_len(nrow(rss_plan)), function(i) {
    nm <- as.character(rss_plan$n[i])
    measure_subprocess_rss(
      method = rss_plan$method[i],
      dat    = parent_inputs[[nm]]$dat,
      v_mod  = v_exp,
      f_k    = parent_inputs[[nm]]$f_k
    )
  })
  print(do.call(rbind, rss_rows), row.names = FALSE)
}
