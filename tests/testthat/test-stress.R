# Unusual geometries, scales, trends and fold vectors, against krige.cv().

stress_field <- function(coords, seed, range = 0.3, mult = 1, extra = NULL) {
  withr::with_seed(seed, {
    n <- nrow(coords)
    C <- exp(-as.matrix(dist(coords)) / range) + diag(0.1, n)
    eps <- drop(crossprod(chol(C), stats::rnorm(n)))
  })
  xs <- as.numeric(scale(coords[, 1]))
  ys <- as.numeric(scale(coords[, 2]))
  d <- data.frame(x = coords[, 1], y = coords[, 2], xs = xs, ys = ys,
    z = (2 + 1.2 * xs - 0.7 * ys + eps) * mult)
  if (!is.null(extra)) d <- cbind(d, extra)
  d
}

v_stress <- gstat::vgm(1, "Exp", 0.3, 0.1)

expect_matches_gstat <- function(formula, d, model, folds = NULL,
                                 tolerance = tol) {
  fk <- if (is.null(folds)) {
    kcv(formula, ~ x + y, d, model = model)
  } else {
    kcv(formula, ~ x + y, d, model = model, folds = folds)
  }
  gs <- krige_cv(formula, ~ x + y, d, model = model,
    nfold = if (is.null(folds)) nrow(d) else folds)
  expect_cv_equal(fk, gs, tolerance = tolerance)
  invisible(fk)
}

base_xy <- function() {
  withr::with_seed(106, cbind(stats::runif(150), stats::runif(150)))
}

test_that("a regular lattice and tight clusters match krige.cv", {
  skip_on_cran()
  # Many exactly equal distances on the lattice.
  grid <- as.matrix(expand.grid(seq(0, 1, length.out = 12),
    seq(0, 1, length.out = 12)))
  expect_matches_gstat(z ~ 1, stress_field(grid, 101), v_stress)

  # Six clusters of 25 points with sd 0.002 (nearly coincident points),
  # leave-one-out and leave-one-cluster-out with a linear trend.
  centers <- cbind(c(0.1, 0.5, 0.9, 0.2, 0.7, 0.85),
    c(0.15, 0.2, 0.1, 0.8, 0.75, 0.9))
  id <- rep(1:6, each = 25)
  xy <- withr::with_seed(102, centers[id, ] +
    matrix(stats::rnorm(300, 0, 0.002), ncol = 2))
  d <- stress_field(xy, 103)
  expect_matches_gstat(z ~ 1, d, v_stress)
  expect_matches_gstat(z ~ xs + ys, d, v_stress, folds = id)
})

test_that("a 10,000:1 strip matches krige.cv", {
  skip_on_cran()
  xy <- withr::with_seed(104, cbind(stats::runif(150, 0, 10),
    stats::runif(150, 0, 1e-3)))
  expect_matches_gstat(z ~ xs + ys, stress_field(xy, 105, range = 2),
    gstat::vgm(1, "Exp", 2, 0.1))
})

test_that("results scale with the response and the sill", {
  skip_on_cran()
  # With z and the sill multiplied by s (nugget too), predictions scale by
  # s and variances by s^2, so a fixed absolute tolerance is meaningless:
  # compare relative to the largest value, at the usual 1e-8.
  xy <- base_xy()
  for (s in c(1e-6, 1e6)) {
    d <- stress_field(xy, 107, mult = s)
    v <- gstat::vgm(s^2, "Exp", 0.3, 0.1 * s^2)
    fk <- kcv(z ~ xs + ys, ~ x + y, d, model = v)
    gs <- krige_cv(z ~ xs + ys, ~ x + y, d, model = v)
    for (col in c("var1.pred", "var1.var", "residual")) {
      expect_close(fk[[col]] / max(abs(gs[[col]])),
        gs[[col]] / max(abs(gs[[col]])), label = paste(col, s))
    }
    base <- kcv(z ~ xs + ys, ~ x + y, stress_field(xy, 107), model = v_stress)
    expect_close(fk$var1.pred / s, base$var1.pred, label = paste("scaled", s))
    expect_close(fk$var1.var / s^2, base$var1.var, label = paste("var", s))
  }
})

test_that("universal kriging with 20 covariates matches krige.cv", {
  skip_on_cran()
  xy <- base_xy()
  n <- nrow(xy)
  extra <- withr::with_seed(108, {
    orth <- matrix(stats::rnorm(n * 20), n)
    u <- stats::rnorm(n)
    coll <- sqrt(0.99) * u + sqrt(0.01) * matrix(stats::rnorm(n * 20), n)
    colnames(orth) <- paste0("v", 1:20)
    colnames(coll) <- paste0("c", 1:20)
    as.data.frame(cbind(orth, coll))
  })
  d <- stress_field(xy, 109, extra = extra)
  folds <- withr::with_seed(42, stats::kmeans(xy, centers = 5)$cluster)
  expect_matches_gstat(reformulate(paste0("v", 1:20), "z"), d, v_stress,
    folds = folds)
  # Pairwise correlation 0.99 between the 20 covariates.
  expect_matches_gstat(reformulate(paste0("c", 1:20), "z"), d, v_stress)
})

test_that("unbalanced, factor and character fold vectors match krige.cv", {
  skip_on_cran()
  d <- stress_field(base_xy(), 107)
  # One fold of 120 points (80%), one of 25 and five singletons.
  imbal <- withr::with_seed(42, sample(c(rep(1L, 120), rep(2L, 25), 3:7)))
  expect_matches_gstat(z ~ xs + ys, d, v_stress, folds = imbal)

  # Labels of other types give the same predictions as integer labels in
  # the same sorted order; the fold column keeps the labels as given.
  blk <- withr::with_seed(42, unname(stats::kmeans(d[, c("x", "y")],
    centers = 5)$cluster))
  ref <- kcv(z ~ xs + ys, ~ x + y, d, model = v_stress, folds = blk)
  f_unused <- factor(blk, levels = c(1:5, 99))
  out_f <- kcv(z ~ xs + ys, ~ x + y, d, model = v_stress, folds = f_unused)
  expect_identical(out_f$var1.pred, ref$var1.pred)
  expect_identical(out_f$var1.var, ref$var1.var)
  expect_identical(out_f$fold, f_unused)
  f_chr <- letters[blk]
  out_c <- kcv(z ~ xs + ys, ~ x + y, d, model = v_stress, folds = f_chr)
  expect_identical(out_c$var1.pred, ref$var1.pred)
  expect_identical(out_c$fold, f_chr)
})

test_that("a three-level nested model with mixed anisotropy matches krige.cv", {
  skip_on_cran()
  v3 <- gstat::vgm(0.4, "Exp", 0.15, 0.05, anis = c(30, 0.5),
    add.to = gstat::vgm(0.35, "Sph", 0.5, anis = c(120, 0.4),
      add.to = gstat::vgm(0.25, "Gau", 0.8)))
  expect_matches_gstat(z ~ xs + ys, stress_field(base_xy(), 107), v3)
})
