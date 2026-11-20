m_ok <- gstat::vgm(0.5, "Sph", 900, 0.1)

test_that("local-neighborhood and block arguments are rejected", {
  meuse <- meuse_df()
  for (arg in list(list(nmax = 20), list(nmin = 5), list(maxdist = 500),
                   list(block = c(40, 40)))) {
    expect_error(
      do.call(kcv, c(list(log(zinc) ~ 1, ~ x + y, meuse, model = m_ok), arg)),
      "global kriging"
    )
  }
  expect_error(kcv(log(zinc) ~ 1, ~ x + y, meuse, model = m_ok, foo = 1),
    "Unknown argument")
  expect_error(check_dots(1), "unnamed")
})

test_that("unbounded and unsupported models are rejected", {
  meuse <- meuse_df()
  expect_error(kcv(log(zinc) ~ 1, ~ x + y, meuse, model = gstat::vgm(1, "Lin", 0)),
    "Unbounded.*Lin")
  expect_error(kcv(log(zinc) ~ 1, ~ x + y, meuse, model = gstat::vgm(1, "Pow", 1)),
    "Unbounded.*Pow")
  expect_error(kcv(log(zinc) ~ 1, ~ x + y, meuse, model = gstat::vgm(1, "Cir", 300)),
    "Cir.*not supported")
  expect_error(kcv(log(zinc) ~ 1, ~ x + y, meuse, model = list(1)),
    "variogramModel")
  expect_error(kcv(log(zinc) ~ 1, ~ x + y, meuse), "`model`")
})

test_that("invalid model parameters are rejected", {
  m3d <- gstat::vgm(1, "Exp", 300, anis = c(30, 10, 0, 0.5, 0.5))
  expect_error(validate_model(m3d), "3D anisotropy")
  bad <- m_ok
  bad$psill[2] <- -1
  expect_error(validate_model(bad), "negative partial sill")
  bad <- m_ok
  bad$psill[] <- 0
  expect_error(validate_model(bad), "total sill of zero")
  bad <- m_ok
  bad$range[2] <- 0
  expect_error(validate_model(bad), "positive range")
  bad <- m_ok
  bad$anis1[2] <- 2
  expect_error(validate_model(bad), "ratios")
  expect_error(validate_model(m_ok[0, ]), "no components")
})

test_that("geographic coordinates are rejected", {
  ll <- sf::st_transform(meuse_sf(), 4326)
  expect_error(kcv(log(zinc) ~ 1, ll, model = m_ok), "longitude/latitude")
  ll_sp <- meuse_df()
  sp::coordinates(ll_sp) <- ~ x + y
  ll_sp <- sp::spTransform(meuse_sp(), sp::CRS("EPSG:4326"))
  expect_error(kcv(log(zinc) ~ 1, ll_sp, model = m_ok), "longitude/latitude")
})

test_that("bad inputs are rejected", {
  meuse <- meuse_df()
  na <- meuse
  na$zinc[3] <- NA
  expect_error(kcv(log(zinc) ~ 1, ~ x + y, na, model = m_ok), "missing value")
  na_cov <- meuse
  na_cov$dist[3] <- NA
  expect_error(kcv(log(zinc) ~ sqrt(dist), ~ x + y, na_cov, model = m_ok),
    "covariates")
  dup <- rbind(meuse, meuse[1, ])
  expect_error(kcv(log(zinc) ~ 1, ~ x + y, dup, model = m_ok), "Duplicate locations")
  expect_error(kcv(log(zinc) ~ 1, ~ x + y, meuse, model = m_ok, folds = 1:3),
    "one non-missing value")
  expect_error(kcv(log(zinc) ~ 1, ~ x + y, meuse, model = m_ok, nfold = 0),
    "positive integer")
  expect_error(kcv(log(zinc) ~ 1, ~ x + y, model = m_ok), "`data` must be supplied")
  expect_error(kcv(log(zinc) ~ 1, ~ x + y, as.list(meuse), model = m_ok),
    "must be a data.frame")
  expect_error(kcv(log(zinc) ~ 1, as.matrix(meuse[, 1:2]), model = m_ok),
    "must be a formula")
  expect_error(kcv(~ zinc, ~ x + y, meuse, model = m_ok), "two-sided")
  expect_error(kcv(landuse ~ 1, ~ x + y, meuse, model = m_ok), "numeric vector")
  expect_error(kcv(log(zinc) ~ 1, ~ x + y + dist, meuse, model = m_ok), "2D")
  expect_error(kcv(log(zinc) ~ 1, ~ x + y, meuse, model = m_ok, beta = c(1, 2)),
    "`beta` has length")
  lines <- sf::st_sf(z = 1, geometry = sf::st_sfc(
    sf::st_linestring(rbind(c(0, 0), c(1, 1))), crs = 28992))
  expect_error(kcv(z ~ 1, lines, model = m_ok), "POINT")
  expect_error(kcv(log(zinc) ~ 1, sp::SpatialPoints(meuse_sp()), model = m_ok),
    "Spatial\\*DataFrame")
  na_xy <- meuse
  na_xy$x[2] <- NA
  expect_error(kcv(log(zinc) ~ 1, ~ x + y, na_xy, model = m_ok), "Coordinates")
})

test_that("folds that leave the trend unidentifiable are rejected", {
  meuse <- meuse_df()
  # Factor trend: removing every observation of one level leaves its
  # coefficient unidentifiable.
  folds <- ifelse(meuse$ffreq == "3", 1L, 2L)
  expect_error(
    kcv(log(zinc) ~ ffreq, ~ x + y, meuse, model = m_ok, folds = folds),
    "not identifiable"
  )
  expect_error(
    kcv(log(zinc) ~ 1, ~ x + y, meuse, model = m_ok, nfold = 1),
    "every observation"
  )
  expect_error(
    kcv(log(zinc) ~ 1, ~ x + y, meuse, model = m_ok, beta = 5, folds = rep(1, 155)),
    "every observation"
  )
})

test_that("engine errors are informative", {
  xy <- cbind(c(0, 1, 2), c(0, 0, 0))
  C <- covariance_matrix(gstat::vgm(1, "Exp", 1), xy)
  expect_error(kcv_engine(C[1:2, 1:2], matrix(1, 3, 1), 1:3, 1:3), "n x n")
  expect_error(kcv_engine(C, matrix(1, 2, 1), 1:3, 1:3), "rows")
  expect_error(kcv_engine(C, matrix(1, 3, 1), 1:3, c(1, NA, 2)), "`folds`")
  singular <- matrix(1, 3, 3)
  expect_error(kcv_engine(singular, matrix(1, 3, 1), 1:3, 1:3), "not positive definite")
  expect_error(kcv_engine(C, cbind(rep(1, 3), rep(1, 3)), 1:3, 1:3), "rank deficient")
  # Q[I, I] is positive definite whenever the trend is identifiable, so the
  # numerical-singularity guard is tested directly.
  expect_error(fold_inverse(matrix(c(1, 1, 1, 1), 2), 4), "Fold 4.*singular")
  expect_equal(fold_inverse(matrix(4), 1), matrix(0.25))
})

test_that("C++ engine raises the R engine's fold errors", {
  s_xy <- withr::with_seed(5, matrix(stats::runif(40, 0, 100), ncol = 2))
  C <- covariance_matrix(gstat::vgm(1, "Exp", 30, 0.1), s_xy)
  F <- cbind(1, rep(c(0, 1), c(15, 5)))
  y <- withr::with_seed(6, stats::rnorm(20))
  # Fold 2 holds every row with F[, 2] = 1.
  folds <- rep(c(1, 3, 2), c(8, 7, 5))
  for (eng in list(kcv_engine, kcv_engine_r)) {
    expect_error(eng(C, F, y, folds), "Fold 2: the trend is not identifiable")
    expect_error(eng(C, F, y, rep(7, 20)), "Fold 7 contains every observation")
    expect_error(eng(C, F, y, rep(7, 20), beta = c(0, 1)),
      "Fold 7 contains every observation")
  }
  # Same check through the complement branch (fold larger than n / 2).
  # Fold 2 (rows 4 to 20) leaves only rows with F[, 2] = 1.
  folds_big <- rep(c(1, 2), c(3, 17))
  F_big <- cbind(1, rep(c(1, 0), c(5, 15)))
  expect_error(kcv_engine(C, F_big, y, folds_big), "Fold 2: the trend")
  expect_error(kcv_engine_r(C, F_big, y, folds_big), "Fold 2: the trend")
  # LOO path: a column nonzero in a single row.
  F_loo <- cbind(1, c(1, rep(0, 19)))
  expect_error(kcv_engine(C, F_loo, y, 20:1), "Fold 20: the trend")
  expect_error(kcv_engine_r(C, F_loo, y, 20:1), "Fold 20: the trend")
  # Simple kriging needs no identifiable trend.
  expect_no_error(kcv_engine(C, F, y, folds, beta = c(0, 1)))
  # Single observation, LOO.
  expect_error(kcv_engine(matrix(1), matrix(1), 1, 1), "every observation")
})

test_that("C++ fold inverse detects numerically singular Q[I, I]", {
  expect_null(fold_inverse_cpp(matrix(c(1, 1, 1, 1), 2)))
  expect_null(fold_inverse_cpp(matrix(-1)))
  expect_equal(fold_inverse_cpp(matrix(4)), matrix(0.25))
  A <- matrix(c(2, 1, 1, 3), 2)
  expect_equal(fold_inverse_cpp(A), solve(A), tolerance = tol)
})

test_that("inputs above the memory budget are rejected before allocation", {
  # Budget for n = 155: kcv() needs one 155 x 155 matrix (192,200 bytes),
  # kcv_engine() two.
  withr::local_options(fastkcv.max_gb = 3e-4)
  meuse <- meuse_df()
  expect_no_error(kcv(log(zinc) ~ 1, ~ x + y, meuse, model = m_ok))
  C <- covariance_matrix(m_ok, as.matrix(meuse[, c("x", "y")]))
  expect_error(
    kcv_engine(C, matrix(1, 155, 1), log(meuse$zinc), seq_len(155)),
    "n = 155 .*2 dense n x n matrices .*fastkcv.max_gb = 0.00039\\)"
  )
  expect_error(kcv_engine_prepare(C, matrix(1, 155, 1)), "fastkcv.max_gb")
  withr::local_options(fastkcv.max_gb = 1e-4)
  expect_error(kcv(log(zinc) ~ 1, ~ x + y, meuse, model = m_ok),
    "n = 155 .*1 dense n x n matrix .*memory budget")
  expect_error(kcv_prepare(log(zinc) ~ 1, ~ x + y, meuse, model = m_ok),
    "fastkcv.max_gb")
  expect_error(err_memory(1e5, 2), "about 160 GB.*could not be allocated")
  expect_error(err_memory(28285, 1), "about 6.5 GB")
})

test_that("non-finite model parameters and non-positive kappa are rejected", {
  meuse <- meuse_df()
  # kappa <= 0 used to give a silent pure-nugget model (lgamma(0) = Inf
  # zeroes every correlation); krige.cv returns NA there.
  for (k in c(0, -1)) {
    expect_error(
      kcv(log(zinc) ~ 1, ~ x + y, meuse,
        model = gstat::vgm(0.6, "Mat", 300, 0.05, kappa = k)
      ),
      "Mat or Ste component of `model` needs a positive kappa"
    )
  }
  ste <- gstat::vgm(0.6, "Ste", 300, kappa = 1)
  ste$kappa <- 0
  expect_error(validate_model(ste), "positive kappa")
  for (col in c("psill", "range", "anis1", "ang1")) {
    bad <- m_ok
    bad[[col]][2] <- if (col == "range") Inf else NA
    expect_error(validate_model(bad), "missing or non-finite parameters")
  }
})

test_that("a trend without terms is rejected", {
  meuse <- meuse_df()
  for (f in list(log(zinc) ~ 0, log(zinc) ~ -1)) {
    expect_error(kcv(f, ~ x + y, meuse, model = m_ok), "trend has no terms")
  }
  expect_error(kcv_engine_prepare(diag(3), matrix(0, 3, 0)),
    "trend has no terms")
})

test_that("the engine rejects missing values in y, F and beta", {
  set.seed(1)
  xy <- matrix(runif(40), ncol = 2)
  C <- exp(-as.matrix(dist(xy)) / 0.3)
  F <- cbind(1, rnorm(20))
  y <- rnorm(20)
  y_na <- replace(y, 3, NA)
  # One NA used to turn every prediction into NA.
  expect_error(kcv_engine(C, F, y_na, rep(1:4, 5)), "`y` contains missing")
  p <- kcv_engine_prepare(C, F)
  expect_error(kcv_engine_run(p, replace(y, 2, Inf), 1:20), "non-finite")
  F_na <- F
  F_na[2, 2] <- NA
  expect_error(kcv_engine_prepare(C, F_na), "`F` must be numeric")
  expect_error(kcv_engine_prepare(C, F, beta = c(1, NA)), "`beta` must be")
  expect_error(kcv(log(zinc) ~ 1, ~ x + y, meuse_df(), model = m_ok,
    beta = "a"), "`beta` must be numeric")
})

test_that("an asymmetric covariance matrix is rejected", {
  set.seed(2)
  xy <- matrix(runif(40), ncol = 2)
  C <- exp(-as.matrix(dist(xy)) / 0.3)
  Ca <- C
  Ca[2, 1] <- 5
  expect_error(kcv_engine(Ca, matrix(1, 20, 1), rnorm(20), 1:20),
    "`C` is not symmetric")
  expect_error(kcv_engine_prepare(Ca, matrix(1, 20, 1)), "not symmetric")
  Cn <- C
  Cn[20, 3] <- NaN
  expect_error(kcv_engine_prepare(Cn, matrix(1, 20, 1)), "not symmetric")
  # Rounding-level asymmetry (far below 1e-8 times the diagonal) is fine.
  Cr <- C
  Cr[2, 1] <- C[2, 1] + 1e-12
  expect_no_error(kcv_engine_prepare(Cr, matrix(1, 20, 1)))
})

test_that("n above the 32-bit matrix limit is rejected with a clear error", {
  withr::local_options(fastkcv.max_gb = 1e6)
  expect_error(check_memory(65536, 1), "largest supported n = 65535")
  expect_no_error(check_memory(65535, 1))
})

test_that("running out of memory after the factorization is reported", {
  expect_error(raise_status(6, NA, 100, 0),
    "Out of memory while forming the blocks.*smaller folds")
  expect_error(raise_buffer_status(6, 0, NULL, 50, 100),
    "Out of memory while forming the blocks.*`targets` or a smaller")
  # The prepare stage keeps the n x n message.
  expect_error(raise_status(6, NA, 100, 1), "n = 100 .*could not be allocated")
})
