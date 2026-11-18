m_ok <- gstat::vgm(0.5, "Sph", 900, 0.1)

test_that("local-neighbourhood and block arguments are rejected", {
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

test_that("inputs above the size limit are rejected before allocation", {
  withr::local_options(fastkcv.max_n = 10)
  xy <- cbind(1:11, 0)
  C <- covariance_matrix(gstat::vgm(1, "Exp", 3, 0.1), xy)
  expect_error(kcv_engine(C, matrix(1, 11, 1), rnorm(11), 1:11),
    "n = 11 .*fastkcv.max_n")
  expect_error(kcv(log(zinc) ~ 1, ~ x + y, meuse_df(), model = m_ok),
    "n = 155 .*options\\(fastkcv.max_n = 155\\)")
  expect_error(err_memory(1e5), "about 160 GB.*could not be allocated")
  withr::local_options(fastkcv.max_n = 11)
  expect_no_error(kcv_engine(C, matrix(1, 11, 1), rnorm(11), 1:11))
})
