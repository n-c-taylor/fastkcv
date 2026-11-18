test_that("k-fold with n singleton folds equals LOO", {
  meuse <- meuse_df()
  m <- oracle_models()$NugSph
  f <- oracle_formulas()$UK
  loo <- kcv(f, ~ x + y, meuse, model = m)
  perm_ids <- withr::with_seed(3, sample(nrow(meuse)))
  kf <- kcv(f, ~ x + y, meuse, model = m, folds = perm_ids)
  for (col in cv_cols) expect_equal(kf[[col]], loo[[col]], tolerance = tol)
})

test_that("permuting rows permutes the output", {
  meuse <- meuse_df()
  m <- oracle_models()$Mat15
  f <- oracle_formulas()$UK
  folds <- rep(1:5, length.out = nrow(meuse))
  base <- kcv(f, ~ x + y, meuse, model = m, folds = folds)
  p <- withr::with_seed(5, sample(nrow(meuse)))
  perm <- kcv(f, ~ x + y, meuse[p, ], model = m, folds = folds[p])
  for (col in c(cv_cols, "fold", "x", "y")) {
    expect_equal(perm[[col]], base[[col]][p], tolerance = tol, label = col)
  }
})

test_that("residual and zscore identities hold", {
  cv <- kcv(log(zinc) ~ 1, ~ x + y, meuse_df(), model = oracle_models()$Sph)
  expect_equal(cv$residual, cv$observed - cv$var1.pred)
  expect_equal(cv$zscore, cv$residual / sqrt(cv$var1.var))
})

test_that("data.frame, sf and sp inputs agree and keep their class", {
  m <- oracle_models()$NugSph
  f <- oracle_formulas()$UK
  a <- kcv(f, ~ x + y, meuse_df(), model = m)
  b <- kcv(f, meuse_sf(), model = m)
  c <- kcv(f, meuse_sp(), model = m)
  d <- kcv(f, meuse_sp(), m) # model in third position, as krige.cv allows
  expect_s3_class(a, "data.frame")
  expect_false(inherits(a, "sf"))
  expect_identical(names(a), c(cv_cols, "fold", "x", "y"))
  expect_s3_class(b, "sf")
  expect_identical(sf::st_crs(b), sf::st_crs(meuse_sf()))
  expect_s4_class(c, "SpatialPointsDataFrame")
  expect_identical(sp::proj4string(c), sp::proj4string(meuse_sp()))
  expect_cv_equal(a, b)
  expect_cv_equal(a, c)
  expect_cv_equal(a, d)
})

test_that("verbose and debug.level are accepted and ignored", {
  expect_no_error(
    kcv(log(zinc) ~ 1, ~ x + y, meuse_df(), model = oracle_models()$Exp,
      verbose = FALSE, debug.level = 0)
  )
})
