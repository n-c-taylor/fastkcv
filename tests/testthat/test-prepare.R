m_prep <- gstat::vgm(0.5, "Sph", 900, 0.1)

test_that("kcv_run(kcv_prepare()) equals kcv() for OK, UK and SK", {
  meuse <- meuse_df()
  k7 <- withr::with_seed(11, sample(rep(1:7, length.out = nrow(meuse))))
  cases <- list(
    list(f = log(zinc) ~ 1, beta = NULL),
    list(f = log(zinc) ~ sqrt(dist), beta = NULL),
    list(f = log(zinc) ~ 1, beta = 5.9)
  )
  for (cs in cases) {
    p <- kcv_prepare(cs$f, ~ x + y, meuse, model = m_prep, beta = cs$beta)
    expect_s3_class(p, "kcv_prepared")
    expect_cv_equal(
      kcv_run(p),
      kcv(cs$f, ~ x + y, meuse, model = m_prep, beta = cs$beta)
    )
    expect_cv_equal(
      kcv_run(p, folds = k7),
      kcv(cs$f, ~ x + y, meuse, model = m_prep, beta = cs$beta, folds = k7)
    )
  }
})

test_that("one prepared object serves several fold schemes and input types", {
  for (loc in list(meuse_sf(), meuse_sp())) {
    p <- kcv_prepare(log(zinc) ~ sqrt(dist), loc, model = m_prep)
    schemes <- list(
      seq_len(155),
      rep(1:5, length.out = 155),
      withr::with_seed(4, sample(rep(1:25, length.out = 155)))
    )
    for (folds in schemes) {
      a <- kcv_run(p, nfold = folds)
      b <- kcv(log(zinc) ~ sqrt(dist), loc, model = m_prep, nfold = folds)
      expect_identical(class(a), class(b))
      expect_cv_equal(a, b)
    }
  }
})

test_that("random nfold draws the same folds as kcv() and krige.cv()", {
  meuse <- meuse_df()
  p <- kcv_prepare(log(zinc) ~ 1, ~ x + y, meuse, model = m_prep)
  a <- withr::with_seed(9, kcv_run(p, nfold = 5))
  b <- withr::with_seed(9, kcv(log(zinc) ~ 1, ~ x + y, meuse, model = m_prep,
    nfold = 5))
  g <- withr::with_seed(9, krige_cv(log(zinc) ~ 1, ~ x + y, meuse,
    model = m_prep, nfold = 5))
  expect_cv_equal(a, b)
  expect_cv_equal(a, g)
})

test_that("kcv_engine_prepare/run match the R engine, including T > 1", {
  d <- simulate_field(60, gstat::vgm(1, "Exp", 250, 0.2), seed = 21)
  xy <- as.matrix(d[, c("x", "y")])
  C <- covariance_matrix(gstat::vgm(1, "Exp", 250, 0.2), xy)
  F <- cbind(1, d$x, d$cov1)
  Y <- cbind(d$z, d$z^2, withr::with_seed(5, stats::rnorm(60)))
  C0 <- C
  p <- kcv_engine_prepare(C, F)
  ps <- kcv_engine_prepare(C, F, beta = c(2, 0, 0.3))
  expect_identical(C, C0)
  for (folds in list(seq_len(60), rep(1:6, 10), rep(c(1, 2), c(40, 20)))) {
    a <- kcv_engine_run(p, Y, folds)
    b <- kcv_engine_r(C, F, Y, folds)
    expect_equal(a$pred, b$pred, tolerance = tol)
    expect_equal(a$var, b$var, tolerance = tol)
    expect_identical(a, kcv_engine(C, F, Y, folds))
    a <- kcv_engine_run(ps, Y, folds)
    b <- kcv_engine_r(C, F, Y, folds, beta = c(2, 0, 0.3))
    expect_equal(a$pred, b$pred, tolerance = tol)
    expect_equal(a$var, b$var, tolerance = tol)
  }
})

test_that("prepared objects survive saveRDS()", {
  p <- kcv_prepare(log(zinc) ~ sqrt(dist), ~ x + y, meuse_df(), model = m_prep)
  f <- withr::local_tempfile(fileext = ".rds")
  saveRDS(p, f)
  q <- readRDS(f)
  folds <- rep(1:4, length.out = 155)
  expect_identical(kcv_run(q, folds = folds), kcv_run(p, folds = folds))
})

test_that("print method summarizes the prepared object", {
  meuse <- meuse_df()
  expect_output(print(kcv_prepare(log(zinc) ~ 1, ~ x + y, meuse,
    model = m_prep)), "ordinary kriging, n = 155, p = 1")
  expect_output(print(kcv_prepare(log(zinc) ~ dist, ~ x + y, meuse,
    model = m_prep)), "universal kriging, n = 155, p = 2")
  expect_output(print(kcv_prepare(log(zinc) ~ 1, ~ x + y, meuse,
    model = m_prep, beta = 6)), "simple kriging")
  expect_output(print(kcv_engine_prepare(diag(3), matrix(1, 3, 1))),
    "from a covariance matrix")
})

test_that("prepare/run argument errors", {
  meuse <- meuse_df()
  pe <- kcv_engine_prepare(diag(3), matrix(1, 3, 1))
  expect_error(kcv_run(pe), "kcv_engine_run")
  expect_error(kcv_run(list()), "kcv_prepared")
  expect_error(kcv_engine_run(list(), 1:3, 1:3), "kcv_prepared")
  expect_error(kcv_engine_run(pe, 1:4, 1:3), "`y` must")
  expect_error(kcv_engine_run(pe, 1:3, c(1, NA, 2)), "`folds` must")
  expect_error(kcv_engine_prepare(diag(3)[, 1:2], matrix(1, 3, 1)), "square")
  expect_error(kcv_engine_prepare(diag(3), matrix(1, 2, 1)), "rows")
  expect_error(kcv_engine_prepare(diag(3), matrix(1, 3, 1), beta = 1:2),
    "`beta` has length")
  expect_error(kcv_engine_prepare(matrix(1, 3, 3), matrix(1, 3, 1)),
    "not positive definite")
  expect_error(kcv_engine_prepare(diag(3), cbind(1, rep(1, 3))),
    "rank deficient")
  expect_error(kcv_prepare(log(zinc) ~ 1, ~ x + y, meuse, model = m_prep,
    beta = 1:2), "`beta` has length")
  p <- kcv_prepare(log(zinc) ~ 1, ~ x + y, meuse, model = m_prep)
  expect_error(kcv_run(p, folds = 1:3), "one non-missing value")
  expect_error(kcv_run(p, nfold = 1), "every observation")
})

test_that("prepared objects in the earlier storage layout are refused", {
  p <- kcv_prepare(log(zinc) ~ 1, ~ x + y, meuse_df(), model = m_prep)
  old <- p
  names(old)[names(old) == "Wt"] <- "W"
  expect_error(kcv_run(old), "earlier version of fastkcv")
  expect_error(kcv_engine_run(old, meuse_df()$zinc, 1:155),
    "earlier version of fastkcv")
})
