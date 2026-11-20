# Most comparisons with gstat are skipped on CRAN: nearly all of their time
# is spent in gstat itself. One comparison per scheme runs everywhere.

test_that("simple kriging with a known trend matches krige.cv", {
  meuse <- meuse_df()
  m <- oracle_models()$NugSph
  expect_cv_equal(
    kcv(log(zinc) ~ sqrt(dist), ~ x + y, meuse, model = m, beta = c(6.9, -2.2)),
    krige_cv(log(zinc) ~ sqrt(dist), ~ x + y, meuse, model = m,
      beta = c(6.9, -2.2))
  )
})

test_that("simple kriging with known beta matches krige.cv", {
  skip_on_cran()
  meuse <- meuse_df()
  cases <- list(
    list(f = log(zinc) ~ 1, beta = 5.9),
    list(f = log(zinc) ~ sqrt(dist), beta = c(6.9, -2.2))
  )
  for (cs in cases) {
    for (m in oracle_models()) {
      expect_cv_equal(
        kcv(cs$f, ~ x + y, meuse, model = m, beta = cs$beta),
        krige_cv(cs$f, ~ x + y, meuse, model = m, beta = cs$beta)
      )
    }
  }
  fv <- rep(1:6, length.out = nrow(meuse))
  m <- oracle_models()$NugSph
  expect_cv_equal(
    kcv(log(zinc) ~ 1, ~ x + y, meuse, model = m, beta = 5.9, nfold = fv),
    krige_cv(log(zinc) ~ 1, ~ x + y, meuse, model = m, beta = 5.9, nfold = fv)
  )
})
