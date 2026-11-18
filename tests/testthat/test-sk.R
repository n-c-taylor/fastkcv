test_that("simple kriging with known beta matches krige.cv", {
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
