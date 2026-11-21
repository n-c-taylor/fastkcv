# Most comparisons with gstat are skipped on CRAN: nearly all of their time
# is spent in gstat itself. One comparison per scheme runs everywhere.

test_that("LOO matches krige.cv on meuse (universal kriging)", {
  meuse <- meuse_df()
  m <- oracle_models()$NugSph
  f <- oracle_formulas()$UK
  expect_cv_equal(
    kcv(f, ~ x + y, meuse, model = m),
    krige_cv(f, ~ x + y, meuse, model = m)
  )
})

test_that("LOO matches krige.cv on meuse for each model, OK and UK", {
  skip_on_cran()
  meuse <- meuse_df()
  for (mn in names(oracle_models())) {
    for (fn in names(oracle_formulas())) {
      m <- oracle_models()[[mn]]
      f <- oracle_formulas()[[fn]]
      ours <- kcv(f, ~ x + y, meuse, model = m)
      gs <- krige_cv(f, ~ x + y, meuse, model = m)
      expect_cv_equal(ours, gs)
      expect_identical(names(ours), names(gs))
    }
  }
})

test_that("LOO matches krige.cv for Cir, Pen, Bes, Wav and Hol", {
  skip_on_cran()
  meuse <- meuse_df()
  for (m in oracle_models_extra()) {
    for (f in oracle_formulas()) {
      expect_cv_equal(
        kcv(f, ~ x + y, meuse, model = m),
        krige_cv(f, ~ x + y, meuse, model = m)
      )
    }
  }
})

test_that("LOO matches krige.cv for sp and sf inputs", {
  skip_on_cran()
  m <- oracle_models()$NugSph
  f <- oracle_formulas()$UK
  expect_cv_equal(kcv(f, meuse_sp(), model = m), krige_cv(f, meuse_sp(), model = m))
  expect_cv_equal(kcv(f, meuse_sf(), model = m), krige_cv(f, meuse_sf(), model = m))
})

test_that("LOO matches krige.cv with anisotropy and a Gaussian model", {
  skip_on_cran()
  meuse <- meuse_df()
  m <- gstat::vgm(0.5, "Sph", 900, 0.1, anis = c(40, 0.5))
  expect_cv_equal(
    kcv(log(zinc) ~ 1, ~ x + y, meuse, model = m),
    krige_cv(log(zinc) ~ 1, ~ x + y, meuse, model = m)
  )
  g <- gstat::vgm(0.5, "Gau", 300, 0.1)
  expect_cv_equal(
    kcv(log(zinc) ~ 1, ~ x + y, meuse, model = g),
    krige_cv(log(zinc) ~ 1, ~ x + y, meuse, model = g)
  )
})
