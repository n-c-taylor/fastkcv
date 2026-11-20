test_that("k-fold with a fixed fold vector matches krige.cv", {
  meuse <- meuse_df()
  n <- nrow(meuse)
  folds <- withr::with_seed(7, sample(rep(1:5, length.out = n)))
  # Unequal fold sizes, including a fold of size 1.
  uneven <- folds
  uneven[uneven == 5] <- 4
  uneven[1] <- 9
  for (fv in list(folds, uneven)) {
    for (fn in names(oracle_formulas())) {
      f <- oracle_formulas()[[fn]]
      m <- oracle_models()$NugSph
      expect_cv_equal(
        kcv(f, ~ x + y, meuse, model = m, nfold = fv),
        krige_cv(f, ~ x + y, meuse, model = m, nfold = fv)
      )
      expect_cv_equal(
        kcv(f, ~ x + y, meuse, model = m, folds = fv),
        krige_cv(f, ~ x + y, meuse, model = m, nfold = fv)
      )
    }
  }
})

test_that("random nfold draws the same folds as krige.cv for a given seed", {
  meuse <- meuse_df()
  m <- oracle_models()$Exp
  ours <- withr::with_seed(11, kcv(log(zinc) ~ 1, ~ x + y, meuse, model = m, nfold = 5))
  gs <- withr::with_seed(11, krige_cv(log(zinc) ~ 1, ~ x + y, meuse, model = m, nfold = 5))
  expect_cv_equal(ours, gs)
})

test_that("k-fold matches krige.cv for every oracle model and anisotropy", {
  meuse <- meuse_df()
  fv <- withr::with_seed(13, sample(rep(1:6, length.out = nrow(meuse))))
  models <- c(oracle_models(),
    Aniso = list(gstat::vgm(0.4, "Exp", 500, 0.1, anis = c(30, 0.5),
      add.to = gstat::vgm(0.2, "Gau", 200, anis = c(120, 0.7))))
  )
  for (m in models) {
    for (f in oracle_formulas()) {
      expect_cv_equal(
        kcv(f, ~ x + y, meuse, model = m, folds = fv),
        krige_cv(f, ~ x + y, meuse, model = m, nfold = fv)
      )
    }
    expect_cv_equal(
      kcv(log(zinc) ~ sqrt(dist), ~ x + y, meuse, model = m,
        beta = c(6.9, -2.2), folds = fv),
      krige_cv(log(zinc) ~ sqrt(dist), ~ x + y, meuse, model = m,
        beta = c(6.9, -2.2), nfold = fv)
    )
  }
})
