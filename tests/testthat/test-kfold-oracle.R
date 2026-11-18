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
