# Ill-conditioned covariance matrices: a Gaussian model without a nugget.
# fastkcv must agree with krige.cv as long as the matrix is factorizable;
# both degrade together as the conditioning worsens (near-duplicate points
# give meaningless predictions in both, so those are not compared here).

test_that("Gaussian model without nugget matches krige.cv", {
  meuse <- meuse_df()
  for (rg in c(150, 300)) {
    g <- gstat::vgm(0.6, "Gau", rg)
    for (f in oracle_formulas()) {
      expect_cv_equal(
        kcv(f, ~ x + y, meuse, model = g),
        krige_cv(f, ~ x + y, meuse, model = g)
      )
    }
  }
  # rcond(C) is about 7e-8 at range 300, so these cases do exercise a
  # poorly conditioned factorization.
  C <- covariance_matrix(gstat::vgm(0.6, "Gau", 300),
    as.matrix(meuse[, c("x", "y")]))
  expect_lt(rcond(C), 1e-6)
})

test_that("a close pair without nugget still matches krige.cv", {
  d <- meuse_df()[1:30, ]
  d$x[2] <- d$x[1] + 1
  d$y[2] <- d$y[1]
  g <- gstat::vgm(0.6, "Gau", 300)
  expect_cv_equal(
    kcv(log(zinc) ~ 1, ~ x + y, d, model = g),
    krige_cv(log(zinc) ~ 1, ~ x + y, d, model = g)
  )
})
