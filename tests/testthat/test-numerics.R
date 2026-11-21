# Ill-conditioned covariance matrices: a Gaussian model without a nugget.
# fastkcv must agree with krige.cv as long as the matrix is factorizable;
# both degrade together as the conditioning worsens (near-duplicate points
# give meaningless predictions in both, so those are not compared here).

test_that("Gaussian model without nugget matches krige.cv", {
  skip_on_cran()
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
  skip_on_cran()
  d <- meuse_df()[1:30, ]
  d$x[2] <- d$x[1] + 1
  d$y[2] <- d$y[1]
  g <- gstat::vgm(0.6, "Gau", 300)
  expect_cv_equal(
    kcv(log(zinc) ~ 1, ~ x + y, d, model = g),
    krige_cv(log(zinc) ~ 1, ~ x + y, d, model = g)
  )
})

# A quadratic trend in raw UTM-scale coordinates: the columns of F range
# from 1 to about 2e13, and F'C^-1 F had rcond 1e-35, so the C++ engine used
# to solve with it through Armadillo's approximate (SVD) fallback and was
# off by 0.2. OK/UK depend on F only through its column space, so the
# oracle is the same fit with centered and scaled coordinates (well
# conditioned). krige.cv() on the raw coordinates is itself about 1e-7 off
# that, and kcv_reference() with the raw F fails in solve(), so neither is
# used with the raw F.
test_that("UTM-scale coordinates with a quadratic trend match the centered fit", {
  set.seed(45)
  n <- 150
  u <- cbind(runif(n), runif(n))
  L <- t(chol(exp(-as.matrix(dist(u)) / 0.3) + diag(0.1, n)))
  d <- data.frame(x = 5e5 + 1e5 * u[, 1], y = 4.5e6 + 1e5 * u[, 2])
  d$z <- 2 + 1.5 * u[, 1] - 0.8 * u[, 2] + 0.5 * u[, 1]^2 +
    drop(L %*% rnorm(n))
  d$xs <- (d$x - mean(d$x)) / 1e5
  d$ys <- (d$y - mean(d$y)) / 1e5
  v <- gstat::vgm(1, "Exp", 30000, 0.1)
  f_raw <- z ~ x + y + I(x^2) + I(y^2) + I(x * y)
  f_cen <- z ~ xs + ys + I(xs^2) + I(ys^2) + I(xs * ys)
  F_raw <- model.matrix(f_raw, d)
  F_cen <- model.matrix(f_cen, d)
  C <- covariance_matrix(v, as.matrix(d[, c("x", "y")]))

  # Leave-one-out against gstat and the refit oracle on the centered trend
  # (measured: about 5e-12).
  loo <- kcv(f_raw, ~ x + y, d, model = v)
  g <- krige_cv(f_cen, ~ x + y, d, model = v)
  expect_close(loo$var1.pred, g$var1.pred, 1e-8, "LOO predictions")
  expect_close(loo$var1.var, g$var1.var, 1e-8, "LOO variances")
  ref <- kcv_reference(C, F_cen, d$z, seq_len(n))
  expect_close(loo$var1.pred, ref$pred, 1e-8, "LOO predictions (refit)")
  expect_close(loo$var1.var, ref$var, 1e-8, "LOO variances (refit)")

  # k-fold and buffered leave-one-out from the same prepared object.
  p <- kcv_prepare(f_raw, ~ x + y, d, model = v)
  folds <- rep(1:10, length.out = n)
  k10 <- kcv_run(p, folds = folds)
  ref10 <- kcv_reference(C, F_cen, d$z, folds)
  expect_close(k10$var1.pred, ref10$pred, 1e-8, "10-fold predictions")
  expect_close(k10$var1.var, ref10$var, 1e-8, "10-fold variances")
  buf <- kcv_run(p, buffer = 15000)
  sets <- buffer_sets_r(as.matrix(d[, c("x", "y")]), 15000)
  refb <- kcv_buffer_reference(C, F_cen, d$z, sets)
  expect_gt(mean(buf$n_heldout), 2)
  expect_close(buf$var1.pred, refb$pred, 1e-8, "buffered predictions")
  expect_close(buf$var1.var, refb$var, 1e-8, "buffered variances")

  # The matrix engine takes the same path.
  eng <- kcv_engine(C, F_raw, d$z, seq_len(n))
  expect_close(eng$pred, ref$pred, 1e-8, "engine predictions")
  expect_close(eng$var, ref$var, 1e-8, "engine variances")

  # The C++ backstop: given the raw F directly, F'C^-1 F is numerically
  # singular (status 2) instead of a silent approximate solution; the
  # orthonormal basis is accepted.
  expect_identical(prepare_matrix_cpp(C, F_raw, FALSE)$status, 2L)
  expect_identical(prepare_matrix_cpp(C, qr.Q(qr(F_raw)), FALSE)$status, 0L)
  expect_error(raise_status(2, NA, n, 1), "numerically singular")

  # Simple kriging keeps the raw F (beta is in its parametrization).
  beta <- c(2, rep(0, 5))
  sk <- kcv(f_raw, ~ x + y, d, model = v, beta = beta)
  refsk <- kcv_reference(C, F_raw, d$z, seq_len(n), beta = beta)
  expect_close(sk$var1.pred, refsk$pred, 1e-8, "SK predictions")
  expect_close(sk$var1.var, refsk$var, 1e-8, "SK variances")
})
