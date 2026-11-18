sim_setup <- function(seed = 123, n = 80) {
  m <- gstat::vgm(1, "Exp", 250, 0.2)
  d <- simulate_field(n, m, seed)
  xy <- as.matrix(d[, c("x", "y")])
  list(
    C = covariance_matrix(m, xy),
    F_ok = matrix(1, n, 1),
    F_uk = cbind(1, d$x, d$cov1),
    y = d$z,
    n = n
  )
}

test_that("kcv_engine matches the naive kcv_reference (OK, UK, SK)", {
  s <- sim_setup()
  folds_list <- list(
    loo = seq_len(s$n),
    k7 = withr::with_seed(1, sample(rep(1:7, length.out = s$n)))
  )
  for (folds in folds_list) {
    for (F in list(s$F_ok, s$F_uk)) {
      a <- kcv_engine(s$C, F, s$y, folds)
      b <- kcv_reference(s$C, F, s$y, folds)
      expect_equal(a$pred, b$pred, tolerance = tol)
      expect_equal(a$var, b$var, tolerance = tol)
    }
    a <- kcv_engine(s$C, s$F_uk, s$y, folds, beta = c(2, 0.001, 0.5))
    b <- kcv_reference(s$C, s$F_uk, s$y, folds, beta = c(2, 0.001, 0.5))
    expect_equal(a$pred, b$pred, tolerance = tol)
    expect_equal(a$var, b$var, tolerance = tol)
  }
})

test_that("matrix-valued y equals column-by-column runs", {
  s <- sim_setup()
  Y <- cbind(s$y, withr::with_seed(2, stats::rnorm(s$n)), s$y^2)
  folds <- rep(1:4, length.out = s$n)
  all <- kcv_engine(s$C, s$F_uk, Y, folds)
  for (j in seq_len(ncol(Y))) {
    one <- kcv_engine(s$C, s$F_uk, Y[, j], folds)
    expect_equal(all$pred[, j], drop(one$pred), tolerance = tol)
    expect_equal(all$var, one$var)
  }
  ref <- kcv_reference(s$C, s$F_uk, Y, folds)
  expect_equal(all$pred, ref$pred, tolerance = tol)
})
