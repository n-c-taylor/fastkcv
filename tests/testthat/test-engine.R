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

engine_folds <- function(n) {
  list(
    loo = seq_len(n),
    k7 = withr::with_seed(1, sample(rep(1:7, length.out = n))),
    # uneven sizes, one fold of size 1, labels not 1..k
    uneven = c(rep(10L, 30), 99L, rep(3L, n - 41), rep(-2L, 10)),
    # character labels
    chr = rep(c("b", "a", "c"), length.out = n),
    # one fold larger than n / 2
    big = c(rep(1L, 50), rep(2L, n - 50))
  )
}

expect_engines_agree <- function(C, F, y, folds, beta = NULL) {
  a <- kcv_engine(C, F, y, folds, beta = beta)
  b <- kcv_engine_r(C, F, y, folds, beta = beta)
  r <- kcv_reference(C, F, y, folds, beta = beta)
  expect_equal(a$pred, b$pred, tolerance = tol)
  expect_close(a$pred, b$pred)
  expect_equal(a$var, b$var, tolerance = tol)
  expect_close(a$var, b$var)
  expect_equal(a$pred, r$pred, tolerance = tol, ignore_attr = TRUE)
  expect_close(a$pred, r$pred)
  expect_equal(a$var, r$var, tolerance = tol)
  expect_close(a$var, r$var)
  expect_identical(a$folds, folds)
}

test_that("C++ kcv_engine matches kcv_engine_r and kcv_reference (OK, UK, SK)", {
  s <- sim_setup()
  for (nm in names(engine_folds(s$n))) {
    folds <- engine_folds(s$n)[[nm]]
    expect_engines_agree(s$C, s$F_ok, s$y, folds)
    expect_engines_agree(s$C, s$F_uk, s$y, folds)
    expect_engines_agree(s$C, s$F_uk, s$y, folds, beta = c(2, 0.001, 0.5))
    expect_engines_agree(s$C, s$F_ok, s$y, folds, beta = 2)
  }
})

test_that("matrix-valued y equals column-by-column runs", {
  s <- sim_setup()
  Y <- cbind(s$y, withr::with_seed(2, stats::rnorm(s$n)), s$y^2)
  for (folds in list(rep(1:4, length.out = s$n), seq_len(s$n))) {
    all <- kcv_engine(s$C, s$F_uk, Y, folds)
    for (j in seq_len(ncol(Y))) {
      one <- kcv_engine(s$C, s$F_uk, Y[, j], folds)
      expect_equal(all$pred[, j], drop(one$pred), tolerance = tol)
      expect_close(all$pred[, j], drop(one$pred))
      expect_equal(all$var, one$var)
    }
    expect_engines_agree(s$C, s$F_uk, Y, folds)
    expect_engines_agree(s$C, s$F_uk, Y, folds, beta = c(2, 0.001, 0.5))
  }
})

test_that("LOO fast path equals the general fold loop", {
  s <- sim_setup()
  loo <- kcv_engine(s$C, s$F_uk, s$y, seq_len(s$n))
  # Same singleton folds, but with one extra observation paired up so the
  # general loop runs; compare the untouched singletons.
  folds <- seq_len(s$n)
  folds[2] <- 1L
  gen <- kcv_engine(s$C, s$F_uk, s$y, folds)
  keep <- -(1:2)
  expect_equal(gen$pred[keep, ], loo$pred[keep, ], tolerance = tol)
  expect_close(gen$pred[keep, ], loo$pred[keep, ])
  expect_equal(gen$var[keep], loo$var[keep], tolerance = tol)
  expect_close(gen$var[keep], loo$var[keep])
  # Permuted labels give the same LOO result.
  perm <- withr::with_seed(3, sample(s$n))
  expect_equal(kcv_engine(s$C, s$F_uk, s$y, perm)$pred, loo$pred, tolerance = tol)
  expect_close(kcv_engine(s$C, s$F_uk, s$y, perm)$pred, loo$pred)
})

test_that("kcv_engine keeps the dimnames of a matrix y", {
  s <- sim_setup(n = 20)
  Y <- matrix(c(s$y, -s$y), ncol = 2,
    dimnames = list(paste0("s", 1:20), c("a", "b")))
  res <- kcv_engine(s$C, s$F_ok, Y, rep(1:4, 5))
  expect_identical(dimnames(res$pred), dimnames(Y))
})
