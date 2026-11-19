metrics_example <- function() {
  # z-scores 1, -1, 0.5, 1.5, -2, 0.
  data.frame(
    var1.pred = 0,
    var1.var = c(1, 4, 1, 4, 0.25, 1),
    observed = 0,
    residual = c(1, -2, 0.5, 3, -1, 0),
    zscore = 99, # ignored: z is recomputed from residual and var1.var
    fold = c(1, 1, 2, 2, 2, 3)
  )
}

test_that("kcv_metrics() matches hand-computed values", {
  m <- kcv_metrics(metrics_example())
  expect_s3_class(m, "data.frame")
  expect_identical(nrow(m), 1L)
  expect_named(m, c(
    "n", "rmse", "mae", "bias", "mean_z", "var_z",
    "coverage_50", "coverage_80", "coverage_90", "coverage_95"
  ))
  expect_identical(m$n, 6L)
  expect_equal(m$rmse, sqrt(15.25 / 6), tolerance = tol)
  expect_equal(m$mae, 1.25, tolerance = tol)
  expect_equal(m$bias, 0.25, tolerance = tol)
  expect_equal(m$mean_z, 0, tolerance = tol)
  expect_equal(m$var_z, 8.5 / 5, tolerance = tol)
  # Half-widths 0.674, 1.282, 1.645 and 1.960 standard errors.
  expect_equal(m$coverage_50, 2 / 6, tolerance = tol)
  expect_equal(m$coverage_80, 4 / 6, tolerance = tol)
  expect_equal(m$coverage_90, 5 / 6, tolerance = tol)
  expect_equal(m$coverage_95, 5 / 6, tolerance = tol)
})

test_that("interval bounds are inclusive and levels name the columns", {
  # qnorm((1 + L) / 2) is exactly 1 and 2 for these levels, so |z| = 1 and
  # |z| = 2 lie exactly on the bounds.
  lv <- 2 * stats::pnorm(c(1, 2)) - 1
  expect_identical(stats::qnorm((1 + lv) / 2), c(1, 2))
  m <- kcv_metrics(metrics_example(), level = lv)
  expect_identical(
    names(m)[7:8],
    paste0("coverage_", format(100 * lv, trim = TRUE, drop0trailing = TRUE))
  )
  expect_equal(m[[7]], 4 / 6, tolerance = tol) # 1, -1, 0.5, 0
  expect_equal(m[[8]], 1, tolerance = tol) # -2 is on the bound
  expect_named(
    kcv_metrics(metrics_example(), level = c(0.5, 0.99)),
    c("n", "rmse", "mae", "bias", "mean_z", "var_z", "coverage_50",
      "coverage_99")
  )
})

test_that("by_fold gives one row per fold equal to the subset metrics", {
  x <- metrics_example()
  m <- kcv_metrics(x, by_fold = TRUE)
  expect_identical(m$fold, c(1, 2, 3))
  expect_identical(m$n, c(2L, 3L, 1L))
  for (f in 1:3) {
    expect_equal(
      m[f, -1],
      kcv_metrics(x[x$fold == f, ]),
      tolerance = tol,
      ignore_attr = TRUE
    )
  }
  # A single-observation fold has no z-score variance.
  expect_true(is.na(m$var_z[3]))
  expect_equal(m$bias[2], 2.5 / 3, tolerance = tol)
})

test_that("kcv_metrics() works on data.frame, sf and sp output", {
  v <- gstat::vgm(0.5, "Sph", 900, 0.1)
  folds <- rep(1:5, length.out = 155)
  a <- kcv(log(zinc) ~ 1, ~ x + y, meuse_df(), model = v, folds = folds)
  b <- kcv(log(zinc) ~ 1, meuse_sf(), model = v, folds = folds)
  s <- kcv(log(zinc) ~ 1, meuse_sp(), model = v, folds = folds)
  ma <- kcv_metrics(a)
  expect_equal(ma$rmse, sqrt(mean(a$residual^2)), tolerance = tol)
  expect_equal(ma$mean_z, mean(a$zscore), tolerance = tol)
  expect_equal(ma$var_z, stats::var(a$zscore), tolerance = tol)
  expect_equal(kcv_metrics(b), ma, tolerance = tol)
  expect_equal(kcv_metrics(s), ma, tolerance = tol)
  expect_equal(kcv_metrics(s, by_fold = TRUE), kcv_metrics(a, by_fold = TRUE),
    tolerance = tol
  )
  # Buffered output has one fold per point.
  bu <- kcv(log(zinc) ~ 1, ~ x + y, meuse_df(), model = v, buffer = 200)
  expect_identical(nrow(kcv_metrics(bu, by_fold = TRUE)), 155L)
})

test_that("kcv_metrics() rejects invalid input", {
  x <- metrics_example()
  expect_error(kcv_metrics(as.matrix(x)), "output of kcv")
  expect_error(kcv_metrics(x[, -4]), "lacks column\\(s\\) `residual`")
  expect_error(kcv_metrics(x[, -6], by_fold = TRUE), "`fold`")
  expect_no_error(kcv_metrics(x[, -6]))
  expect_error(kcv_metrics(x, level = 1), "`level`")
  expect_error(kcv_metrics(x, level = c(0.5, NA)), "`level`")
  expect_error(kcv_metrics(x, level = numeric(0)), "`level`")
  expect_error(kcv_metrics(x, by_fold = NA), "`by_fold`")
  bad <- x
  bad$residual <- as.character(bad$residual)
  expect_error(kcv_metrics(bad), "numeric")
  expect_error(kcv_metrics(sp::SpatialPoints(meuse_sp())), "attribute table")
})
