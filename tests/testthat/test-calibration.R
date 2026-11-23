# Statistical validity: when the data are simulated from the kriging model,
# each cross-validation residual divided by its kriging standard error is
# exactly standard normal, for every fold design and for ordinary and
# universal kriging (the variance includes the trend-estimation term). So
# the mean squared z-score is 1 and the 90% interval covers 90% of the
# observations, on average over fields. Agreement with krige.cv() cannot
# catch a convention that gstat shares; this test can.

test_that("z-scores and coverage are calibrated under the model", {
  n <- 120
  n_fields <- 500
  model <- gstat::vgm(1, "Exp", 200, 0.1)
  xy <- withr::with_seed(11, matrix(stats::runif(2 * n, 0, 1000), ncol = 2))
  d <- data.frame(x = xy[, 1], y = xy[, 2])
  L <- t(chol(covariance_matrix(model, xy)))
  E <- withr::with_seed(12, L %*% matrix(stats::rnorm(n * n_fields), n))
  trends <- list(OK = rep(5, n), UK = 5 + d$x / 500 - d$y / 1000)
  formulas <- list(OK = z ~ 1, UK = z ~ x + y)
  blocks <- withr::with_seed(13, fold_blocks(xy, blocks = 4, k = 5))
  q90 <- stats::qnorm(0.95)
  for (tr in names(trends)) {
    Y <- trends[[tr]] + E
    d$z <- Y[, 1]
    p <- kcv_prepare(formulas[[tr]], ~ x + y, d, model = model)
    runs <- list(
      LOO = kcv_run(p, y = Y),
      blocks = kcv_run(p, folds = blocks, y = Y),
      buffer = kcv_run(p, buffer = 100, y = Y)
    )
    expect_gt(mean(runs$buffer$n_heldout), 2)
    for (s in names(runs)) {
      out <- runs[[s]]
      z <- matrix(out$zscore, n)
      msdr <- colMeans(z^2)
      cover <- colMeans(abs(z) <= q90)
      lab <- paste(tr, s)
      # Within 4 Monte Carlo standard errors over the independent fields.
      expect_lt(abs(mean(msdr) - 1), 4 * stats::sd(msdr) / sqrt(n_fields),
        label = paste(lab, "mean squared z-score"))
      expect_lt(abs(mean(cover) - 0.9), 4 * stats::sd(cover) / sqrt(n_fields),
        label = paste(lab, "90% coverage"))
    }
  }
})
