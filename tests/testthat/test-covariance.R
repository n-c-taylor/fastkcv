# Oracle: gstat::variogramLine() with the lag direction passed via `dir`, so
# gstat itself applies the (possibly anisotropic) model to each pair.
gstat_pairwise <- function(model, xy) {
  n <- nrow(xy)
  out <- matrix(NA_real_, n, n)
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      h <- xy[i, ] - xy[j, ]
      d <- sqrt(sum(h^2))
      dir <- if (d > 0) c(h / d, 0) else c(1, 0, 0)
      out[i, j] <- gstat::variogramLine(
        model, dist_vector = d, dir = dir, covariance = TRUE
      )$gamma
    }
  }
  out
}

xy <- withr::with_seed(42, matrix(stats::runif(30, 0, 600), ncol = 2))

test_that("covariance matches variogramLine for each supported model", {
  models <- list(
    Nug = gstat::vgm(0.3, "Nug", 0),
    Exp = gstat::vgm(1, "Exp", 200),
    Sph = gstat::vgm(1, "Sph", 400),
    Gau = gstat::vgm(1, "Gau", 150, 0.1),
    Mat05 = gstat::vgm(1, "Mat", 200, kappa = 0.5),
    Mat15 = gstat::vgm(1, "Mat", 200, kappa = 1.5),
    Mat25 = gstat::vgm(1, "Mat", 200, kappa = 2.5),
    Ste = gstat::vgm(1, "Ste", 200, kappa = 1.2),
    NugExpSph = gstat::vgm(0.5, "Sph", 500, add.to = gstat::vgm(0.4, "Exp", 100, 0.1))
  )
  for (nm in names(models)) {
    expect_equal(
      covariance_matrix(models[[nm]], xy),
      gstat_pairwise(models[[nm]], xy),
      tolerance = tol, label = nm
    )
  }
})

test_that("diagonal is the total sill including the nugget", {
  m <- gstat::vgm(0.5, "Sph", 900, 0.1)
  expect_equal(diag(covariance_matrix(m, xy)), rep(0.6, nrow(xy)))
})

test_that("covariance matches variogramLine under 2D geometric anisotropy", {
  for (ang in c(0, 45, 90)) {
    for (ratio in c(1, 0.5, 0.2)) {
      m <- gstat::vgm(1, "Exp", 300, 0.1, anis = c(ang, ratio))
      expect_equal(
        covariance_matrix(m, xy), gstat_pairwise(m, xy),
        tolerance = tol, label = paste("anis", ang, ratio)
      )
    }
  }
})

test_that("nested components with different anisotropies are handled", {
  m <- gstat::vgm(0.6, "Sph", 400, anis = c(30, 0.5),
    add.to = gstat::vgm(0.4, "Exp", 200, 0.1, anis = c(120, 0.3)))
  expect_equal(covariance_matrix(m, xy), gstat_pairwise(m, xy), tolerance = tol)
})

test_that("cross-covariance between two coordinate sets", {
  m <- gstat::vgm(1, "Exp", 200, 0.1)
  C <- covariance_matrix(m, xy)
  expect_equal(covariance_matrix(m, xy[1:4, ], xy[5:15, ]), C[1:4, 5:15])
})
