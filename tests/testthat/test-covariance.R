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
    Cir = gstat::vgm(1, "Cir", 400, 0.1),
    Pen = gstat::vgm(1, "Pen", 400),
    Bes = gstat::vgm(1, "Bes", 100),
    Wav = gstat::vgm(1, "Wav", 150, 0.1),
    Hol = gstat::vgm(1, "Hol", 80),
    NugExpSph = gstat::vgm(0.5, "Sph", 500, add.to = gstat::vgm(0.4, "Exp", 100, 0.1))
  )
  for (nm in names(models)) {
    expect_equal(
      covariance_matrix(models[[nm]], xy),
      gstat_pairwise(models[[nm]], xy),
      tolerance = tol, label = nm
    )
    expect_close(
      covariance_matrix(models[[nm]], xy),
      gstat_pairwise(models[[nm]], xy),
      label = nm
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
      expect_close(covariance_matrix(m, xy), gstat_pairwise(m, xy),
        label = paste("anis", ang, ratio)
      )
    }
  }
})

test_that("nested components with different anisotropies are handled", {
  m <- gstat::vgm(0.6, "Sph", 400, anis = c(30, 0.5),
    add.to = gstat::vgm(0.4, "Exp", 200, 0.1, anis = c(120, 0.3)))
  expect_equal(covariance_matrix(m, xy), gstat_pairwise(m, xy), tolerance = tol)
  expect_close(covariance_matrix(m, xy), gstat_pairwise(m, xy))
})

# The C++ covariance against the variogramLine() version (covariance_matrix_r)
# for every model, nested and anisotropic case, on more points.
xy2 <- withr::with_seed(7, matrix(stats::runif(120, 0, 800), ncol = 2))

cpp_cases <- function() {
  list(
    Nug = gstat::vgm(0.3, "Nug", 0),
    Exp = gstat::vgm(1, "Exp", 200),
    Sph = gstat::vgm(1, "Sph", 400),
    Gau = gstat::vgm(1, "Gau", 150, 0.1),
    Mat01 = gstat::vgm(1, "Mat", 200, kappa = 0.1),
    Mat05 = gstat::vgm(1, "Mat", 200, kappa = 0.5),
    Mat15 = gstat::vgm(1, "Mat", 200, kappa = 1.5),
    Mat25 = gstat::vgm(1, "Mat", 200, kappa = 2.5),
    Mat5 = gstat::vgm(1, "Mat", 100, kappa = 5),
    Mat10 = gstat::vgm(1, "Mat", 50, kappa = 10),
    Ste01 = gstat::vgm(1, "Ste", 200, kappa = 0.1),
    Ste12 = gstat::vgm(1, "Ste", 200, kappa = 1.2),
    Ste5 = gstat::vgm(1, "Ste", 300, kappa = 5),
    Ste10 = gstat::vgm(1, "Ste", 300, kappa = 10),
    NugExpSph = gstat::vgm(0.5, "Sph", 500,
      add.to = gstat::vgm(0.4, "Exp", 100, 0.1)),
    AnisNested = gstat::vgm(0.6, "Sph", 400, anis = c(30, 0.5),
      add.to = gstat::vgm(0.4, "Mat", 200, 0.1, kappa = 1.5,
        anis = c(120, 0.3))),
    AnisSte = gstat::vgm(1, "Ste", 300, 0.2, kappa = 0.8, anis = c(135, 0.4)),
    Cir = gstat::vgm(1, "Cir", 400, 0.1),
    Pen = gstat::vgm(1, "Pen", 400),
    Bes = gstat::vgm(1, "Bes", 100),
    Wav = gstat::vgm(1, "Wav", 150, 0.1),
    Hol = gstat::vgm(1, "Hol", 80),
    AnisPen = gstat::vgm(1, "Pen", 400, 0.1, anis = c(60, 0.3)),
    NestedBesWav = gstat::vgm(0.6, "Bes", 120,
      add.to = gstat::vgm(0.4, "Wav", 200, 0.1, anis = c(150, 0.5)))
  )
}

test_that("C++ covariance matches the variogramLine() version", {
  for (nm in names(cpp_cases())) {
    m <- cpp_cases()[[nm]]
    expect_equal(covariance_matrix(m, xy2), covariance_matrix_r(m, xy2),
      tolerance = tol, label = nm
    )
    expect_close(covariance_matrix(m, xy2), covariance_matrix_r(m, xy2),
      label = nm
    )
  }
  for (ang in c(0, 45, 90)) {
    for (ratio in c(1, 0.5, 0.2)) {
      for (type in c("Exp", "Sph", "Gau")) {
        m <- gstat::vgm(1, type, 300, 0.1, anis = c(ang, ratio))
        expect_equal(covariance_matrix(m, xy2), covariance_matrix_r(m, xy2),
          tolerance = tol, label = paste(type, ang, ratio)
        )
        expect_close(covariance_matrix(m, xy2), covariance_matrix_r(m, xy2),
          label = paste(type, ang, ratio)
        )
      }
    }
  }
})

test_that("C++ covariance is exactly symmetric with the total sill on the diagonal", {
  m <- cpp_cases()$AnisNested
  C <- covariance_matrix(m, xy2)
  expect_identical(C, t(C))
  expect_identical(diag(C), rep(sum(m$psill), nrow(xy2)))
})

test_that("the nugget applies only to coincident points", {
  # A tiny lag already drops the nugget, as in gstat's variogramLine().
  m <- gstat::vgm(1, "Exp", 200, 0.25)
  C <- covariance_matrix(m, rbind(c(0, 0), c(1e-9, 0), c(5, 5)))
  expect_equal(diag(C), rep(1.25, 3))
  expect_lt(C[1, 2], 1)
  expect_gt(C[1, 2], 1 - 1e-10)
})

test_that("Matern with large kappa stays finite where gstat returns NaN", {
  # gstat evaluates x^kappa * K_kappa(x) directly; for kappa = 50 and small
  # lags x^kappa underflows while K_kappa(x) overflows, giving NaN. The C++
  # version works on the log scale.
  m <- gstat::vgm(1, "Mat", 200, kappa = 50)
  xy3 <- cbind(c(0, 1e-3, 0.5, 10, 150), 0)
  C <- covariance_matrix(m, xy3)
  expect_true(all(is.finite(C)))
  expect_true(anyNA(covariance_matrix_r(m, xy3)))
  x <- abs(outer(xy3[, 1], xy3[, 1], "-")) / 200
  k <- 50
  ref <- exp(k * log(x) - (k - 1) * log(2) - lgamma(k) +
    log(besselK(x, k, expon.scaled = TRUE)) - x)
  ref[x == 0] <- 1
  # At x = 5e-6 even the scaled K_50 overflows; there C++ uses the small-x
  # series 1 - x^2 / (4 (kappa - 1)) + ...
  over <- !is.finite(ref)
  expect_true(any(over))
  ref[over] <- 1 - x[over]^2 / (4 * (k - 1))
  expect_equal(C, ref, tolerance = tol)
  expect_close(C, ref)
  expect_true(all(C[over] < 1))
  # Where gstat is finite it agrees.
  ok <- is.finite(covariance_matrix_r(m, xy3))
  expect_equal(C[ok], covariance_matrix_r(m, xy3)[ok], tolerance = tol)
  expect_close(C[ok], covariance_matrix_r(m, xy3)[ok])
})

test_that("Ste with large kappa approaches the Gaussian model", {
  # An independent check of the large-kappa evaluation (the test above
  # re-derives the same log-scale formula): Ste(kappa, a) uses
  # x = 2 sqrt(kappa) h / a and converges to Gau(a) as kappa grows, with an
  # error of order 1 / kappa (hence the loose tolerance). Mat with range
  # a / (2 sqrt(kappa)) is the same function.
  xy3 <- cbind(c(0, 1e-3, 0.5, 10, 50, 150, 300), 0)
  gau <- covariance_matrix(gstat::vgm(1, "Gau", 200), xy3)
  err <- vapply(c(50, 200), function(k) {
    ste <- covariance_matrix(gstat::vgm(1, "Ste", 200, kappa = k), xy3)
    mat <- covariance_matrix(
      gstat::vgm(1, "Mat", 200 / (2 * sqrt(k)), kappa = k), xy3
    )
    expect_equal(mat, ste, tolerance = tol)
    expect_close(mat, ste)
    max(abs(ste - gau))
  }, numeric(1))
  expect_lt(err[1], 0.25 / 50)
  expect_lt(err[2], 0.25 / 200)
  expect_lt(err[2], err[1] / 3)
})
