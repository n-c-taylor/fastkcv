tol <- 1e-8

meuse_df <- function() {
  e <- new.env()
  utils::data("meuse", package = "sp", envir = e)
  e$meuse
}

meuse_sp <- function() {
  m <- meuse_df()
  sp::coordinates(m) <- ~ x + y
  sp::proj4string(m) <- sp::CRS("EPSG:28992")
  m
}

meuse_sf <- function() {
  sf::st_as_sf(meuse_df(), coords = c("x", "y"), crs = 28992, remove = FALSE)
}

oracle_models <- function() {
  list(
    Exp = gstat::vgm(0.7, "Exp", 400),
    Sph = gstat::vgm(0.59, "Sph", 897),
    Mat15 = gstat::vgm(0.6, "Mat", 300, kappa = 1.5),
    NugSph = gstat::vgm(0.5, "Sph", 900, 0.1)
  )
}

oracle_formulas <- function() {
  list(OK = log(zinc) ~ 1, UK = log(zinc) ~ sqrt(dist))
}

cv_cols <- c("var1.pred", "var1.var", "observed", "residual", "zscore")

expect_cv_equal <- function(a, b, tolerance = tol) {
  a <- as.data.frame(a)
  b <- as.data.frame(b)
  for (col in cv_cols) {
    expect_equal(a[[col]], b[[col]], tolerance = tolerance, label = col)
  }
  expect_equal(as.vector(a$fold), as.vector(b$fold))
}

krige_cv <- function(...) gstat::krige.cv(..., verbose = FALSE)

# Buffered leave-one-out by brute force with gstat: drop the hold-out set of
# point i and krige point i from the remaining data. Used by test-buffer.R
# and test-multi.R.
krige_buffer <- function(formula, data, model, sets, targets, beta = NULL) {
  out <- lapply(targets, function(i) {
    k <- gstat::krige(formula, ~ x + y, data[-sets[[i]], ], data[i, ],
      model = model, beta = beta, debug.level = 0
    )
    c(k$var1.pred, k$var1.var)
  })
  out <- do.call(rbind, out)
  list(pred = out[, 1], var = out[, 2])
}

simulate_field <- function(n, model, seed) {
  withr::with_seed(seed, {
    xy <- matrix(stats::runif(2 * n, 0, 1000), ncol = 2)
    C <- covariance_matrix(model, xy)
    z <- drop(crossprod(chol(C), stats::rnorm(n)))
    data.frame(x = xy[, 1], y = xy[, 2], cov1 = stats::rnorm(n), z = z + 2)
  })
}
