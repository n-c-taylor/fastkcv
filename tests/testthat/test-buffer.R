m_buf <- gstat::vgm(0.5, "Sph", 900, 0.1)

# Buffered leave-one-out by brute force with gstat: drop the hold-out set of
# point i and krige point i from the remaining data.
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

buffer_cases <- function() {
  list(
    OK = list(f = log(zinc) ~ 1, beta = NULL),
    UK = list(f = log(zinc) ~ sqrt(dist), beta = NULL),
    SK = list(f = log(zinc) ~ sqrt(dist), beta = c(6.6, -2.4))
  )
}

test_that("buffered LOO matches gstat::krige on the remaining data (meuse)", {
  meuse <- meuse_df()
  xy <- as.matrix(meuse[, c("x", "y")])
  targets <- seq(1, nrow(meuse), by = 6)
  for (radius in c(0, 150, 400)) {
    sets <- buffer_sets_r(xy, radius)
    for (cs in buffer_cases()) {
      out <- kcv(cs$f, ~ x + y, meuse, model = m_buf, beta = cs$beta,
        buffer = radius
      )
      ref <- krige_buffer(cs$f, meuse, m_buf, sets, targets, cs$beta)
      expect_equal(out$var1.pred[targets], ref$pred, tolerance = tol)
      expect_equal(out$var1.var[targets], ref$var, tolerance = tol)
      expect_identical(out$n_heldout, lengths(sets))
    }
  }
  # The radii give hold-out sets of several sizes.
  expect_gt(max(lengths(buffer_sets_r(xy, 400))), 20)
})

test_that("buffered LOO matches gstat::krige on a simulated field", {
  v <- gstat::vgm(1, "Mat", 150, 0.2, kappa = 1.5)
  d <- simulate_field(120, v, seed = 31)
  xy <- as.matrix(d[, c("x", "y")])
  targets <- seq(1, 120, by = 5)
  cases <- list(
    list(f = z ~ 1, beta = NULL),
    list(f = z ~ x + y + cov1, beta = NULL),
    list(f = z ~ 1, beta = 2)
  )
  for (radius in c(60, 200)) {
    sets <- buffer_sets_r(xy, radius)
    for (cs in cases) {
      out <- kcv(cs$f, ~ x + y, d, model = v, beta = cs$beta, buffer = radius)
      ref <- krige_buffer(cs$f, d, v, sets, targets, cs$beta)
      expect_equal(out$var1.pred[targets], ref$pred, tolerance = tol)
      expect_equal(out$var1.var[targets], ref$var, tolerance = tol)
    }
  }
})

test_that("buffer = 0 is exactly leave-one-out", {
  meuse <- meuse_df()
  cols <- c(cv_cols, "fold", "x", "y")
  for (cs in buffer_cases()) {
    loo <- kcv(cs$f, ~ x + y, meuse, model = m_buf, beta = cs$beta)
    b0 <- kcv(cs$f, ~ x + y, meuse, model = m_buf, beta = cs$beta, buffer = 0)
    expect_identical(b0[cols], loo[cols])
    expect_identical(b0$n_heldout, rep(1L, nrow(meuse)))
    p <- kcv_prepare(cs$f, ~ x + y, meuse, model = m_buf, beta = cs$beta)
    expect_identical(kcv_run(p, buffer = 0)[cols], kcv_run(p)[cols])
  }
  # Also for sf output.
  a <- kcv(log(zinc) ~ 1, meuse_sf(), model = m_buf, buffer = 0)
  b <- kcv(log(zinc) ~ 1, meuse_sf(), model = m_buf)
  expect_identical(a$var1.pred, b$var1.pred)
  expect_identical(a$var1.var, b$var1.var)
})

test_that("the C++ buffer engine matches both R oracles", {
  set.seed(5)
  n <- 60
  xy <- matrix(runif(2 * n, 0, 10), ncol = 2)
  C <- 0.8 * exp(-as.matrix(dist(xy)) / 3) + diag(0.2, n)
  Y <- matrix(rnorm(3 * n), n, 3)
  Fs <- list(
    OK = matrix(1, n, 1),
    UK = cbind(1, xy[, 1], rnorm(n))
  )
  for (radius in c(0, 0.8, 2, 3.5)) {
    sets <- buffer_sets_r(xy, radius)
    for (nm in c("OK", "UK", "SK")) {
      F <- if (nm == "SK") Fs$UK else Fs[[nm]]
      beta <- if (nm == "SK") c(0.3, -0.1, 0.2) else NULL
      p <- kcv_engine_prepare(C, F, beta)
      p$inp <- list(coords = xy)
      res <- run_buffer_prepared(p, Y, radius)
      r1 <- kcv_buffer_r(C, F, Y, sets, beta)
      r2 <- kcv_buffer_reference(C, F, Y, sets, beta)
      for (r in list(r1, r2)) {
        expect_equal(unname(res$pred), unname(r$pred), tolerance = tol)
        expect_equal(res$var, r$var, tolerance = tol)
      }
      expect_identical(res$n_heldout, lengths(sets))
      expect_identical(res$folds, seq_len(n))
      # Each column of a matrix-valued response equals a separate run.
      expect_equal(
        res$pred[, 2],
        drop(run_buffer_prepared(p, Y[, 2, drop = FALSE], radius)$pred),
        tolerance = tol
      )
    }
  }
})

test_that("one prepared object serves several radii and input types", {
  meuse <- meuse_df()
  p <- kcv_prepare(log(zinc) ~ sqrt(dist), ~ x + y, meuse, model = m_buf)
  for (radius in c(100, 300)) {
    a <- kcv_run(p, buffer = radius)
    b <- kcv(log(zinc) ~ sqrt(dist), ~ x + y, meuse, model = m_buf,
      buffer = radius
    )
    expect_identical(a, b)
    s_f <- kcv(log(zinc) ~ sqrt(dist), meuse_sf(), model = m_buf,
      buffer = radius
    )
    s_p <- kcv(log(zinc) ~ sqrt(dist), meuse_sp(), model = m_buf,
      buffer = radius
    )
    expect_s3_class(s_f, "sf")
    expect_s4_class(s_p, "SpatialPointsDataFrame")
    expect_cv_equal(s_f, a)
    expect_cv_equal(s_p, a)
    expect_identical(s_f$n_heldout, a$n_heldout)
    expect_identical(s_p$n_heldout, a$n_heldout)
  }
  # Larger buffers hold out more points.
  expect_true(all(kcv_run(p, buffer = 300)$n_heldout >=
    kcv_run(p, buffer = 100)$n_heldout))
})

test_that("buffer arguments are validated", {
  meuse <- meuse_df()
  for (bad in list(-1, NA_real_, c(1, 2), "a", Inf)) {
    expect_error(
      kcv(log(zinc) ~ 1, ~ x + y, meuse, model = m_buf, buffer = bad),
      "`buffer` must be a single non-negative number"
    )
  }
  expect_error(
    kcv(log(zinc) ~ 1, ~ x + y, meuse, model = m_buf, buffer = 100,
      folds = rep(1:5, 31)
    ),
    "cannot be combined"
  )
  expect_error(
    kcv(log(zinc) ~ 1, ~ x + y, meuse, model = m_buf, buffer = 100, nfold = 5),
    "cannot be combined"
  )
  p <- kcv_prepare(log(zinc) ~ 1, ~ x + y, meuse, model = m_buf)
  expect_error(kcv_run(p, nfold = 5, buffer = 100), "cannot be combined")
  expect_error(kcv_run(p, buffer = -2), "non-negative")
  pe <- kcv_engine_prepare(diag(3), matrix(1, 3, 1))
  expect_error(kcv_run(pe, buffer = 1), "kcv_engine_run")
})

test_that("a buffer holding out every observation names the point", {
  meuse <- meuse_df()
  msg <- paste0(
    "Point 1 \\(the buffer of radius 1e\\+05 holds out 155 of 155 ",
    "observations\\): every observation is held out"
  )
  expect_error(
    kcv(log(zinc) ~ 1, ~ x + y, meuse, model = m_buf, buffer = 1e5),
    msg
  )
  xy <- as.matrix(meuse[, c("x", "y")])
  C <- covariance_matrix(m_buf, xy)
  expect_error(
    kcv_buffer_r(C, matrix(1, 155, 1), meuse$zinc, buffer_sets_r(xy, 1e5),
      radius = 1e5
    ),
    msg
  )
  # SK too: no observation is left.
  expect_error(
    kcv(log(zinc) ~ 1, ~ x + y, meuse, model = m_buf, beta = 6, buffer = 1e5),
    msg
  )
})

test_that("a buffer that leaves the trend unidentifiable names the point", {
  meuse <- meuse_df()
  xy <- as.matrix(meuse[, c("x", "y")])
  d1 <- sqrt(colSums((t(xy) - xy[1, ])^2))
  # An indicator covariate for the three points closest to point 1.
  meuse$near <- as.numeric(rank(d1) <= 3)
  r <- max(d1[meuse$near == 1])
  m <- sum(d1 <= r)
  # Plain LOO is fine; a buffer around point 1 removes all three points.
  expect_no_error(kcv(log(zinc) ~ near, ~ x + y, meuse, model = m_buf))
  msg <- paste0(
    "Point 1 \\(the buffer of radius ", format(r), " holds out ", m,
    " of 155 observations\\): the trend is not identifiable from the data ",
    "outside the buffer.*smaller radius"
  )
  expect_error(
    kcv(log(zinc) ~ near, ~ x + y, meuse, model = m_buf, buffer = r),
    msg
  )
  C <- covariance_matrix(m_buf, xy)
  F <- cbind(1, meuse$near)
  for (oracle in list(kcv_buffer_r, kcv_buffer_reference)) {
    expect_error(
      oracle(C, F, log(meuse$zinc), buffer_sets_r(xy, r), radius = r),
      msg
    )
  }
  # With a known trend (SK) nothing needs identifying.
  expect_no_error(kcv(log(zinc) ~ near, ~ x + y, meuse,
    model = m_buf,
    beta = c(6, 0.5), buffer = r
  ))
})

test_that("a numerically singular Q[I, I] names the point", {
  expect_null(buffer_solve_cpp(matrix(c(1, 1, 1, 1), 2), 0L))
  Qii <- matrix(c(2, 1, 1, 2), 2)
  expect_equal(drop(buffer_solve_cpp(Qii, 1L)), solve(Qii)[, 2], tolerance = tol)
  expect_error(
    raise_buffer_status(5, 4, 2, 100, 10),
    paste0(
      "Point 4 \\(the buffer of radius 100 holds out 2 of 10 ",
      "observations\\): the matrix Q\\[I, I\\] is numerically singular, so ",
      "this point cannot be predicted"
    )
  )
})

test_that("the pair-cache and per-point strategies agree with the oracles", {
  set.seed(8)
  n <- 80
  xy <- matrix(runif(2 * n, 0, 10), ncol = 2)
  C <- 0.7 * exp(-as.matrix(dist(xy)) / 2.5) + diag(0.3, n)
  Y <- matrix(rnorm(2 * n), n, 2)
  F <- cbind(1, xy[, 2], rnorm(n))
  for (radius in c(0.5, 1.5, 3)) {
    sets <- buffer_sets_r(xy, radius)
    for (beta in list(NULL, c(0.1, 0.05, -0.2))) {
      for (FF in list(F[, 1, drop = FALSE], F)) {
        if (!is.null(beta) && ncol(FF) != length(beta)) next
        p <- kcv_engine_prepare(C, FF, beta)
        p$inp <- list(coords = xy)
        pr <- run_buffer_prepared(p, Y, radius, method = "pairs")
        po <- run_buffer_prepared(p, Y, radius, method = "points")
        expect_identical(pr$method, "pairs")
        expect_identical(po$method, "points")
        ref <- kcv_buffer_reference(C, FF, Y, sets, beta)
        for (r in list(pr, po)) {
          expect_equal(unname(r$pred), unname(ref$pred), tolerance = tol)
          expect_equal(r$var, ref$var, tolerance = tol)
          expect_identical(r$n_heldout, lengths(sets))
        }
      }
    }
  }
})

test_that("auto uses the pair cache only when it fits the memory budget", {
  meuse <- meuse_df()
  p <- kcv_prepare(log(zinc) ~ 1, ~ x + y, meuse, model = m_buf)
  y <- p$inp$y
  expect_identical(run_buffer_prepared(p, y, 400)$method, "pairs")
  # A budget that only just holds the factorization leaves no room.
  withr::local_options(fastkcv.max_gb = 8 * 155^2 / 1e9 + 1e-9)
  auto <- run_buffer_prepared(p, y, 400)
  expect_identical(auto$method, "points")
  # Forcing the cache ignores the budget (internal, for testing).
  forced <- run_buffer_prepared(p, y, 400, method = "pairs")
  expect_identical(forced$method, "pairs")
  expect_equal(forced$pred, auto$pred, tolerance = tol)
  expect_equal(forced$var, auto$var, tolerance = tol)
})

test_that("targets returns the target rows of the full run", {
  meuse <- meuse_df()
  tg <- c(3L, 40L, 41L, 150L)
  full <- kcv(log(zinc) ~ sqrt(dist), ~ x + y, meuse, model = m_buf,
    buffer = 300
  )
  a <- kcv(log(zinc) ~ sqrt(dist), ~ x + y, meuse, model = m_buf,
    buffer = 300, targets = rev(tg)
  )
  expect_identical(nrow(a), length(tg))
  expect_identical(a$fold, tg)
  expect_equal(a, full[tg, ], tolerance = tol, ignore_attr = TRUE)
  lg <- seq_len(155) %in% tg
  p <- kcv_prepare(log(zinc) ~ sqrt(dist), ~ x + y, meuse, model = m_buf)
  expect_equal(kcv_run(p, buffer = 300, targets = lg), a, tolerance = tol)
  # Spatial outputs keep the target geometries.
  s_f <- kcv(log(zinc) ~ sqrt(dist), meuse_sf(), model = m_buf, buffer = 300,
    targets = tg
  )
  s_p <- kcv(log(zinc) ~ sqrt(dist), meuse_sp(), model = m_buf, buffer = 300,
    targets = tg
  )
  expect_cv_equal(s_f, a)
  expect_cv_equal(s_p, a)
  expect_equal(unname(sf::st_coordinates(s_f)),
    unname(as.matrix(meuse[tg, c("x", "y")])))
  expect_equal(unname(sp::coordinates(s_p)),
    unname(as.matrix(meuse[tg, c("x", "y")])))
  expect_identical(kcv_metrics(a)$n, 4L)
  # buffer = 0 with targets is plain LOO at those rows.
  loo <- kcv(log(zinc) ~ 1, ~ x + y, meuse, model = m_buf)
  b0 <- kcv(log(zinc) ~ 1, ~ x + y, meuse, model = m_buf, buffer = 0,
    targets = tg
  )
  expect_identical(b0$var1.pred, loo$var1.pred[tg])
})

test_that("targets errors name the original row and invalid targets fail", {
  meuse <- meuse_df()
  xy <- as.matrix(meuse[, c("x", "y")])
  d1 <- sqrt(colSums((t(xy) - xy[1, ])^2))
  meuse$near <- as.numeric(rank(d1) <= 3)
  r <- max(d1[meuse$near == 1])
  # Point 1's buffer removes the covariate; scoring rows 2 and 1 reports 1,
  # and scoring only rows far away succeeds.
  expect_error(
    kcv(log(zinc) ~ near, ~ x + y, meuse, model = m_buf, buffer = r,
      targets = c(150, 1)
    ),
    "Point 1 \\("
  )
  far <- which(d1 > 3 * r)[1:5]
  expect_no_error(kcv(log(zinc) ~ near, ~ x + y, meuse, model = m_buf,
    buffer = r, targets = far
  ))
  m <- function(...) kcv(log(zinc) ~ 1, ~ x + y, meuse, model = m_buf, ...)
  expect_error(m(targets = 1:3), "only be used with `buffer`")
  expect_error(m(buffer = 10, targets = c(1, 1)), "duplicate")
  expect_error(m(buffer = 10, targets = c(0, 2)), "row numbers between 1")
  expect_error(m(buffer = 10, targets = 156), "row numbers between 1")
  expect_error(m(buffer = 10, targets = c(1, NA)), "row numbers between 1")
  expect_error(m(buffer = 10, targets = 1.5), "row numbers between 1")
  expect_error(m(buffer = 10, targets = "a"), "row numbers between 1")
  expect_error(m(buffer = 10, targets = c(TRUE, FALSE)), "length n = 155")
  expect_error(m(buffer = 10, targets = rep(FALSE, 155)), "no observations")
  expect_error(m(buffer = 10, targets = integer(0)), "no observations")
})

test_that("a point at distance about r is held out from both sides", {
  # x[i] + r rounds below x[j] although the distance test accepts j, so the
  # x prefilter used to drop j from i's set (but not i from j's).
  xi <- -12.917385902255774
  r <- 18.62176014110446
  xj <- 5.7043742388486871
  expect_lte(sqrt((xj - xi)^2), r)
  expect_lt(xi + r, xj)
  d <- data.frame(
    x = c(xi, xj, 200, 400, 600, 300), y = c(0, 0, 100, 50, 300, 500),
    z = c(1, 2, 3, 2, 1, 4)
  )
  v <- gstat::vgm(1, "Exp", 100, 0.1)
  out <- kcv(z ~ 1, ~ x + y, d, model = v, buffer = r)
  sets <- buffer_sets_r(as.matrix(d[, c("x", "y")]), r)
  expect_identical(out$n_heldout, lengths(sets))
  expect_identical(out$n_heldout[1:2], c(2L, 2L))
  ref <- krige_buffer(z ~ 1, d, v, sets, 1:6)
  expect_equal(out$var1.pred, ref$pred, tolerance = tol)
  expect_equal(out$var1.var, ref$var, tolerance = tol)
})

test_that("buffered LOO with a nested anisotropic model matches gstat::krige", {
  # The buffer is isotropic in the coordinates while the covariance is not.
  meuse <- meuse_df()
  v <- gstat::vgm(0.4, "Sph", 900, 0.1, anis = c(40, 0.4),
    add.to = gstat::vgm(0.2, "Exp", 150, anis = c(130, 0.6)))
  targets <- seq(2, nrow(meuse), by = 7)
  sets <- buffer_sets_r(as.matrix(meuse[, c("x", "y")]), 250)
  for (cs in buffer_cases()) {
    out <- kcv(cs$f, ~ x + y, meuse, model = v, beta = cs$beta,
      buffer = 250, targets = targets
    )
    ref <- krige_buffer(cs$f, meuse, v, sets, targets, cs$beta)
    expect_equal(out$var1.pred, ref$pred, tolerance = tol)
    expect_equal(out$var1.var, ref$var, tolerance = tol)
  }
})
