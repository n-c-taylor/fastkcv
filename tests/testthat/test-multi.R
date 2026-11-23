# Several responses cross-validated with one factorization: each response's
# block of the output must equal the single-response kcv() call.
#
# Predictions are compared at 1e-12, not with expect_identical(): with T
# columns the engine forms Q y with one matrix-matrix product (dgemm), with
# one column a matrix-vector product (dgemv), and an optimized BLAS such as
# OpenBLAS sums in a different order in the two kernels, so the last bits
# can differ (about 1e-15 here). Variances, folds and hold-out sizes do not
# depend on y and are compared with expect_identical().

m_multi <- gstat::vgm(0.5, "Sph", 900, 0.1)

multi_formulas <- function() {
  list(
    OK = list(
      multi = cbind(log(zinc), log(copper), lz = log(lead)) ~ 1,
      single = list(log(zinc) ~ 1, log(copper) ~ 1, log(lead) ~ 1)
    ),
    UK = list(
      multi = cbind(log(zinc), log(copper), lz = log(lead)) ~ sqrt(dist),
      single = list(
        log(zinc) ~ sqrt(dist), log(copper) ~ sqrt(dist),
        log(lead) ~ sqrt(dist)
      )
    )
  )
}
multi_names <- c("log(zinc)", "log(copper)", "lz")

# Block of response `r` as a plain data.frame without the response column.
response_block <- function(cv, r) {
  cv <- as.data.frame(cv)
  b <- cv[cv$response == r, setdiff(names(cv), c("response", "geometry"))]
  rownames(b) <- NULL
  b
}

expect_block_matches <- function(multi, single, r) {
  b <- response_block(multi, r)
  s <- as.data.frame(single)
  s <- s[, setdiff(names(s), "geometry")]
  expect_identical(names(b), names(s))
  for (col in setdiff(names(s), c("var1.pred", "residual", "zscore"))) {
    expect_identical(b[[col]], s[[col]], label = paste(r, col))
  }
  for (col in c("var1.pred", "residual", "zscore")) {
    expect_equal(b[[col]], s[[col]], tolerance = 1e-12, label = paste(r, col))
  }
}

expect_all_blocks <- function(multi, singles) {
  expect_identical(levels(multi$response), multi_names)
  for (j in seq_along(singles)) {
    expect_block_matches(multi, singles[[j]], multi_names[j])
  }
}

test_that("LOO of several responses equals single-response kcv()", {
  meuse <- meuse_df()
  for (fm in multi_formulas()) {
    multi <- kcv(fm$multi, ~ x + y, meuse, model = m_multi)
    expect_s3_class(multi$response, "factor")
    expect_identical(nrow(multi), 3L * nrow(meuse))
    expect_identical(
      names(multi),
      c(cv_cols, "fold", "response", "x", "y")
    )
    singles <- lapply(fm$single, kcv, ~ x + y, meuse, model = m_multi)
    expect_all_blocks(multi, singles)
  }
})

test_that("k-fold of several responses equals single-response kcv()", {
  meuse <- meuse_df()
  fv <- rep(1:7, length.out = nrow(meuse))
  for (fm in multi_formulas()) {
    multi <- withr::with_seed(3, kcv(fm$multi, ~ x + y, meuse,
      model = m_multi, nfold = 5
    ))
    singles <- lapply(fm$single, function(f) {
      withr::with_seed(3, kcv(f, ~ x + y, meuse, model = m_multi, nfold = 5))
    })
    expect_all_blocks(multi, singles)

    multi <- kcv(fm$multi, ~ x + y, meuse, model = m_multi, folds = fv)
    singles <- lapply(fm$single, kcv, ~ x + y, meuse,
      model = m_multi,
      folds = fv
    )
    expect_all_blocks(multi, singles)
  }
})

test_that("buffered LOO and targets of several responses equal single runs", {
  meuse <- meuse_df()
  targets <- c(150, 3, 77, 20)
  for (fm in multi_formulas()) {
    for (tg in list(NULL, targets)) {
      multi <- kcv(fm$multi, ~ x + y, meuse,
        model = m_multi, buffer = 300,
        targets = tg
      )
      expect_identical(
        names(multi),
        c(cv_cols, "fold", "n_heldout", "response", "x", "y")
      )
      singles <- lapply(fm$single, kcv, ~ x + y, meuse,
        model = m_multi,
        buffer = 300, targets = tg
      )
      expect_all_blocks(multi, singles)
    }
  }
  # Both strategies for Q[I, I] give the same blocks.
  p <- kcv_prepare(log(zinc) ~ 1, ~ x + y, meuse, model = m_multi)
  Y <- cbind(log(meuse$zinc), log(meuse$copper))
  for (method in c("pairs", "points")) {
    all <- run_buffer_prepared(p, Y, 300, method = method)
    one <- run_buffer_prepared(p, Y[, 2, drop = FALSE], 300, method = method)
    expect_identical(all$var, one$var)
    expect_equal(all$pred[, 2], one$pred[, 1], tolerance = 1e-12)
  }
})

test_that("simple kriging with a shared or per-response beta", {
  meuse <- meuse_df()
  f <- cbind(log(zinc), log(copper), lz = log(lead)) ~ sqrt(dist)
  singles_f <- multi_formulas()$UK$single
  B <- cbind(c(6.6, -2.4), c(3.8, -1.5), c(5.1, -1.9))
  for (scheme in list(list(), list(folds = rep(1:4, length.out = 155)),
    list(buffer = 250))) {
    shared <- do.call(kcv, c(
      list(f, ~ x + y, meuse, model = m_multi, beta = c(6, -2)),
      scheme
    ))
    singles <- lapply(singles_f, function(g) {
      do.call(kcv, c(
        list(g, ~ x + y, meuse, model = m_multi, beta = c(6, -2)),
        scheme
      ))
    })
    expect_all_blocks(shared, singles)

    per <- do.call(kcv, c(
      list(f, ~ x + y, meuse, model = m_multi, beta = B),
      scheme
    ))
    singles <- lapply(seq_along(singles_f), function(j) {
      do.call(kcv, c(
        list(singles_f[[j]], ~ x + y, meuse, model = m_multi, beta = B[, j]),
        scheme
      ))
    })
    expect_all_blocks(per, singles)
  }
  # A p x 1 matrix is the same as the shared vector.
  expect_identical(
    kcv(f, ~ x + y, meuse, model = m_multi, beta = matrix(c(6, -2))),
    kcv(f, ~ x + y, meuse, model = m_multi, beta = c(6, -2))
  )
})

test_that("kcv_run(y = ) reuses the factorization for new responses", {
  meuse <- meuse_df()
  fv <- rep(1:5, length.out = nrow(meuse))
  Y <- cbind(log(meuse$zinc), log(meuse$copper), log(meuse$lead))
  colnames(Y) <- multi_names
  for (fm in multi_formulas()) {
    rhs <- fm$single[[1]]
    p <- kcv_prepare(rhs, ~ x + y, meuse, model = m_multi)
    pm <- kcv_prepare(fm$multi, ~ x + y, meuse, model = m_multi)
    # The same computation as from a matrix left-hand side.
    expect_identical(kcv_run(p, y = Y), kcv(fm$multi, ~ x + y, meuse,
      model = m_multi
    ))
    expect_identical(kcv_run(p, folds = fv, y = Y), kcv_run(pm, folds = fv))
    expect_identical(
      kcv_run(p, buffer = 200, targets = 1:9, y = Y),
      kcv_run(pm, buffer = 200, targets = 1:9)
    )
    # A vector y gives the single-response output.
    expect_identical(
      kcv_run(p, folds = fv, y = log(meuse$copper)),
      kcv(fm$single[[2]], ~ x + y, meuse, model = m_multi, folds = fv)
    )
    # A prepared object for several responses runs one new response.
    expect_identical(
      kcv_run(pm, y = log(meuse$lead)),
      kcv(fm$single[[3]], ~ x + y, meuse, model = m_multi)
    )
  }
  # Unnamed columns are numbered; a one-column matrix keeps the response
  # column.
  p <- kcv_prepare(log(zinc) ~ 1, ~ x + y, meuse, model = m_multi)
  out <- kcv_run(p, y = unname(Y))
  expect_identical(levels(out$response), c("1", "2", "3"))
  out1 <- kcv_run(p, y = Y[, 2, drop = FALSE])
  expect_identical(levels(out1$response), "log(copper)")
  expect_block_matches(out1, kcv_run(p, y = Y[, 2]), "log(copper)")
  expect_false("response" %in% names(kcv_run(p)))
})

test_that("kcv_run(beta = ) replaces the simple-kriging coefficients", {
  meuse <- meuse_df()
  p <- kcv_prepare(log(zinc) ~ sqrt(dist), ~ x + y, meuse,
    model = m_multi,
    beta = c(6, -2)
  )
  B <- cbind(c(6.6, -2.4), c(3.8, -1.5))
  Y <- cbind(a = log(meuse$zinc), b = log(meuse$copper))
  out <- kcv_run(p, buffer = 150, y = Y, beta = B)
  for (j in 1:2) {
    single <- kcv_run(p, buffer = 150, y = Y[, j], beta = B[, j])
    expect_block_matches(out, single, colnames(Y)[j])
    ref <- kcv(
      stats::as.formula(paste(
        c("log(zinc)", "log(copper)")[j],
        "~ sqrt(dist)"
      )), ~ x + y, meuse,
      model = m_multi, beta = B[, j], buffer = 150
    )
    expect_identical(single, ref)
  }
  # A prepared per-response beta needs as many responses.
  pB <- kcv_prepare(cbind(log(zinc), log(copper)) ~ sqrt(dist), ~ x + y,
    meuse,
    model = m_multi, beta = B
  )
  expect_identical(
    kcv_run(pB, y = Y),
    kcv_run(p, y = Y, beta = B)
  )
  expect_error(kcv_run(pB, y = Y[, 1]), "has 2 columns \\(one per response\\)")
  expect_no_error(kcv_run(pB, y = Y[, 1], beta = B[, 1]))
})

test_that("sf and sp output of several responses", {
  f <- cbind(log(zinc), log(copper)) ~ 1
  m_sf <- meuse_sf()
  m_sp <- meuse_sp()
  out_sf <- kcv(f, m_sf, model = m_multi, buffer = 200, targets = 5:9)
  out_sp <- kcv(f, m_sp, model = m_multi, buffer = 200, targets = 5:9)
  expect_s3_class(out_sf, "sf")
  expect_s4_class(out_sp, "SpatialPointsDataFrame")
  expect_identical(nrow(out_sf), 10L)
  expect_identical(
    unname(sf::st_coordinates(out_sf)[6:10, ]),
    unname(sf::st_coordinates(m_sf)[5:9, ])
  )
  expect_identical(unname(sp::coordinates(out_sp)[6:10, ]),
    unname(sp::coordinates(m_sp)[5:9, ]))
  for (r in c("log(zinc)", "log(copper)")) {
    single_sf <- kcv(
      stats::as.formula(paste(r, "~ 1")), m_sf,
      model = m_multi,
      buffer = 200, targets = 5:9
    )
    expect_block_matches(out_sf, single_sf, r)
    expect_block_matches(out_sp@data, single_sf, r)
    b <- out_sf[out_sf$response == r, ]
    expect_identical(sf::st_geometry(b), sf::st_geometry(single_sf),
      ignore_attr = TRUE)
  }
})

test_that("kcv_metrics() summarizes each response", {
  meuse <- meuse_df()
  fm <- multi_formulas()$OK
  fv <- rep(1:4, length.out = nrow(meuse))
  multi <- kcv(fm$multi, ~ x + y, meuse, model = m_multi, folds = fv)
  singles <- lapply(fm$single, kcv, ~ x + y, meuse,
    model = m_multi,
    folds = fv
  )
  met <- kcv_metrics(multi)
  expect_identical(met$response, factor(multi_names, levels = multi_names))
  ref <- do.call(rbind, lapply(singles, kcv_metrics))
  rownames(ref) <- NULL
  expect_equal(met[, -1], ref, tolerance = 1e-12)

  met_f <- kcv_metrics(multi, by_fold = TRUE)
  ref_f <- do.call(rbind, lapply(singles, kcv_metrics, by_fold = TRUE))
  rownames(ref_f) <- NULL
  expect_identical(nrow(met_f), 12L)
  expect_identical(as.character(met_f$response), rep(multi_names, each = 4))
  expect_equal(met_f[, -1], ref_f, tolerance = 1e-12)

  pooled <- kcv_metrics(multi, by_response = FALSE)
  expect_identical(pooled$n, 3L * nrow(meuse))
  expect_equal(pooled$rmse, sqrt(mean(multi$residual^2)))
  expect_error(kcv_metrics(multi, by_response = NA), "`by_response`")
  # sf output works the same way.
  out_sf <- kcv(cbind(log(zinc), log(copper)) ~ 1, meuse_sf(),
    model = m_multi
  )
  expect_identical(nrow(kcv_metrics(out_sf)), 2L)
})

test_that("a matrix column as the response", {
  meuse <- meuse_df()
  S <- withr::with_seed(9, matrix(stats::rnorm(155 * 4), 155, 4))
  meuse$sims <- S
  out <- kcv(sims ~ 1, ~ x + y, meuse, model = m_multi, nfold = meuse$ffreq)
  expect_identical(levels(out$response), c("1", "2", "3", "4"))
  p <- kcv_prepare(zinc ~ 1, ~ x + y, meuse, model = m_multi)
  expect_identical(out, kcv_run(p, nfold = meuse$ffreq, y = S))
  colnames(S) <- paste0("s", 1:4)
  meuse$sims <- S
  out <- kcv(sims ~ 1, ~ x + y, meuse, model = m_multi)
  expect_identical(levels(out$response), paste0("s", 1:4))
  # One column: a single response (model.response() drops it to a vector,
  # as for scale(z) ~ 1), unless the left-hand side is an explicit cbind().
  meuse$one <- S[, 2, drop = FALSE]
  meuse$s2 <- S[, 2]
  single <- kcv(s2 ~ 1, ~ x + y, meuse, model = m_multi)
  expect_identical(kcv(one ~ 1, ~ x + y, meuse, model = m_multi), single)
  sc <- kcv(scale(zinc) ~ 1, ~ x + y, meuse, model = m_multi)
  expect_false("response" %in% names(sc))
  expect_equal(sc$observed, as.vector(scale(meuse$zinc)))
  cb <- kcv(cbind(s2) ~ 1, ~ x + y, meuse, model = m_multi)
  expect_identical(levels(cb$response), "s2")
  expect_block_matches(cb, single, "s2")
  cb <- kcv(cbind(log(zinc)) ~ 1, ~ x + y, meuse, model = m_multi)
  expect_identical(levels(cb$response), "log(zinc)")
})

test_that("several responses match gstat::krige.cv() (LOO, k-fold, SK)", {
  skip_on_cran()
  meuse <- meuse_df()
  fv <- rep(1:5, length.out = nrow(meuse))
  B <- cbind(c(6.6, -2.4), c(3.8, -1.5), c(5.1, -1.9))
  for (fm in multi_formulas()) {
    loo <- kcv(fm$multi, ~ x + y, meuse, model = m_multi)
    k5 <- kcv(fm$multi, ~ x + y, meuse, model = m_multi, nfold = fv)
    for (j in seq_along(fm$single)) {
      ref <- krige_cv(fm$single[[j]], ~ x + y, meuse, model = m_multi)
      expect_cv_equal(response_block(loo, multi_names[j]), ref)
      ref <- krige_cv(fm$single[[j]], ~ x + y, meuse,
        model = m_multi,
        nfold = fv
      )
      expect_cv_equal(response_block(k5, multi_names[j]), ref)
    }
  }
  fm <- multi_formulas()$UK
  sk <- kcv(fm$multi, ~ x + y, meuse, model = m_multi, beta = B, nfold = fv)
  for (j in seq_along(fm$single)) {
    ref <- krige_cv(fm$single[[j]], ~ x + y, meuse,
      model = m_multi,
      beta = B[, j], nfold = fv
    )
    expect_cv_equal(response_block(sk, multi_names[j]), ref)
  }
})

test_that("several responses match gstat::krige for buffered LOO", {
  skip_on_cran()
  meuse <- meuse_df()
  xy <- as.matrix(meuse[, c("x", "y")])
  sets <- buffer_sets_r(xy, 300)
  targets <- seq(2, nrow(meuse), by = 9)
  B <- cbind(c(6.6, -2.4), c(3.8, -1.5), c(5.1, -1.9))
  cases <- list(
    OK = list(fm = multi_formulas()$OK, beta = NULL),
    UK = list(fm = multi_formulas()$UK, beta = NULL),
    SK = list(fm = multi_formulas()$UK, beta = B)
  )
  for (cs in cases) {
    all <- kcv(cs$fm$multi, ~ x + y, meuse, model = m_multi, beta = cs$beta,
      buffer = 300)
    sub <- kcv(cs$fm$multi, ~ x + y, meuse, model = m_multi, beta = cs$beta,
      buffer = 300, targets = targets)
    for (j in seq_along(cs$fm$single)) {
      ref <- krige_buffer(cs$fm$single[[j]], meuse, m_multi, sets, targets,
        beta = if (!is.null(cs$beta)) cs$beta[, j]
      )
      b <- response_block(all, multi_names[j])
      expect_equal(b$var1.pred[targets], ref$pred, tolerance = tol)
      expect_close(b$var1.pred[targets], ref$pred)
      expect_equal(b$var1.var[targets], ref$var, tolerance = tol)
      expect_close(b$var1.var[targets], ref$var)
      b <- response_block(sub, multi_names[j])
      expect_equal(b$var1.pred, ref$pred, tolerance = tol)
      expect_close(b$var1.pred, ref$pred)
      expect_equal(b$var1.var, ref$var, tolerance = tol)
      expect_close(b$var1.var, ref$var)
      expect_identical(b$n_heldout, lengths(sets)[targets])
    }
  }
})

test_that("several responses: kcv_engine_run() with a p x T beta", {
  s <- withr::with_seed(4, {
    xy <- matrix(stats::runif(80, 0, 100), ncol = 2)
    list(
      C = exp(-as.matrix(stats::dist(xy)) / 30) + diag(0.1, 40),
      F = cbind(1, xy[, 1]), Y = matrix(stats::rnorm(120), 40, 3)
    )
  })
  B <- rbind(c(1, 2, 3), c(0.01, -0.02, 0))
  folds <- rep(1:4, 10)
  p <- kcv_engine_prepare(s$C, s$F, beta = c(0, 0))
  out <- kcv_engine_run(p, s$Y, folds, beta = B)
  ref <- kcv_engine_r(s$C, s$F, s$Y, folds, beta = B)
  expect_equal(out$pred, ref$pred, tolerance = tol)
  expect_close(out$pred, ref$pred)
  ref <- kcv_reference(s$C, s$F, s$Y, folds, beta = B)
  expect_equal(out$pred, ref$pred, tolerance = tol)
  expect_close(out$pred, ref$pred)
  expect_identical(kcv_engine(s$C, s$F, s$Y, folds, beta = B), out)
  for (j in 1:3) {
    one <- kcv_engine_run(p, s$Y[, j], folds, beta = B[, j])
    expect_equal(out$pred[, j], drop(one$pred), tolerance = 1e-12)
  }
  expect_error(kcv_engine_run(kcv_engine_prepare(s$C, s$F), s$Y, folds,
    beta = B), "only be given for a prepared object made for simple")
  expect_error(kcv_engine(s$C, s$F, s$Y, folds, beta = B[, 1:2]),
    "has 2 columns but there are 3 response")
})

test_that("errors for several responses", {
  meuse <- meuse_df()
  meuse$cu <- log(meuse$copper)
  meuse$cu[c(4, 9)] <- NA
  expect_error(
    kcv(cbind(log(zinc), cu) ~ 1, ~ x + y, meuse, model = m_multi),
    "Response `cu` has 2 missing or non-finite value\\(s\\)\\. Every response"
  )
  meuse$cu[c(4, 9)] <- c(Inf, 1)
  meuse$cd <- meuse$cu
  expect_error(
    kcv(cbind(cu, log(zinc), cd) ~ 1, ~ x + y, meuse, model = m_multi),
    "Response `cu` has 1 missing .*\\(and 1 other response\\(s\\)\\)"
  )
  expect_error(
    kcv(cbind(a = zinc, a = lead) ~ 1, ~ x + y, meuse, model = m_multi),
    "duplicate column names \\(`a`\\)"
  )

  p <- kcv_prepare(log(zinc) ~ 1, ~ x + y, meuse, model = m_multi)
  Y <- cbind(log(meuse$zinc), log(meuse$lead))
  expect_error(kcv_run(p, y = Y[-1, ]), "`y` must be a numeric vector of length n = 155")
  expect_error(kcv_run(p, y = 1:3), "`y` must be a numeric vector")
  expect_error(kcv_run(p, y = as.character(Y[, 1])), "`y` must be a numeric")
  expect_error(kcv_run(p, y = array(0, c(155, 2, 2))), "`y` must be a numeric")
  expect_error(kcv_run(p, y = matrix(0, 155, 0)), "`y` has no columns")
  Yb <- Y
  Yb[3, 2] <- NaN
  expect_error(kcv_run(p, y = Yb), "Response `2` has 1 missing or non-finite")
  yb <- Y[, 1]
  yb[7] <- NA
  expect_error(kcv_run(p, y = yb), "`y` has 1 missing value\\(s\\)")
  expect_error(kcv_run(p, y = replace(Y[, 1], 2, Inf)),
    "`y` has 1 infinite value\\(s\\)")
  expect_error(
    kcv_run(p, y = `colnames<-`(Y, c("a", "a"))),
    "`y` has duplicate column names"
  )
  # NA column names used to fail with "missing value where TRUE/FALSE
  # needed"; they are treated like empty names.
  Yna <- `colnames<-`(Y, c("a", NA))
  expect_identical(levels(kcv_run(p, y = Yna)$response), c("1", "2"))
  meuse$Yna <- Yna
  expect_identical(
    levels(kcv(Yna ~ 1, ~ x + y, meuse, model = m_multi)$response),
    c("1", "2")
  )
  expect_identical(
    levels(kcv(cbind(a = log(zinc), log(lead)) ~ 1, ~ x + y, meuse,
      model = m_multi)$response),
    c("a", "log(lead)")
  )
  expect_error(kcv_run(p, y = Y, beta = 6), "only be given for a prepared")

  ps <- kcv_prepare(log(zinc) ~ sqrt(dist), ~ x + y, meuse,
    model = m_multi,
    beta = c(6, -2)
  )
  expect_error(kcv_run(ps, y = Y, beta = matrix(1, 3, 2)),
    "`beta` has 3 row\\(s\\) but the trend has 2 coefficient\\(s\\)")
  expect_error(kcv_run(ps, y = Y, beta = matrix(1, 2, 3)),
    "`beta` has 3 columns but there are 2 response\\(s\\)")
  expect_error(kcv_run(ps, y = Y, beta = matrix(NA_real_, 2, 2)),
    "`beta` must be numeric, with no missing")
  expect_error(kcv(cbind(log(zinc), log(lead)) ~ 1, ~ x + y, meuse,
    model = m_multi, beta = matrix(1, 1, 3)), "`beta` has 3 columns")

  pe <- kcv_engine_prepare(diag(155), matrix(1, 155, 1))
  expect_error(kcv_run(pe, y = Y), "use kcv_engine_run\\(\\) with it")
})

test_that("the memory budget limits the number of responses", {
  meuse <- meuse_df()
  p <- kcv_prepare(log(zinc) ~ 1, ~ x + y, meuse, model = m_multi)
  Y <- matrix(rep(log(meuse$zinc), 10), 155, 10)
  # 8 n^2 = 192,200 bytes for the factor, 8 (4 n + 10 n) = 17,360 bytes per
  # response: a budget of 2.5e-4 GB leaves room for 3 responses.
  withr::with_options(list(fastkcv.max_gb = 2.5e-4), {
    expect_error(kcv_run(p, y = Y),
      "10 responses at n = 155 need about 0.00017 GB.*batches of at most 3")
    expect_no_error(kcv_run(p, y = Y[, 1:3]))
    expect_error(kcv_run(p, buffer = 100, targets = 1:100, y = Y),
      "batches of at most 4 columns")
    # One response is never checked.
    expect_no_error(kcv_run(p, y = Y[, 1]))
    expect_error(
      kcv(cbind(log(zinc), log(lead), log(copper), cadmium) ~ 1, ~ x + y,
        meuse, model = m_multi),
      "4 responses at n = 155"
    )
  })
  withr::with_options(list(fastkcv.max_gb = 2e-4), {
    expect_error(kcv_run(p, y = Y), "budget.*\\. Raise the budget")
  })
  # kcv_engine_run() returns only the n x T predictions, so it counts one
  # output column per response, not ten: 8 (4 n + n) = 6,200 bytes per
  # response at n = 155. A budget of 2.5e-4 GB leaves 57,800 bytes beside
  # the factor: room for 9 responses (it used to refuse more than 3).
  pe <- kcv_engine_prepare(diag(155), matrix(1, 155, 1))
  withr::with_options(list(fastkcv.max_gb = 2.5e-4), {
    expect_no_error(kcv_engine_run(pe, Y[, 1:9], seq_len(155)))
    expect_error(kcv_engine_run(pe, Y, seq_len(155)),
      "10 responses at n = 155.*batches of at most 9")
  })
})
