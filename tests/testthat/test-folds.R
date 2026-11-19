is_partition <- function(f, n) {
  expect_type(f, "integer")
  expect_length(f, n)
  expect_false(anyNA(f))
  # Fold ids are 1..K with no gaps.
  expect_identical(sort(unique(f)), seq_len(max(f)))
}

location_inputs <- function() {
  m <- meuse_df()
  list(
    matrix = as.matrix(m[, c("x", "y")]),
    df = list(~ x + y, m),
    sf = list(meuse_sf(), NULL),
    sp = list(meuse_sp(), NULL),
    sp_points = list(sp::SpatialPoints(meuse_sp()), NULL)
  )
}

call_builder <- function(fun, loc, ...) {
  if (is.list(loc)) fun(loc[[1]], loc[[2]], ...) else fun(loc, ...)
}

test_that("fold builders give one fold per point for every input type", {
  n <- nrow(meuse_df())
  builders <- list(
    random = function(l, d = NULL) fold_random(l, d, k = 7),
    blocks = function(l, d = NULL) fold_blocks(l, d, blocks = c(4, 3)),
    blocks_k = function(l, d = NULL) fold_blocks(l, d, size = 600, k = 3),
    kmeans = function(l, d = NULL) fold_kmeans(l, d, k = 8)
  )
  for (b in builders) {
    res <- lapply(location_inputs(), function(loc) {
      withr::with_seed(42, call_builder(b, loc))
    })
    for (f in res) {
      is_partition(f, n)
      expect_identical(f, res[[1]])
    }
  }
})

test_that("fold builders are reproducible under a seed", {
  m <- meuse_df()
  draws <- list(
    function() fold_random(~ x + y, m, k = 10),
    function() fold_kmeans(~ x + y, m, k = 10),
    function() fold_blocks(~ x + y, m, blocks = 5, k = 4)
  )
  for (d in draws) {
    a <- withr::with_seed(1, d())
    expect_identical(a, withr::with_seed(1, d()))
    expect_false(identical(a, withr::with_seed(2, d())))
  }
})

test_that("fold_random() draws the same folds as krige.cv()", {
  m <- meuse_df()
  for (k in c(2, 5, 10)) {
    ref <- withr::with_seed(7, krige_cv(log(zinc) ~ 1, ~ x + y, m,
      model = gstat::vgm(0.5, "Sph", 900, 0.1), nfold = k
    ))
    expect_identical(
      withr::with_seed(7, fold_random(~ x + y, m, k = k)),
      as.integer(ref$fold)
    )
    # And kcv(nfold = k) uses the same draw.
    expect_identical(
      withr::with_seed(7, fold_random(~ x + y, m, k = k)),
      as.integer(withr::with_seed(7, kcv(log(zinc) ~ 1, ~ x + y, m,
        model = gstat::vgm(0.5, "Sph", 900, 0.1), nfold = k
      ))$fold)
    )
  }
  expect_identical(fold_random(~ x + y, m, k = 155), seq_len(155))
  expect_identical(fold_random(~ x + y, m, k = 1000), seq_len(155))
})

test_that("fold_blocks() gives the requested number of blocks", {
  # A 4 x 3 lattice with two points per cell, away from the cell edges.
  g <- expand.grid(x = seq(0.5, 3.5, by = 1), y = seq(0.5, 2.5, by = 1))
  xy <- as.matrix(rbind(g, g + 0.1))
  f <- fold_blocks(xy, blocks = c(4, 3))
  is_partition(f, nrow(xy))
  expect_identical(max(f), 12L)
  expect_true(all(table(f) == 2))
  # Cells are numbered x fastest, then y, from the minimum corner.
  expect_identical(f[1:12], 1:12)
  expect_identical(f[13:24], 1:12)
  # A scalar is used for both axes; size gives the matching grid.
  expect_identical(max(fold_blocks(xy, blocks = 2)), 4L)
  # The bounding box is [0.5, 3.6] x [0.5, 2.6]: cells of side 0.8 give
  # ceiling(3.1 / 0.8) = 4 by ceiling(2.1 / 0.8) = 3 cells.
  expect_identical(fold_blocks(xy, size = 0.8), f)
  expect_identical(max(fold_blocks(xy, size = c(1.6, 1.5))), 4L)
})

test_that("fold_blocks() drops empty blocks and handles edges", {
  # Points in three corners of a 2 x 2 grid only: the empty cell is dropped
  # and the ids stay consecutive.
  xy <- rbind(c(0, 0), c(0.1, 0.2), c(10, 0), c(0, 10), c(1, 9))
  f <- fold_blocks(xy, blocks = 2)
  expect_identical(f, c(1L, 1L, 2L, 3L, 3L))
  # The point on the maximum edge belongs to the last cell.
  xy <- rbind(c(0, 0), c(5, 0), c(10, 0))
  expect_identical(fold_blocks(xy, blocks = c(2, 1)), c(1L, 2L, 2L))
  expect_identical(fold_blocks(xy, size = 5), c(1L, 2L, 2L))
  # A zero-extent axis gives a single row of cells.
  expect_identical(fold_blocks(xy, blocks = c(1, 5)), c(1L, 1L, 1L))
})

test_that("fold_blocks(k =) groups blocks into k balanced folds", {
  m <- meuse_df()
  blocks <- fold_blocks(~ x + y, m, blocks = 4)
  f <- withr::with_seed(3, fold_blocks(~ x + y, m, blocks = 4, k = 3))
  is_partition(f, nrow(m))
  expect_identical(max(f), 3L)
  # Every block lies within one fold; blocks per fold differ by at most one.
  per_block <- tapply(f, blocks, function(v) length(unique(v)))
  expect_true(all(per_block == 1))
  n_blocks <- table(tapply(f, blocks, `[`, 1))
  expect_lte(max(n_blocks) - min(n_blocks), 1)
})

test_that("fold_kmeans() gives k clusters, as stats::kmeans()", {
  m <- meuse_df()
  for (k in c(1, 5, 12)) {
    f <- withr::with_seed(9, fold_kmeans(~ x + y, m, k = k))
    is_partition(f, nrow(m))
    expect_identical(max(f), as.integer(k))
  }
  ref <- withr::with_seed(9, stats::kmeans(as.matrix(m[, c("x", "y")]), 5,
    nstart = 10, iter.max = 100
  )$cluster)
  expect_identical(
    withr::with_seed(9, fold_kmeans(~ x + y, m, k = 5)),
    as.integer(ref)
  )
})

test_that("kcv() with spatial folds matches krige.cv()", {
  m <- meuse_df()
  v <- gstat::vgm(0.5, "Sph", 900, 0.1)
  for (f in list(
    fold_blocks(~ x + y, m, blocks = 3),
    withr::with_seed(5, fold_kmeans(~ x + y, m, k = 6))
  )) {
    expect_cv_equal(
      kcv(log(zinc) ~ sqrt(dist), ~ x + y, m, model = v, folds = f),
      krige_cv(log(zinc) ~ sqrt(dist), ~ x + y, m, model = v, nfold = f)
    )
  }
})

test_that("fold builders reject invalid arguments", {
  m <- meuse_df()
  xy <- as.matrix(m[, c("x", "y")])
  expect_error(fold_blocks(xy), "exactly one of")
  expect_error(fold_blocks(xy, blocks = 2, size = 100), "exactly one of")
  expect_error(fold_blocks(xy, blocks = 0), "`blocks`")
  expect_error(fold_blocks(xy, blocks = c(1, 2, 3)), "`blocks`")
  expect_error(fold_blocks(xy, blocks = 1.5), "`blocks`")
  expect_error(fold_blocks(xy, size = -1), "`size`")
  expect_error(fold_blocks(xy, size = NA_real_), "`size`")
  expect_error(fold_blocks(xy, blocks = 2, k = 5), "non-empty blocks \\(4\\)")
  expect_error(fold_blocks(xy, blocks = 2, k = 0), "`k`")
  expect_error(fold_random(xy, k = 2.5), "`k`")
  expect_error(fold_random(xy, k = c(2, 3)), "`k`")
  expect_error(fold_random(xy), "`k`.*required")
  expect_error(fold_kmeans(xy), "`k`.*required")
  expect_error(fold_kmeans(xy, k = NA), "`k`")
  expect_error(
    fold_kmeans(rbind(xy[1:3, ], xy[1:3, ]), k = 4),
    "distinct locations \\(3\\)"
  )
  expect_error(fold_random(xy[, 1, drop = FALSE], k = 2), "2D")
  expect_error(fold_random(cbind(xy[, 1], NA), k = 2), "Coordinates")
  expect_error(fold_random(m, k = 2), "numeric two-column coordinate matrix")
  expect_error(fold_random(~ x + y, as.list(m), k = 2), "must be a data.frame")
  ll <- sf::st_as_sf(data.frame(x = c(5, 6), y = c(52, 53)),
    coords = c("x", "y"), crs = 4326
  )
  expect_error(fold_blocks(ll, blocks = 2), "longitude/latitude")
})
