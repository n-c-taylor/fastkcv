#' Fold builders for cross-validation
#'
#' Assign observations to cross-validation folds. Each function returns a
#' plain integer vector with one fold identifier per observation, to be passed
#' as `folds` to [kcv()] or [kcv_run()] (or as `nfold` to
#' [gstat::krige.cv()]).
#'
#' * `fold_random()` draws folds exactly as [gstat::krige.cv()] does for a
#'   scalar `nfold`: `sample(k, n, replace = TRUE)` when `k < n`, so the same
#'   seed gives the same folds in both packages. Fold sizes are random and a
#'   fold may be empty (then it is simply absent). With `k >= n` it returns
#'   `1:n` (leave-one-out).
#' * `fold_blocks()` lays a regular grid over the bounding box of the
#'   locations and makes each non-empty grid cell a fold. Specify the grid
#'   either by the number of cells (`blocks`) or by the cell size (`size`).
#'   Empty cells are dropped and the remaining cells are numbered
#'   consecutively (x fastest, then y, starting at the minimum corner), so the
#'   number of folds can be smaller than the number of cells. With `k`, the
#'   blocks are grouped at random into `k` folds, each receiving the same
#'   number of blocks (up to one).
#' * `fold_kmeans()` clusters the coordinates with [stats::kmeans()] and makes
#'   each cluster a fold. The result depends on the random starts, so set a
#'   seed for reproducibility. On large, strongly clustered data the default
#'   Hartigan-Wong algorithm can reach its internal step limit
#'   ("Quick-TRANSfer stage steps exceeded"); the clustering is then still a
#'   valid set of folds, possibly not a locally optimal one, so this warning
#'   is not shown unless `algorithm` is given explicitly.
#'
#' Coordinates are used as given (projected units); all three functions reject
#' geographic (longitude/latitude) coordinates, as [kcv()] does.
#'
#' @param locations A numeric two-column coordinate matrix, a formula giving
#'   the coordinates (e.g. `~x + y`) with a data.frame in `data`, an `sf`
#'   object with POINT geometries, or an `sp` `Spatial*` object.
#' @param data A data.frame, when `locations` is a formula.
#' @param k Number of folds. For `fold_blocks()`, optional: the number of
#'   folds into which the blocks are grouped at random.
#' @param blocks Number of grid cells along x and y: a single integer (used
#'   for both) or a vector `c(nx, ny)`. The cells divide the bounding box of
#'   the locations into equal parts.
#' @param size Cell side length along x and y, in coordinate units: a single
#'   number (square cells) or `c(dx, dy)`. The grid starts at the minimum
#'   corner of the bounding box.
#' @param nstart,iter.max Passed to [stats::kmeans()].
#' @param ... Further arguments to [stats::kmeans()], such as `algorithm`.
#'
#' @return An integer vector of fold identifiers, one per observation.
#' @seealso [kcv()], [kcv_run()], [kcv_metrics()]; for buffered
#'   leave-one-out, the `buffer` argument of [kcv()].
#' @name folds
#' @examples
#' data(meuse, package = "sp")
#' v <- gstat::vgm(0.59, "Sph", 897, 0.05)
#' set.seed(1)
#' f_rand <- fold_random(~ x + y, meuse, k = 10)
#' f_block <- fold_blocks(~ x + y, meuse, blocks = 4)
#' f_km <- fold_kmeans(~ x + y, meuse, k = 10)
#' table(f_block)
#'
#' p <- kcv_prepare(log(zinc) ~ 1, ~ x + y, meuse, model = v)
#' rbind(
#'   random = kcv_metrics(kcv_run(p, folds = f_rand)),
#'   blocks = kcv_metrics(kcv_run(p, folds = f_block)),
#'   kmeans = kcv_metrics(kcv_run(p, folds = f_km))
#' )
NULL

#' @rdname folds
#' @export
fold_random <- function(locations, data = NULL, k) {
  n <- nrow(fold_coords(locations, data))
  if (missing(k)) stop("`k` (the number of folds) is required.", call. = FALSE)
  check_count(k, "k")
  if (k < n) {
    sample(k, n, replace = TRUE)
  } else {
    seq_len(n)
  }
}

#' @rdname folds
#' @export
fold_blocks <- function(locations, data = NULL, blocks = NULL, size = NULL,
                        k = NULL) {
  xy <- fold_coords(locations, data)
  if (is.null(blocks) == is.null(size)) {
    stop("Give exactly one of `blocks` (number of cells) and `size` ",
      "(cell side length).",
      call. = FALSE
    )
  }
  lo <- c(min(xy[, 1]), min(xy[, 2]))
  ext <- c(max(xy[, 1]), max(xy[, 2])) - lo
  if (!is.null(blocks)) {
    if (!is.numeric(blocks) || !length(blocks) %in% 1:2 || anyNA(blocks) ||
      any(blocks < 1) || any(blocks != round(blocks))) {
      stop("`blocks` must be one or two positive integers.", call. = FALSE)
    }
    nb <- rep_len(as.double(blocks), 2)
    width <- ext / nb
  } else {
    if (!is.numeric(size) || !length(size) %in% 1:2 || anyNA(size) ||
      any(!is.finite(size)) || any(size <= 0)) {
      stop("`size` must be one or two positive numbers.", call. = FALSE)
    }
    width <- rep_len(as.double(size), 2)
    nb <- pmax(1, ceiling(ext / width))
  }
  if (prod(nb) > 2^53) {
    stop("The grid has too many cells; use larger `size` or fewer `blocks`.",
      call. = FALSE
    )
  }
  # Cell numbers in double arithmetic (exact up to 2^53): with small cells
  # the counts and their product overflow integers.
  cell_index <- function(v, lo, width, nb) {
    if (nb == 1 || width == 0) {
      return(rep(1, length(v)))
    }
    pmin(pmax(floor((v - lo) / width) + 1, 1), nb)
  }
  ix <- cell_index(xy[, 1], lo[1], width[1], nb[1])
  iy <- cell_index(xy[, 2], lo[2], width[2], nb[2])
  cell <- (iy - 1) * nb[1] + ix
  block <- match(cell, sort(unique(cell)))
  if (is.null(k)) {
    return(block)
  }
  check_count(k, "k")
  n_blocks <- max(block)
  if (k > n_blocks) {
    stop("`k` = ", k, " is larger than the number of non-empty blocks (",
      n_blocks, ").",
      call. = FALSE
    )
  }
  sample(rep_len(seq_len(k), n_blocks))[block]
}

#' @rdname folds
#' @export
fold_kmeans <- function(locations, data = NULL, k, nstart = 10,
                        iter.max = 100, ...) {
  xy <- fold_coords(locations, data)
  if (missing(k)) stop("`k` (the number of folds) is required.", call. = FALSE)
  check_count(k, "k")
  n_distinct <- nrow(unique(xy))
  if (k > n_distinct) {
    stop("`k` = ", k, " is larger than the number of distinct locations (",
      n_distinct, ").",
      call. = FALSE
    )
  }
  if (k == 1) {
    return(rep(1L, nrow(xy)))
  }
  kmeans_folds(xy, k, nstart = nstart, iter.max = iter.max, ...)
}

# stats::kmeans() cluster labels. The default Hartigan-Wong algorithm can
# stop at its Quick-TRANSfer step limit on large clustered data and warn;
# the clustering it returns is still a valid partition (only possibly not a
# local optimum), which is all fold assignment needs, so that one warning
# is muffled unless the user chose the algorithm.
kmeans_folds <- function(xy, k, ..., .kmeans = stats::kmeans) {
  muffle <- !"algorithm" %in% names(list(...))
  res <- withCallingHandlers(
    .kmeans(xy, centers = k, ...),
    warning = function(w) {
      if (muffle && grepl("Quick-TRANSfer", conditionMessage(w))) {
        invokeRestart("muffleWarning")
      }
    }
  )
  as.integer(res$cluster)
}

fold_coords <- function(locations, data) {
  parse_coords(locations, data, matrix_ok = TRUE, need_data = FALSE)$coords
}

check_count <- function(k, name) {
  if (!is.numeric(k) || length(k) != 1 || is.na(k) || k < 1 ||
    k != round(k)) {
    stop("`", name, "` must be a single positive integer.", call. = FALSE)
  }
  invisible(TRUE)
}
