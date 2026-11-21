#' Summary statistics of a cross-validation
#'
#' Summarizes the output of [kcv()] or [kcv_run()] (or of
#' [gstat::krige.cv()], which has the same columns), overall or per fold.
#'
#' With residuals \eqn{e_i} (observed minus predicted), kriging variances
#' \eqn{\sigma_i^2} and z-scores \eqn{z_i = e_i / \sigma_i}:
#'
#' * `rmse`: \eqn{\sqrt{\mathrm{mean}(e_i^2)}};
#' * `mae`: \eqn{\mathrm{mean}(|e_i|)};
#' * `bias`: \eqn{\mathrm{mean}(e_i)} (positive when the predictions are too
#'   low on average);
#' * `mean_z`, `var_z`: mean and variance of the z-scores, ideally close to 0
#'   and 1. `var_z` uses the \eqn{n - 1} denominator of [stats::var()], so it
#'   is `NA` for a fold with a single observation;
#' * `coverage_<100 L>` for each `level` \eqn{L}: the proportion of
#'   observations inside the Gaussian prediction interval
#'   \eqn{\hat y_i \pm q \sigma_i}, \eqn{q} = `qnorm((1 + L) / 2)`, i.e.
#'   with \eqn{|z_i| \le q}. Well-calibrated variances give coverage close to
#'   \eqn{L}. The default levels 50, 80, 90 and 95% probe both the center and
#'   the tails of the predictive distribution.
#'
#' The z-scores are recomputed from the `residual` and `var1.var` columns.
#'
#' @param x Output of [kcv()] or [kcv_run()]: a data.frame, an `sf` object or
#'   an `sp` `SpatialPointsDataFrame`, with columns `residual` and `var1.var`
#'   (and `fold` when `by_fold = TRUE`).
#' @param by_fold If `TRUE`, one row per fold (in `sort(unique(fold))` order)
#'   instead of one overall row.
#' @param level Nominal coverage levels of the prediction intervals, in
#'   (0, 1).
#' @param by_response If `TRUE` (the default) and `x` has a `response`
#'   column (the output of a cross-validation of several responses), one row
#'   per response (in the order of its levels), or per response and fold
#'   with `by_fold = TRUE`. `FALSE` pools the responses. Ignored without a
#'   `response` column.
#' @return A data.frame with columns `n`, `rmse`, `mae`, `bias`, `mean_z`,
#'   `var_z` and one `coverage_*` column per level, preceded by `fold` when
#'   `by_fold = TRUE`, and by `response` for several responses.
#' @seealso [kcv()], [fold_blocks()] and the other fold builders.
#' @export
#' @examples
#' data(meuse, package = "sp")
#' v <- gstat::vgm(0.59, "Sph", 897, 0.05)
#' cv <- kcv(log(zinc) ~ 1, ~ x + y, meuse, model = v)
#' kcv_metrics(cv)
#'
#' cv5 <- kcv(log(zinc) ~ 1, ~ x + y, meuse, model = v,
#'   folds = fold_blocks(~ x + y, meuse, blocks = 2))
#' kcv_metrics(cv5, by_fold = TRUE)
kcv_metrics <- function(x, by_fold = FALSE,
                        level = c(0.5, 0.8, 0.9, 0.95), by_response = TRUE) {
  if (inherits(x, "sf")) {
    x <- sf::st_drop_geometry(x)
  } else if (inherits(x, "Spatial")) {
    if (!methods::.hasSlot(x, "data")) {
      stop("`x` must have an attribute table (a Spatial*DataFrame).",
        call. = FALSE
      )
    }
    x <- x@data
  }
  if (!is.data.frame(x)) {
    stop("`x` must be the output of kcv() or kcv_run(): a data.frame, sf or ",
      "sp object.",
      call. = FALSE
    )
  }
  if (!is.logical(by_fold) || length(by_fold) != 1 || is.na(by_fold)) {
    stop("`by_fold` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(by_response) || length(by_response) != 1 ||
    is.na(by_response)) {
    stop("`by_response` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.numeric(level) || length(level) == 0 || anyNA(level) ||
    any(level <= 0 | level >= 1) || anyDuplicated(level)) {
    stop("`level` must be a numeric vector of distinct levels strictly ",
      "between 0 and 1.",
      call. = FALSE
    )
  }
  need <- c("residual", "var1.var", if (by_fold) "fold")
  missing_cols <- setdiff(need, names(x))
  if (length(missing_cols) > 0) {
    stop("`x` lacks column(s) ", paste0("`", missing_cols, "`", collapse = ", "),
      "; it should be the output of kcv() or kcv_run().",
      call. = FALSE
    )
  }
  e <- x$residual
  v <- x$var1.var
  if (!is.numeric(e) || !is.numeric(v)) {
    stop("Columns `residual` and `var1.var` must be numeric.", call. = FALSE)
  }
  q <- stats::qnorm((1 + level) / 2)
  cov_names <- paste0("coverage_", format(100 * level, trim = TRUE,
    drop0trailing = TRUE))

  summarize <- function(e, v) {
    z <- e / sqrt(v)
    out <- data.frame(
      n = length(e),
      rmse = sqrt(mean(e^2)),
      mae = mean(abs(e)),
      bias = mean(e),
      mean_z = mean(z),
      var_z = if (length(z) > 1) stats::var(z) else NA_real_
    )
    for (j in seq_along(q)) {
      out[[cov_names[j]]] <- mean(abs(e) <= q[j] * sqrt(v))
    }
    out
  }

  one <- function(i) {
    if (!by_fold) {
      return(summarize(e[i], v[i]))
    }
    fold <- x$fold[i]
    folds <- sort(unique(fold))
    rows <- lapply(folds, function(f) {
      j <- i[fold == f]
      summarize(e[j], v[j])
    })
    cbind(fold = folds, do.call(rbind, rows))
  }

  if (!by_response || !("response" %in% names(x))) {
    return(one(seq_along(e)))
  }
  resp <- x$response
  lev <- if (is.factor(resp)) levels(resp) else sort(unique(resp))
  lev <- lev[lev %in% resp]
  parts <- lapply(lev, function(r) {
    out <- one(which(resp == r))
    cbind(response = rep(r, nrow(out)), out)
  })
  out <- do.call(rbind, parts)
  if (is.factor(resp)) out$response <- factor(out$response, levels = lev)
  rownames(out) <- NULL
  out
}
