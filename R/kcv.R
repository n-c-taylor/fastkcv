#' Fast exact cross-validation for global kriging
#'
#' A drop-in companion to [gstat::krige.cv()] for global kriging (every
#' observation in the kriging neighbourhood). Instead of refitting the kriging
#' system once per fold, `kcv()` factorises the covariance matrix once and
#' obtains every fold's predictions and variances from it. Results match
#' `krige.cv()` to numerical precision.
#'
#' Ordinary kriging (`formula = z ~ 1`) and universal kriging (`z ~ x + y`,
#' external drift, ...) re-estimate the trend by generalised least squares in
#' each fold, as gstat does. Simple kriging is used when `beta` is supplied.
#'
#' @param formula Formula defining the response and the trend, e.g.
#'   `log(zinc) ~ 1` or `log(zinc) ~ sqrt(dist)`.
#' @param locations A formula giving the coordinates (e.g. `~x + y`) when
#'   `data` is a data.frame, or an `sf` object or `sp` `Spatial*DataFrame`
#'   holding both the locations and the variables (then `data` is not used).
#' @param data A data.frame, when `locations` is a formula.
#' @param model A `gstat` variogram model ([gstat::vgm()]) with fixed
#'   parameters. Supported: `Nug`, `Exp`, `Sph`, `Gau`, `Mat`, `Ste` and nested
#'   sums of them, with optional 2D geometric anisotropy.
#' @param nfold As in [gstat::krige.cv()]: missing or `nrow(data)` for
#'   leave-one-out; a single integer `k < n` for random k-fold assignment
#'   (drawn exactly as `krige.cv()` does, so the same seed gives the same
#'   folds); or a vector of length `n` of fold identifiers.
#' @param folds Optional vector of length `n` of fold identifiers. Takes
#'   precedence over `nfold`.
#' @param beta Optional trend coefficients; if given, simple kriging with this
#'   known trend is used.
#' @param ... Accepted for compatibility with `krige.cv()`: `verbose` and
#'   `debug.level` are ignored. Local-neighbourhood or block arguments
#'   (`nmax`, `nmin`, `maxdist`, `block`) are not supported and raise an error.
#'
#' @return An object of the same kind as the input (a data.frame with the
#'   coordinate columns appended, an `sf` object, or an `sp`
#'   `SpatialPointsDataFrame`) with columns `var1.pred`, `var1.var`,
#'   `observed`, `residual` (observed minus predicted), `zscore` (residual
#'   divided by the kriging standard error) and `fold`, as in `krige.cv()`.
#' @export
#' @examples
#' library(sp)
#' data(meuse, package = "sp")
#' v <- gstat::vgm(0.59, "Sph", 897, 0.05)
#' cv <- kcv(log(zinc) ~ 1, ~ x + y, meuse, model = v)
#' head(cv)
#' sqrt(mean(cv$residual^2))
kcv <- function(formula, locations, data, model, nfold, folds = NULL,
                beta = NULL, ...) {
  check_dots(...)

  # krige.cv(formula, <sf/sp>, model) passes the model in third position;
  # accept that too.
  if (!inherits(locations, "formula") && missing(model) && !missing(data) &&
    inherits(data, "variogramModel")) {
    model <- data
  }
  if (missing(model)) {
    stop("`model` (a gstat variogram model) is required.", call. = FALSE)
  }
  if (inherits(locations, "formula") && missing(data)) {
    stop("When `locations` is a formula, `data` must be supplied.",
      call. = FALSE
    )
  }
  validate_model(model)

  inp <- parse_inputs(formula, locations, if (missing(data)) NULL else data)
  n <- nrow(inp$y)

  fold <- if (!is.null(folds)) {
    if (length(folds) != n || anyNA(folds)) {
      stop("`folds` must have one non-missing value per observation (n = ", n,
        ").",
        call. = FALSE
      )
    }
    folds
  } else if (missing(nfold)) {
    seq_len(n)
  } else {
    resolve_nfold(nfold, n)
  }

  C <- covariance_matrix(model, inp$coords)
  res <- kcv_engine(C, inp$F, inp$y, fold, beta = beta)
  build_output(inp, res)
}

# Mirrors krige.cv()'s interpretation of `nfold`.
resolve_nfold <- function(nfold, n) {
  if (length(nfold) == n) {
    if (anyNA(nfold)) {
      stop("`nfold` must not contain missing values.", call. = FALSE)
    }
    return(nfold)
  }
  if (length(nfold) != 1 || !is.numeric(nfold) || is.na(nfold) || nfold < 1 ||
    nfold != round(nfold)) {
    stop(
      "`nfold` must be a single positive integer or a vector of fold ids of ",
      "length n = ", n, ".",
      call. = FALSE
    )
  }
  if (nfold < n) {
    sample(nfold, n, replace = TRUE)
  } else {
    seq_len(n)
  }
}

check_dots <- function(...) {
  dots <- names(list(...))
  if (length(list(...)) > 0 && (is.null(dots) || any(dots == ""))) {
    stop("Unexpected unnamed argument(s) in `...`.", call. = FALSE)
  }
  local_args <- intersect(dots, c("nmax", "nmin", "maxdist", "block"))
  if (length(local_args) > 0) {
    stop(
      "Argument(s) ", paste0("`", local_args, "`", collapse = ", "),
      " are not supported. fastkcv performs global kriging (every ",
      "observation in the neighbourhood) with point support; local ",
      "neighbourhoods and block kriging change the kriging system for each ",
      "prediction, so the single-factorisation shortcut does not apply. ",
      "Use gstat::krige.cv() for those.",
      call. = FALSE
    )
  }
  unknown <- setdiff(dots, c("verbose", "debug.level"))
  if (length(unknown) > 0) {
    stop("Unknown argument(s): ", paste0("`", unknown, "`", collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

build_output <- function(inp, res) {
  observed <- drop(inp$y)
  pred <- drop(res$pred)
  residual <- observed - pred
  out <- data.frame(
    var1.pred = pred,
    var1.var = res$var,
    observed = observed,
    residual = residual,
    zscore = residual / sqrt(res$var),
    fold = res$folds
  )
  switch(inp$kind,
    data.frame = {
      cc <- as.data.frame(inp$coords)
      names(cc) <- inp$coord_names
      cbind(out, cc)
    },
    sf = sf::st_sf(out, geometry = sf::st_geometry(inp$locations)),
    sp = sp::SpatialPointsDataFrame(
      inp$coords, out,
      proj4string = inp$locations@proj4string
    )
  )
}
