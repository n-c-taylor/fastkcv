#' Parse kcv() inputs into the engine's matrix form
#'
#' Accepts the same shapes as [gstat::krige.cv()]: a `data.frame` with
#' `locations = ~x + y`, an `sf` object, or an `sp` `Spatial*DataFrame`.
#'
#' @return A list with `y` (n x 1 matrix), `F` (n x p trend matrix), `coords`
#'   (n x 2 matrix), `kind` ("data.frame", "sf" or "sp"), `locations` (the
#'   spatial input, for rebuilding output), and `coord_names`.
#' @noRd
parse_inputs <- function(formula, locations, data) {
  if (!inherits(formula, "formula") || length(formula) != 3) {
    stop("`formula` must be a two-sided formula such as log(zinc) ~ 1.",
      call. = FALSE
    )
  }

  cc <- parse_coords(locations, data, matrix_ok = FALSE, need_data = TRUE)
  coords <- cc$coords
  df <- cc$df
  if (anyDuplicated(coords) > 0) {
    dup <- which(duplicated(coords))
    stop(
      "Duplicate locations (rows ", paste(utils::head(dup, 5), collapse = ", "),
      if (length(dup) > 5) ", ..." else "", ") make the covariance matrix ",
      "singular for kriging. Remove or average duplicates first ",
      "(see sp::zerodist()).",
      call. = FALSE
    )
  }

  mf <- stats::model.frame(formula, df, na.action = stats::na.pass)
  y <- stats::model.response(mf)
  if (!is.numeric(y) || !is.null(dim(y))) {
    stop("The response in `formula` must be a numeric vector.", call. = FALSE)
  }
  if (anyNA(y)) {
    stop(
      "The response has ", sum(is.na(y)), " missing value(s). ",
      "Remove them before cross-validation.",
      call. = FALSE
    )
  }
  F <- stats::model.matrix(stats::delete.response(stats::terms(formula)), mf)
  if (anyNA(F)) {
    stop("Trend covariates in `formula` contain missing values.", call. = FALSE)
  }
  check_trend(F)
  attr(F, "assign") <- NULL
  attr(F, "contrasts") <- NULL

  list(
    y = matrix(as.numeric(y), ncol = 1),
    F = F,
    coords = coords,
    coord_names = cc$coord_names,
    kind = cc$kind,
    locations = if (cc$kind == "data.frame") NULL else locations
  )
}

lonlat_message <- function() {
  paste0(
    "Geographic (longitude/latitude) coordinates are not supported in ",
    "fastkcv: gstat uses great-circle distances for them. Project the ",
    "data first, e.g. with sf::st_transform()."
  )
}

#' Coordinates of the observations
#'
#' Shared by `parse_inputs()` and the fold builders.
#'
#' @param matrix_ok Accept a numeric two-column coordinate matrix.
#' @param need_data For sp input, require a `Spatial*DataFrame`.
#' @return A list with `coords` (unnamed n x 2 double matrix), `coord_names`,
#'   `kind` ("data.frame", "sf", "sp" or "matrix") and `df` (the attribute
#'   table, `NULL` for a matrix).
#' @noRd
parse_coords <- function(locations, data, matrix_ok = FALSE,
                         need_data = TRUE) {
  df <- NULL
  if (inherits(locations, "formula")) {
    if (!is.data.frame(data)) {
      stop("When `locations` is a formula, `data` must be a data.frame.",
        call. = FALSE
      )
    }
    kind <- "data.frame"
    coords <- as.matrix(stats::model.frame(locations, data,
      na.action = stats::na.pass
    ))
    df <- as.data.frame(data)
  } else if (inherits(locations, "sf")) {
    kind <- "sf"
    geom_type <- unique(as.character(sf::st_geometry_type(locations)))
    if (!identical(geom_type, "POINT")) {
      stop("`locations` (sf) must contain POINT geometries only.", call. = FALSE)
    }
    if (isTRUE(sf::st_is_longlat(locations))) {
      stop(lonlat_message(), call. = FALSE)
    }
    coords <- sf::st_coordinates(locations)
    df <- sf::st_drop_geometry(locations)
  } else if (inherits(locations, "Spatial")) {
    kind <- "sp"
    if (identical(sp::is.projected(locations), FALSE)) {
      stop(lonlat_message(), call. = FALSE)
    }
    has_data <- methods::.hasSlot(locations, "data")
    if (need_data && !has_data) {
      stop("`locations` (sp) must be a Spatial*DataFrame holding the variables.",
        call. = FALSE
      )
    }
    coords <- sp::coordinates(locations)
    # as.data.frame() of an sp object includes the coordinate columns, as in
    # gstat, so formulas such as z ~ x + y work.
    if (has_data) df <- as.data.frame(locations)
  } else if (matrix_ok && is.matrix(locations) && is.numeric(locations)) {
    kind <- "matrix"
    coords <- locations
  } else {
    stop(
      "`locations` must be a formula (with a data.frame in `data`), ",
      if (matrix_ok) "a numeric two-column coordinate matrix, " else "",
      "an sf object, or an sp Spatial*DataFrame.",
      call. = FALSE
    )
  }

  if (ncol(coords) != 2) {
    stop("fastkcv supports 2D coordinates only; got ", ncol(coords),
      " coordinate columns.",
      call. = FALSE
    )
  }
  coord_names <- colnames(coords)
  coords <- unname(as.matrix(coords))
  storage.mode(coords) <- "double"
  if (any(!is.finite(coords))) {
    stop("Coordinates contain missing or non-finite values.", call. = FALSE)
  }
  list(coords = coords, coord_names = coord_names, kind = kind, df = df)
}
