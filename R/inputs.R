#' Parse kcv() inputs into the engine's matrix form
#'
#' Accepts the same shapes as [gstat::krige.cv()]: a `data.frame` with
#' `locations = ~x + y`, an `sf` object, or an `sp` `Spatial*DataFrame`.
#'
#' @return A list with `y` (n x T matrix; T = 1 for a vector response),
#'   `multi` (whether the response was a matrix), `F` (n x p trend matrix),
#'   `coords` (n x 2 matrix), `kind` ("data.frame", "sf" or "sp"),
#'   `locations` (the spatial input, for rebuilding output), and
#'   `coord_names`.
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
  # A matrix response (cbind(z1, z2) ~ ..., or a matrix column) gives
  # several responses cross-validated with one factorization.
  # model.response() drops a one-column matrix to a vector (as for
  # scale(z) ~ 1, which stays a single response); an explicit cbind() keeps
  # it as one of "several" responses.
  multi <- is.matrix(y) || (is_cbind(formula[[2]]) && is.matrix(mf[[1L]]))
  if (multi) {
    if (!is.matrix(y)) y <- mf[[1L]]
    y <- check_responses(y, response_names(formula, y), "The response")
  } else {
    if (!is.null(dim(y))) {
      stop("The response in `formula` must be a numeric vector.",
        call. = FALSE
      )
    }
    y <- check_responses(y, what = "The response")
  }
  F <- stats::model.matrix(stats::delete.response(stats::terms(formula)), mf)
  if (any(!is.finite(F))) {
    stop("Trend covariates in `formula` contain missing or non-finite ",
      "values (e.g. log(0)).",
      call. = FALSE
    )
  }
  check_trend(F)
  attr(F, "assign") <- NULL
  attr(F, "contrasts") <- NULL

  list(
    y = y,
    multi = multi,
    F = F,
    coords = coords,
    coord_names = cc$coord_names,
    kind = cc$kind,
    locations = if (cc$kind == "data.frame") NULL else locations
  )
}

# Names of the columns of a matrix response: its column names, else the
# deparsed arguments of a cbind() left-hand side (cbind() names only bare
# symbols), else "1", ..., "T".
response_names <- function(formula, y) {
  nm <- colnames(y)
  if (is.null(nm)) nm <- rep("", ncol(y))
  lhs <- formula[[2]]
  if (is_cbind(lhs) && length(lhs) - 1 == ncol(y)) {
    args <- vapply(as.list(lhs)[-1], function(a) {
      paste(deparse(a, width.cutoff = 500L), collapse = " ")
    }, "")
    unnamed <- is.na(nm) | nm == ""
    nm[unnamed] <- args[unnamed]
  }
  if (any(is.na(nm) | nm == "")) nm <- as.character(seq_len(ncol(y)))
  nm
}

is_cbind <- function(lhs) {
  is.call(lhs) && identical(lhs[[1]], as.name("cbind"))
}

# Checks the responses and returns them as a double matrix without row
# names: a vector (one response) as an n x 1 matrix, an n x T matrix with
# column names `nm`. `what` starts the messages. Every response must be
# observed at every location; a single response reports missing and
# infinite values separately (log(0) is the usual source of the latter).
check_responses <- function(y, nm = NULL, what) {
  if (is.null(dim(y))) {
    if (!is.numeric(y)) {
      stop(what, " must be a numeric vector or matrix.", call. = FALSE)
    }
    if (anyNA(y)) {
      stop(
        what, " has ", sum(is.na(y)), " missing value(s). ",
        "Remove them before cross-validation.",
        call. = FALSE
      )
    }
    if (any(is.infinite(y))) {
      stop(
        what, " has ", sum(is.infinite(y)), " infinite value(s) ",
        "(e.g. log(0)). Remove them before cross-validation.",
        call. = FALSE
      )
    }
    return(matrix(as.double(y), ncol = 1))
  }
  if (!is.numeric(y) || length(dim(y)) != 2) {
    stop(what, " must be a numeric vector or matrix.", call. = FALSE)
  }
  if (ncol(y) == 0) {
    stop(what, " has no columns.", call. = FALSE)
  }
  if (anyDuplicated(nm)) {
    stop(what, " has duplicate column names (",
      paste0("`", unique(nm[duplicated(nm)]), "`", collapse = ", "),
      "); each response needs a distinct name.",
      call. = FALSE
    )
  }
  y <- matrix(as.double(y), nrow(y), ncol(y), dimnames = list(NULL, nm))
  bad <- colSums(!is.finite(y))
  if (any(bad > 0)) {
    j <- which(bad > 0)[1]
    stop(
      "Response `", nm[j], "` has ", bad[j], " missing or non-finite ",
      "value(s)",
      if (sum(bad > 0) > 1) {
        paste0(" (and ", sum(bad > 0) - 1, " other response(s))")
      },
      ". Every response must be observed at every location (a response ",
      "cannot drop its own rows); remove those locations for all responses ",
      "or cross-validate that response separately.",
      call. = FALSE
    )
  }
  y
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
