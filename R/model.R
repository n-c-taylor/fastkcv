# Supported variogram models: all are bounded, so C(h) = sill - gamma(h)
# is a valid covariance.
supported_models <- c("Nug", "Exp", "Sph", "Gau", "Mat", "Ste")

#' Validate a gstat variogram model for use in fastkcv
#'
#' Checks that `model` is a `gstat` `variogramModel` built only from bounded
#' components (`Nug`, `Exp`, `Sph`, `Gau`, `Mat`, `Ste`, or nested sums of
#' them) with at most 2D geometric anisotropy.
#'
#' @param model A `variogramModel`, as returned by [gstat::vgm()].
#' @return `model`, invisibly.
#' @noRd
validate_model <- function(model) {
  if (!inherits(model, "variogramModel")) {
    stop(
      "`model` must be a gstat variogram model (class 'variogramModel'), ",
      "e.g. gstat::vgm(1, \"Exp\", 300).",
      call. = FALSE
    )
  }
  if (nrow(model) == 0) {
    stop("`model` has no components.", call. = FALSE)
  }
  types <- as.character(model$model)
  bad <- setdiff(types, supported_models)
  if (length(bad) > 0) {
    unbounded <- intersect(bad, c("Lin", "Pow", "Log"))
    why <- if (length(unbounded) > 0) {
      paste0(
        "Unbounded variogram model(s) ", paste(unbounded, collapse = ", "),
        " have no covariance function C(h) = sill - gamma(h), which the ",
        "fast cross-validation formulas require."
      )
    } else {
      paste0(
        "Variogram model(s) ", paste(bad, collapse = ", "),
        " are not supported in fastkcv."
      )
    }
    stop(
      why, " Supported models: ", paste(supported_models, collapse = ", "),
      " (and nested sums of them).",
      call. = FALSE
    )
  }
  pars <- c("psill", "range", "kappa", "ang1", "anis1")
  if (!all(vapply(model[pars], function(v) all(is.finite(v)), logical(1)))) {
    stop("`model` has missing or non-finite parameters (",
      paste(pars, collapse = ", "), ").",
      call. = FALSE
    )
  }
  if (any(model$psill < 0)) {
    stop("`model` has a negative partial sill.", call. = FALSE)
  }
  if (sum(model$psill) <= 0) {
    stop("`model` has a total sill of zero.", call. = FALSE)
  }
  if (any(model$range[types != "Nug"] <= 0)) {
    stop("Every non-nugget component of `model` needs a positive range.",
      call. = FALSE
    )
  }
  if (any(model$kappa[types %in% c("Mat", "Ste")] <= 0)) {
    stop("Every Mat or Ste component of `model` needs a positive kappa.",
      call. = FALSE
    )
  }
  if (any(model$ang2 != 0) || any(model$ang3 != 0) || any(model$anis2 != 1)) {
    stop(
      "3D anisotropy (anis2, ang2, ang3) is not supported; fastkcv handles ",
      "2D coordinates with 2D geometric anisotropy only.",
      call. = FALSE
    )
  }
  if (any(model$anis1 <= 0 | model$anis1 > 1)) {
    stop("Anisotropy ratios (anis1) must lie in (0, 1].", call. = FALSE)
  }
  invisible(model)
}

#' Covariance matrix from a gstat variogram model
#'
#' Builds `C[i, j] = C(s_i - s_j)` where `C(h) = total sill - gamma(h)`, in
#' C++ (`src/covariance.cpp`) with gstat's model parameterizations. The nugget
#' contributes only where the lag vector is exactly zero. Each component keeps
#' its own 2D geometric anisotropy. When `coords2` is omitted the result is
#' exactly symmetric.
#'
#' @param model A `variogramModel`.
#' @param coords1 n1 x 2 coordinate matrix.
#' @param coords2 n2 x 2 coordinate matrix; defaults to `coords1`.
#' @return An n1 x n2 covariance matrix.
#' @noRd
covariance_matrix <- function(model, coords1, coords2 = NULL) {
  validate_model(model)
  symmetric <- is.null(coords2)
  coords1 <- as.matrix(coords1)
  coords2 <- if (symmetric) coords1 else as.matrix(coords2)
  if (ncol(coords1) != 2 || ncol(coords2) != 2) {
    stop("Coordinates must have exactly 2 columns.", call. = FALSE)
  }
  storage.mode(coords1) <- "double"
  storage.mode(coords2) <- "double"
  do.call(cov_matrix_cpp, c(
    list(coords1, coords2), model_cpp_args(model),
    list(symmetric = symmetric)
  ))
}

# The columns of a validated variogram model as passed to the C++ kernels
# (type code = position in `supported_models` minus 1).
model_cpp_args <- function(model) {
  list(
    type = match(as.character(model$model), supported_models) - 1L,
    psill = as.double(model$psill),
    range = as.double(model$range),
    kappa = as.double(model$kappa),
    ang1 = as.double(model$ang1),
    anis1 = as.double(model$anis1)
  )
}

#' Covariance matrix via gstat::variogramLine() (internal oracle)
#'
#' Builds `C[i, j] = C(s_i - s_j)` where `C(h) = total sill - gamma(h)`,
#' evaluated with [gstat::variogramLine()] (`covariance = TRUE`) so that the
#' model parameterization is exactly gstat's. Used to test
#' `covariance_matrix()`. The nugget contributes only at
#' distance exactly zero. Components with 2D geometric anisotropy are handled
#' by mapping lag vectors to gstat's reduced (isotropic) distance.
#'
#' @param model A `variogramModel`.
#' @param coords1 n1 x 2 coordinate matrix.
#' @param coords2 n2 x 2 coordinate matrix; defaults to `coords1`.
#' @return An n1 x n2 covariance matrix.
#' @noRd
covariance_matrix_r <- function(model, coords1, coords2 = coords1) {
  validate_model(model)
  coords1 <- as.matrix(coords1)
  coords2 <- as.matrix(coords2)
  if (ncol(coords1) != 2 || ncol(coords2) != 2) {
    stop("Coordinates must have exactly 2 columns.", call. = FALSE)
  }
  dx <- outer(coords1[, 1], coords2[, 1], "-")
  dy <- outer(coords1[, 2], coords2[, 2], "-")

  # Group components sharing the same anisotropy so that each group needs one
  # reduced-distance matrix and one variogramLine() call.
  key <- paste(model$ang1, model$anis1)
  out <- matrix(0, nrow(coords1), nrow(coords2))
  for (k in unique(key)) {
    rows <- model[key == k, , drop = FALSE]
    d <- reduced_distance(dx, dy, rows$ang1[1], rows$anis1[1])
    rows$ang1 <- 0
    rows$anis1 <- 1
    cv <- gstat::variogramLine(rows, dist_vector = as.vector(d), covariance = TRUE)
    out <- out + cv$gamma
  }
  out
}

#' Reduced distance under gstat's 2D geometric anisotropy
#'
#' gstat's `anis = c(p, s)`: `p` is the azimuth of the major axis (degrees,
#' clockwise from north, i.e. from the positive y axis) and `s` the ratio of
#' the minor to the major range. Lags along the minor axis are stretched by
#' `1 / s` before evaluating the isotropic model.
#'
#' @noRd
reduced_distance <- function(dx, dy, ang1, anis1) {
  if (ang1 == 0 && anis1 == 1) {
    return(sqrt(dx^2 + dy^2))
  }
  a <- ang1 * pi / 180
  h_major <- dx * sin(a) + dy * cos(a)
  h_minor <- dx * cos(a) - dy * sin(a)
  sqrt(h_major^2 + (h_minor / anis1)^2)
}
