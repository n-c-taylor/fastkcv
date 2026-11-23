# Supported variogram models: all are bounded and positive (semi)definite in
# two dimensions, so C(h) = sill - gamma(h) is a valid covariance. The C++
# type code is the position minus 1; append new models at the end.
supported_models <- c(
  "Nug", "Exp", "Sph", "Gau", "Mat", "Ste",
  "Cir", "Pen", "Bes", "Wav", "Hol"
)

#' Validate a gstat variogram model for use in fastkcv
#'
#' Checks that `model` is a `gstat` `variogramModel` built only from bounded
#' components (`Nug`, `Exp`, `Sph`, `Gau`, `Mat`, `Ste`, `Cir`, `Pen`,
#' `Bes`, `Wav`, `Hol`, or nested sums of them) with at most 2D geometric
#' anisotropy.
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
    # gstat's Lin with a positive range is bounded (sill at the range), but
    # its covariance, the triangular function, is positive definite only in
    # one dimension; with range 0 it is unbounded.
    lin_bounded <- types == "Lin" & model$range > 0
    unbounded <- intersect(types[!lin_bounded], c("Lin", "Pow", "Log"))
    why <- if (any(lin_bounded)) {
      paste0(
        "The bounded linear model (Lin with a positive range) does not ",
        "give a valid covariance in two dimensions: its covariance ",
        "function, the triangular function, is positive definite only in ",
        "one dimension."
      )
    } else if (length(unbounded) > 0) {
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
#' its own 2D geometric anisotropy. The result is exactly symmetric. Used by
#' the tests and benchmarks; kcv() builds the matrix in place.
#'
#' @param model A `variogramModel`.
#' @param coords n x 2 coordinate matrix.
#' @return An n x n covariance matrix.
#' @noRd
covariance_matrix <- function(model, coords) {
  validate_model(model)
  coords <- as.matrix(coords)
  if (ncol(coords) != 2) {
    stop("Coordinates must have exactly 2 columns.", call. = FALSE)
  }
  storage.mode(coords) <- "double"
  do.call(cov_matrix_cpp, c(list(coords), model_cpp_args(model)))
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
