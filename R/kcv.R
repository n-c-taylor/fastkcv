#' Fast exact cross-validation for global kriging
#'
#' A drop-in companion to [gstat::krige.cv()] for global kriging (every
#' observation in the kriging neighborhood). Instead of refitting the kriging
#' system once per fold, `kcv()` factorizes the covariance matrix once and
#' obtains every fold's predictions and variances from it. Results match
#' `krige.cv()` to numerical precision.
#'
#' Ordinary kriging (`formula = z ~ 1`) and universal kriging (`z ~ x + y`,
#' external drift, ...) re-estimate the trend by generalized least squares in
#' each fold, as gstat does. Simple kriging is used when `beta` is supplied.
#'
#' @param formula Formula defining the response and the trend, e.g.
#'   `log(zinc) ~ 1` or `log(zinc) ~ sqrt(dist)`. A matrix response, such
#'   as `cbind(log(zinc), log(copper)) ~ 1` or a matrix column of `data`,
#'   cross-validates each column with the same model, trend and folds; see
#'   "Several responses" below.
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
#'   known trend is used. A vector of length p (the number of trend
#'   coefficients) is shared by all responses; with several responses, a
#'   p x T matrix gives one column of coefficients per response.
#' @param buffer Optional radius for buffered leave-one-out, in coordinate
#'   units (\eqn{\ge 0}). Each observation is predicted from the observations
#'   farther than `buffer` from it; see "Buffered leave-one-out" below.
#'   Cannot be combined with `nfold` or `folds`.
#' @param targets Optional, with `buffer`: the observations to predict, as
#'   row numbers or a logical vector of length `n`. The hold-out sets are
#'   still formed from all observations. Useful to estimate buffered
#'   leave-one-out statistics from a sample of points when the radius is
#'   large. By default every observation is predicted.
#' @param ... Accepted for compatibility with `krige.cv()`: `verbose` and
#'   `debug.level` are ignored. Local-neighborhood or block arguments
#'   (`nmax`, `nmin`, `maxdist`, `block`) are not supported and raise an error.
#'
#' @return An object of the same kind as the input (a data.frame with the
#'   coordinate columns appended, an `sf` object, or an `sp`
#'   `SpatialPointsDataFrame`) with columns `var1.pred`, `var1.var`,
#'   `observed`, `residual` (observed minus predicted), `zscore` (residual
#'   divided by the kriging standard error) and `fold`, as in `krige.cv()`.
#'   With `buffer`, `fold` is the row number (each observation is its own
#'   hold-out set) and a column `n_heldout` gives the number of observations
#'   held out for it, itself included; with `targets`, only the target rows
#'   are returned, in increasing row order. With a matrix response, the
#'   results of the responses are stacked (all rows of the first response,
#'   then the second, ...) and a factor column `response`, after `fold`
#'   (and `n_heldout`), names the response of each row.
#'
#' @section Several responses:
#' The kriging weights depend only on the locations, the variogram model,
#' the trend and the folds, not on the observed values. A matrix response
#' with T columns is therefore cross-validated from one factorization: each
#' additional response costs a few matrix products (no further
#' factorization), and in buffered leave-one-out the hold-out sets and their
#' factorizations are shared too. Typical uses are many simulated fields,
#' several variables or time slices observed at the same locations, or
#' bootstrap replicates. [kcv_run()] with `y` cross-validates new responses
#' from an existing prepared object.
#'
#' A response with one column (such as `scale(z) ~ 1`, or a one-column
#' matrix column) is a single response, as in [stats::model.response()],
#' unless the left-hand side is an explicit `cbind()`.
#'
#' Every response is predicted at every location, so a missing or
#' non-finite value in any column is an error. Response names come from
#' the column names, else from the arguments of `cbind()`, else they are
#' `"1"`, `"2"`, .... Each response's rows are the output of the
#' single-response call, up to rounding in the last bits (with several
#' columns the BLAS uses matrix-matrix products). [kcv_metrics()]
#' summarizes each response separately. The output holds about
#' \eqn{10 n} numbers per response; several responses are checked against
#' the `fastkcv.max_gb` budget (see [kcv_engine()]) left after the
#' factorization, and the error message gives a batch size that fits.
#'
#' @section Buffered leave-one-out:
#' With `buffer = r`, observation \eqn{i} is held out together with every
#' observation within Euclidean distance \eqn{r} of it (distance
#' \eqn{\le r}, computed on the coordinates as given, without any variogram
#' anisotropy), and only \eqn{i} is predicted and reported. This removes the
#' spatially correlated neighbors that make ordinary leave-one-out
#' optimistic. `buffer = 0` gives ordinary leave-one-out exactly.
#'
#' Hold-out sets of neighboring points overlap, so each point is a separate
#' kriging problem, solved from the same factorization with a hold-out set
#' \eqn{I} of size \eqn{m}, which needs the block \eqn{Q_{II}} of the
#' matrix \eqn{Q} described in [kcv_engine()]. Two strategies are used,
#' whichever needs fewer operations:
#'
#' * per point: form \eqn{Q_{II}} directly, about \eqn{m^2 n} operations
#'   per point and \eqn{8 m n} bytes of temporary memory (160 MB for
#'   \eqn{m = 1000}, \eqn{n = 20000});
#' * shared pairs: two points in a common hold-out set are at most
#'   \eqn{2r} apart, so the entries of \eqn{Q} for all pairs of points
#'   within \eqn{2r} are computed once (about \eqn{n} operations each,
#'   roughly \eqn{2 \bar m n^2} in total) and each \eqn{Q_{II}} is
#'   assembled from them. The cache takes 12 bytes per pair (about
#'   \eqn{24 \bar m n} bytes) and is used only if it fits, together with the
#'   factorization, in the `fastkcv.max_gb` budget (see [kcv_engine()]).
#'
#' Each point also needs a Cholesky factorization of \eqn{Q_{II}}
#' (\eqn{m^3 / 3} operations). Since \eqn{m} grows with the area of the
#' buffer, large radii are expensive at large \eqn{n}; `targets` restricts
#' the predictions to a sample of points, at a proportional saving. These
#' computations consist of many small matrix operations, which a
#' multithreaded BLAS does not speed up and can slow down: with OpenBLAS,
#' setting the environment variable `OPENBLAS_NUM_THREADS=1` before starting
#' R is usually as fast or faster. Use [kcv_prepare()] and [kcv_run()] to try
#' several radii from one factorization.
#'
#' An error names the first point that cannot be predicted: a buffer that
#' holds out every observation, or (universal kriging) one that removes all
#' observations that identify the trend, for example every observation of
#' one level of a factor.
#'
#' @section Performance:
#' Almost all of the time goes into one dense Cholesky factorization and
#' triangular inverse, which R delegates to its BLAS/LAPACK libraries. The
#' reference BLAS shipped with many R builds is single-threaded and slow; an
#' optimized BLAS is typically many times faster on these operations. Check
#' which one R uses with `extSoftVersion()["BLAS"]` or `sessionInfo()`.
#'
#' * Linux (Debian/Ubuntu): install an optimized BLAS such as OpenBLAS
#'   (`sudo apt install libopenblas0-pthread`); the system's alternatives
#'   mechanism switches R to it without reinstalling anything.
#' * macOS: the CRAN build of R can use Apple's Accelerate framework; see the
#'   R for macOS FAQ for how to switch.
#' * Windows: see "R Installation and Administration" for replacing R's
#'   BLAS.
#'
#' Multithreaded BLAS libraries use all cores by default; when running
#' several R processes in parallel, limit the threads (for OpenBLAS, the
#' environment variable `OPENBLAS_NUM_THREADS`).
#'
#' The factorization depends only on the locations, the model and the trend,
#' not on the folds. To cross-validate several fold schemes for the same
#' data, call [kcv_prepare()] once and [kcv_run()] for each scheme.
#'
#' See [kcv_engine()] for memory requirements.
#' @export
#' @examples
#' library(sp)
#' data(meuse, package = "sp")
#' v <- gstat::vgm(0.59, "Sph", 897, 0.05)
#' cv <- kcv(log(zinc) ~ 1, ~ x + y, meuse, model = v)
#' head(cv)
#' sqrt(mean(cv$residual^2))
#'
#' # Buffered leave-one-out: neighbors within 300 m are held out too.
#' cv_buf <- kcv(log(zinc) ~ 1, ~ x + y, meuse, model = v, buffer = 300)
#' summary(cv_buf$n_heldout)
#' kcv_metrics(cv_buf)
#'
#' # The same statistics from a sample of 50 points.
#' set.seed(1)
#' kcv_metrics(kcv(log(zinc) ~ 1, ~ x + y, meuse, model = v, buffer = 300,
#'   targets = sample(nrow(meuse), 50)))
#'
#' # Several responses with the same model, trend and folds.
#' cv3 <- kcv(cbind(log(zinc), log(copper), log(lead)) ~ 1, ~ x + y, meuse,
#'   model = v, nfold = 5)
#' kcv_metrics(cv3)
kcv <- function(formula, locations, data, model, nfold, folds = NULL,
                beta = NULL, buffer = NULL, targets = NULL, ...) {
  inp <- kcv_inputs(formula, locations, data, model, beta, ...)
  # `p` is a promise: the factorization runs only after the fold and buffer
  # arguments have been checked.
  run_scheme(kcv_factor(inp$inp, inp$model, beta), inp$inp, nfold, folds,
    buffer, targets)
}

#' Factorize once, cross-validate many times
#'
#' `kcv_prepare()` does the expensive part of [kcv()] (building and
#' factorizing the covariance matrix, which depends only on the locations,
#' the variogram model and the trend) and returns it as an object.
#' `kcv_run()` then computes the cross-validation for one fold scheme from
#' it, at a fraction of the cost. `kcv(...)` is equivalent to
#' `kcv_run(kcv_prepare(...), ...)`.
#'
#' The prepared object holds an \eqn{n \times n} matrix (\eqn{8 n^2} bytes)
#' and can be saved with [saveRDS()].
#'
#' `kcv_run(prepared, y = )` cross-validates new responses at the same
#' locations, with the trend matrix and model stored in `prepared`, without
#' another factorization: for example, many fields simulated at the
#' observation locations. See "Several responses" in [kcv()].
#'
#' @inheritParams kcv
#' @param ... As in [kcv()].
#' @return `kcv_prepare()`: a `kcv_prepared` object. `kcv_run()`: the same
#'   output as [kcv()].
#' @seealso [kcv_engine_prepare()] for a precomputed covariance matrix; the
#'   "Performance" section of [kcv()].
#' @export
#' @examples
#' data(meuse, package = "sp")
#' v <- gstat::vgm(0.59, "Sph", 897, 0.05)
#' p <- kcv_prepare(log(zinc) ~ sqrt(dist), ~ x + y, meuse, model = v)
#' p
#' loo <- kcv_run(p)
#' set.seed(1)
#' rmse_10fold <- replicate(5, sqrt(mean(kcv_run(p, nfold = 10)$residual^2)))
#' rmse_10fold
#' # Buffered leave-one-out at several radii from the same factorization.
#' sapply(c(0, 200, 400), function(r) kcv_metrics(kcv_run(p, buffer = r))$rmse)
#'
#' # New responses at the same locations: here 20 random permutations of
#' # the observations, cross-validated from the same factorization.
#' set.seed(2)
#' perm <- replicate(20, sample(log(meuse$zinc)))
#' cv_perm <- kcv_run(p, nfold = 10, y = perm)
#' summary(kcv_metrics(cv_perm)$rmse)
kcv_prepare <-function(formula, locations, data, model, beta = NULL, ...) {
  inp <- kcv_inputs(formula, locations, data, model, beta, ...)
  kcv_factor(inp$inp, inp$model, beta)
}

#' @rdname kcv_prepare
#' @param prepared A `kcv_prepared` object from `kcv_prepare()`.
#' @param y Optional new responses, replacing those of the data given to
#'   `kcv_prepare()`: a numeric vector of length n (output as from a single
#'   response) or an n x T matrix (T responses, output with a `response`
#'   column), with rows in the order of that data. The locations, trend
#'   covariates and model are those of `prepared`. Missing or non-finite
#'   values are an error.
#' @param beta Optional, for a `prepared` object made for simple kriging
#'   (with `beta`): new trend coefficients, a vector of length p or a p x T
#'   matrix as in [kcv()]. The factorization of simple kriging does not
#'   depend on them.
#' @export
kcv_run <- function(prepared, nfold, folds = NULL, buffer = NULL,
                    targets = NULL, y = NULL, beta = NULL) {
  check_prepared(prepared)
  if (is.null(prepared$inp)) {
    stop(
      "`prepared` was made by kcv_engine_prepare() from a covariance matrix; ",
      "use kcv_engine_run() with it.",
      call. = FALSE
    )
  }
  inp <- prepared$inp
  if (!is.null(y)) {
    inp <- new_responses(inp, y)
  }
  prepared <- use_beta(prepared, beta, ncol(inp$y))
  run_scheme(prepared, inp, nfold, folds, buffer, targets)
}

# Replaces the responses of the parsed inputs `inp` by `y` (kcv_run(y = )):
# a vector of length n (one response) or an n x T matrix (T responses).
new_responses <- function(inp, y) {
  n <- nrow(inp$y)
  if (!is.numeric(y) || (!is.null(dim(y)) && length(dim(y)) != 2) ||
    NROW(y) != n) {
    stop(
      "`y` must be a numeric vector of length n = ", n, " or a matrix with ",
      "n rows (one per observation of the prepared data, in the same order).",
      call. = FALSE
    )
  }
  inp$multi <- is.matrix(y)
  if (inp$multi) {
    nm <- colnames(y)
    if (is.null(nm) || any(nm == "")) nm <- as.character(seq_len(ncol(y)))
    inp$y <- check_responses(y, nm, "`y`")
  } else {
    if (any(!is.finite(y))) {
      stop(
        "`y` has ", sum(!is.finite(y)), " missing or non-finite value(s); ",
        "every observation of the prepared data needs a response value.",
        call. = FALSE
      )
    }
    inp$y <- matrix(as.double(y), ncol = 1)
  }
  inp
}

# One fold scheme (or buffered leave-one-out) from a prepared factorization
# `p`, which is first used after the arguments are checked. `nfold` may be
# missing, as in kcv().
run_scheme <- function(p, inp, nfold, folds, buffer, targets) {
  n <- nrow(inp$y)
  if (!is.null(buffer) || !is.null(targets)) {
    targets <- check_buffer_args(buffer, missing(nfold), folds, targets, n)
    check_run_memory(n, ncol(inp$y), if (is.null(targets)) n else
      length(targets))
    return(build_output(inp, run_buffer_prepared(p, inp$y, buffer, targets)))
  }
  fold <- resolve_folds(nfold, folds, n)
  check_run_memory(n, ncol(inp$y), n)
  build_output(inp, run_prepared(p, inp$y, fold))
}

# Checks the buffered leave-one-out arguments; returns the target rows.
check_buffer_args <- function(buffer, nfold_missing, folds, targets, n) {
  if (is.null(buffer)) {
    stop("`targets` can only be used with `buffer` (buffered ",
      "leave-one-out); use `buffer = 0` for plain leave-one-out.",
      call. = FALSE
    )
  }
  check_buffer(buffer)
  if (!nfold_missing || !is.null(folds)) {
    stop("`buffer` (buffered leave-one-out) cannot be combined with ",
      "`nfold` or `folds`.",
      call. = FALSE
    )
  }
  resolve_targets(targets, n)
}

# Argument checks and input parsing shared by kcv() and kcv_prepare().
kcv_inputs <- function(formula, locations, data, model, beta, ...) {
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
  check_beta(beta, inp$F, ncol(inp$y))
  list(inp = inp, model = model)
}

# Builds the covariance directly in the factorization buffer (one n x n
# matrix) and factorizes it.
kcv_factor <- function(inp, model, beta) {
  n <- nrow(inp$y)
  check_memory(n, matrices = 1)
  F <- inp$F
  storage.mode(F) <- "double"
  res <- do.call(prepare_model_cpp, c(
    list(inp$coords), model_cpp_args(model),
    list(F = F, simple = !is.null(beta))
  ))
  new_prepared(res, F, beta, inp = inp, matrices = 1)
}

# Fold vector from kcv()'s `nfold` / `folds` arguments.
resolve_folds <- function(nfold, folds, n) {
  if (!is.null(folds)) {
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
      "observation in the neighborhood) with point support; local ",
      "neighborhoods and block kriging change the kriging system for each ",
      "prediction, so the single-factorization shortcut does not apply. ",
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

# Output in the format of krige.cv(). With a matrix response (`inp$multi`),
# one block per response, stacked in column order, with a factor column
# `response`.
build_output <- function(inp, res) {
  rows <- res$rows
  if (!is.null(rows)) {
    inp$y <- inp$y[rows, , drop = FALSE]
    inp$coords <- inp$coords[rows, , drop = FALSE]
    if (!is.null(inp$locations)) inp$locations <- inp$locations[rows, ]
  }
  if (isTRUE(inp$multi)) {
    return(build_output_multi(inp, res))
  }
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
  if (!is.null(res$n_heldout)) {
    out$n_heldout <- res$n_heldout
  }
  attach_locations(out, inp)
}

build_output_multi <- function(inp, res) {
  m <- nrow(inp$y)
  n_resp <- ncol(inp$y)
  names <- colnames(inp$y)
  observed <- as.vector(inp$y)
  pred <- as.vector(res$pred)
  var <- rep(res$var, n_resp)
  residual <- observed - pred
  out <- data.frame(
    var1.pred = pred,
    var1.var = var,
    observed = observed,
    residual = residual,
    zscore = residual / sqrt(var),
    fold = rep(res$folds, n_resp)
  )
  if (!is.null(res$n_heldout)) {
    out$n_heldout <- rep(res$n_heldout, n_resp)
  }
  out$response <- factor(rep(names, each = m), levels = names)
  attach_locations(out, inp, rep(seq_len(m), n_resp))
}

# Appends the locations of rows `idx` of the (subset) input to `out`; all
# rows, unchanged, when `idx` is NULL.
attach_locations <- function(out, inp, idx = NULL) {
  if (!is.null(idx)) {
    inp$coords <- inp$coords[idx, , drop = FALSE]
    if (!is.null(inp$locations)) inp$locations <- inp$locations[idx, ]
  }
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
