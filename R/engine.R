#' Exact kriging cross-validation from a covariance matrix
#'
#' Low-level engine behind [kcv()], implemented in C++. Given the covariance
#' matrix `C` of the observations, the trend matrix `F` and fold
#' assignments, it computes cross-validation predictions and variances for
#' every fold from a single factorization of `C`:
#'
#' * ordinary/universal kriging (`beta = NULL`, trend re-estimated by GLS in
#'   each fold): \eqn{Q = C^{-1} - C^{-1}F(F'C^{-1}F)^{-1}F'C^{-1}};
#' * simple kriging (known `beta`): \eqn{Q = C^{-1}}, applied to
#'   \eqn{y - F\beta}.
#'
#' For a fold \eqn{I}, the residual is \eqn{e_I = Q_{II}^{-1}(Qy)_I}, the
#' prediction \eqn{y_I - e_I}, and the prediction-error variance
#' \eqn{\mathrm{diag}(Q_{II}^{-1})}.
#'
#' Only the diagonal blocks \eqn{Q_{II}} are formed, never the full
#' \eqn{n \times n} matrix \eqn{Q}. Leave-one-out (every fold of size one)
#' uses \eqn{e_i = (Qy)_i / Q_{ii}} and \eqn{1 / Q_{ii}} directly.
#'
#' @section Memory:
#' Factorizing `C` needs one \eqn{n \times n} workspace in addition to `C`,
#' so the peak memory of `kcv_engine()` is about \eqn{16 n^2} bytes. (The
#' model-based [kcv()] builds the covariance directly in the workspace and
#' needs half of that.) Because operating systems often over-commit memory,
#' an allocation that is too large may not fail cleanly, so the estimate is
#' checked up front against a budget of `getOption("fastkcv.max_gb", 6.4)`
#' GB, which allows \eqn{n} up to about 20,000 here and 28,000 in [kcv()].
#' Raise it with `options(fastkcv.max_gb = )` on machines with enough
#' memory. A prepared object keeps its \eqn{n \times n} factor (\eqn{8 n^2}
#' bytes) for as long as it exists.
#'
#' @param C n x n covariance matrix of the observations (symmetric positive
#'   definite; an asymmetric matrix or a non-finite entry is an error).
#' @param F n x p trend matrix (a column of ones for ordinary kriging).
#' @param y Response: a numeric vector or an n x T matrix. Each column is
#'   cross-validated with the same kriging weights.
#' @param folds Vector of length n of fold identifiers.
#' @param beta Optional trend coefficients for simple kriging: a length-p
#'   vector shared by every column of `y`, or a p x T matrix with one column
#'   per column of `y`. In `kcv_engine_run()`, replaces the coefficients
#'   given to `kcv_engine_prepare()`, which must have been made with `beta`
#'   (the factorization of simple kriging does not depend on `beta`).
#' @return `kcv_engine()` and `kcv_engine_run()` return a list with `pred`
#'   (n x T matrix), `var` (length-n vector; the same for every column of
#'   `y`) and `folds`. `kcv_engine_prepare()` returns a `kcv_prepared` object
#'   holding the factorization, to be passed to `kcv_engine_run()` any number
#'   of times with different `y` or `folds`.
#' @seealso [kcv_prepare()] and [kcv_run()] for the model-based equivalent.
#' @export
#' @examples
#' set.seed(1)
#' xy <- matrix(runif(100), ncol = 2)
#' C <- exp(-as.matrix(dist(xy)) / 0.3)
#' y <- rnorm(50)
#' res <- kcv_engine(C, matrix(1, 50, 1), y, folds = seq_len(50))
#' head(res$pred)
#'
#' # Factorize once, then cross-validate several fold schemes.
#' p <- kcv_engine_prepare(C, matrix(1, 50, 1))
#' loo <- kcv_engine_run(p, y, seq_len(50))
#' k5 <- kcv_engine_run(p, y, rep(1:5, 10))
kcv_engine <- function(C, F, y, folds, beta = NULL) {
  args <- check_engine_args(C, F, y, folds, beta)
  p <- kcv_engine_prepare(C, args$F, beta)
  kcv_engine_run(p, args$y, folds)
}

#' @rdname kcv_engine
#' @export
kcv_engine_prepare <- function(C, F, beta = NULL) {
  if (!is.matrix(C) || nrow(C) != ncol(C)) {
    stop("`C` must be a square n x n matrix.", call. = FALSE)
  }
  n <- nrow(C)
  F <- as.matrix(F)
  if (nrow(F) != n) {
    stop("`F` must have n = nrow(C) rows.", call. = FALSE)
  }
  check_trend(F)
  if (!is.numeric(F) || any(!is.finite(F))) {
    stop("`F` must be numeric, with no missing or non-finite values.",
      call. = FALSE
    )
  }
  check_beta(beta, F)
  if (is.null(beta)) check_trend_rank(F)
  check_memory(n, matrices = 2)
  # Assigning a storage mode copies C even when it is already double, which
  # would add a third n x n matrix to the peak.
  if (!is.double(C)) storage.mode(C) <- "double"
  storage.mode(F) <- "double"
  U <- trend_basis(F, beta)
  res <- prepare_matrix_cpp(C, if (is.null(beta)) U else F, !is.null(beta))
  new_prepared(res, F, U, beta, inp = NULL, matrices = 2)
}

#' @rdname kcv_engine
#' @param prepared A `kcv_prepared` object from `kcv_engine_prepare()` or
#'   [kcv_prepare()].
#' @export
kcv_engine_run <- function(prepared, y, folds, beta = NULL) {
  check_prepared(prepared)
  n <- length(prepared$qdiag)
  if (!is.matrix(y)) {
    y <- matrix(y, ncol = 1)
  }
  if (nrow(y) != n || !is.numeric(y)) {
    stop("`y` must be a numeric vector of length n or a matrix with n = ", n,
      " rows.",
      call. = FALSE
    )
  }
  if (any(!is.finite(y))) {
    stop("`y` contains missing or non-finite values.", call. = FALSE)
  }
  if (length(folds) != n || anyNA(folds)) {
    stop("`folds` must have length n = ", n, " and no missing values.",
      call. = FALSE
    )
  }
  prepared <- use_beta(prepared, beta, ncol(y))
  # The engine returns only the n x T predictions (and one variance vector).
  check_run_memory(n, ncol(y), n, out_cols = 1)
  run_prepared(prepared, y, folds)
}

#' @export
print.kcv_prepared <- function(x, ...) {
  n <- length(x$qdiag)
  type <- if (!is.null(x$beta)) {
    "simple kriging"
  } else if (ncol(x$F) == 1 && all(x$F == 1)) {
    "ordinary kriging"
  } else {
    "universal kriging"
  }
  cat(
    "<kcv_prepared> ", type, ", n = ", n, ", p = ", ncol(x$F), ", ",
    format(structure(8 * as.double(n)^2, class = "object_size"),
      units = "auto"
    ),
    if (is.null(x$inp)) ", from a covariance matrix" else "", "\n",
    sep = ""
  )
  invisible(x)
}

# Orthonormal basis U of the column space of F for ordinary/universal
# kriging (empty for simple kriging). Q depends on F only through its column
# space, so the factorization uses U in place of F: F'C^-1 F can be
# numerically singular for a full-rank F whose columns differ in scale by
# many orders of magnitude (a quadratic trend in UTM coordinates: F'C^-1 F
# had rcond 1e-35), while cond(U'C^-1 U) <= cond(C). U also serves the
# per-fold identifiability checks.
trend_basis <- function(F, beta) {
  if (is.null(beta)) qr.Q(qr(F)) else matrix(0, 0, 0)
}

# Wraps the C++ factorization into a kcv_prepared object. `F` is kept as
# given (simple kriging uses F beta, in the parametrization of `beta`); B, S
# and G of ordinary/universal kriging are in the basis `U`.
new_prepared <- function(res, F, U, beta, inp, matrices) {
  if (res$status != 0) {
    raise_status(res$status, NA, nrow(F), matrices)
  }
  simple <- !is.null(beta)
  structure(
    list(
      Wt = res$Wt, B = res$B, S = res$S, G = res$G, U = U,
      qdiag = res$qdiag, F = F,
      beta = if (simple) store_beta(beta) else NULL,
      inp = inp
    ),
    class = "kcv_prepared"
  )
}

# A vector beta is kept as a double vector (as before); a p x T matrix as a
# double matrix without dimnames.
store_beta <- function(beta) {
  if (is.matrix(beta)) {
    matrix(as.double(beta), nrow(beta), ncol(beta))
  } else {
    as.double(beta)
  }
}

# beta for the C++ engine: p x 1 (shared) or p x T; empty for OK/UK.
beta_cpp <- function(p) {
  if (is.null(p$beta)) matrix(0, 0, 0) else as.matrix(p$beta)
}

# Returns `p` with `beta` replacing its simple-kriging coefficients (when
# given), after checking it against the `n_resp` responses.
use_beta <- function(p, beta, n_resp) {
  if (!is.null(beta)) {
    if (is.null(p$beta)) {
      stop(
        "`beta` can only be given for a prepared object made for simple ",
        "kriging (with `beta`); ordinary and universal kriging estimate ",
        "the trend in every fold.",
        call. = FALSE
      )
    }
    check_beta(beta, p$F, n_resp)
    p$beta <- store_beta(beta)
  } else if (is.matrix(p$beta) && ncol(p$beta) > 1 &&
    ncol(p$beta) != n_resp) {
    stop(
      "The prepared object's `beta` has ", ncol(p$beta), " columns (one ",
      "per response) but there are ", n_resp, " response(s); pass a ",
      "matching `beta`.",
      call. = FALSE
    )
  }
  p
}

# Up-front memory check of the run stage for several responses: the C++
# work matrices (y, Wt y, Q y, predictions: about 4 n doubles per response)
# and the output (`out_cols` columns of n_out doubles per response: about
# 10 for the long kcv() output, 1 for the prediction matrix of
# kcv_engine_run()) must fit in the budget left after the n x n factor. One
# response is never checked, so single-response behavior is unchanged.
check_run_memory <- function(n, n_resp, n_out, out_cols = 10) {
  if (n_resp <= 1) {
    return(invisible(TRUE))
  }
  max_gb <- getOption("fastkcv.max_gb", 6.4)
  per_col <- 8 * (4 * as.double(n) + out_cols * as.double(n_out))
  avail <- max_gb * 1e9 - 8 * as.double(n)^2
  if (per_col * n_resp > avail) {
    err_responses_memory(n, n_resp, per_col * n_resp / 1e9, max_gb,
      max(0, floor(avail / per_col)))
  }
  invisible(TRUE)
}

err_responses_memory <- function(n, n_resp, gb, max_gb, batch) {
  stop(
    n_resp, " responses at n = ", n, " need about ",
    format(signif(gb, 2)), " GB for the cross-validation work matrices ",
    "and output, more than the memory budget of ", format(max_gb),
    " GB leaves after the n x n factor. ",
    if (batch >= 1) {
      paste0("Run the responses in batches of at most ", batch,
        " columns (kcv_run(prepared, y = )), or raise ")
    } else {
      "Raise "
    },
    "the budget with options(fastkcv.max_gb = ).",
    call. = FALSE
  )
}

# Buffered leave-one-out from a prepared factorization: each target point i
# (row numbers, all rows by default) is predicted from the observations
# farther than `radius` from it. `rows` is set only when `targets` is given,
# so that build_output() subsets the input only then. `method` ("auto",
# "pairs" or "points") selects how the blocks Q[I, I] are formed; see
# run_buffer_cpp().
run_buffer_prepared <- function(p, y, radius, targets = NULL,
                                method = c("auto", "pairs", "points")) {
  method <- match.arg(method)
  storage.mode(y) <- "double"
  n <- nrow(y)
  rows <- targets
  if (is.null(targets)) targets <- seq_len(n)
  simple <- !is.null(p$beta)
  max_gb <- getOption("fastkcv.max_gb", 6.4)
  res <- run_buffer_cpp(
    p$Wt, p$B, p$S, p$G, p$U, p$qdiag, p$F, y, p$inp$coords,
    as.double(radius), as.integer(targets) - 1L, simple, beta_cpp(p),
    match(method, c("auto", "pairs", "points")) - 1L,
    max(0, max_gb * 1e9 - 8 * as.double(n)^2)
  )
  if (res$status != 0) {
    raise_buffer_status(res$status, res$fold, res$size, radius, n)
  }
  pred <- res$pred
  dimnames(pred) <- list(rownames(y)[targets], colnames(y))
  list(pred = pred, var = res$var, folds = targets,
    n_heldout = res$size, rows = rows, method = res$method)
}

# Resolves `targets` to sorted row numbers.
resolve_targets <- function(targets, n) {
  if (is.null(targets)) {
    return(NULL)
  }
  if (is.logical(targets)) {
    if (length(targets) != n || anyNA(targets)) {
      stop("A logical `targets` must have length n = ", n, " and no ",
        "missing values.",
        call. = FALSE
      )
    }
    targets <- which(targets)
  } else if (!is.numeric(targets) || anyNA(targets) ||
    any(targets != round(targets)) || any(targets < 1 | targets > n)) {
    stop("`targets` must be row numbers between 1 and n = ", n,
      " or a logical vector of length n.",
      call. = FALSE
    )
  }
  if (anyDuplicated(targets)) {
    stop("`targets` contains duplicate row numbers.", call. = FALSE)
  }
  if (length(targets) == 0) {
    stop("`targets` selects no observations.", call. = FALSE)
  }
  sort(as.integer(targets))
}

check_buffer <- function(buffer) {
  if (!is.numeric(buffer) || length(buffer) != 1 || !is.finite(buffer) ||
    buffer < 0) {
    stop("`buffer` must be a single non-negative number (a radius in ",
      "coordinate units).",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

# Cross-validates one fold scheme from a prepared factorization.
run_prepared <- function(p, y, folds) {
  storage.mode(y) <- "double"
  simple <- !is.null(p$beta)
  labels <- sort(unique(folds))
  # One pass, in `labels` order (`which(folds == f)` per fold is O(n^2) for
  # leave-one-out and left about one n x n matrix of garbage).
  fold_idx <- unname(split(seq_along(folds) - 1L, match(folds, labels)))
  res <- run_cpp(
    p$Wt, p$B, p$S, p$G, p$U, p$qdiag, p$F, y, fold_idx, simple,
    beta_cpp(p)
  )
  if (res$status != 0) {
    f <- if (res$fold > 0) labels[res$fold] else NA
    raise_status(res$status, f, nrow(y), 0)
  }
  pred <- res$pred
  dimnames(pred) <- dimnames(y)
  list(pred = pred, var = res$var, folds = folds)
}

raise_status <- function(status, f, n, matrices) {
  switch(status,
    err_not_pd(),
    err_trend_singular(),
    err_every_obs(f),
    err_unidentifiable(f),
    err_singular(f),
    if (matrices == 0) {
      err_run_memory("more, smaller folds")
    } else {
      err_memory(n, matrices)
    },
    err_not_symmetric(),
    err_not_finite_c()
  )
}

check_prepared <- function(x) {
  if (!inherits(x, "kcv_prepared")) {
    stop("`prepared` must be a kcv_prepared object from kcv_prepare() or ",
      "kcv_engine_prepare().",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

# `beta` is a length-p vector shared by every response or a p x T matrix
# (column t for response t); `n_resp` is T when known.
check_beta <- function(beta, F, n_resp = NULL) {
  if (is.null(beta)) {
    return(invisible(TRUE))
  }
  if (is.matrix(beta)) {
    if (nrow(beta) != ncol(F)) {
      stop(
        "`beta` has ", nrow(beta), " row(s) but the trend has ", ncol(F),
        " coefficient(s); a matrix `beta` must be p x T, one column per ",
        "response.",
        call. = FALSE
      )
    }
    if (!is.null(n_resp) && ncol(beta) != 1 && ncol(beta) != n_resp) {
      stop(
        "`beta` has ", ncol(beta), " columns but there are ", n_resp,
        " response(s); give one column per response, or a vector shared ",
        "by all.",
        call. = FALSE
      )
    }
  } else if (length(beta) != ncol(F)) {
    stop(
      "`beta` has length ", length(beta), " but the trend has ", ncol(F),
      " coefficient(s).",
      call. = FALSE
    )
  }
  if (!is.numeric(beta) || any(!is.finite(beta))) {
    stop("`beta` must be numeric, with no missing or non-finite values.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

# A trend matrix without columns (formula z ~ 0) is neither ordinary nor
# universal kriging; gstat rejects it too.
check_trend <- function(F) {
  if (ncol(F) == 0) {
    stop("The trend has no terms (e.g. formula z ~ 0); use z ~ 1 for ",
      "ordinary kriging, or a known mean with `beta`.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

# OK/UK need a full-rank trend matrix. Checked with a pivoted QR (R's
# default tolerance 1e-7, as in check_fold_identifiable()) before the
# factorization. The C++ check of F'C^-1 F (status 2) cannot replace it:
# with F itself, whether its Cholesky failed for a rank-deficient F depended
# on the BLAS's rounding of Wt F (MKL and ATLAS accepted F = cbind(1, 1)),
# and it now sees only the orthonormal basis of trend_basis().
# `uses_coords`: the trend formula uses the coordinate names (kcv() and
# kcv_prepare() only). A polynomial in raw coordinates with a large offset
# (e.g. UTM) is the common way to get here with a full-rank but
# ill-conditioned F (unit square + 1e3, full quadratic: condition number of
# the column-scaled F 1e8), so the message then suggests centering.
check_trend_rank <- function(F, uses_coords = FALSE) {
  q <- qr(F)
  p <- ncol(F)
  if (q$rank < p) {
    dropped <- q$pivot[(q$rank + 1):p]
    nm <- colnames(F)
    cols <- if (is.null(nm)) dropped else paste0("`", nm[dropped], "`")
    norms <- sqrt(colSums(F^2))
    norms[norms == 0] <- 1
    cond <- kappa(sweep(F, 2, norms, "/"), exact = TRUE)
    err_rank_deficient(cols, cond, uses_coords)
  }
  invisible(TRUE)
}

#' Exact kriging cross-validation in pure R (internal oracle)
#'
#' The same algebra as [kcv_engine()], written with base R matrix functions
#' and an explicit \eqn{Q}. Used to test the C++ engine.
#'
#' @inheritParams kcv_engine
#' @noRd
kcv_engine_r <- function(C, F, y, folds, beta = NULL) {
  args <- check_engine_args(C, F, y, folds, beta)
  y <- args$y
  F <- args$F

  qp <- q_parts_r(C, F, y, beta)
  Q <- qp$Q
  Qy <- qp$Qy

  n <- nrow(y)
  pred <- matrix(NA_real_, n, ncol(y), dimnames = dimnames(y))
  var <- rep(NA_real_, n)
  for (f in sort(unique(folds))) {
    idx <- which(folds == f)
    check_fold_identifiable(F, idx, f, simple = !is.null(beta))
    V <- fold_inverse(Q[idx, idx, drop = FALSE], f)
    pred[idx, ] <- y[idx, , drop = FALSE] - V %*% Qy[idx, , drop = FALSE]
    var[idx] <- diag(V)
  }
  list(pred = pred, var = var, folds = folds)
}

# Explicit Q and Q y (Q (y - F beta) for simple kriging) for the R oracles.
q_parts_r <- function(C, F, y, beta) {
  R <- tryCatch(chol(C), error = function(e) err_not_pd())
  Ci <- chol2inv(R)
  if (is.null(beta)) {
    CiF <- Ci %*% F
    RA <- tryCatch(chol(crossprod(F, CiF)), error = function(e) {
      check_trend_rank(F)
      err_trend_singular()
    })
    G <- backsolve(RA, t(CiF), transpose = TRUE)
    Q <- Ci - crossprod(G)
    Qy <- Q %*% y
  } else {
    Q <- Ci
    Qy <- Q %*% (y - drop(F %*% beta))
  }
  list(Q = Q, Qy = Qy)
}

# Hold-out sets of buffered leave-one-out: element i holds the (sorted)
# indices of the points within Euclidean distance `radius` of point i,
# including i. Pure R, independent of the C++ neighbor search.
buffer_sets_r <- function(coords, radius) {
  D <- as.matrix(stats::dist(coords))
  lapply(seq_len(nrow(D)), function(i) which(D[i, ] <= radius))
}

# Checks shared by the buffered-LOO oracles, with the same statuses and
# messages as the C++ engine.
check_buffer_set_r <- function(F, I, i, radius, simple) {
  n <- nrow(F)
  if (length(I) == n) {
    raise_buffer_status(3, i, length(I), radius, n)
  }
  if (!simple && qr(F[-I, , drop = FALSE])$rank < ncol(F)) {
    raise_buffer_status(4, i, length(I), radius, n)
  }
  invisible(TRUE)
}

#' Buffered leave-one-out in pure R (internal oracles)
#'
#' Point `i` is predicted from the observations outside `sets[[i]]` (which
#' contains `i`). `kcv_buffer_r()` uses the explicit \eqn{Q} as
#' `kcv_engine_r()` does; `kcv_buffer_reference()` refits kriging on the
#' remaining data as `kcv_reference()` does. Used to test the C++ engine.
#'
#' @inheritParams kcv_engine
#' @param sets List of hold-out index sets, e.g. from `buffer_sets_r()`.
#' @param radius Used only in error messages.
#' @noRd
kcv_buffer_r <- function(C, F, y, sets, beta = NULL, radius = NA) {
  args <- check_engine_args(C, F, y, seq_along(sets), beta)
  y <- args$y
  F <- args$F
  qp <- q_parts_r(C, F, y, beta)
  n <- nrow(y)
  pred <- matrix(NA_real_, n, ncol(y), dimnames = dimnames(y))
  var <- rep(NA_real_, n)
  for (i in seq_len(n)) {
    I <- sets[[i]]
    check_buffer_set_r(F, I, i, radius, simple = !is.null(beta))
    Rq <- tryCatch(chol(qp$Q[I, I, drop = FALSE]), error = function(e) {
      raise_buffer_status(5, i, length(I), radius, n)
    })
    V <- chol2inv(Rq)
    k <- match(i, I)
    pred[i, ] <- y[i, ] - V[k, ] %*% qp$Qy[I, , drop = FALSE]
    var[i] <- V[k, k]
  }
  list(pred = pred, var = var, folds = seq_len(n),
    n_heldout = lengths(sets))
}

#' @rdname kcv_buffer_r
#' @noRd
kcv_buffer_reference <- function(C, F, y, sets, beta = NULL, radius = NA) {
  args <- check_engine_args(C, F, y, seq_along(sets), beta)
  y <- args$y
  F <- args$F
  n <- nrow(y)
  pred <- matrix(NA_real_, n, ncol(y))
  var <- rep(NA_real_, n)
  for (i in seq_len(n)) {
    I <- sets[[i]]
    check_buffer_set_r(F, I, i, radius, simple = !is.null(beta))
    r <- setdiff(seq_len(n), I)
    k <- krige_subset_r(C, F, y, r, i, beta)
    pred[i, ] <- k$pred
    var[i] <- k$var
  }
  list(pred = pred, var = var, folds = seq_len(n),
    n_heldout = lengths(sets))
}

# Kriging prediction at rows `I` from rows `r`, with textbook formulas and
# solve(): the body of kcv_reference() for one fold.
krige_subset_r <- function(C, F, y, r, I, beta) {
  Crr <- C[r, r, drop = FALSE]
  Cri <- C[r, I, drop = FALSE]
  Fr <- F[r, , drop = FALSE]
  Fi <- F[I, , drop = FALSE]
  W <- solve(Crr, Cri) # simple kriging weights, |r| x |I|
  if (is.null(beta)) {
    CrrF <- solve(Crr, Fr)
    A <- crossprod(Fr, CrrF)
    b <- solve(A, crossprod(CrrF, y[r, , drop = FALSE]))
  } else {
    b <- matrix(beta, ncol(F), ncol(y))
  }
  pred <- Fi %*% b + crossprod(W, y[r, , drop = FALSE] - Fr %*% b)
  V <- C[I, I, drop = FALSE] - crossprod(Cri, W)
  if (is.null(beta)) {
    D <- Fi - crossprod(W, Fr)
    V <- V + D %*% solve(A, t(D))
  }
  list(pred = pred, var = diag(V))
}

#' Naive kriging cross-validation (internal oracle)
#'
#' Refits simple or universal kriging on the data outside each fold with
#' textbook formulas and `solve()`. Slow, used only to test [kcv_engine()].
#'
#' @inheritParams kcv_engine
#' @noRd
kcv_reference <- function(C, F, y, folds, beta = NULL) {
  args <- check_engine_args(C, F, y, folds, beta)
  y <- args$y
  F <- args$F
  n <- nrow(y)
  pred <- matrix(NA_real_, n, ncol(y))
  var <- rep(NA_real_, n)
  for (f in sort(unique(folds))) {
    I <- which(folds == f)
    r <- setdiff(seq_len(n), I)
    check_fold_identifiable(F, I, f, simple = !is.null(beta))
    k <- krige_subset_r(C, F, y, r, I, beta)
    pred[I, ] <- k$pred
    var[I] <- k$var
  }
  list(pred = pred, var = var, folds = folds)
}

# Inverse of Q[I, I] via its Cholesky factor. Q[I, I] is positive definite
# whenever the trend is identifiable without fold I (checked beforehand), so
# failure here signals numerical singularity.
fold_inverse <- function(Qii, f) {
  Rq <- tryCatch(chol(Qii), error = function(e) err_singular(f))
  chol2inv(Rq)
}

check_engine_args <- function(C, F, y, folds, beta) {
  if (!is.matrix(y)) {
    y <- matrix(y, ncol = 1)
  }
  n <- nrow(y)
  if (!is.matrix(C) || nrow(C) != n || ncol(C) != n) {
    stop("`C` must be an n x n matrix with n = nrow(y).", call. = FALSE)
  }
  F <- as.matrix(F)
  if (nrow(F) != n) {
    stop("`F` must have n = nrow(y) rows.", call. = FALSE)
  }
  if (length(folds) != n || anyNA(folds)) {
    stop("`folds` must have length n = nrow(y) and no missing values.",
      call. = FALSE
    )
  }
  check_beta(beta, F, ncol(y))
  list(y = y, F = F)
}

check_fold_identifiable <- function(F, idx, f, simple) {
  n_left <- nrow(F) - length(idx)
  if (n_left == 0) {
    err_every_obs(f)
  }
  if (!simple && qr(F[-idx, , drop = FALSE])$rank < ncol(F)) {
    err_unidentifiable(f)
  }
  invisible(TRUE)
}

# Guard against inputs whose dense n x n matrices cannot reasonably be held
# in memory: `matrices` n x n double matrices must fit in the budget
# getOption("fastkcv.max_gb"). See the "Memory" section of ?kcv_engine.
check_memory <- function(n, matrices) {
  # Armadillo's 32-bit indices (arma::uword) cap a matrix at 2^32 - 1
  # elements.
  if (n > 65535) {
    stop("n = ", n, " observations exceed the largest supported n = 65535 ",
      "(the n x n matrix would have more than 2^32 - 1 elements).",
      call. = FALSE
    )
  }
  max_gb <- getOption("fastkcv.max_gb", 6.4)
  if (matrices * 8 * as.double(n)^2 / 1e9 > max_gb) {
    err_memory(n, matrices, max_gb)
  }
  invisible(TRUE)
}

# Error conditions shared by the C++ and R engines.
err_not_pd <- function() {
  stop(
    "The covariance matrix is not positive definite (Cholesky failed). ",
    "This usually means duplicate or near-duplicate locations, or a ",
    "Gaussian model without nugget; consider adding a small nugget.",
    call. = FALSE
  )
}

err_rank_deficient <- function(cols, cond, uses_coords) {
  stop(
    "The trend matrix F is rank deficient or nearly so: ",
    if (length(cols) == 1) {
      paste0("column ", cols, " is (nearly) a linear combination")
    } else {
      paste0(
        "columns ", paste(cols, collapse = ", "),
        " are (nearly) linear combinations"
      )
    },
    " of the others (condition number of F with unit-length columns: ",
    format(signif(cond, 2)), "). ",
    if (uses_coords) {
      paste0(
        "The trend uses the coordinates: with large coordinate values ",
        "(e.g. UTM), powers and products of them are nearly collinear. ",
        "Center the coordinates first (e.g. x - mean(x), y - mean(y)); for ",
        "a full polynomial trend this does not change the predictions."
      )
    } else {
      "Drop or combine those columns."
    },
    call. = FALSE
  )
}

# C++ status 2: the Cholesky factor of U'C^-1 U (U an orthonormal basis of
# the trend) failed or has reciprocal condition number below machine
# precision. Since cond(U'C^-1 U) <= cond(C), this means C is too
# ill-conditioned, not that the trend is badly scaled.
err_trend_singular <- function() {
  stop(
    "The generalized least squares matrix F'C^-1 F of the trend is ",
    "numerically singular, so the trend cannot be estimated: the ",
    "covariance matrix is too ill-conditioned. Consider adding a nugget or ",
    "using a simpler trend.",
    call. = FALSE
  )
}

err_every_obs <- function(f, label = NULL) {
  if (!is.null(label)) {
    stop(label, ": every observation is held out; nothing is left to ",
      "predict from.",
      call. = FALSE
    )
  }
  stop("Fold ", f, " contains every observation; nothing is left to ",
    "predict from.",
    call. = FALSE
  )
}

err_unidentifiable <- function(f, label = paste("Fold", f),
                               region = "this fold",
                               advice = "different folds") {
  stop(
    label, ": the trend is not identifiable from the data outside ",
    region, " (the remaining rows of the trend matrix are rank ",
    "deficient). Use ", advice, " or a simpler trend.",
    call. = FALSE
  )
}

err_singular <- function(f, label = paste("Fold", f), region = "fold") {
  stop(
    label, ": the matrix Q[I, I] is numerically singular, so this ",
    region, " cannot be predicted from the remaining data.",
    call. = FALSE
  )
}

# Errors of buffered leave-one-out, naming the point (row number) and its
# hold-out set.
raise_buffer_status <- function(status, i, m, radius, n) {
  label <- paste0(
    "Point ", i, " (the buffer of radius ", format(radius), " holds out ",
    m, " of ", n, " observations)"
  )
  # Status 1 belongs to the prepare stage and never occurs here; 2 only for a
  # modified prepared object whose S fails the solve in compute_Qy().
  unexpected <- function() {
    stop("Unexpected status ", status, ".", call. = FALSE)
  }
  switch(status,
    unexpected(),
    err_trend_singular(),
    err_every_obs(i, label),
    err_unidentifiable(i, label,
      region = "the buffer", advice = "a smaller radius"
    ),
    err_singular(i, label, region = "point"),
    err_run_memory("`targets` or a smaller `buffer`")
  )
}

err_not_finite_c <- function() {
  stop("`C` contains missing or non-finite values.", call. = FALSE)
}

err_not_symmetric <- function() {
  stop("`C` is not symmetric (C[i, j] and C[j, i] differ by more than ",
    "1e-8 times the largest diagonal entry).",
    call. = FALSE
  )
}

# Out of memory after the factorization, while forming the blocks Q[I, I].
err_run_memory <- function(advice) {
  stop("Out of memory while forming the blocks Q[I, I] of the ",
    "cross-validation (after the factorization). Use ", advice, ".",
    call. = FALSE
  )
}

err_memory <- function(n, matrices, max_gb = NULL) {
  gb <- matrices * 8 * as.double(n)^2 / 1e9
  # Rounded up to two significant digits, so the suggested budget suffices.
  step <- 10^(floor(log10(gb)) - 1)
  need <- ceiling(gb / step) * step
  stop(
    "n = ", n, " observations need ", matrices, " dense n x n matri",
    if (matrices == 1) "x" else "ces", " (about ", format(need), " GB)",
    if (!is.null(max_gb)) {
      paste0(
        ", above the memory budget of ", format(max_gb), " GB. If this ",
        "machine has enough memory, raise it with options(fastkcv.max_gb = ",
        format(need), ")"
      )
    } else {
      ", which could not be allocated"
    },
    ".",
    call. = FALSE
  )
}
