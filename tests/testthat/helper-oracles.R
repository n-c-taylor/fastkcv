# Test oracles: independent implementations that the package is checked
# against. Sourced by testthat (in an environment whose parent is the
# package namespace, so internal helpers such as the err_*() functions
# resolve) and by inst/bench/common.R.

# Per-fold kriging refit through a Cholesky factorization
#
# The competent way to refit kriging for each fold: one Cholesky
# factorization of `C[r, r]` (r the rows outside the fold) and one
# triangular solve with all right-hand sides (the covariances with the
# fold, the trend columns and the responses) at once. About
# (n - m)^3 / 3 operations per fold of size m, so it beats the single
# factorization of kcv() when there are only a few large folds.
kcv_refit_chol <- function(C, F, y, folds, beta = NULL) {
  y <- as.matrix(y)
  F <- as.matrix(F)
  n <- nrow(y)
  n_resp <- ncol(y)
  p <- ncol(F)
  pred <- matrix(NA_real_, n, n_resp)
  var <- rep(NA_real_, n)
  for (f in sort(unique(folds))) {
    I <- which(folds == f)
    r <- which(folds != f)
    m <- length(I)
    yr <- y[r, , drop = FALSE]
    if (!is.null(beta)) yr <- yr - drop(F[r, , drop = FALSE] %*% beta)
    R <- chol(C[r, r, drop = FALSE])
    B <- cbind(C[r, I, drop = FALSE], F[r, , drop = FALSE], yr)
    X <- backsolve(R, backsolve(R, B, transpose = TRUE))
    K <- X[, seq_len(m), drop = FALSE] # C_rr^-1 C_rI
    CiF <- X[, m + seq_len(p), drop = FALSE]
    Ciy <- X[, m + p + seq_len(n_resp), drop = FALSE]
    Cri <- C[r, I, drop = FALSE]
    v <- diag(C)[I] - colSums(Cri * K)
    if (is.null(beta)) {
      A <- crossprod(F[r, , drop = FALSE], CiF)
      b <- solve(A, crossprod(F[r, , drop = FALSE], Ciy))
      pred[I, ] <- F[I, , drop = FALSE] %*% b + crossprod(Cri, Ciy - CiF %*% b)
      D <- F[I, , drop = FALSE] - crossprod(K, F[r, , drop = FALSE])
      v <- v + rowSums((D %*% solve(A)) * D)
    } else {
      pred[I, ] <- drop(F[I, , drop = FALSE] %*% beta) + crossprod(Cri, Ciy)
    }
    var[I] <- v
  }
  list(pred = pred, var = var, folds = folds)
}

# Exact kriging cross-validation in pure R (internal oracle)
#
# The same algebra as kcv_engine(), written with base R matrix functions
# and an explicit Q. Used to test the C++ engine.
#
# Arguments as in kcv_engine().
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
  if (!simple && !identifiable_r(F, I)) {
    raise_buffer_status(4, i, length(I), radius, n)
  }
  invisible(TRUE)
}

# Buffered leave-one-out in pure R (internal oracles)
#
# Point `i` is predicted from the observations outside `sets[[i]]` (which
# contains `i`). `kcv_buffer_r()` uses the explicit \eqn{Q} as
# `kcv_engine_r()` does; `kcv_buffer_reference()` refits kriging on the
# remaining data as `kcv_reference()` does. Used to test the C++ engine.
#
# Arguments as in kcv_engine().
# `sets`: List of hold-out index sets, e.g. from `buffer_sets_r()`.
# `radius`: Used only in error messages.
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

# @rdname kcv_buffer_r
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

# Naive kriging cross-validation (internal oracle)
#
# Refits simple or universal kriging on the data outside each fold with
# textbook formulas and `solve()`. Slow, used only to test [kcv_engine()].
#
# Arguments as in kcv_engine().
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
  if (!simple && !identifiable_r(F, idx)) {
    err_unidentifiable(f)
  }
  invisible(TRUE)
}

# The criterion of the C++ engine (identifiable() in src/engine.cpp), which
# is authoritative: with U an orthonormal basis of the columns of F, the
# trend is identifiable without rows `idx` when the smallest eigenvalue of
# U[-idx, ]'U[-idx, ] is at least 1e-12 (smallest singular value 1e-6).
identifiable_r <- function(F, idx) {
  U <- qr.Q(qr(F))
  ev <- eigen(crossprod(U[-idx, , drop = FALSE]), symmetric = TRUE,
    only.values = TRUE)$values
  min(ev) >= 1e-12
}

# Covariance matrix via gstat::variogramLine() (internal oracle)
#
# Builds `C[i, j] = C(s_i - s_j)` where `C(h) = total sill - gamma(h)`,
# evaluated with [gstat::variogramLine()] (`covariance = TRUE`) so that the
# model parameterization is exactly gstat's. Used to test
# `covariance_matrix()`. The nugget contributes only at
# distance exactly zero. Components with 2D geometric anisotropy are handled
# by mapping lag vectors to gstat's reduced (isotropic) distance.
#
# `model`: A `variogramModel`.
# `coords1`: n1 x 2 coordinate matrix.
# `coords2`: n2 x 2 coordinate matrix; defaults to `coords1`.
# Returns An n1 x n2 covariance matrix.
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

# Reduced distance under gstat's 2D geometric anisotropy
#
# gstat's `anis = c(p, s)`: `p` is the azimuth of the major axis (degrees,
# clockwise from north, i.e. from the positive y axis) and `s` the ratio of
# the minor to the major range. Lags along the minor axis are stretched by
# `1 / s` before evaluating the isotropic model.
#
reduced_distance <- function(dx, dy, ang1, anis1) {
  if (ang1 == 0 && anis1 == 1) {
    return(sqrt(dx^2 + dy^2))
  }
  a <- ang1 * pi / 180
  h_major <- dx * sin(a) + dy * cos(a)
  h_minor <- dx * cos(a) - dy * sin(a)
  sqrt(h_major^2 + (h_minor / anis1)^2)
}
