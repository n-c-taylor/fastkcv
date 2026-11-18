#' Exact kriging cross-validation from a covariance matrix
#'
#' Low-level engine behind [kcv()]. Given the covariance matrix `C` of the
#' observations, the trend matrix `F` and fold assignments, it computes
#' cross-validation predictions and variances for every fold from a single
#' factorisation of `C`:
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
#' @param C n x n covariance matrix of the observations (positive definite).
#' @param F n x p trend matrix (a column of ones for ordinary kriging).
#' @param y Response: a numeric vector or an n x T matrix. Each column is
#'   cross-validated with the same kriging weights.
#' @param folds Vector of length n of fold identifiers.
#' @param beta Optional trend coefficients (length p) for simple kriging.
#' @return A list with `pred` (n x T matrix), `var` (length-n vector; the same
#'   for every column of `y`) and `folds`.
#' @export
#' @examples
#' set.seed(1)
#' xy <- matrix(runif(100), ncol = 2)
#' C <- exp(-as.matrix(dist(xy)) / 0.3)
#' y <- rnorm(50)
#' res <- kcv_engine(C, matrix(1, 50, 1), y, folds = seq_len(50))
#' head(res$pred)
kcv_engine <- function(C, F, y, folds, beta = NULL) {
  args <- check_engine_args(C, F, y, folds, beta)
  y <- args$y
  F <- args$F

  R <- tryCatch(chol(C), error = function(e) {
    stop(
      "The covariance matrix is not positive definite (Cholesky failed). ",
      "This usually means duplicate or near-duplicate locations, or a ",
      "Gaussian model without nugget; consider adding a small nugget.",
      call. = FALSE
    )
  })
  Ci <- chol2inv(R)

  if (is.null(beta)) {
    CiF <- Ci %*% F
    RA <- tryCatch(chol(crossprod(F, CiF)), error = function(e) {
      stop("The trend matrix F is rank deficient.", call. = FALSE)
    })
    G <- backsolve(RA, t(CiF), transpose = TRUE)
    Q <- Ci - crossprod(G)
  } else {
    Q <- Ci
    y_tilde <- y - drop(F %*% beta)
  }
  Qy <- Q %*% (if (is.null(beta)) y else y_tilde)

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
    pred[I, ] <- Fi %*% b + crossprod(W, y[r, , drop = FALSE] - Fr %*% b)
    V <- C[I, I, drop = FALSE] - crossprod(Cri, W)
    if (is.null(beta)) {
      D <- Fi - crossprod(W, Fr)
      V <- V + D %*% solve(A, t(D))
    }
    var[I] <- diag(V)
  }
  list(pred = pred, var = var, folds = folds)
}

# Inverse of Q[I, I] via its Cholesky factor. Q[I, I] is positive definite
# whenever the trend is identifiable without fold I (checked beforehand), so
# failure here signals numerical singularity.
fold_inverse <- function(Qii, f) {
  Rq <- tryCatch(chol(Qii), error = function(e) {
    stop(
      "Fold ", f, ": the matrix Q[I, I] is numerically singular, so this ",
      "fold cannot be predicted from the remaining data.",
      call. = FALSE
    )
  })
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
  if (!is.null(beta) && length(beta) != ncol(F)) {
    stop(
      "`beta` has length ", length(beta), " but the trend has ", ncol(F),
      " coefficient(s).",
      call. = FALSE
    )
  }
  list(y = y, F = F)
}

check_fold_identifiable <- function(F, idx, f, simple) {
  n_left <- nrow(F) - length(idx)
  if (n_left == 0) {
    stop("Fold ", f, " contains every observation; nothing is left to ",
      "predict from.",
      call. = FALSE
    )
  }
  if (!simple && qr(F[-idx, , drop = FALSE])$rank < ncol(F)) {
    stop(
      "Fold ", f, ": the trend is not identifiable from the data outside ",
      "this fold (the remaining rows of the trend matrix are rank ",
      "deficient). Use different folds or a simpler trend.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}
