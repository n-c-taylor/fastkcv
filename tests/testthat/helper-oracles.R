# Test oracles: independent implementations that the package is checked
# against. Sourced by testthat (in an environment whose parent is the
# package namespace, so internal helpers such as the err_*() functions
# resolve) and by inst/bench/common.R.

#' Per-fold kriging refit through a Cholesky factorization
#'
#' The competent way to refit kriging for each fold: one Cholesky
#' factorization of `C[r, r]` (r the rows outside the fold) and one
#' triangular solve with all right-hand sides (the covariances with the
#' fold, the trend columns and the responses) at once. About
#' (n - m)^3 / 3 operations per fold of size m, so it beats the single
#' factorization of kcv() when there are only a few large folds.
#' @noRd
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
