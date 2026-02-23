.select_lambda_cv <- function(y, M, weights = NULL, grid = 10^seq(-4, 3, length.out = 80), Kfold = 10) {
  G <- length(y)
  folds <- sample(rep(seq_len(Kfold), length.out = G))
  mse <- numeric(length(grid))

  for (i in seq_along(grid)) {
    lam <- grid[i]
    err <- numeric(Kfold)

    for (f in seq_len(Kfold)) {
      te <- which(folds == f)
      tr <- setdiff(seq_len(G), te)

      fit <- .fit_multiset(
        y[tr], M[tr, , drop = FALSE],
        lambda = lam, weights = if (!is.null(weights)) weights[tr]
      )

      yhat <- fit$intercept + as.numeric(M[te, , drop = FALSE] %*% fit$beta)
      err[f] <- mean((y[te] - yhat)^2)
    }
    mse[i] <- mean(err)
  }
  grid[which.min(mse)]
}
