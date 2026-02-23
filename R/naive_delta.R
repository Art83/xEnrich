.naive_delta <- function(y, M, weights = NULL) {
  if (!is.null(weights)) {
    w <- weights / sum(weights)
    yin  <- Matrix::t(M) %*% (w * y)
    nin  <- Matrix::t(M) %*% w
    yout <- sum(w * y) - yin
    nout <- 1 - nin
    as.numeric(yin / nin - yout / nout)
  } else {
    yin  <- Matrix::t(M) %*% y
    nin  <- Matrix::colSums(M)
    yout <- sum(y) - yin
    nout <- length(y) - nin
    as.numeric(yin / nin - yout / nout)
  }
}
