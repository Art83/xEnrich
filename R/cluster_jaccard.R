.cluster_jaccard_knn <- function(M, k = 5, threshold = 0.2) {
  inter <- as.matrix(Matrix::crossprod(M))
  size  <- diag(inter)
  union <- outer(size, size, "+") - inter
  J <- inter / union
  J[is.na(J)] <- 0
  diag(J) <- 0

  K <- ncol(M)
  adj <- matrix(FALSE, K, K)

  for (i in 1:K) {
    # take top-k neighbors by Jaccard, but only those above threshold
    ord <- order(J[i, ], decreasing = TRUE)
    ord <- ord[J[i, ord] >= threshold]
    if (length(ord) > k) ord <- ord[1:k]
    if (length(ord) > 0) adj[i, ord] <- TRUE
  }

  # symmetrize (undirected)
  adj <- adj | t(adj)

  # connected components
  comp <- rep(NA_integer_, K)
  cid <- 0
  for (i in 1:K) {
    if (is.na(comp[i])) {
      cid <- cid + 1
      stack <- i
      while (length(stack)) {
        v <- stack[1]; stack <- stack[-1]
        if (is.na(comp[v])) {
          comp[v] <- cid
          stack <- c(stack, which(adj[v, ]))
        }
      }
    }
  }
  comp
}
