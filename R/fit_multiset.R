.fit_multiset <- function(y, M, lambda = 0, weights = NULL, intercept = TRUE,
                          jitter = 1e-8, max_tries = 3) {
  G <- length(y)
  K <- ncol(M)

  if (!is.null(weights)) {
    sw <- sqrt(weights)
    y  <- sw * y
    M  <- Matrix::Diagonal(x = sw) %*% M
  }

  if (intercept) {
    y_mean <- mean(y)
    y <- y - y_mean
    M_means <- Matrix::colMeans(M)
    M <- sweep(M, 2, M_means)
  } else {
    y_mean <- 0
    M_means <- rep(0, K)
  }

  XtX <- Matrix::crossprod(M)        # K x K (sparse/dense Matrix)
  Xty <- Matrix::crossprod(M, y)     # K x 1

  # Ensure lambda is non-negative and allow vector penalty later
  lam <- lambda
  if (length(lam) == 1L) lam_vec <- rep(lam, K) else lam_vec <- lam
  if (length(lam_vec) != K) stop("lambda must be scalar or length ncol(M)")

  # Add penalty + jitter (jitter is crucial for near-singular cases)
  add_diag <- function(A, d) {
    # d numeric vector length K
    A + Matrix::Diagonal(x = d)
  }

  # Try solve with increasing jitter if needed
  cur_jitter <- jitter
  for (attempt in 0:max_tries) {
    A <- add_diag(XtX, lam_vec + cur_jitter)

    beta <- tryCatch(
      {
        # Prefer Cholesky for SPD matrices
        cholA <- Matrix::Cholesky(A, LDL = FALSE, Imult = 0)
        as.numeric(Matrix::solve(cholA, Xty))
      },
      error = function(e) {
        # Fallback generic solve
        tryCatch(as.numeric(Matrix::solve(A, Xty)), error = function(e2) NULL)
      }
    )

    if (!is.null(beta)) {
      intercept_hat <- if (intercept) y_mean - sum(M_means * beta) else 0
      return(list(beta = beta, intercept = intercept_hat, jitter_used = cur_jitter))
    }
    cur_jitter <- cur_jitter * 10
  }

  stop("Failed to solve ridge system even after jitter escalation. Consider increasing jitter/max_tries.")
}
