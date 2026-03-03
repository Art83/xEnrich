#' Conditional enrichment via ridge regression
#'
#' Estimates pathway-level effect sizes by fitting a ridge regression of a
#' gene-level response \code{y} (e.g. log fold-change, z-score) on a binary
#' membership matrix \code{M}. Unlike marginal enrichment, ridge regression
#' accounts for pathway overlap: the coefficient for pathway \eqn{k}
#' reflects its unique contribution after conditioning on all other pathways.
#'
#' The \code{drop} column (\code{delta_naive - beta}) quantifies how much of
#' a pathway's apparent marginal signal is explained by co-membership with
#' other pathways — a direct measure of redundancy.
#'
#' @section Clustering:
#' When \code{cluster = TRUE}, pathways are grouped into clusters by Jaccard
#' overlap using a kNN graph (see \code{cluster_k} and
#' \code{cluster_threshold}). The \code{axis} table summarises the total
#' ridge weight per cluster, identifying the dominant pathway axes.
#'
#' @section Lambda selection:
#' \describe{
#'   \item{\code{"cv"}}{10-fold cross-validation on prediction MSE.
#'     Computationally expensive but data-driven.}
#'   \item{\code{"fixed"}}{Heuristic: \code{0.1 * mean(colSums(M))}.
#'     Fast; reasonable default when CV is too slow.}
#' }
#'
#' @param y Named numeric vector of per-gene scores (e.g. log-FC, z-score).
#'   Names must match \code{rownames(M)} or \code{universe}.
#' @param gene_sets Named list of character vectors. Required if \code{M} is
#'   not supplied directly.
#' @param universe Character vector of background genes. Defaults to
#'   \code{names(y)}.
#' @param M Optional precomputed membership matrix (genes x pathways) as a
#'   \code{Matrix::sparseMatrix}. Overrides \code{gene_sets}.
#' @param subset Optional character vector of pathway names to retain from
#'   \code{M}.
#' @param method One of \code{"ridge"} (default) or \code{"ols"}.
#'   \code{"ols"} sets \code{lambda = 0}; clustering is disabled.
#' @param lambda Numeric. Penalty parameter. If \code{NULL} (default),
#'   selected automatically according to \code{lambda_select}.
#' @param lambda_select One of \code{"cv"} (default) or \code{"fixed"}.
#'   Ignored when \code{lambda} is supplied explicitly or \code{method = "ols"}.
#' @param weights Optional named numeric vector of per-gene weights.
#' @param min_size Integer. Minimum pathway size after universe intersection.
#'   Default \code{5}.
#' @param max_size Integer. Maximum pathway size. Default \code{Inf}.
#' @param cluster Logical. Cluster pathways by Jaccard overlap?
#'   Default \code{TRUE}. Ignored when \code{method = "ols"}.
#' @param cluster_k Integer. Maximum kNN neighbours in the Jaccard graph.
#'   Default \code{5}.
#' @param cluster_threshold Numeric in \[0, 1\]. Minimum Jaccard similarity
#'   to form an edge. If \code{NULL}, auto-selected as the 95th percentile
#'   of pairwise similarities, floored at \code{0.08}.
#'
#' @return An object of class \code{xenrich_multi}, a list with:
#'   \describe{
#'     \item{\code{table}}{Data frame with one row per pathway:
#'       \code{pathway}, \code{size}, \code{delta_naive} (marginal mean
#'       difference), \code{beta} (ridge coefficient),
#'       \code{beta_std} (standardised beta), \code{drop} (redundancy
#'       measure), and \code{cluster} (integer cluster ID, if clustering
#'       was performed).}
#'     \item{\code{axis}}{Data frame of cluster-level summed ridge weights
#'       (\code{axis_strength}), or \code{NULL} if all clusters are
#'       singletons.}
#'     \item{\code{lambda}}{Lambda value used.}
#'     \item{\code{method}}{Method used.}
#'     \item{\code{M}}{The membership matrix used in fitting.}
#'   }
#'
#' @seealso [run_info_enrichment()], [run_info_assoc()], [run_enrichment()]
#'
#' @examples
#' set.seed(1)
#' genes     <- paste0("G", 1:300)
#' y         <- setNames(rnorm(300), genes)
#' gene_sets <- list(
#'   path_A = genes[1:50],
#'   path_B = genes[30:80],   # overlaps path_A
#'   path_C = genes[150:200]
#' )
#'
#' res <- run_conditional_enrichment(
#'   y         = y,
#'   gene_sets = gene_sets,
#'   universe  = genes,
#'   seed      = 42
#' )
#' res$table
#' res$axis
#'
#' @export
run_conditional_enrichment <- function(
    y,
    gene_sets         = NULL,
    universe          = names(y),
    M                 = NULL,
    subset            = NULL,
    method            = c("ridge", "ols"),
    lambda            = NULL,
    lambda_select     = c("cv", "fixed"),
    weights           = NULL,
    min_size          = 5L,
    max_size          = Inf,
    cluster           = TRUE,
    cluster_k         = 5L,
    cluster_threshold = NULL,
    seed              = NULL
) {
  # --- Validation -------------------------------------------------------------
  if (is.null(M) && is.null(gene_sets))
    stop("Supply either `gene_sets` or a precomputed membership matrix `M`.")
  if (!is.numeric(y))
    stop("`y` must be a named numeric vector.")

  method        <- match.arg(method)
  lambda_select <- match.arg(lambda_select)
  if (!is.null(seed)) set.seed(seed)

  # --- Build / validate M -----------------------------------------------------
  if (is.null(M)) {
    M <- .make_membership_matrix(gene_sets, universe)
  }

  if (!is.null(subset)) M <- M[, intersect(subset, colnames(M)), drop = FALSE]

  # Size filtering
  set_sizes <- Matrix::colSums(M)
  keep      <- set_sizes >= min_size & set_sizes <= max_size
  M         <- M[, keep, drop = FALSE]
  n_dropped <- sum(!keep)
  if (n_dropped > 0L) message(n_dropped, " pathway(s) dropped by size filter.")
  if (ncol(M) == 0L)  stop("No pathways remain after size filtering.")

  # Align y and weights to M rows
  y <- y[rownames(M)]
  if (!is.null(weights)) weights <- weights[rownames(M)]

  # --- Marginal delta ---------------------------------------------------------
  delta <- .naive_delta(y, M, weights)

  # --- Lambda -----------------------------------------------------------------
  if (method == "ols") {
    lambda <- 0
  } else if (is.null(lambda)) {
    lambda <- if (lambda_select == "cv") {
      message("Selecting lambda by ", Kfold <- 10, "-fold CV...")
      .select_lambda_cv(y, M, weights)
    } else {
      0.1 * mean(Matrix::colSums(M))
    }
    message("Lambda = ", signif(lambda, 3))
  }

  # --- Fit --------------------------------------------------------------------
  fit  <- .fit_multiset(y, M, lambda, weights)
  sy   <- stats::sd(y)
  pk   <- as.numeric(Matrix::colMeans(M))
  sdMk <- sqrt(pk * (1 - pk))

  res <- data.frame(
    pathway     = colnames(M),
    size        = as.integer(Matrix::colSums(M)),
    delta_naive = delta,
    beta        = fit$beta,
    beta_std    = fit$beta * sdMk / sy,
    drop        = delta - fit$beta,
    stringsAsFactors = FALSE
  )

  # --- Clustering (ridge only) ------------------------------------------------
  axis <- NULL
  if (cluster && method == "ridge" && ncol(M) >= 2L) {
    J   <- .jaccard_matrix(M)                        # computed once, passed down
    thr <- if (is.null(cluster_threshold)) {
      max(as.numeric(stats::quantile(J[upper.tri(J)], 0.95)), 0.08)
    } else {
      cluster_threshold
    }

    cl         <- .cluster_jaccard_knn(J, k = cluster_k, threshold = thr)
    res$cluster <- cl

    cs <- table(cl)
    if (max(cs) >= 2L) {
      axis <- stats::aggregate(beta ~ cluster, data = res, FUN = sum)
      names(axis)[2] <- "axis_strength"
    }
  }

  structure(
    list(
      table  = res,
      axis   = axis,
      lambda = lambda,
      method = method,
      M      = M
    ),
    class = "xenrich_multi"
  )
}

#' Union-Find (disjoint set) connected components
#'
#' O(K * alpha(K)) — effectively linear. Used internally by
#' [run_conditional_enrichment()] to cluster overlapping pathways.
#'
#' @param adj Logical symmetric adjacency matrix (K x K).
#' @return Integer vector of component IDs, length K.
#' @keywords internal
#' @noRd
.connected_components <- function(adj) {
  K      <- nrow(adj)
  parent <- seq_len(K)

  find <- function(x) {
    while (parent[x] != x) {
      parent[x] <<- parent[parent[x]]   # path compression
      x <- parent[x]
    }
    x
  }

  union <- function(a, b) {
    ra <- find(a); rb <- find(b)
    if (ra != rb) parent[ra] <<- rb
  }

  for (i in seq_len(K)) {
    neighbours <- which(adj[i, ])
    for (j in neighbours[neighbours > i]) union(i, j)
  }

  # Canonicalise labels to 1..n_components
  roots <- vapply(seq_len(K), find, integer(1L))
  as.integer(factor(roots))
}


#' Sparse-aware column centring for a Matrix object
#'
#' Subtracting column means from a sparse 0/1 matrix would densify it.
#' Instead we store the means separately and apply them only when needed
#' (i.e. when forming X'X and X'y), keeping M sparse throughout.
#'
#' Returns a list with \code{M} (unchanged) and \code{col_means}.
#' Callers use the identity:
#'   (M - 1*mu') ' (M - 1*mu') = M'M  - n * mu mu'
#'
#' @keywords internal
#' @noRd
.center_sparse <- function(M) {
  col_means <- Matrix::colMeans(M)
  list(M = M, col_means = col_means)
}


#' Compute Jaccard similarity between columns of a membership matrix
#'
#' @param M Sparse membership matrix (genes x pathways).
#' @return Dense numeric matrix (pathways x pathways).
#' @keywords internal
#' @noRd
.jaccard_matrix <- function(M) {
  inter <- as.matrix(Matrix::crossprod(M))
  size  <- diag(inter)
  union <- outer(size, size, "+") - inter
  J     <- inter / union
  J[is.na(J)] <- 0
  diag(J)      <- 0
  J
}


#' Build a binary membership matrix from a named list of gene sets
#'
#' @param sets Named list of character vectors.
#' @param universe Character vector of background genes.
#'
#' @return A sparse \code{Matrix::sparseMatrix} (genes x pathways).
#' @keywords internal
#' @noRd
.make_membership_matrix <- function(sets, universe) {
  if (!is.list(sets) || is.null(names(sets)))
    stop("`sets` must be a named list.")

  universe <- unique(universe)
  G <- length(universe)
  K <- length(sets)

  idx <- lapply(sets, function(s) match(intersect(s, universe), universe))
  i   <- unlist(idx, use.names = FALSE)
  j   <- rep(seq_len(K), lengths(idx))

  Matrix::sparseMatrix(
    i        = i,
    j        = j,
    x        = 1L,
    dims     = c(G, K),
    dimnames = list(universe, names(sets))
  )
}


#' Naive (marginal) mean difference per pathway
#'
#' Computes the simple in-vs-out mean difference for each pathway,
#' optionally weighted. Used as a baseline for comparison with ridge betas.
#'
#' @keywords internal
#' @noRd
.naive_delta <- function(y, M, weights = NULL) {
  if (!is.null(weights)) {
    w    <- weights / sum(weights)
    yin  <- as.numeric(Matrix::t(M) %*% (w * y))
    nin  <- as.numeric(Matrix::t(M) %*% w)
    yout <- sum(w * y) - yin
    nout <- 1 - nin
  } else {
    yin  <- as.numeric(Matrix::t(M) %*% y)
    nin  <- Matrix::colSums(M)
    yout <- sum(y) - yin
    nout <- length(y) - nin
  }
  as.numeric(yin / nin - yout / nout)
}


#' Cross-validated lambda selection for ridge pathway regression
#'
#' @keywords internal
#' @noRd
.select_lambda_cv <- function(y, M, weights = NULL,
                              grid   = 10^seq(-4, 3, length.out = 80),
                              Kfold  = 10) {
  G     <- length(y)
  folds <- sample(rep(seq_len(Kfold), length.out = G))
  mse   <- numeric(length(grid))

  for (i in seq_along(grid)) {
    err <- numeric(Kfold)
    for (f in seq_len(Kfold)) {
      te  <- which(folds == f)
      tr  <- setdiff(seq_len(G), te)
      fit <- .fit_multiset(
        y[tr], M[tr, , drop = FALSE],
        lambda  = grid[i],
        weights = if (!is.null(weights)) weights[tr] else NULL
      )
      yhat   <- fit$intercept + as.numeric(M[te, , drop = FALSE] %*% fit$beta)
      err[f] <- mean((y[te] - yhat)^2)
    }
    mse[i] <- mean(err)
  }
  grid[which.min(mse)]
}


#' Fit ridge (or OLS) regression with a sparse membership matrix
#'
#' Uses the centred normal equations. Column-centring is handled analytically
#' to keep M sparse: X'X and X'y are adjusted using stored column means
#' rather than materialising a dense centred matrix.
#'
#' @keywords internal
#' @noRd
.fit_multiset <- function(y, M, lambda = 0, weights = NULL,
                          intercept = TRUE,
                          jitter = 1e-8, max_tries = 3) {
  K <- ncol(M)
  n <- length(y)

  if (!is.null(weights)) {
    sw <- sqrt(weights)
    y  <- sw * y
    M  <- Matrix::Diagonal(x = sw) %*% M
  }

  # Centre y and analytically centre M (avoids densification)
  if (intercept) {
    y_mean    <- mean(y)
    y_c       <- y - y_mean
    col_means <- Matrix::colMeans(M)
    # (M - 1*mu')' (M - 1*mu') = M'M - n * mu mu'
    XtX <- Matrix::crossprod(M) - n * Matrix::tcrossprod(col_means)
    Xty <- as.numeric(Matrix::crossprod(M, y_c)) -
      as.numeric(col_means * sum(y_c))
  } else {
    y_mean    <- 0
    y_c       <- y
    col_means <- rep(0, K)
    XtX       <- Matrix::crossprod(M)
    Xty       <- as.numeric(Matrix::crossprod(M, y_c))
  }

  lam_vec <- if (length(lambda) == 1L) rep(lambda, K) else lambda
  if (length(lam_vec) != K)
    stop("`lambda` must be a scalar or a vector of length ncol(M).")

  add_diag <- function(A, d) A + Matrix::Diagonal(x = d)

  cur_jitter <- jitter
  for (attempt in seq(0, max_tries)) {
    A    <- add_diag(XtX, lam_vec + cur_jitter)
    beta <- tryCatch({
      ch <- Matrix::Cholesky(A, LDL = FALSE, Imult = 0)
      as.numeric(Matrix::solve(ch, Xty))
    }, error = function(e)
      tryCatch(as.numeric(Matrix::solve(A, Xty)),
               error = function(e2) NULL)
    )
    if (!is.null(beta)) {
      intercept_hat <- if (intercept) y_mean - sum(col_means * beta) else 0
      return(list(beta = beta, intercept = intercept_hat,
                  jitter_used = cur_jitter))
    }
    cur_jitter <- cur_jitter * 10
  }
  stop("Ridge system failed after jitter escalation. ",
       "Try increasing `jitter` or `max_tries`.")
}




.make_membership_matrix <- function(sets, universe) {
  stopifnot(is.list(sets))
  if (is.null(names(sets))) stop("sets must be a named list")

  universe <- unique(universe)
  G <- length(universe)
  K <- length(sets)

  idx <- lapply(sets, function(s) match(intersect(s, universe), universe))
  i <- unlist(idx)
  j <- rep(seq_len(K), lengths(idx))

  Matrix::sparseMatrix(
    i = i,
    j = j,
    x = 1L,
    dims = c(G, K),
    dimnames = list(universe, names(sets))
  )
}


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
