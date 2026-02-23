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
