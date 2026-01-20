#' Run enrichment analysis using information theory. Which pathways explain the phenotype
#'
#' @param expr samples x proteins/genes
#' @param phenotype outcome vector
#' @param universe Character vector of all genes in the background set.
#' @param gene_sets list of characters with relevant pathways.
#' @param n_perm Number of permutations for null model.
#' @param seed Optional random seed for reproducibility.
#'
#' @return A dataframe with results
#' @export
run_info_assoc <- function(
    expr,
    gene_sets,
    phenotype,
    universe = NULL,
    score = c("mean_z", "pc1"),
    nbins = 10,
    n_perm = 1000,
    min_set_size = 10,
    seed = NULL
) {
  if (!is.null(seed)) set.seed(seed)

  score <- match.arg(score)

  # -----------------------------
  # Checks
  # -----------------------------
  if (!is.matrix(expr))
    stop("expr must be a numeric matrix (samples x genes)")

  if (length(phenotype) != nrow(expr))
    stop("phenotype length must equal number of samples")

  if (is.null(colnames(expr)))
    stop("expr must have gene names as colnames")

  if (is.null(universe))
    universe <- colnames(expr)

  universe <- intersect(universe, colnames(expr))
  expr <- expr[, universe, drop = FALSE]

  # intersect gene sets and filter
  gene_sets <- lapply(gene_sets, intersect, universe)
  gene_sets <- gene_sets[lengths(gene_sets) >= min_set_size]

  if (length(gene_sets) == 0)
    stop("No gene sets left after filtering by universe/min_set_size")

  # -----------------------------
  # Pathway scores (computed once)
  # -----------------------------
  Z <- scale(expr)

  score_fun <- switch(
    score,
    mean_z = function(G) rowMeans(Z[, G, drop = FALSE]),
    pc1 = function(G) {
      x <- Z[, G, drop = FALSE]
      prcomp(x, center = FALSE, scale. = FALSE)$x[, 1]
    }
  )

  scores <- lapply(gene_sets, score_fun)

  # -----------------------------
  # Discretize scores
  # -----------------------------
  library(infotheo)

  scores_disc <- lapply(
    scores,
    function(S)
      discretize(S, disc = "equalfreq", nbins = nbins)$X
  )

  y_fac <- as.factor(phenotype)

  # -----------------------------
  # Observed MI (effect size)
  # -----------------------------
  MI_obs <- vapply(
    scores_disc,
    function(Sd)
      mutinformation(Sd, y_fac, method = "emp"),
    numeric(1)
  )

  # -----------------------------
  # Permutation null (label-based)
  # -----------------------------
  MI_null <- matrix(
    NA_real_,
    nrow = length(scores_disc),
    ncol = n_perm
  )

  for (b in seq_len(n_perm)) {
    y_perm <- sample(y_fac)
    MI_null[, b] <- vapply(
      scores_disc,
      function(Sd)
        mutinformation(Sd, y_perm, method = "emp"),
      numeric(1)
    )
  }

  # -----------------------------
  # p-values + FDR
  # -----------------------------
  p_val <- (1 + rowSums(MI_null >= MI_obs)) / (n_perm + 1)
  padj  <- p.adjust(p_val, method = "BH")

  # -----------------------------
  # Output
  # -----------------------------
  data.frame(
    set = names(gene_sets),
    set_size = lengths(gene_sets),
    MI_bits = MI_obs,
    p_value = p_val,
    padj = padj,
    null_mean = rowMeans(MI_null),
    null_sd = apply(MI_null, 1, sd),
    row.names = NULL
  )
}
