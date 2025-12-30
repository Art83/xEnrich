#' Run enrichment analysis using information theory
#'
#' @param gene_list Character vector of user-supplied gene symbols.
#' @param universe Character vector of all genes in the background set.
#' @param gene_sets list of characters with relevant pathways.
#' @param n_perm Number of permutations for null model.
#' @param seed Optional random seed for reproducibility.
#'
#' @return A dataframe with results
#' @export
run_info_enrichment <- function(
    gene_list,
    gene_sets,
    universe = NULL,
    n_perm = 1000,
    seed = NULL
) {
  set.seed(seed)

  # universe handling
  if(is.null(universe)){
    stop("No universe detected")
  }

  gene_list <- intersect(gene_list, universe)
  gene_sets <- lapply(gene_sets, intersect, universe)

  N <- length(universe)
  n <- length(gene_list)

  results <- lapply(names(gene_sets), function(set_name) {

    G <- gene_sets[[set_name]]
    m <- length(G)
    k <- length(intersect(gene_list, G))

    obs_bits <- .ite_overlap_bits(N, n, m, k)

    expected <- n * m / N
    ratio <- ifelse(expected > 0, k / expected, NA_real_)

    # exact one-sided enrichment p-value
    p_enrich <- phyper(k - 1, m, N - m, n, lower.tail = FALSE)

    # optional: depletion p-value
    p_deplete <- phyper(k, m, N - m, n, lower.tail = TRUE)

    k_null <- rhyper(n_perm, m, N - m, n)

    bits_null <- vapply(
      k_null,
      function(x) .ite_overlap_bits(N, n, m, x),
      numeric(1)
    )

    data.frame(
      set = set_name,
      overlap = k,
      set_size = m,
      list_size = n,
      universe_size = N,
      expected = expected,
      ratio = ratio,
      info_bits = obs_bits,
      signed_bits = sign(k - expected) * obs_bits,
      p_enrich = p_enrich,
      p_deplete = p_deplete,
      null_mean = mean(bits_null),
      null_sd = sd(bits_null)
    )
  })

  res <- do.call(rbind, results)
  res$padj_enrich <- p.adjust(res$p_enrich, method = "BH")
  res$padj_deplete <- p.adjust(res$p_deplete, method = "BH")
  res
}




