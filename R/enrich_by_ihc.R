#' Enrichment Analysis Using IHC Expression Data
#'
#' Perform hypergeometric enrichment of a gene list across cell types or tissues
#' using categorical immunohistochemistry (IHC) data from the Human Protein Atlas.
#'
#' @param gene_list Character vector of target genes (e.g., DEGs) to test for enrichment.
#' @param ihc_reference A data frame of IHC expression with columns for gene names,
#'   expression level (e.g., "Not detected", "Low", "Medium", "High"),
#'   and groupings such as cell type or tissue.
#' @param expression_levels Character vector indicating which IHC levels to treat as "expressed".
#'   Defaults to \code{c("Medium", "High")}.
#' @param group_col Column in \code{ihc_reference} used for grouping (e.g., "CellType" or "Tissue").
#' @param universe Optional character vector of background genes. If NULL, all genes
#'   found in \code{ihc_reference} are used.
#' @param seed seed for downstream gsea
#' @param ... Reserved for future options.
#'
#' @return A data frame with enrichment statistics per group:
#' \itemize{
#'   \item group (e.g. cell type)
#'   \item p_value
#'   \item overlap (genes shared with expressed set)
#'   \item input_set_size
#'   \item universe_size
#'   \item leading_edge (string of overlapping genes)
#' }
#'
#' @export
#'
#' @examples
#' # ihc_df <- readr::read_tsv("ihc_data.tsv")
#' # degs <- c("TP53", "BCL2", "CDKN1A")
#' # enrich_by_ihc(degs, ihc_df, expression_levels = c("High"))
enrich_by_ihc <- function(gene_list,
                          ihc_reference,
                          expression_levels = c("Medium", "High"),
                          group_col = "Tissue",
                          method = "hypergeometric",
                          universe = NULL,
                          gene_stats = NULL,
                          n_perm = NULL,
                          min_genes = 10,
                          seed=NULL,
                          alternative="greater",
                          expr_col = "Level") {

  # Validate
  if (!group_col %in% names(ihc_reference)) {
    stop("'", group_col, "' not found in IHC reference.")
  }

  gene_col <- if ("Gene.name" %in% names(ihc_reference)) "Gene.name" else "Gene"

  if (!expr_col %in% names(ihc_reference)) {
    stop("'Expression' column not found in IHC reference.")
  }

  remove_levels <- c("N/A", "Not representative","Ascending","Descending")

  ihc_reference <- ihc_reference %>%
    dplyr::filter(!Level %in% remove_levels)

  if (is.null(universe)) {
    universe <- unique(ihc_reference[[gene_col]])
  }

  results <- ihc_reference %>%
    dplyr::group_by(.data[[group_col]]) %>%
    dplyr::group_map(~{
      ref_group <- .x
      group_name <- unique(ref_group[[group_col]])

      # Expressed genes for this group
      expressed_genes <- ref_group %>%
        dplyr::filter(.data[[expr_col]] %in% expression_levels) %>%
        dplyr::pull(!!gene_col) %>%
        unique()

      if (length(expressed_genes) < min_genes) return(NULL)

      # GSEA logic
      if (method == "gsea") {
        if (is.null(n_perm)) stop("Number of permutations is not defined")

        local_gene_stats <- gene_stats
        if (is.null(local_gene_stats)) {
          ihc_levels <- c("Not detected" = 0, "Low" = 1, "Medium" = 2, "High" = 3)
          ref_group$expr_numeric <- ihc_levels[ref_group[[expr_col]]]
          local_gene_stats <- tapply(ref_group$expr_numeric, ref_group[[gene_col]], mean, na.rm = TRUE)
        } else {
          local_gene_stats <- local_gene_stats[names(local_gene_stats) %in% universe]
          if (length(local_gene_stats) == 0) return(NULL)
        }

        res <- run_enrichment(
          gene_list = gene_list,
          gene_stats = local_gene_stats,
          method = "gsea",
          n_perm = n_perm,
          seed = seed,
          alternative = alternative,
          universe = universe
        )
      } else {
        # Hypergeometric
        res <- run_enrichment(
          gene_list = gene_list,
          enriched_genes = expressed_genes,
          method = "hypergeometric",
          universe = universe
        )
      }

      if (!is.null(res)) {
        res[[group_col]] <- group_name
      }

      return(res)
    }, .keep = TRUE) %>%
    dplyr::bind_rows()

  return(results)
}
