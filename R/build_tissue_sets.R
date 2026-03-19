#' Build the expressed-gene sets used internally by enrich_by_hpa
#'
#' Materialises the per-group gene sets that \code{\link{enrich_by_hpa}}
#' and \code{\link{enrich_by_hpa_ihc}} construct and test internally but
#' do not expose. The returned named list can be passed directly to
#' \code{\link{run_redundancy_selection}} to assess whether significant
#' tissue or cell-type signals are independent of one another.
#'
#' @section Why this closes the Track B redundancy loop:
#' \code{\link{enrich_by_hpa}} asks "is my gene list enriched in each
#' tissue's expressed genes?" \code{\link{run_redundancy_selection}} asks
#' "which of these enriched tissues explain non-redundant biology?" The gene
#' sets connecting both questions are exactly the expressed-gene sets per
#' tissue — what this function returns.
#'
#' Because \code{\link{run_redundancy_selection}} uses the same CMI engine
#' as the rest of Track A, tissue redundancy is handled by an existing
#' function, not a new one. \code{build_tissue_sets} is only the thin
#' adapter that makes Track B output compatible with the Track A redundancy
#' API.
#'
#' @section cutoff_type = "specificity":
#' When \code{cutoff_type = "specificity"}, gene sets are built using the
#' same cross-tissue fold-change definition as TissueEnrich: a gene is
#' included if its mean expression in that tissue is at least
#' \code{cutoff_value}-fold above the mean across all other tissues. This
#' produces small, precise gene sets (typically 1-100 genes per tissue)
#' suitable for identifying tissue-exclusive markers. The \code{reference_df}
#' must contain all tissues for the fold computation to be meaningful.
#'
#' Typical usage:
#' \preformatted{
#' # Match enrich_by_hpa settings exactly
#' hpa_res <- enrich_by_hpa(gene_list, hpa_data,
#'                           cutoff_type = "specificity", cutoff_value = 5)
#' tissue_sets <- build_tissue_sets(hpa_data,
#'                                   cutoff_type = "specificity", cutoff_value = 5)
#'
#' sig_tissues <- hpa_res$group[hpa_res$padj < 0.05]
#' tissue_sel  <- run_redundancy_selection(
#'   gene_list = gene_list,
#'   gene_sets = tissue_sets[sig_tissues],
#'   universe  = unique(hpa_data$Gene)
#' )
#' }
#'
#' @param reference_df Data frame of HPA data. Same object passed to
#'   \code{\link{enrich_by_hpa}} or \code{\link{enrich_by_hpa_ihc}}.
#' @param source Character. \code{"rna"} (default) or \code{"ihc"}.
#' @param group_col Character scalar. Column to group by. Auto-detected
#'   from the same candidate list as \code{\link{enrich_by_hpa}}.
#' @param expression_col Character scalar. Expression value column
#'   (RNA source only). Auto-detected from the same candidate list as
#'   \code{\link{enrich_by_hpa}}.
#' @param cutoff_type Character. \code{"quantile"} (default),
#'   \code{"threshold"}, \code{"top_k"}, or \code{"specificity"}.
#'   Must match the value used in \code{\link{enrich_by_hpa}}.
#'   See \code{\link{enrich_by_hpa}} for definitions.
#' @param cutoff_value Numeric. Cutoff value matching the value used in
#'   \code{\link{enrich_by_hpa}}. Default \code{0.9}. Use \code{5} for
#'   \code{"specificity"} (TissueEnrich default).
#' @param expression_levels Character vector. IHC levels treated as
#'   "expressed" (IHC source only). Default \code{c("Low", "Medium", "High")}.
#' @param remove_levels_of_reliability Character vector. Reliability
#'   labels to exclude (IHC source only). Default \code{"Uncertain"}.
#' @param min_genes Integer. Minimum number of genes a tissue set must
#'   contain to be included. Default \code{10L}.
#' @param universe Optional character vector. If supplied, gene sets are
#'   restricted to these genes before returning.
#'
#' @return A named list of character vectors, one entry per group that
#'   passed \code{min_genes}. Names match \code{group_col} values and are
#'   directly compatible with \code{\link{run_redundancy_selection}}.
#'
#' @seealso \code{\link{enrich_by_hpa}}, \code{\link{enrich_by_hpa_ihc}},
#'   \code{\link{run_redundancy_selection}}.
#'
#' @examples
#' \dontrun{
#' hpa_rna <- load_reference("rna_hpa", source = "hpa")[[1]]
#' hpa_rna$Gene <- hpa_rna[["Gene name"]]   # use HGNC symbols
#'
#' gene_list <- c("UMOD", "SLC12A1", "NPHS1", "PODXL")
#'
#' # --- Quantile (broad high-expression) ---
#' hpa_res    <- enrich_by_hpa(gene_list, hpa_rna,
#'                              cutoff_type = "quantile", cutoff_value = 0.9)
#' tissue_sets <- build_tissue_sets(hpa_rna,
#'                                   cutoff_type = "quantile", cutoff_value = 0.9)
#'
#' # --- Specificity (tissue-exclusive, matches TissueEnrich) ---
#' hpa_res_sp  <- enrich_by_hpa(gene_list, hpa_rna,
#'                               cutoff_type = "specificity", cutoff_value = 5)
#' tissue_sets_sp <- build_tissue_sets(hpa_rna,
#'                                      cutoff_type = "specificity",
#'                                      cutoff_value = 5)
#'
#' # Redundancy selection works the same way for both
#' sig <- hpa_res$group[hpa_res$padj < 0.05]
#' tissue_sel <- run_redundancy_selection(
#'   gene_list = gene_list,
#'   gene_sets = tissue_sets[sig],
#'   universe  = unique(hpa_rna$Gene)
#' )
#' plot_gains(tissue_sel)
#' }
#'
#' @export
build_tissue_sets <- function(
    reference_df,
    source                       = c("rna", "ihc"),
    group_col                    = c("Brain.region", "Subregion",
                                     "Immune.cell", "Tissue",
                                     "Cell.type", "IHC.tissue.name"),
    expression_col               = c("pTPM", "nTPM", "TPM",
                                     "Tags.per.million",
                                     "Scaled.tags.per.million",
                                     "Normalized.tags.per.million"),
    cutoff_type                  = c("quantile", "threshold", "top_k",
                                     "specificity"),
    cutoff_value                 = 0.9,
    expression_levels            = c("Low", "Medium", "High"),
    remove_levels_of_reliability = "Uncertain",
    min_genes                    = 10L,
    universe                     = NULL
) {
  ## --- Argument matching -----------------------------------------------------
  source      <- match.arg(source)
  cutoff_type <- match.arg(cutoff_type)

  if (!is.data.frame(reference_df))
    stop("`reference_df` must be a data frame.")

  ## --- Resolve group_col -----------------------------------------------------
  group_col <- if (length(group_col) > 1L) {
    .resolve_col(group_col, reference_df, "group_col")
  } else {
    group_col
  }
  if (!group_col %in% colnames(reference_df))
    stop("Column '", group_col, "' not found in `reference_df`.")
  if (!"Gene" %in% colnames(reference_df))
    stop("Column 'Gene' not found in `reference_df`.")

  ## --- Source-specific prep --------------------------------------------------
  if (source == "rna") {

    expression_col <- if (length(expression_col) > 1L) {
      .resolve_col(expression_col, reference_df, "expression_col")
    } else {
      expression_col
    }
    if (!expression_col %in% colnames(reference_df))
      stop("Column '", expression_col, "' not found in `reference_df`.")

    ## Validate cutoff
    if (cutoff_type == "quantile" && (cutoff_value <= 0 || cutoff_value >= 1))
      stop("For `cutoff_type = 'quantile'`, `cutoff_value` must be in (0, 1).")
    if (cutoff_type == "top_k" &&
        (cutoff_value %% 1 != 0 || cutoff_value < 1))
      stop("For `cutoff_type = 'top_k'`, `cutoff_value` must be a positive integer.")
    if (cutoff_type == "specificity" && cutoff_value <= 1)
      stop("For `cutoff_type = 'specificity'`, `cutoff_value` must be > 1 ",
           "(fold-change over mean of other tissues). ",
           "Typical values: 4, 5 (TissueEnrich default), 10.")

    ## Specificity: pre-compute cross-tissue folds before group split
    if (cutoff_type == "specificity") {
      n_groups <- length(unique(reference_df[[group_col]]))
      if (n_groups < 2L)
        stop("`cutoff_type = 'specificity'` requires >= 2 groups.")
      message("Computing cross-tissue specificity folds across ",
              n_groups, " groups...")
      specific_genes <- .compute_specificity_sets(
        reference_df, group_col, expression_col, cutoff_value
      )
    }

    group_list <- split(reference_df, reference_df[[group_col]])

    sets <- lapply(names(group_list), function(grp) {
      df <- group_list[[grp]]

      if (cutoff_type == "specificity") {
        return(specific_genes[[grp]])
      }

      ## Within-tissue cutoffs: compute mean expression per gene
      genes       <- unique(df$Gene)
      local_stats <- vapply(genes, function(g) {
        mean(df[[expression_col]][df$Gene == g], na.rm = TRUE)
      }, numeric(1L))

      switch(
        cutoff_type,
        quantile  = {
          q <- quantile(local_stats, cutoff_value, na.rm = TRUE)
          names(local_stats)[local_stats >= q]
        },
        threshold = names(local_stats)[local_stats > cutoff_value],
        top_k     = {
          k <- min(as.integer(cutoff_value), length(local_stats))
          names(sort(local_stats, decreasing = TRUE))[seq_len(k)]
        }
      )
    })
    names(sets) <- names(group_list)

  } else {
    ## IHC source ---------------------------------------------------------------
    if (!"Level" %in% colnames(reference_df))
      stop("Column 'Level' not found in `reference_df`. ",
           "Provide an HPA IHC data frame for source = 'ihc'.")
    if (cutoff_type == "specificity")
      stop("`cutoff_type = 'specificity'` is not supported for IHC data. ",
           "Use source = 'rna' for specificity-based sets.")

    if ("Reliability" %in% colnames(reference_df)) {
      reference_df <- reference_df[
        !is.na(reference_df$Reliability) &
          !reference_df$Reliability %in% remove_levels_of_reliability, ]
    }

    group_list <- split(reference_df, reference_df[[group_col]])
    sets <- lapply(group_list, function(df) {
      unique(df$Gene[df$Level %in% expression_levels])
    })
  }

  ## --- Apply universe and min_genes filter -----------------------------------
  if (!is.null(universe))
    sets <- lapply(sets, intersect, universe)

  n_genes <- lengths(sets)
  keep    <- n_genes >= min_genes
  n_drop  <- sum(!keep)
  if (n_drop > 0L)
    message(n_drop, " group(s) dropped: fewer than ", min_genes,
            " genes in expressed set.")

  sets[keep]
}
