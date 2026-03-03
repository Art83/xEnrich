#' Group-wise RNA Expression Enrichment Analysis
#'
#' Applies enrichment testing (GSEA or hypergeometric) across multiple groups
#' (e.g., tissues, brain regions, subregions or immune cell types) using Human Protein
#' Atlas (HPA) RNA expression data.
#'
#' @param gene_list A character vector of gene symbols (the "test set").
#' @param reference_df A data frame containing HPA expression data. Must include
#'   a \code{Gene} column, the column specified in \code{group_col}, and the
#'   column specified in \code{expression_col}.
#' @param group_col Character string specifying the column to group by
#'   (e.g., \code{"Tissue"}, \code{"Brain.region"}). Must be length 1.
#'   Defaults to the first match among \code{c("Brain.region", "Subregion",
#'   "Immune.cell", "Tissue")} found in \code{reference_df}.
#' @param expression_col Character string specifying the expression value column
#'   to use. Must be length 1. Defaults to the first match among
#'   \code{c("pTPM", "nTPM", "TPM", "Tags.per.million",
#'   "Scaled.tags.per.million", "Normalized.tags.per.million")} found in
#'   \code{reference_df}.
#' @param method Character string, either \code{"hypergeometric"} or
#'   \code{"gsea"}.
#' @param cutoff_type Character string specifying how to discretize expression
#'   for the hypergeometric test. One of:
#'   \itemize{
#'     \item \code{"threshold"}: genes with expression strictly above
#'       \code{cutoff_value}.
#'     \item \code{"quantile"}: genes at or above the \code{cutoff_value}
#'       quantile (must be between 0 and 1).
#'     \item \code{"top_k"}: the top \code{cutoff_value} genes by expression
#'       (must be a positive integer).
#'   }
#'   Ignored when \code{method = "gsea"}.
#' @param cutoff_value Numeric value controlling the cutoff. Interpretation
#'   depends on \code{cutoff_type} (see above). Default \code{0.9}.
#' @param universe Optional character vector defining the gene universe.
#'   Only used when \code{universe_type = "provided"}.
#' @param universe_type Character string defining how the background gene
#'   universe is constructed:
#'   \itemize{
#'     \item \code{"global"}: all unique genes in \code{reference_df}.
#'     \item \code{"group"}: only genes measured within the specific group
#'       being tested.
#'     \item \code{"provided"}: uses the \code{universe} argument directly
#'       (intersected with genes present in \code{reference_df}).
#'   }
#' @param min_overlap Minimum number of query genes (\code{gene_list}) that
#'   must overlap the group universe for the group to be tested. Groups below
#'   this threshold emit a warning but still run; set higher to skip
#'   low-overlap groups. Default \code{1}.
#' @param min_background Minimum size of the background gene set required to
#'   attempt enrichment. For \code{"hypergeometric"}, this is the number of
#'   genes in \code{local_stats} after universe filtering; for \code{"gsea"},
#'   the length of the ranking vector. Groups below this threshold are silently
#'   skipped. Default \code{10}.
#' @param gene_stats Optional named numeric vector of pre-calculated gene-level
#'   statistics. Required when \code{method = "gsea"} and
#'   \code{universe_type = "provided"}; ignored otherwise.
#' @param seed Optional integer seed for reproducibility of GSEA permutations.
#' @param alternative Character string specifying the alternative hypothesis
#'   for the enrichment test. One of \code{"greater"}, \code{"less"}, or
#'   \code{"two.sided"}.
#' @param ... Additional named arguments passed to \code{run_enrichment} for
#'   GSEA (e.g., \code{n_perm}, \code{gsea_weight}, \code{adaptive},
#'   \code{adaptive_mode}, \code{alpha}, \code{eps}).
#'
#' @return A data frame of enrichment results, one row per group that passed
#'   filtering, with a column added for the group name.
#'
#' @export
enrich_by_hpa <- function(gene_list,
                          reference_df,
                          group_col      = c("Brain.region", "Subregion", "Immune.cell", "Tissue"),
                          expression_col = c("pTPM", "nTPM", "TPM", "Tags.per.million",
                                             "Scaled.tags.per.million", "Normalized.tags.per.million"),
                          method         = c("hypergeometric", "gsea"),
                          cutoff_type    = c("threshold", "quantile", "top_k"),
                          cutoff_value   = 0.9,
                          universe       = NULL,
                          universe_type  = c("global", "group", "provided"),
                          min_overlap    = 1,
                          min_background = 10,
                          gene_stats     = NULL,
                          seed           = NULL,
                          alternative    = c("greater", "less", "two.sided"),
                          ...) {

  # --- Argument standardisation ---
  stopifnot(is.data.frame(reference_df))
  method        <- match.arg(method)
  universe_type <- match.arg(universe_type)
  alternative   <- match.arg(alternative)
  cutoff_type   <- match.arg(cutoff_type)

  # match.arg picks the first valid element from vector defaults,
  # or validates a user-supplied scalar against the allowed set
  group_col      <- match.arg(group_col,      choices = c("Brain.region", "Subregion", "Immune.cell", "Tissue"))
  expression_col <- match.arg(expression_col, choices = c("pTPM", "nTPM", "TPM", "Tags.per.million",
                                                          "Scaled.tags.per.million", "Normalized.tags.per.million"))

  # --- Column presence checks ---
  if (!group_col %in% colnames(reference_df))
    stop(sprintf("Grouping variable '%s' is absent in reference_df.", group_col))
  if (!expression_col %in% colnames(reference_df))
    stop(sprintf("Expression column '%s' is absent in reference_df.", expression_col))
  if (!'Gene' %in% colnames(reference_df))
    stop("'Gene' column is absent in reference_df.")

  # --- Scalar argument checks ---
  if (!is.null(universe) && (!is.character(universe) || length(universe) == 0))
    stop("'universe' must be a non-empty character vector.")
  if (universe_type == "provided" && is.null(universe))
    stop("universe_type = 'provided' but universe is NULL.")
  if (!is.numeric(min_overlap)    || length(min_overlap)    != 1 || min_overlap    < 1)
    stop("'min_overlap' must be a single integer >= 1.")
  if (!is.numeric(min_background) || length(min_background) != 1 || min_background < 1)
    stop("'min_background' must be a single integer >= 1.")

  # --- Cutoff checks ---
  if (cutoff_type == "quantile" && (cutoff_value <= 0 || cutoff_value >= 1))
    stop("For cutoff_type = 'quantile', cutoff_value must be strictly between 0 and 1.")
  if (cutoff_type == "top_k" && (!is.numeric(cutoff_value) || cutoff_value %% 1 != 0 || cutoff_value < 1))
    stop("For cutoff_type = 'top_k', cutoff_value must be a positive integer.")

  # --- gene_stats relevance warning ---
  if (!is.null(gene_stats) && !(method == "gsea" && universe_type == "provided"))
    warning("`gene_stats` is only used when method='gsea' and universe_type='provided'. It will be ignored.")

  # --- GSEA dots validation ---
  dots <- list(...)
  if (method == "gsea" && length(dots) > 0) {
    if (is.null(names(dots)) || any(names(dots) == ""))
      stop("All '...' arguments must be named.")
    allowed <- c("n_perm", "gsea_weight", "adaptive", "adaptive_mode", "alpha", "eps")
    bad <- setdiff(names(dots), allowed)
    if (length(bad))
      stop("Unsupported GSEA arguments: ", paste(bad, collapse = ", "))
  }

  # --- gene_stats validation for GSEA + provided ---
  if (method == "gsea" && universe_type == "provided") {
    if (is.null(gene_stats))
      stop("method='gsea' and universe_type='provided' require a `gene_stats` named numeric vector.")
    if (!is.numeric(gene_stats) || is.null(names(gene_stats)))
      stop("`gene_stats` must be a named numeric vector.")
  }

  # --- Universe construction ---
  global_universe <- if (universe_type == "global") {
    unique(reference_df$Gene)
  } else if (universe_type == "provided") {
    original_len    <- length(universe)
    trimmed         <- intersect(universe, unique(reference_df$Gene))
    dropped         <- original_len - length(trimmed)
    if (dropped > 0)
      message(dropped, " of ", original_len,
              " provided universe genes are absent from reference_df and will be excluded.")
    if (length(trimmed) == 0)
      stop("Provided universe has no overlap with reference_df$Gene.")
    trimmed
  } else {
    NULL  # universe_type == "group": determined per-group below
  }

  # --- Per-group enrichment ---
  group_list <- split(reference_df, reference_df[[group_col]])

  results_list <- lapply(names(group_list), function(group_name) {

    ref_group        <- group_list[[group_name]]
    current_universe <- if (universe_type == "group") unique(ref_group$Gene) else global_universe
    ref_group        <- ref_group[ref_group$Gene %in% current_universe, ]

    # Compute expression-derived stats (used by hypergeometric always,
    # and by gsea when universe_type != "provided")
    local_stats <- tapply(ref_group[[expression_col]], ref_group$Gene, mean, na.rm = TRUE)
    local_stats <- setNames(as.vector(local_stats), names(local_stats))
    local_stats <- local_stats[names(local_stats) %in% current_universe]
    current_universe <- intersect(current_universe, names(local_stats))

    if (length(local_stats) < min_background) return(NULL)
    if (length(current_universe) == 0)        return(NULL)

    # Override with user stats for gsea + provided
    if (method == "gsea" && universe_type == "provided") {
      local_stats      <- gene_stats[names(gene_stats) %in% current_universe]
      current_universe <- intersect(current_universe, names(local_stats))
      if (length(local_stats) < min_background) return(NULL)
      if (length(current_universe) == 0)        return(NULL)
    }

    # Overlap check
    gene_list_local <- intersect(gene_list, current_universe)
    n_overlap       <- length(gene_list_local)
    if (n_overlap == 0) return(NULL)
    if (n_overlap < min_overlap)
      message("Group '", group_name, "': only ", n_overlap, " of ", length(gene_list),
              " query genes found in universe. Results may be unreliable.")

    # Enrichment
    if (method == "gsea") {
      res <- run_enrichment(
        gene_list  = gene_list_local,
        gene_stats = local_stats,
        method     = "gsea",
        alternative = alternative,
        universe   = current_universe,
        seed       = seed,
        ...
      )
    } else {
      enriched_genes <- switch(cutoff_type,
                               "threshold" = names(local_stats)[local_stats > cutoff_value],
                               "quantile"  = {
                                 q_thresh <- quantile(local_stats, cutoff_value, na.rm = TRUE)
                                 names(local_stats)[local_stats >= q_thresh]
                               },
                               "top_k" = {
                                 names(sort(local_stats, decreasing = TRUE))[seq_len(min(cutoff_value, length(local_stats)))]
                               }
      )
      res <- run_enrichment(
        gene_list      = gene_list_local,
        enriched_genes = enriched_genes,
        method         = "hypergeometric",
        universe       = current_universe,
        alternative    = alternative
      )
    }

    if (!is.null(res)) {
      res[[group_col]] <- group_name
      return(as.data.frame(res))
    }
    return(NULL)
  })

  # --- Merge results ---
  results_list <- results_list[!vapply(results_list, is.null, logical(1))]
  if (length(results_list) == 0) return(data.frame())
  do.call(rbind, results_list)
}
