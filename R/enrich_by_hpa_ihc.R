#' Group-wise IHC Enrichment Analysis
#'
#' Applies enrichment testing (GSEA or hypergeometric) across multiple groups
#' (e.g., tissues, or cell types) using Human Protein Atlas (HPA)
#' Immunohistochemistry (IHC) data.
#'
#' @param gene_list A character vector of gene symbols (the "test set").
#' @param reference_df A data frame containing HPA IHC data. Must include 'Gene'
#'   and 'Level' columns, as well as the column specified in \code{group_col}.
#' @param group_col Character string specifying the column to group by
#'   (e.g., "Tissue", "Cell.type"). Defaults to the first element of
#'   \code{c("Tissue", "Cell.type")}.
#' @param remove_levels_of_reliability Character vector of reliability scores
#'   to exclude (e.g., "Uncertain").
#' @param method Character string, either \code{"hypergeometric"} or \code{"gsea"}.
#' @param universe Optional character vector defining the gene universe.
#' @param universe_type Character string defining how the universe is constructed:
#'   \itemize{
#'     \item \code{"global"}: All unique genes in \code{reference_df}.
#'     \item \code{"group"}: Only genes measured within the specific group being tested.
#'     \item \code{"provided"}: Uses the \code{universe} argument.
#'   }
#' @param min_overlap Minimum number of query genes that must overlap the
#'   universe for a group to be tested. Applied after universe intersection.
#'   Groups with fewer overlapping query genes are silently skipped. Default 1.
#' @param min_background Minimum size of the background gene set required to
#'   attempt enrichment. For \code{"hypergeometric"}, this is the number of
#'   expressed genes in the group (i.e. genes at levels in
#'   \code{expression_levels}); for \code{"gsea"}, this is the length of the
#'   gene ranking vector. Default 10.
#' @param gene_stats Optional named numeric vector of pre-calculated statistics
#'   (required only if \code{method="gsea"} and \code{universe_type="provided"}).
#' @param seed Integer seed for reproducibility in GSEA permutations.
#' @param expression_levels Character vector of IHC levels considered "enriched"
#'   for the hypergeometric test (e.g., \code{c("Low", "Medium", "High")}).
#' @param alternative Character string specifying the alternative hypothesis.
#' @param ... Additional arguments passed to \code{run_enrichment} (e.g.,
#'   \code{n_perm}, \code{gsea_weight}, \code{adaptive}, \code{alpha}, \code{eps}).
#'
#' @return A data frame containing enrichment results for each group.
#' @export
enrich_by_hpa_ihc <- function(gene_list,
                              reference_df,
                              group_col = c("Tissue","Cell.type"),
                              remove_levels_of_reliability = c("Uncertain"),
                              method = c("hypergeometric", "gsea"),
                              universe = NULL,
                              universe_type = c("global", "group", "provided"),
                              min_overlap = 1,
                              min_background = 10,
                              gene_stats = NULL,
                              seed=NULL,
                              expression_levels = c("Low", "Medium", "High"),
                              alternative = c("greater", "less", "two.sided"),
                              ...) {
  # Mandatory args
  method <- match.arg(method)
  universe_type <- match.arg(universe_type)
  alternative <- match.arg(alternative)
  if (missing(group_col) || length(group_col) > 1) {
    if (!missing(group_col) && length(group_col) > 1)
      warning("'group_col' has length > 1; using only the first element: '", group_col[1], "'")
    group_col <- group_col[1]
  }

  stopifnot(is.data.frame(reference_df))
  if( !is.null(universe) && (!is.character(universe) || length(universe)==0) ) stop("Universe must be a non-empty character vector")
  if ( any( !c("Level","Gene") %in% colnames(reference_df)) )stop("Check the reference table. It should have Level and Gene")
  if( !group_col %in% colnames(reference_df) ) stop("Grouping varaible is absent in the reference table")
  if( min_overlap < 1 || !is.numeric(min_overlap) || length(min_overlap) > 1 ) stop("Check 'min_overlap'. Should be one integer > 0")
  if( min_background < 1 || !is.numeric(min_background) || length(min_background) > 1 ) stop("Check 'min_background'. Should be one integer > 0")
  if( universe_type=="provided" && is.null(universe) ) stop("universe_type = provided but universe is NULL")
  if (!is.null(gene_stats) && !(method == "gsea" && universe_type == "provided")) {
    warning("`gene_stats` is only used when method='gsea' and universe_type='provided'. It will be ignored in this run.")
  }
  if ("Reliability" %in% colnames(reference_df)) {
    reference_df <- reference_df[!reference_df$Reliability %in% remove_levels_of_reliability & !is.na(reference_df$Reliability), ]
    if(nrow(reference_df) == 0) stop("0 rows after reliability filtering.")
  } else {
    warning("Column 'Reliability' not found, skipping reliability filtering.")
  }
  if(method == "hypergeometric" && !all(expression_levels %in% unique(reference_df$Level)) ) stop("Check expression_levels.No such levels in reference table")
  ihc_levels <- c("Not detected" = 0, "Low" = 1, "Medium" = 2, "High" = 3)
  reference_df$expr_numeric <- ihc_levels[reference_df$Level]
  if( any(is.na(reference_df$expr_numeric))) {
    warning("NAs in Level detected")
    reference_df <- reference_df[!is.na(reference_df$expr_numeric),]
    if(nrow(reference_df) == 0) stop("No rows left after removing NAs in ihc_levels")
  }
  dots <- list(...)
  if (method == "gsea" && length(dots) ){
    if( is.null(names(dots)) || any(names(dots)=="")){
      stop("Check the names for arguments in ellipsis")
    }
    allowed <- c("n_perm","gsea_weight","adaptive","adaptive_mode","alpha","eps")
    bad <- setdiff(names(dots), allowed)
    if (length(bad)) stop("Unsupported `...` for GSEA: ", paste(bad, collapse=", "))
  }
  # Define global universe if needed
  global_universe <- if (universe_type == "global") {
    unique(reference_df$Gene)
  } else if (universe_type == "provided"){
    universe
  } else {
    NULL
  }
  if (universe_type == "provided") {
    original_len <- length(global_universe)
    global_universe <- intersect(global_universe, unique(reference_df$Gene))
    dropped <- original_len - length(global_universe)
    if (dropped > 0)
      message(dropped, " of ", original_len, " provided universe genes are absent from reference_df and will be excluded.")
    if (length(global_universe) == 0) stop("Provided universe has no overlap with reference_df$Gene.")
  }

  if (method == "gsea" && universe_type == "provided") {
    if (is.null(gene_stats)) stop("For universe_type='provided', `gene_stats` must be provided.")
    if (!is.numeric(gene_stats) || is.null(names(gene_stats))) {
      stop("`gene_stats` must be a *named numeric* vector (names = gene IDs).")
    }
  }


  ihc_list <- split(reference_df, reference_df[[group_col]])


  results_list <- lapply(names(ihc_list), function(group_name) {
    ref_group <- ihc_list[[group_name]]

    current_universe <- if (universe_type == "group") unique(ref_group$Gene) else global_universe
    ref_group <- ref_group[ref_group$Gene %in% current_universe, ]

    if (method == "gsea") {
      if (universe_type == "provided") {


        # restrict stats to the current universe
        local_stats <- gene_stats[names(gene_stats) %in% current_universe]

        # ensure universe matches the stats support (avoid passing genes with no stats)
        current_universe <- intersect(current_universe, names(local_stats))

        if (length(local_stats) < min_background) return(NULL)
        if (length(current_universe) == 0) return(NULL)

      } else {
        # derive stats from IHC levels
        local_stats <- tapply(ref_group$expr_numeric, ref_group$Gene, mean, na.rm = TRUE)
        local_stats <- setNames(as.vector(local_stats), names(local_stats))
        # restrict to current universe (good hygiene even here)
        local_stats <- local_stats[names(local_stats) %in% current_universe]
        current_universe <- intersect(current_universe, names(local_stats))

        if (length(local_stats) < min_background) return(NULL)
        if (length(current_universe) == 0) return(NULL)
      }

      gene_list_local <- intersect(gene_list, current_universe)
      n_overlap <- length(gene_list_local)
      if (n_overlap == 0) return(NULL)
      if (n_overlap < min_overlap)
        message("Group '", group_name, "': only ", n_overlap, " of ", length(gene_list),
                " query genes found in universe. Results may be unreliable.")

      res <- run_enrichment(
        gene_list = gene_list_local,
        gene_stats = local_stats,
        method = "gsea",
        alternative = alternative,
        universe = current_universe,
        seed = seed,
        ...
      )
    } else {

      expressed_genes <- unique(ref_group$Gene[ref_group$Level %in% expression_levels])

      if (length(expressed_genes) < min_background) return(NULL)

      gene_list_local <- intersect(gene_list, current_universe)
      n_overlap <- length(gene_list_local)
      if (n_overlap == 0) return(NULL)
      if (n_overlap < min_overlap)
        message("Group '", group_name, "': only ", n_overlap, " of ", length(gene_list),
                " query genes found in universe. Results may be unreliable.")
      res <- run_enrichment(
        gene_list = gene_list_local,
        enriched_genes = expressed_genes,
        method = "hypergeometric",
        universe = current_universe,
        alternative = alternative
      )
    }
    if (!is.null(res)) {
      res[[group_col]] <- group_name
      return(as.data.frame(res))
    }
    return(NULL)
  })
  results_list <- results_list[!vapply(results_list, is.null, logical(1))]
  if (length(results_list) == 0) return(data.frame())
  do.call(rbind, results_list)
}





