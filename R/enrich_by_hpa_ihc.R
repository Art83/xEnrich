#' Group-wise IHC Enrichment Analysis
#'
#' Applies enrichment testing across multiple groups (e.g. tissues, cell types)
#' using Human Protein Atlas (HPA) Immunohistochemistry (IHC) data. For each
#' group, genes are ranked or categorised by IHC expression level and tested
#' for enrichment of \code{gene_list}.
#'
#' @param gene_list Character vector of gene symbols (the query set).
#' @param reference_df Data frame of HPA IHC data. Must contain columns
#'   \code{Gene}, \code{Level}, and the column named in \code{group_col}.
#'   A \code{Reliability} column is used for filtering if present.
#' @param group_col Character scalar. Column to group by (e.g. \code{"Tissue"},
#'   \code{"Cell.type"}). Default \code{"Tissue"}.
#' @param remove_levels_of_reliability Character vector of reliability labels
#'   to exclude. Default \code{"Uncertain"}.
#' @param method Character. \code{"hypergeometric"} (default) or \code{"gsea"}.
#' @param universe Optional character vector. Required when
#'   \code{universe_type = "provided"}.
#' @param universe_type Character. How to construct the per-group background:
#'   \describe{
#'     \item{\code{"global"}}{All unique genes in \code{reference_df}.}
#'     \item{\code{"group"}}{Only genes measured within the group being tested.}
#'     \item{\code{"provided"}}{Use the \code{universe} argument directly.}
#'   }
#' @param min_overlap Integer. Minimum number of \code{gene_list} genes that
#'   must overlap the group universe. Groups below this threshold are silently
#'   skipped. Default \code{1}.
#' @param min_background Integer. Minimum background size required to attempt
#'   enrichment. For hypergeometric: number of expressed genes in the group.
#'   For GSEA: length of the ranking vector. Groups below this are silently
#'   skipped. Default \code{10}.
#' @param gene_stats Optional named numeric vector of pre-computed gene scores.
#'   Used only when \code{method = "gsea"} and \code{universe_type = "provided"}.
#' @param expression_levels Character vector of IHC levels treated as
#'   "expressed" for the hypergeometric test. Default
#'   \code{c("Low", "Medium", "High")}.
#' @param alternative Character. \code{"greater"} (default), \code{"less"},
#'   or \code{"two.sided"}.
#' @param seed Optional integer seed for GSEA permutations.
#' @param ... Additional arguments passed to [run_enrichment()] for GSEA
#'   (\code{n_perm}, \code{gsea_weight}, \code{adaptive}, \code{adaptive_mode},
#'   \code{alpha}, \code{eps}).
#'
#' @return A data frame with one row per group that passed filters, containing
#'   enrichment results plus a \code{padj} column (BH adjustment across groups).
#'   Returns an empty data frame if no groups pass.
#'
#' @seealso [run_enrichment()], [enrich_by_hpa()]
#' @export
enrich_by_hpa_ihc <- function(
    gene_list,
    reference_df,
    group_col                    = "Tissue",
    remove_levels_of_reliability = "Uncertain",
    method                       = c("hypergeometric", "gsea"),
    universe                     = NULL,
    universe_type                = c("global", "group", "provided"),
    min_overlap                  = 1L,
    min_background               = 10L,
    gene_stats                   = NULL,
    expression_levels            = c("Low", "Medium", "High"),
    alternative                  = c("greater", "less", "two.sided"),
    seed                         = NULL,
    ...
) {
  # --- Argument matching ------------------------------------------------------
  method        <- match.arg(method)
  universe_type <- match.arg(universe_type)
  alternative   <- match.arg(alternative)

  # --- Input validation -------------------------------------------------------
  if (!is.data.frame(reference_df))
    stop("`reference_df` must be a data frame.")

  required_cols <- c("Gene", "Level", group_col)
  missing_cols  <- setdiff(required_cols, colnames(reference_df))
  if (length(missing_cols) > 0)
    stop("Columns missing from `reference_df`: ",
         paste(missing_cols, collapse = ", "), ".")

  if (!is.numeric(min_overlap) || length(min_overlap) != 1L || min_overlap < 1)
    stop("`min_overlap` must be a single integer >= 1.")
  if (!is.numeric(min_background) || length(min_background) != 1L ||
      min_background < 1)
    stop("`min_background` must be a single integer >= 1.")
  if (universe_type == "provided" && is.null(universe))
    stop("`universe_type = 'provided'` but `universe` is NULL.")
  if (!is.null(universe) && (!is.character(universe) || length(universe) == 0))
    stop("`universe` must be a non-empty character vector.")

  # Warn about ignored arguments
  if (!is.null(gene_stats) &&
      !(method == "gsea" && universe_type == "provided"))
    warning("`gene_stats` is only used when method = 'gsea' and ",
            "universe_type = 'provided'. Ignoring.")

  # Validate dots
  if (length(list(...)) > 0) {
    dot_names <- names(list(...))
    if (is.null(dot_names) || any(dot_names == ""))
      stop("All `...` arguments must be named.")
    gsea_only <- c("n_perm", "gsea_weight", "adaptive",
                   "adaptive_mode", "alpha", "eps")
    bad <- setdiff(dot_names, gsea_only)
    if (length(bad))
      stop("Unsupported `...` arguments: ", paste(bad, collapse = ", "), ".")
    if (method == "hypergeometric" && length(dot_names) > 0)
      warning("Extra `...` arguments are ignored for ",
              "method = 'hypergeometric': ",
              paste(dot_names, collapse = ", "), ".")
  }

  # --- Reliability filtering --------------------------------------------------
  if ("Reliability" %in% colnames(reference_df)) {
    reference_df <- reference_df[
      !is.na(reference_df$Reliability) &
        !reference_df$Reliability %in% remove_levels_of_reliability, ]
    if (nrow(reference_df) == 0L)
      stop("No rows remain after reliability filtering.")
  } else {
    warning("Column 'Reliability' not found, skipping reliability filtering.")
  }

  # --- IHC numeric encoding ---------------------------------------------------
  ihc_levels <- c("Not detected" = 0L, "Low" = 1L, "Medium" = 2L, "High" = 3L)
  reference_df$expr_numeric <- ihc_levels[reference_df$Level]

  if (any(is.na(reference_df$expr_numeric))) {
    warning("Unrecognised Level values detected, corresponding rows removed.")
    reference_df <- reference_df[!is.na(reference_df$expr_numeric), ]
    if (nrow(reference_df) == 0L)
      stop("No rows remain after removing unrecognised Level values.")
  }

  if (method == "hypergeometric" &&
      !all(expression_levels %in% unique(reference_df$Level)))
    stop("Some `expression_levels` values are absent from `reference_df$Level`.")

  # --- Global universe --------------------------------------------------------
  global_universe <- switch(
    universe_type,
    global   = unique(reference_df$Gene),
    provided = {
      u       <- intersect(universe, unique(reference_df$Gene))
      dropped <- length(universe) - length(u)
      if (dropped > 0)
        message(dropped, " of ", length(universe),
                " provided universe genes absent from reference_df.Excluded.")
      if (length(u) == 0)
        stop("Provided universe has no overlap with reference_df$Gene.")
      u
    },
    NULL   # "group"-computed per group below
  )

  if (method == "gsea" && universe_type == "provided") {
    if (is.null(gene_stats))
      stop("Provide `gene_stats` when method = 'gsea' and ",
           "universe_type = 'provided'.")
    if (!is.numeric(gene_stats) || is.null(names(gene_stats)))
      stop("`gene_stats` must be a named numeric vector.")
  }

  # --- Per-group enrichment ---------------------------------------------------
  ihc_list <- split(reference_df, reference_df[[group_col]])

  results_list <- lapply(names(ihc_list), function(group_name) {

    ref_group       <- ihc_list[[group_name]]
    current_universe <- if (universe_type == "group") {
      unique(ref_group$Gene)
    } else {
      global_universe
    }
    ref_group <- ref_group[ref_group$Gene %in% current_universe, ]

    # ---- GSEA path -----------------------------------------------------------
    if (method == "gsea") {
      local_stats <- if (universe_type == "provided") {
        gs <- gene_stats[names(gene_stats) %in% current_universe]
        current_universe <- intersect(current_universe, names(gs))
        gs
      } else {
        # Derive stats from mean IHC level per gene.explicit vapply
        genes_in_group <- unique(ref_group$Gene)
        gs <- vapply(genes_in_group, function(g) {
          mean(ref_group$expr_numeric[ref_group$Gene == g], na.rm = TRUE)
        }, numeric(1L))
        gs <- gs[names(gs) %in% current_universe]
        current_universe <- intersect(current_universe, names(gs))
        gs
      }

      if (length(local_stats) < min_background) return(NULL)
      if (length(current_universe) == 0L)       return(NULL)

      gene_list_local <- intersect(gene_list, current_universe)
      if (length(gene_list_local) < min_overlap) return(NULL)

      res <- run_enrichment(
        gene_list   = gene_list_local,
        gene_stats  = local_stats,
        method      = "gsea",
        alternative = alternative,
        universe    = current_universe,
        seed        = seed,
        ...
      )

      # ---- Hypergeometric path -------------------------------------------------
    } else {
      expressed_genes <- unique(
        ref_group$Gene[ref_group$Level %in% expression_levels]
      )
      if (length(expressed_genes) < min_background) return(NULL)

      gene_list_local <- intersect(gene_list, current_universe)
      if (length(gene_list_local) < min_overlap) return(NULL)

      res <- run_enrichment(
        gene_list   = gene_list_local,
        gene_set    = expressed_genes,
        method      = "hypergeometric",
        universe    = current_universe,
        alternative = alternative
      )
    }

    if (is.null(res)) return(NULL)

    # Flatten to one-row data frame,handle leading_edge safely
    row <- data.frame(
      group          = group_name,
      method         = res$method,
      p_value        = res$p_value,
      overlap        = res$overlap,
      input_set_size = res$input_set_size,
      universe_size  = res$universe_size,
      stringsAsFactors = FALSE
    )
    row[[group_col]] <- group_name

    if (res$method == "hypergeometric") {
      row$fold_change <- res$fold_change
      row$group_set   <- res$group_set
    } else {
      row$enrichment_score   <- res$enrichment_score
      row$universe_size_used <- res$universe_size_used
      # Store leading edge as a semicolon-delimited string.
      row$leading_edge <- paste(res$leading_edge, collapse = ";")
    }
    row
  })

  # --- Assemble and adjust ----------------------------------------------------
  results_list <- Filter(Negate(is.null), results_list)
  if (length(results_list) == 0L) {
    message("No groups passed the filters.")
    return(data.frame())
  }

  out      <- do.call(rbind, results_list)
  out$padj <- p.adjust(out$p_value, method = "BH")
  out[order(out$padj), ]
}
