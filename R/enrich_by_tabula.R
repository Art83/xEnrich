#' Cell-type enrichment analysis using Tabula Sapiens scRNA data
#'
#' Applies enrichment testing across cell types within a single Tabula Sapiens
#' organ dataset. Genes are filtered by minimum expression prevalence
#' (\code{PctExpr}) before enrichment, reducing noise from genes detected in
#' only a handful of cells. Expression values (\code{MeanLogNorm}) are used
#' either directly as GSEA rankings or discretised for hypergeometric testing.
#'
#' This function is designed as the second step of the location-aware
#' pipeline: after [enrich_by_hpa_ihc()] or [enrich_by_hpa()] identifies
#' enriched organs, [enrich_by_tabula()] resolves which cell types within
#' those organs drive the signal. Use [run_location_pipeline()] to run both
#' steps automatically.
#'
#' @section PctExpr filtering:
#' \code{PctExpr} is the fraction of cells in a cell type that express a
#' given gene (expression > 0). Genes below \code{min_pct_expr} are removed
#' from the background for that cell type before enrichment. This prevents
#' \code{MeanLogNorm} values inflated by a few outlier cells from entering
#' the analysis. A value of \code{0.1} (10\% of cells) is a standard
#' threshold in scRNA analysis.
#'
#' @param gene_list Character vector of gene symbols (the query set).
#' @param reference_df Data frame of Tabula Sapiens expression data for one
#'   organ. Must contain columns \code{GeneSymbol}, \code{CellType},
#'   \code{MeanLogNorm}, and \code{PctExpr}.
#' @param method Character. \code{"hypergeometric"} (default) or \code{"gsea"}.
#' @param cutoff_type Character. How to discretise expression for the
#'   hypergeometric test:
#'   \describe{
#'     \item{\code{"threshold"}}{Genes with mean expression strictly above
#'       \code{cutoff_value}.}
#'     \item{\code{"quantile"}}{Genes at or above the \code{cutoff_value}
#'       quantile (must be in (0, 1)).}
#'     \item{\code{"top_k"}}{Top \code{cutoff_value} genes by mean expression
#'       (\code{cutoff_value} must be a positive integer).}
#'   }
#'   Ignored when \code{method = "gsea"}.
#' @param cutoff_value Numeric. Threshold interpreted per \code{cutoff_type}.
#'   Default \code{0.3} (log-normalised expression; suitable for
#'   \code{"threshold"}).
#' @param min_pct_expr Numeric in \[0, 1\]. Minimum fraction of cells in a
#'   cell type that must express a gene for it to be included in the
#'   background. Default \code{0.1}. Set to \code{0} to disable.
#' @param universe Optional character vector. Required when
#'   \code{universe_type = "provided"}.
#' @param universe_type Character. Background gene universe construction:
#'   \describe{
#'     \item{\code{"global"}}{All genes passing \code{min_pct_expr} across
#'       all cell types in \code{reference_df}.}
#'     \item{\code{"group"}}{Only genes passing \code{min_pct_expr} within
#'       the cell type being tested.}
#'     \item{\code{"provided"}}{Use \code{universe} directly.}
#'   }
#' @param min_overlap Integer. Minimum overlap between \code{gene_list} and
#'   the cell-type universe. Cell types below this threshold are silently
#'   skipped. Default \code{1}.
#' @param min_background Integer. Minimum background size to attempt
#'   enrichment. Default \code{10}.
#' @param gene_stats Optional named numeric vector. Used only when
#'   \code{method = "gsea"} and \code{universe_type = "provided"}.
#' @param alternative Character. \code{"greater"} (default), \code{"less"},
#'   or \code{"two.sided"}.
#' @param seed Optional integer seed for GSEA permutations.
#' @param ... Additional named arguments passed to [run_enrichment()] for
#'   GSEA (\code{n_perm}, \code{gsea_weight}, \code{adaptive},
#'   \code{adaptive_mode}, \code{alpha}, \code{eps}).
#'
#' @return A data frame with one row per cell type that passed filters,
#'   containing enrichment results and a \code{padj} column (BH adjustment
#'   across cell types). Returns an empty data frame if no cell types pass.
#'
#' @seealso [run_location_pipeline()], [enrich_by_hpa_ihc()],
#'   [enrich_by_hpa()]
#' @export
enrich_by_tabula <- function(
    gene_list,
    reference_df,
    method         = c("hypergeometric", "gsea"),
    cutoff_type    = c("quantile", "threshold", "top_k"),
    cutoff_value   = 0.9,
    min_pct_expr   = 0.1,
    universe       = NULL,
    universe_type  = c("global", "group", "provided"),
    min_overlap    = 1L,
    min_background = 10L,
    gene_stats     = NULL,
    alternative    = c("greater", "less", "two.sided"),
    seed           = NULL,
    ...
) {
  # --- Argument matching ------------------------------------------------------
  method        <- match.arg(method)
  cutoff_type   <- match.arg(cutoff_type)
  universe_type <- match.arg(universe_type)
  alternative   <- match.arg(alternative)

  # --- Input validation -------------------------------------------------------
  if (!is.data.frame(reference_df))
    stop("`reference_df` must be a data frame.")

  required_cols <- c("GeneSymbol", "CellType", "MeanLogNorm", "PctExpr")
  missing_cols  <- setdiff(required_cols, colnames(reference_df))
  if (length(missing_cols) > 0L)
    stop("Columns missing from `reference_df`: ",
         paste(missing_cols, collapse = ", "),
         ". Are you passing a Tabula Sapiens dataset from load_reference()?")

  if (!is.numeric(min_pct_expr) || min_pct_expr < 0 || min_pct_expr > 1)
    stop("`min_pct_expr` must be a number in [0, 1].")
  if (!is.numeric(min_overlap) || length(min_overlap) != 1L || min_overlap < 1)
    stop("`min_overlap` must be a single integer >= 1.")
  if (!is.numeric(min_background) || length(min_background) != 1L ||
      min_background < 1)
    stop("`min_background` must be a single integer >= 1.")
  if (universe_type == "provided" && is.null(universe))
    stop("`universe_type = 'provided'` but `universe` is NULL.")
  if (!is.null(universe) &&
      (!is.character(universe) || length(universe) == 0L))
    stop("`universe` must be a non-empty character vector.")

  # --- Cutoff checks ----------------------------------------------------------
  if (cutoff_type == "quantile" && (cutoff_value <= 0 || cutoff_value >= 1))
    stop("For `cutoff_type = 'quantile'`, `cutoff_value` must be in (0, 1).")
  if (cutoff_type == "top_k" &&
      (!is.numeric(cutoff_value) || cutoff_value %% 1 != 0 || cutoff_value < 1))
    stop("For `cutoff_type = 'top_k'`, `cutoff_value` must be a positive integer.")

  # --- Ignored argument warnings ----------------------------------------------
  if (!is.null(gene_stats) &&
      !(method == "gsea" && universe_type == "provided"))
    warning("`gene_stats` is only used when method = 'gsea' and ",
            "universe_type = 'provided'. Ignoring.")

  # --- dots validation --------------------------------------------------------
  dots <- list(...)
  if (length(dots) > 0L) {
    if (is.null(names(dots)) || any(names(dots) == ""))
      stop("All `...` arguments must be named.")
    allowed <- c("n_perm", "gsea_weight", "adaptive",
                 "adaptive_mode", "alpha", "eps")
    bad <- setdiff(names(dots), allowed)
    if (length(bad))
      stop("Unsupported `...` arguments: ", paste(bad, collapse = ", "), ".")
    if (method == "hypergeometric")
      warning("`...` arguments are ignored for method = 'hypergeometric': ",
              paste(names(dots), collapse = ", "), ".")
  }

  # --- GSEA + provided validation ---------------------------------------------
  if (method == "gsea" && universe_type == "provided") {
    if (is.null(gene_stats))
      stop("Provide `gene_stats` when method = 'gsea' and ",
           "universe_type = 'provided'.")
    if (!is.numeric(gene_stats) || is.null(names(gene_stats)))
      stop("`gene_stats` must be a named numeric vector.")
  }

  # --- Remove NA expression ---------------------------------------------------
  reference_df <- reference_df[!is.na(reference_df$MeanLogNorm), ]
  if (nrow(reference_df) == 0L)
    stop("No rows remain after removing NA values in MeanLogNorm.")

  # --- PctExpr filter ---------------------------------------------------------
  # Applied globally before splitting so the universe reflects
  # genuinely detectable genes only. The cutoff then selects the
  # highly-expressed subset within each cell type's own distribution.
  if (min_pct_expr > 0) {
    reference_df <- reference_df[reference_df$PctExpr >= min_pct_expr, ]
    if (nrow(reference_df) == 0L)
      stop("No rows remain after PctExpr filtering (min_pct_expr = ",
           min_pct_expr, "). Consider lowering this threshold.")
  }

  # --- Global universe --------------------------------------------------------
  global_universe <- switch(
    universe_type,
    global   = unique(reference_df$GeneSymbol),
    provided = {
      u       <- intersect(universe, unique(reference_df$GeneSymbol))
      dropped <- length(universe) - length(u)
      if (dropped > 0L)
        message(dropped, " of ", length(universe),
                " provided universe genes absent from reference_df.Excluded.")
      if (length(u) == 0L)
        stop("Provided universe has no overlap with reference_df$GeneSymbol.")
      u
    },
    NULL   # "group": determined per cell type below
  )

  # --- Per cell-type loop -----------------------------------------------------
  cell_types   <- unique(reference_df$CellType)
  ct_list      <- split(reference_df, reference_df$CellType)

  results_list <- lapply(names(ct_list), function(ct_name) {
    ct_df            <- ct_list[[ct_name]]
    current_universe <- if (universe_type == "group") {
      unique(ct_df$GeneSymbol)
    } else {
      global_universe
    }

    ct_df <- ct_df[ct_df$GeneSymbol %in% current_universe, ]

    # Per-gene mean MeanLogNorm within this cell type
    genes_in_ct  <- unique(ct_df$GeneSymbol)
    local_stats  <- vapply(genes_in_ct, function(g) {
      mean(ct_df$MeanLogNorm[ct_df$GeneSymbol == g], na.rm = TRUE)
    }, numeric(1L))
    local_stats      <- local_stats[names(local_stats) %in% current_universe]
    current_universe <- intersect(current_universe, names(local_stats))

    if (length(local_stats) < min_background) return(NULL)
    if (length(current_universe) == 0L)       return(NULL)

    # Override with user stats for gsea + provided
    if (method == "gsea" && universe_type == "provided") {
      local_stats      <- gene_stats[names(gene_stats) %in% current_universe]
      current_universe <- intersect(current_universe, names(local_stats))
      if (length(local_stats) < min_background) return(NULL)
      if (length(current_universe) == 0L)       return(NULL)
    }

    gene_list_local <- intersect(gene_list, current_universe)
    if (length(gene_list_local) < min_overlap) return(NULL)

    # --- Enrichment -----------------------------------------------------------
    if (method == "gsea") {
      res <- run_enrichment(
        gene_list   = gene_list_local,
        gene_stats  = local_stats,
        method      = "gsea",
        alternative = alternative,
        universe    = current_universe,
        seed        = seed,
        ...
      )
    } else {
      gene_set <- switch(
        cutoff_type,
        threshold = names(local_stats)[local_stats > cutoff_value],
        quantile  = {
          q_thresh <- quantile(local_stats, cutoff_value, na.rm = TRUE)
          names(local_stats)[local_stats >= q_thresh]
        },
        top_k = {
          k <- min(as.integer(cutoff_value), length(local_stats))
          names(sort(local_stats, decreasing = TRUE))[seq_len(k)]
        }
      )
      if (length(gene_set) < min_background) return(NULL)

      res <- run_enrichment(
        gene_list   = gene_list_local,
        gene_set    = gene_set,
        method      = "hypergeometric",
        universe    = current_universe,
        alternative = alternative
      )
    }

    if (is.null(res)) return(NULL)

    # Safe row construction
    row <- data.frame(
      cell_type      = ct_name,
      method         = res$method,
      p_value        = res$p_value,
      overlap        = res$overlap,
      input_set_size = res$input_set_size,
      universe_size  = res$universe_size,
      stringsAsFactors = FALSE
    )

    if (res$method == "hypergeometric") {
      row$fold_change <- res$fold_change
      row$group_set   <- res$group_set
    } else {
      row$enrichment_score   <- res$enrichment_score
      row$universe_size_used <- res$universe_size_used
      row$leading_edge       <- paste(res$leading_edge, collapse = ";")
    }
    row
  })

  # --- Assemble and adjust ----------------------------------------------------
  results_list <- Filter(Negate(is.null), results_list)
  if (length(results_list) == 0L) {
    message("No cell types passed the filters.")
    return(data.frame())
  }

  out      <- do.call(rbind, results_list)
  out$padj <- p.adjust(out$p_value, method = "BH")
  out[order(out$padj), ]
}
