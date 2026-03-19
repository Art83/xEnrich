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
#' @section cutoff_type = "specificity":
#' Unlike the other cutoff types, which rank genes by expression level within
#' a single cell type, \code{"specificity"} computes a cross-cell-type fold
#' change: a gene is included if its mean expression in the cell type is at
#' least \code{cutoff_value}-fold above the mean across all \emph{other} cell
#' types in the same organ. This identifies genes that are distinctive to a
#' cell type rather than merely highly expressed in it.
#'
#' The distinction matters for interpretation: a quantile cutoff answers "which
#' cell types express my query genes at high absolute levels?"; specificity
#' answers "which cell types are my query genes distinctive to?" For most
#' biomarker discovery questions, specificity is the more informative question.
#'
#' This is conceptually equivalent to TissueEnrich's tissue-enriched gene
#' definition applied at the cell-type level within a single organ, using
#' the same cross-group fold-change approach as
#' \code{\link{enrich_by_hpa}(..., cutoff_type = "specificity")}.
#'
#' Typical \code{cutoff_value} for cell-type specificity: \code{2} (moderate,
#' includes co-expressed genes), \code{3} (default, selective), \code{5}
#' (stringent, tissue-enriched equivalent). Lower than HPA tissue specificity
#' because cell types within one organ share more genes than tissues across
#' the whole body.
#'
#' @param gene_list Character vector of gene symbols (the query set).
#' @param reference_df Data frame of Tabula Sapiens expression data for one
#'   organ. Must contain columns \code{GeneSymbol}, \code{CellType},
#'   \code{MeanLogNorm}, and \code{PctExpr}.
#' @param method Character. \code{"hypergeometric"} (default) or \code{"gsea"}.
#' @param cutoff_type Character. How to define the expressed gene set per
#'   cell type for the hypergeometric test. Ignored when \code{method = "gsea"}.
#'   \describe{
#'     \item{\code{"quantile"}}{Genes at or above the \code{cutoff_value}
#'       quantile within the cell type (must be in (0, 1)). Default.}
#'     \item{\code{"threshold"}}{Genes with mean expression strictly above
#'       \code{cutoff_value}.}
#'     \item{\code{"top_k"}}{Top \code{cutoff_value} genes by mean expression
#'       (\code{cutoff_value} must be a positive integer).}
#'     \item{\code{"specificity"}}{Genes whose mean expression in the cell type
#'       is at least \code{cutoff_value}-fold above the mean across all other
#'       cell types in the organ. See the specificity section above.}
#'   }
#' @param cutoff_value Numeric. Threshold interpreted per \code{cutoff_type}.
#'   Default \code{0.9} for quantile. Use \code{3} for specificity.
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
#'   the cell-type universe. Default \code{1}.
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
    cutoff_type    = c("quantile", "threshold", "top_k", "specificity"),
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
  if (cutoff_type == "specificity" && cutoff_value <= 1)
    stop("For `cutoff_type = 'specificity'`, `cutoff_value` must be > 1 ",
         "(fold-change over mean of other cell types). ",
         "Typical values: 2 (permissive), 3 (default), 5 (stringent).")
  if (cutoff_type == "specificity" && method == "gsea")
    stop("`cutoff_type = 'specificity'` is only applicable to ",
         "method = 'hypergeometric'. For GSEA, expression rankings are ",
         "derived from MeanLogNorm directly.")

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
  # Applied globally before splitting so the universe reflects genuinely
  # detectable genes only. Also applied before specificity computation so that
  # sporadically-expressed genes do not inflate cross-cell-type fold changes.
  if (min_pct_expr > 0) {
    reference_df <- reference_df[reference_df$PctExpr >= min_pct_expr, ]
    if (nrow(reference_df) == 0L)
      stop("No rows remain after PctExpr filtering (min_pct_expr = ",
           min_pct_expr, "). Consider lowering this threshold.")
  }

  # --- Pre-compute specificity sets (needs full reference_df) ----------------
  # Specificity requires the cross-cell-type view: fold = mean_in_this_ct /
  # mean_across_all_other_cts. Computed once before the per-cell-type loop.
  specific_genes_ct <- if (cutoff_type == "specificity") {
    n_cts <- length(unique(reference_df$CellType))
    if (n_cts < 2L)
      stop("`cutoff_type = 'specificity'` requires at least 2 cell types in ",
           "`reference_df`. Only ", n_cts, " cell type found.")
    message("Computing cross-cell-type specificity folds across ",
            n_cts, " cell types...")
    .compute_specificity_sets_ct(reference_df, cutoff_value)
  } else NULL

  # --- Global universe --------------------------------------------------------
  global_universe <- switch(
    universe_type,
    global   = unique(reference_df$GeneSymbol),
    provided = {
      u       <- intersect(universe, unique(reference_df$GeneSymbol))
      dropped <- length(universe) - length(u)
      if (dropped > 0L)
        message(dropped, " of ", length(universe),
                " provided universe genes absent from reference_df. Excluded.")
      if (length(u) == 0L)
        stop("Provided universe has no overlap with reference_df$GeneSymbol.")
      u
    },
    NULL   # "group": determined per cell type below
  )

  # --- Per cell-type loop -----------------------------------------------------
  ct_list <- split(reference_df, reference_df$CellType)

  results_list <- lapply(names(ct_list), function(ct_name) {
    ct_df <- ct_list[[ct_name]]
    current_universe <- if (universe_type == "group") {
      unique(ct_df$GeneSymbol)
    } else {
      global_universe
    }

    ct_df <- ct_df[ct_df$GeneSymbol %in% current_universe, ]

    # Per-gene mean MeanLogNorm within this cell type
    genes_in_ct <- unique(ct_df$GeneSymbol)
    local_stats <- vapply(genes_in_ct, function(g) {
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
        threshold   = names(local_stats)[local_stats > cutoff_value],
        quantile    = {
          q_thresh <- quantile(local_stats, cutoff_value, na.rm = TRUE)
          names(local_stats)[local_stats >= q_thresh]
        },
        top_k       = {
          k <- min(as.integer(cutoff_value), length(local_stats))
          names(sort(local_stats, decreasing = TRUE))[seq_len(k)]
        },
        # specificity: pre-computed, intersect with current universe
        specificity = intersect(
          specific_genes_ct[[ct_name]],
          current_universe
        )
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


# =============================================================================
# Internal helpers
# =============================================================================

#' Compute cross-cell-type specificity gene sets from Tabula data
#'
#' For each cell type in a single-organ Tabula Sapiens dataset, identifies
#' genes whose mean expression in that cell type is at least
#' \code{fold_threshold}-fold above the mean across all other cell types.
#' This is the cell-type-level equivalent of \code{.compute_specificity_sets()}
#' for HPA tissue data.
#'
#' The PctExpr filter should be applied to \code{reference_df} before calling
#' this function, so that sporadically-expressed genes do not inflate folds.
#' Uses the same vectorised matrix approach (tapply + row-sum trick) as
#' \code{.compute_specificity_sets()} for efficiency.
#'
#' @param reference_df Data frame with GeneSymbol, CellType, MeanLogNorm.
#'   PctExpr filtering must already be applied.
#' @param fold_threshold Numeric > 1. Minimum fold over mean of other cell types.
#'
#' @return Named list (one entry per cell type) of character vectors of
#'   cell-type-specific gene symbols.
#' @keywords internal
#' @noRd
.compute_specificity_sets_ct <- function(reference_df, fold_threshold) {
  # Build genes x cell types mean-expression matrix via tapply.
  mat <- tapply(
    reference_df$MeanLogNorm,
    list(reference_df$GeneSymbol, reference_df$CellType),
    FUN   = mean,
    na.rm = TRUE
  )
  # mat: rows = GeneSymbol, cols = CellType. Missing combinations are NA.
  all_genes  <- rownames(mat)
  all_cts    <- colnames(mat)
  n_cts      <- ncol(mat)

  if (n_cts < 2L)
    stop("Specificity requires >= 2 cell types.")

  # Row-sum trick:
  # mean_others[g, ct] = (rowSum[g] - mat[g, ct]) / (n_nonNA[g] - 1)
  # Adding a small epsilon (1e-6) avoids division by zero for genes with
  # zero expression in all other cell types (rare but possible after PctExpr
  # filtering reduces some cell types to zero entries).
  row_sums <- rowSums(mat, na.rm = TRUE)
  row_nobs <- rowSums(!is.na(mat))

  lapply(all_cts, function(ct) {
    this_col    <- mat[, ct]
    not_na_this <- !is.na(this_col)

    others_sum <- row_sums - ifelse(not_na_this, this_col, 0)
    others_n   <- row_nobs - as.integer(not_na_this)

    mean_others <- ifelse(others_n > 0L, others_sum / others_n, 0)
    fold        <- this_col / (mean_others + 1e-6)

    all_genes[not_na_this & !is.na(fold) & fold >= fold_threshold]
  }) |> setNames(all_cts)
}
