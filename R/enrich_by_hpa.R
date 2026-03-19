#' Group-wise RNA Expression Enrichment Analysis
#'
#' Applies enrichment testing across multiple groups (e.g. tissues, brain
#' regions, immune cell types) using Human Protein Atlas (HPA) RNA expression
#' data. Expression values are aggregated per gene per group and used either
#' directly as GSEA rankings or discretised for hypergeometric testing.
#'
#' @param gene_list Character vector of gene symbols (the query set).
#' @param reference_df Data frame of HPA RNA expression data. Must contain a
#'   \code{Gene} column, a grouping column, and an expression value column.
#' @param group_col Character scalar. Column to group by. Auto-detected from
#'   \code{c("Brain.region", "Subregion", "Immune.cell", "Tissue")} if not
#'   supplied.
#' @param expression_col Character scalar. Expression value column. Auto-
#'   detected from \code{c("pTPM", "nTPM", "TPM", "Tags.per.million",
#'   "Scaled.tags.per.million", "Normalized.tags.per.million")} if not
#'   supplied.
#' @param method Character. \code{"hypergeometric"} (default) or \code{"gsea"}.
#' @param cutoff_type Character. How to define expressed genes for the
#'   hypergeometric test. Ignored when \code{method = "gsea"}.
#'   \describe{
#'     \item{\code{"threshold"}}{Genes with mean expression strictly above
#'       \code{cutoff_value}.}
#'     \item{\code{"quantile"}}{Genes at or above the \code{cutoff_value}
#'       quantile within the tissue (must be in (0, 1)).}
#'     \item{\code{"top_k"}}{Top \code{cutoff_value} genes by mean expression
#'       (\code{cutoff_value} must be a positive integer).}
#'     \item{\code{"specificity"}}{Genes whose mean expression in the tissue
#'       is at least \code{cutoff_value}-fold above the mean across all
#'       \emph{other} tissues. Equivalent to TissueEnrich's tissue-enriched
#'       gene definition (\code{cutoff_value = 5}). Unlike the other cutoff
#'       types, this is computed from the full \code{reference_df} before
#'       the per-group test, so the data frame must contain all tissues.
#'       Produces much smaller gene sets than \code{"quantile"} and is
#'       more appropriate when you want tissue-exclusive rather than
#'       tissue-high genes.}
#'   }
#' @param cutoff_value Numeric. Threshold value interpreted according to
#'   \code{cutoff_type}. Default \code{0.9} (suitable for
#'   \code{"threshold"} and \code{"quantile"}; set \code{5} for
#'   \code{"specificity"}; set an integer for \code{"top_k"}).
#' @param universe Optional character vector. Required when
#'   \code{universe_type = "provided"}.
#' @param universe_type Character. Background gene universe construction:
#'   \describe{
#'     \item{\code{"global"}}{All unique genes in \code{reference_df}.}
#'     \item{\code{"group"}}{Only genes measured within the group being tested.}
#'     \item{\code{"provided"}}{Use \code{universe} directly.}
#'   }
#' @param min_overlap Integer. Minimum overlap between \code{gene_list} and
#'   the group universe. Groups below this threshold are silently skipped.
#'   Default \code{1}.
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
#' @return A data frame with one row per group that passed filters, containing
#'   enrichment results and a \code{padj} column (BH adjustment across groups).
#'   Returns an empty data frame if no groups pass.
#'
#' @seealso [run_enrichment()], [enrich_by_hpa_ihc()], [build_tissue_sets()]
#' @export
enrich_by_hpa <- function(
    gene_list,
    reference_df,
    group_col      = c("Brain.region", "Subregion", "Immune.cell", "Tissue"),
    expression_col = c("pTPM", "nTPM", "TPM", "Tags.per.million",
                       "Scaled.tags.per.million",
                       "Normalized.tags.per.million"),
    method         = c("hypergeometric", "gsea"),
    cutoff_type    = c("threshold", "quantile", "top_k", "specificity"),
    cutoff_value   = 0.9,
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
  if (!is.data.frame(reference_df))
    stop("`reference_df` must be a data frame.")

  method        <- match.arg(method)
  universe_type <- match.arg(universe_type)
  alternative   <- match.arg(alternative)
  cutoff_type   <- match.arg(cutoff_type)

  # Auto-detect group_col and expression_col from available columns
  group_col <- if (length(group_col) > 1L) {
    .resolve_col(group_col, reference_df, "group_col")
  } else {
    match.arg(group_col,
              choices = c("Brain.region", "Subregion", "Immune.cell", "Tissue"))
  }

  expression_col <- if (length(expression_col) > 1L) {
    .resolve_col(expression_col, reference_df, "expression_col")
  } else {
    match.arg(expression_col,
              choices = c("pTPM", "nTPM", "TPM", "Tags.per.million",
                          "Scaled.tags.per.million",
                          "Normalized.tags.per.million"))
  }

  # --- Column presence checks -------------------------------------------------
  required_cols <- c("Gene", group_col, expression_col)
  missing_cols  <- setdiff(required_cols, colnames(reference_df))
  if (length(missing_cols) > 0L)
    stop("Columns missing from `reference_df`: ",
         paste(missing_cols, collapse = ", "), ".")

  # --- Scalar argument checks -------------------------------------------------
  if (!is.null(universe) &&
      (!is.character(universe) || length(universe) == 0L))
    stop("`universe` must be a non-empty character vector.")
  if (universe_type == "provided" && is.null(universe))
    stop("`universe_type = 'provided'` but `universe` is NULL.")
  if (!is.numeric(min_overlap) || length(min_overlap) != 1L || min_overlap < 1)
    stop("`min_overlap` must be a single integer >= 1.")
  if (!is.numeric(min_background) || length(min_background) != 1L ||
      min_background < 1)
    stop("`min_background` must be a single integer >= 1.")

  # --- Cutoff checks ----------------------------------------------------------
  if (cutoff_type == "quantile" && (cutoff_value <= 0 || cutoff_value >= 1))
    stop("For `cutoff_type = 'quantile'`, `cutoff_value` must be in (0, 1).")
  if (cutoff_type == "top_k" &&
      (!is.numeric(cutoff_value) || cutoff_value %% 1 != 0 || cutoff_value < 1))
    stop("For `cutoff_type = 'top_k'`, `cutoff_value` must be a positive integer.")
  if (cutoff_type == "specificity" && cutoff_value <= 1)
    stop("For `cutoff_type = 'specificity'`, `cutoff_value` must be > 1 ",
         "(fold-change over mean of other tissues). ",
         "Typical values: 4 (permissive), 5 (TissueEnrich default), 10 (stringent).")
  if (cutoff_type == "specificity" && method == "gsea")
    stop("`cutoff_type = 'specificity'` is only applicable to ",
         "method = 'hypergeometric'. For GSEA, absolute expression rankings ",
         "are derived internally from the expression column.")

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

  # --- Pre-compute specificity gene sets (needs full reference_df) -----------
  # All other cutoff types work within a single group slice.
  # Specificity requires the global picture: fold = mean_in_group /
  # mean_across_all_other_groups. Computed once before the per-group loop.
  specific_genes <- if (cutoff_type == "specificity") {
    n_groups <- length(unique(reference_df[[group_col]]))
    if (n_groups < 2L)
      stop("`cutoff_type = 'specificity'` requires at least 2 groups in ",
           "`reference_df`. Only ", n_groups, " group found.")
    message("Computing cross-tissue specificity folds across ",
            n_groups, " groups...")
    .compute_specificity_sets(reference_df, group_col, expression_col,
                              cutoff_value)
  } else NULL

  # --- Universe construction --------------------------------------------------
  global_universe <- switch(
    universe_type,
    global   = unique(reference_df$Gene),
    provided = {
      u       <- intersect(universe, unique(reference_df$Gene))
      dropped <- length(universe) - length(u)
      if (dropped > 0L)
        message(dropped, " of ", length(universe),
                " provided universe genes absent from reference_df. Excluded.")
      if (length(u) == 0L)
        stop("Provided universe has no overlap with reference_df$Gene.")
      u
    },
    NULL   # "group": determined per group below
  )

  # --- Per-group loop ---------------------------------------------------------
  group_list <- split(reference_df, reference_df[[group_col]])

  results_list <- lapply(names(group_list), function(group_name) {
    ref_group <- group_list[[group_name]]
    current_universe <- if (universe_type == "group") {
      unique(ref_group$Gene)
    } else {
      global_universe
    }
    ref_group <- ref_group[ref_group$Gene %in% current_universe, ]

    # Derive per-gene mean expression (used for GSEA ranking or non-specificity cutoffs)
    genes_in_group <- unique(ref_group$Gene)
    local_stats <- vapply(genes_in_group, function(g) {
      mean(ref_group[[expression_col]][ref_group$Gene == g], na.rm = TRUE)
    }, numeric(1L))
    local_stats      <- local_stats[names(local_stats) %in% current_universe]
    current_universe <- intersect(current_universe, names(local_stats))

    if (length(local_stats) < min_background) return(NULL)
    if (length(current_universe) == 0L)       return(NULL)

    # Override with user-supplied stats for gsea + provided
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
        },
        # specificity: use pre-computed sets, intersect with current universe
        specificity = intersect(
          specific_genes[[group_name]],
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

    # Flatten to one-row data frame safely
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
      row$leading_edge       <- paste(res$leading_edge, collapse = ";")
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


# =============================================================================
# Internal helpers
# =============================================================================

#' Compute cross-tissue specificity gene sets
#'
#' For each group (tissue), identifies genes whose mean expression in that
#' group is at least \code{fold_threshold}-fold above the mean across all
#' other groups. This is the cross-tissue relative enrichment definition
#' used by TissueEnrich (default threshold: 5).
#'
#' Uses vectorised matrix operations: builds a genes x groups mean-expression
#' matrix via \code{tapply}, then for each column computes fold =
#' col_value / mean(all_other_columns). The row-sum trick avoids an
#' O(n_genes x n_groups^2) loop.
#'
#' @param reference_df Data frame with Gene, group_col, and expression_col.
#' @param group_col Character. Group column name.
#' @param expression_col Character. Expression value column name.
#' @param fold_threshold Numeric > 1. Minimum fold over other groups.
#'
#' @return Named list (one entry per group) of character vectors of
#'   tissue-specific gene symbols.
#' @keywords internal
#' @noRd
.compute_specificity_sets <- function(reference_df, group_col,
                                      expression_col, fold_threshold) {
  # Build genes x groups mean-expression matrix via tapply.
  # tapply with two index vectors returns a 2D array.
  mat <- tapply(
    reference_df[[expression_col]],
    list(reference_df$Gene, reference_df[[group_col]]),
    FUN   = mean,
    na.rm = TRUE
  )
  # mat: rows = Gene, cols = group. Missing combinations are NA.
  all_genes  <- rownames(mat)
  all_groups <- colnames(mat)
  n_groups   <- ncol(mat)

  if (n_groups < 2L)
    stop("Specificity requires >= 2 groups.")

  # Row-sum trick:
  # mean_others[g, t] = (rowSum[g] - mat[g,t]) / (n_nonNA[g] - 1)
  row_sums <- rowSums(mat, na.rm = TRUE)
  row_nobs <- rowSums(!is.na(mat))

  lapply(all_groups, function(grp) {
    this_col    <- mat[, grp]
    not_na_this <- !is.na(this_col)

    # Sum of other groups for each gene
    others_sum <- row_sums - ifelse(not_na_this, this_col, 0)
    # Count of non-NA other groups
    others_n   <- row_nobs - as.integer(not_na_this)

    # Mean of other groups; 0 when no other groups have data
    mean_others <- ifelse(others_n > 0L,
                          others_sum / others_n,
                          0)

    # Fold change; use small epsilon to avoid division by exactly zero
    fold <- this_col / (mean_others + 1e-6)

    # Return gene names passing the threshold
    all_genes[not_na_this & !is.na(fold) & fold >= fold_threshold]
  }) |> setNames(all_groups)
}


#' Auto-detect the first available column from a priority list
#'
#' Used to resolve vector-defaulted \code{group_col} and
#' \code{expression_col} against the actual columns present in
#' \code{reference_df}. Fails informatively if none are found.
#'
#' @param candidates Character vector of column names in priority order.
#' @param df Data frame to search.
#' @param arg_name Character scalar. Argument name for error messages.
#'
#' @return The first element of \code{candidates} present in
#'   \code{colnames(df)}.
#' @keywords internal
#' @noRd
.resolve_col <- function(candidates, df, arg_name) {
  found <- candidates[candidates %in% colnames(df)]
  if (length(found) == 0L)
    stop("Could not auto-detect `", arg_name, "`. ",
         "None of the expected columns (",
         paste(candidates, collapse = ", "),
         ") are present in `reference_df`. ",
         "Supply `", arg_name, "` explicitly.")
  found[1L]
}
