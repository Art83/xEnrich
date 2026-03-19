# =============================================================================
# run_batch_gsea — batch GSEA with adaptive permutations
# =============================================================================

#' Run GSEA across a collection of gene sets
#'
#' Applies \code{\link{run_enrichment}(..., method = "gsea")} to every
#' pathway in \code{gene_sets} and returns both a named list of per-pathway
#' results (compatible with \code{\link{run_gsea_redundancy_selection}}) and
#' a summary data frame with BH-adjusted p-values across all pathways.
#'
#' This function replaces the boilerplate \code{lapply} loop that users would
#' otherwise write to produce input for \code{\link{run_gsea_redundancy_selection}}.
#' It adds size filtering, optional parallelism, progress reporting,
#' per-pathway error recovery, and BH adjustment.
#'
#' @section Adaptive permutations for batch runs:
#' With \code{adaptive = TRUE} and \code{adaptive_mode = "refine"}
#' (default), each pathway first runs \code{n_perm} fixed permutations.
#' Only pathways whose p-value lands in the borderline zone near
#' \code{alpha} receive additional permutations. In a typical Reactome
#' analysis most pathways are either clearly non-significant (stop early)
#' or very significant (also stop early) -- only a minority require the
#' full refinement budget. Across 1 500 pathways this can halve the total
#' permutation count compared to fixed \code{n_perm} for every pathway.
#'
#' @section Parallelism:
#' Uses \code{parallel::mclapply} (fork-based, Linux/macOS). Each worker
#' receives a deterministic seed derived from \code{seed + pathway_index},
#' so results are reproducible regardless of \code{nthreads}. Silently
#' falls back to sequential on Windows where forking is unavailable.
#' Memory cost per worker is negligible because \code{gene_stats} (a
#' named numeric vector) is small relative to the expression matrix used
#' in \code{\link{run_info_assoc}}.
#'
#' @param gene_sets Named list of character vectors (pathways to test).
#' @param gene_stats Named numeric vector of per-gene scores (e.g. log2FC,
#'   t-statistic, signed -log10 p). Names are gene symbols and define the
#'   universe. Identical to the input of \code{\link{run_enrichment}}.
#' @param n_perm Integer. Initial permutations per pathway. Default
#'   \code{1000L}. For borderline pathways under \code{adaptive_mode =
#'   "refine"}, additional permutations are drawn until the Wilson CI
#'   resolves.
#' @param adaptive Logical. Use adaptive permutation stopping? Default
#'   \code{TRUE}. Strongly recommended for batch runs: saves substantial
#'   computation on clear results.
#' @param adaptive_mode Character. One of \code{"refine"} (default),
#'   \code{"pvalue"}, or \code{"decision"}. See
#'   \code{\link{run_enrichment}} for details. \code{"refine"} is best
#'   for batch use.
#' @param alpha Numeric. Significance level for BH adjustment and
#'   adaptive stopping decisions. Default \code{0.05}.
#' @param eps Numeric. RSE threshold for \code{adaptive_mode = "pvalue"}.
#'   Default \code{0.01}.
#' @param min_size Integer. Pathways with fewer than \code{min_size} genes
#'   (after intersection with \code{names(gene_stats)}) are silently
#'   skipped. Default \code{10L}.
#' @param max_size Integer. Pathways with more than \code{max_size} genes
#'   are silently skipped. Default \code{500L}.
#' @param nthreads Positive integer. Parallel workers via
#'   \code{parallel::mclapply}. Default \code{1L}. Silently reduced to 1
#'   on Windows.
#' @param seed Optional integer. Sets the global RNG seed and derives
#'   per-pathway seeds to ensure reproducible results with any
#'   \code{nthreads} value.
#' @param ... Additional arguments passed to \code{\link{run_enrichment}}
#'   (e.g. \code{gsea_weight}, \code{alternative}).
#'
#' @return A list with two named elements:
#'   \describe{
#'     \item{\code{results}}{Named list of \code{\link{run_enrichment}}
#'       outputs, one per pathway that passed size filtering and ran
#'       successfully. Directly compatible with
#'       \code{\link{run_gsea_redundancy_selection}}.}
#'     \item{\code{summary}}{Data frame (sorted by \code{padj}) with one
#'       row per pathway. Columns: \code{pathway}, \code{set_size},
#'       \code{enrichment_score}, \code{p_value}, \code{padj} (BH),
#'       \code{overlap}, \code{leading_edge_size}. When
#'       \code{adaptive = TRUE}: also \code{perms_used} and
#'       \code{converged}.}
#'   }
#'
#' @seealso \code{\link{run_enrichment}} for single-pathway GSEA,
#'   \code{\link{run_gsea_redundancy_selection}} for the downstream
#'   redundancy selection step.
#'
#' @examples
#' set.seed(1)
#' universe   <- paste0("G", 1:1000)
#' gene_stats <- setNames(rnorm(1000), universe)
#' gene_stats[universe[1:50]] <- gene_stats[universe[1:50]] + 3
#'
#' gene_sets <- list(
#'   pathway_A  = universe[1:50],
#'   pathway_A2 = c(universe[1:45], universe[51:55]),
#'   pathway_B  = universe[500:540],
#'   pathway_N  = sample(universe[200:400], 40)
#' )
#'
#' batch <- run_batch_gsea(
#'   gene_sets  = gene_sets,
#'   gene_stats = gene_stats,
#'   n_perm     = 500L,
#'   adaptive   = TRUE,
#'   seed       = 1L
#' )
#'
#' print(batch$summary)
#'
#' # Pipe directly into redundancy selection
#' sel <- run_gsea_redundancy_selection(
#'   gsea_results = batch$results,
#'   gene_stats   = gene_stats,
#'   alpha        = 0.05,
#'   seed         = 1L
#' )
#' print(sel)
#'
#' @export
run_batch_gsea <- function(
    gene_sets,
    gene_stats,
    n_perm        = 1000L,
    adaptive      = TRUE,
    adaptive_mode = c("refine", "pvalue", "decision"),
    alpha         = 0.05,
    eps           = 0.01,
    min_size      = 10L,
    max_size      = 500L,
    nthreads      = 1L,
    seed          = NULL,
    ...
) {
  ## --- Validation ------------------------------------------------------------
  if (!is.list(gene_sets) || is.null(names(gene_sets)))
    stop("`gene_sets` must be a named list.")
  if (anyDuplicated(names(gene_sets)))
    stop("`gene_sets` names must be unique.")
  if (!is.numeric(gene_stats) || is.null(names(gene_stats)))
    stop("`gene_stats` must be a named numeric vector.")
  if (anyDuplicated(names(gene_stats)))
    stop("`gene_stats` must have unique names.")
  if (n_perm < 1L)
    stop("`n_perm` must be a positive integer.")
  if (alpha <= 0 || alpha >= 1)
    stop("`alpha` must be in (0, 1).")

  adaptive_mode <- match.arg(adaptive_mode)

  .check_pathway_names(names(gene_sets))

  ## --- Size filtering --------------------------------------------------------
  universe <- names(gene_stats)[!is.na(names(gene_stats))]
  sizes    <- vapply(
    gene_sets,
    function(gs) length(intersect(gs, universe)),
    integer(1L)
  )
  keep      <- sizes >= min_size & sizes <= max_size
  n_dropped <- sum(!keep)

  if (n_dropped > 0L)
    message(n_dropped, " pathway(s) dropped by size filter (min_size = ",
            min_size, ", max_size = ", max_size, ").")

  gene_sets <- gene_sets[keep]

  if (length(gene_sets) == 0L)
    stop("No pathways remain after size filtering.")

  n_sets <- length(gene_sets)

  ## --- Parallelism setup -----------------------------------------------------
  nthreads <- if (.Platform$OS.type == "windows") 1L else as.integer(nthreads)

  message(sprintf(
    "Running GSEA on %d pathways [adaptive = %s, mode = '%s', nthreads = %d].",
    n_sets, adaptive, adaptive_mode, nthreads
  ))

  if (!is.null(seed)) set.seed(seed)

  ## --- Run GSEA per pathway --------------------------------------------------
  # Per-pathway seed derived from global seed ensures reproducibility
  # regardless of nthreads value (mclapply workers cannot share RNG state).
  run_one <- function(idx) {
    nm      <- names(gene_sets)[idx]
    gs      <- gene_sets[[idx]]
    pw_seed <- if (!is.null(seed)) seed + idx else NULL

    tryCatch(
      run_enrichment(
        gene_list     = gs,
        gene_stats    = gene_stats,
        method        = "gsea",
        universe      = universe,
        n_perm        = as.integer(n_perm),
        adaptive      = adaptive,
        adaptive_mode = adaptive_mode,
        alpha         = alpha,
        eps           = eps,
        seed          = pw_seed,
        ...
      ),
      error = function(e) {
        message("  Error for pathway '", nm, "': ", conditionMessage(e))
        NULL
      }
    )
  }

  results_list <- if (nthreads > 1L) {
    parallel::mclapply(seq_len(n_sets), run_one, mc.cores = nthreads)
  } else {
    lapply(seq_len(n_sets), run_one)
  }
  names(results_list) <- names(gene_sets)

  ## --- Remove failed pathways ------------------------------------------------
  failed <- vapply(results_list, is.null, logical(1L))
  if (any(failed)) {
    message(sum(failed), " pathway(s) failed during GSEA and were removed.")
    results_list <- results_list[!failed]
  }

  if (length(results_list) == 0L)
    stop("All GSEA runs failed. Check `gene_stats` and `gene_sets` inputs.")

  ## --- Summary data frame ----------------------------------------------------
  summary_rows <- lapply(names(results_list), function(nm) {
    r  <- results_list[[nm]]
    gs <- gene_sets[[nm]]

    row <- data.frame(
      pathway           = nm,
      set_size          = length(intersect(gs, universe)),
      enrichment_score  = if (!is.null(r$enrichment_score)) r$enrichment_score
      else NA_real_,
      p_value           = r$p_value,
      overlap           = r$overlap,
      leading_edge_size = length(r$leading_edge),
      stringsAsFactors  = FALSE
    )

    if (adaptive && !is.null(r$inference)) {
      # Extract total permutations used. The structure varies by mode:
      #   adaptive_pvalue/decision: inference$res$B
      #   fixed_not_refined:        inference$res$n_perm (not B)
      #   fixed_refined:            inference$refined$B
      m <- r$inference$method
      perms <- if (m == "fixed_not_refined") {
        r$inference$res$n_perm
      } else if (m == "fixed_refined") {
        r$inference$refined$B
      } else if (!is.null(r$inference$res$B)) {
        r$inference$res$B
      } else NA_integer_

      row$perms_used <- as.integer(perms)
      row$converged  <- if (m == "fixed_not_refined") {
        TRUE   # resolved without refinement
      } else if (!is.null(r$inference$res$converged)) {
        r$inference$res$converged
      } else if (!is.null(r$inference$refined$converged)) {
        r$inference$refined$converged
      } else NA
    }

    row
  })

  summary_df       <- do.call(rbind, summary_rows)
  summary_df$padj  <- p.adjust(summary_df$p_value, method = "BH")
  summary_df       <- summary_df[order(summary_df$padj), ]
  rownames(summary_df) <- NULL

  n_sig <- sum(!is.na(summary_df$padj) & summary_df$padj <= alpha)
  message(sprintf(
    "Done. %d / %d pathways significant at padj <= %s.",
    n_sig, nrow(summary_df), alpha
  ))

  if (adaptive && "perms_used" %in% colnames(summary_df)) {
    total_perms <- sum(summary_df$perms_used, na.rm = TRUE)
    fixed_equiv <- n_sets * as.integer(n_perm)

    if (adaptive_mode == "refine") {
      n_refined <- sum(summary_df$perms_used > n_perm, na.rm = TRUE)
      n_at_base <- sum(summary_df$perms_used == n_perm, na.rm = TRUE)
      extra     <- total_perms - fixed_equiv
      message(sprintf(
        "Refine mode: %d/%d pathways resolved at %d perms; %d borderline pathways refined (%s extra perms).",
        n_at_base, n_sets, n_perm, n_refined,
        format(extra, big.mark = ",")
      ))
    } else {
      saving_pct <- round(100 * (1 - total_perms / fixed_equiv))
      message(sprintf(
        "Adaptive saving: %s permutations used vs %s fixed equivalent (%d%% reduction).",
        format(total_perms, big.mark = ","),
        format(fixed_equiv, big.mark = ","),
        saving_pct
      ))
    }
  }

  list(
    results = results_list,
    summary = summary_df
  )
}
