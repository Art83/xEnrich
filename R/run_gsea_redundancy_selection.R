#' Non-redundant pathway selection after GSEA via leading-edge CMI
#'
#' Bridges adaptive GSEA (\code{\link{run_enrichment}}) and information-
#' theoretic redundancy selection (\code{\link{run_redundancy_selection}})
#' by using each pathway's \strong{leading edge},the genes driving the
#' enrichment score peak,as the unit of comparison rather than full
#' annotated pathway membership.
#'
#' This addresses a limitation of applying redundancy selection directly to
#' annotated gene sets after GSEA: two pathways may share 60\% of their
#' annotated members yet have completely non-overlapping leading edges,
#' meaning they are driven by independent biology despite superficial
#' overlap. Conversely, pathways with modest annotation overlap may share
#' their entire leading edge, making them functionally redundant.
#'
#' @section Algorithm:
#' \enumerate{
#'   \item Filter \code{gsea_results} to pathways with \code{p_value <=
#'     alpha}.
#'   \item Extract the \code{leading_edge} from each significant result.
#'   \item Construct a \strong{background gene list} as the union of all
#'     significant leading edges. This represents the set of genes that
#'     drove enrichment in at least one significant pathway. A warning is
#'     emitted if this union exceeds \code{max_union_fraction} of the
#'     universe (default 0.30), as MI becomes less discriminating when the
#'     gene list covers a large fraction of the background.
#'   \item Score each leading edge by mutual information with the union
#'     gene list via \code{\link{run_info_enrichment}}. MI naturally
#'     penalises large leading edges that cover much of the union (low
#'     surprise) and rewards small, concentrated ones (high surprise).
#'   \item Run greedy CMI selection via
#'     \code{\link{run_redundancy_selection}} using the leading edges as
#'     gene sets and the union as the query list.
#'   \item Attach GSEA diagnostics (p-value, enrichment score,
#'     leading edge size) to the selection output.
#' }
#'
#' @section Why leading edges and not full pathway membership:
#' The leading edge is the subset of pathway genes that actually drove the
#' enrichment score, genes ranked before the running-sum peak. Full pathway
#' membership includes many genes with weak or absent signal. Using leading
#' edges means redundancy is assessed on the biologically active subset,
#' not on database annotation boundaries which are often arbitrary.
#'
#' @section Connection to run_redundancy_selection:
#' This function is a wrapper that constructs appropriate inputs for
#' \code{\link{run_redundancy_selection}}. All \code{min_gain},
#' \code{fraction}, and \code{max_pathways} arguments are passed through
#' unchanged. The key difference from calling
#' \code{\link{run_redundancy_selection}} directly is the automatic
#' construction of \code{gene_list} from leading edges and the attachment
#' of GSEA diagnostics to the output.
#'
#' @param gsea_results Named list of results from
#'   \code{\link{run_enrichment}(..., method = "gsea")}. Names are used as
#'   pathway identifiers and must be unique.
#' @param gene_stats Named numeric vector of per-gene scores used in the
#'   original GSEA run (e.g. log fold-change, t-statistic). Names are gene
#'   symbols. Used to define the universe.
#' @param alpha Numeric. Significance threshold for filtering
#'   \code{gsea_results}. Default \code{0.05}.
#' @param min_gain \code{NULL} (default, auto-select) or positive numeric.
#'   Passed to \code{\link{run_redundancy_selection}}.
#' @param fraction Numeric in (0, 1). Relevance fraction for auto
#'   \code{min_gain}. Default \code{0.15}. Passed to
#'   \code{\link{run_redundancy_selection}}.
#' @param n_perm_it Integer. Permutations for IT enrichment null.
#'   Default \code{500L}.
#' @param n_perm_gain Integer. Permutations for auto \code{min_gain} noise
#'   floor. Default \code{200L}.
#' @param max_pathways Integer. Hard cap on selected pathways. Default
#'   \code{20L}.
#' @param max_union_fraction Numeric in (0, 1). Warn if the leading-edge
#'   union exceeds this fraction of the universe. Default \code{0.30}.
#' @param min_leading_edge Integer. Minimum leading edge size to include a
#'   pathway. Pathways with very small leading edges (e.g. 1-2 genes) are
#'   dropped as MI estimation is unreliable. Default \code{3L}.
#' @param seed Optional integer passed to \code{\link{set.seed}}.
#'
#' @return A data frame with one row per selected pathway in selection
#'   order. Inherits all columns from
#'   \code{\link{run_redundancy_selection}} output, plus:
#'   \describe{
#'     \item{\code{p_value}}{GSEA p-value from \code{gsea_results}.}
#'     \item{\code{enrichment_score}}{Signed enrichment score.}
#'     \item{\code{leading_edge_size}}{Number of genes in the leading edge.}
#'     \item{\code{full_set_size}}{Full annotated pathway size (overlap with
#'       universe), for comparison with \code{leading_edge_size}.}
#'   }
#'   Attributes:
#'   \describe{
#'     \item{\code{min_gain}}{Threshold used (auto or supplied).}
#'     \item{\code{union_size}}{Size of the leading-edge union gene list.}
#'     \item{\code{universe_size}}{Size of the universe.}
#'     \item{\code{union_fraction}}{Union as fraction of universe.}
#'     \item{\code{n_sig}}{Number of significant pathways entering selection.}
#'     \item{\code{n_dropped_small_le}}{Pathways dropped for small leading
#'       edge.}
#'   }
#'
#' @seealso
#'   \code{\link{run_enrichment}} for GSEA with adaptive permutation testing,
#'   \code{\link{run_redundancy_selection}} for the underlying CMI selection,
#'   \code{\link{run_info_enrichment}} for marginal IT enrichment scores,
#'   \code{\link{plot_gains}} to visualise the selection profile.
#'
#' @examples
#' set.seed(1)
#' universe   <- paste0("G", 1:1000)
#' gene_stats <- setNames(rnorm(1000), universe)
#'
#' # Three pathway groups:
#' #   A1/A2/A3 - redundant, driven by top-ranked genes
#' #   B1 - independent signal, different gene region
#' #   N1 - noise
#' gene_sets <- list(
#'   pathway_A1 = universe[1:60],
#'   pathway_A2 = c(universe[1:55], universe[61:65]),
#'   pathway_A3 = c(universe[5:60], universe[66:70]),
#'   pathway_B1 = universe[500:530],
#'   pathway_N1 = sample(universe[200:400], 40)
#' )
#'
#' # Run GSEA on each pathway
#' gsea_results <- lapply(gene_sets, function(gs) {
#'   run_enrichment(
#'     gene_list  = gs,
#'     gene_stats = gene_stats,
#'     method     = "gsea",
#'     n_perm     = 500L,
#'     seed       = 1
#'   )
#' })
#' names(gsea_results) <- names(gene_sets)
#'
#' # Select non-redundant pathways based on leading edges
#' sel <- run_gsea_redundancy_selection(
#'   gsea_results = gsea_results,
#'   gene_stats   = gene_stats,
#'   alpha        = 0.05,
#'   seed         = 1
#' )
#' print(sel)
#'
#' @export
run_gsea_redundancy_selection <- function(
    gsea_results,
    gene_stats,
    alpha              = 0.05,
    min_gain           = NULL,
    fraction           = 0.15,
    n_perm_it          = 500L,
    n_perm_gain        = 200L,
    max_pathways       = 20L,
    max_union_fraction = 0.30,
    min_leading_edge   = 3L,
    seed               = NULL
) {
  ## --- Validation ------------------------------------------------------------
  if (!is.list(gsea_results) || is.null(names(gsea_results)))
    stop("`gsea_results` must be a named list of run_enrichment() outputs.")
  if (!is.numeric(gene_stats) || is.null(names(gene_stats)))
    stop("`gene_stats` must be a named numeric vector.")
  if (anyDuplicated(names(gsea_results)))
    stop("`gsea_results` names must be unique.")
  if (alpha <= 0 || alpha >= 1)
    stop("`alpha` must be in (0, 1).")
  if (max_union_fraction <= 0 || max_union_fraction >= 1)
    stop("`max_union_fraction` must be in (0, 1).")

  if (!is.null(seed)) set.seed(seed)

  universe <- names(gene_stats)
  N        <- length(universe)

  ## --- Check all results are GSEA --------------------------------------------
  methods_used <- sapply(gsea_results, `[[`, "method")
  if (!all(methods_used == "gsea"))
    stop(
      "All entries in `gsea_results` must use method = 'gsea'. ",
      "Found non-GSEA results: ",
      paste(names(gsea_results)[methods_used != "gsea"], collapse = ", ")
    )

  ## --- Filter to significant -------------------------------------------------
  pvals   <- sapply(gsea_results, `[[`, "p_value")
  is_sig  <- !is.na(pvals) & pvals <= alpha
  n_sig   <- sum(is_sig)

  if (n_sig == 0L)
    stop(
      "No significant pathways at alpha = ", alpha, ". ",
      "Consider relaxing alpha or checking your GSEA results."
    )

  message(sprintf("%d / %d pathways significant at alpha = %s.",
                  n_sig, length(gsea_results), alpha))

  sig_results <- gsea_results[is_sig]

  ## --- Extract leading edges -------------------------------------------------
  leading_edges <- lapply(sig_results, `[[`, "leading_edge")

  # Drop pathways with very small leading edges
  le_sizes   <- lengths(leading_edges)
  keep_le    <- le_sizes >= min_leading_edge
  n_dropped  <- sum(!keep_le)

  if (n_dropped > 0L) {
    message(sprintf(
      "%d pathway(s) dropped: leading edge < %d genes.",
      n_dropped, min_leading_edge
    ))
    leading_edges <- leading_edges[keep_le]
    sig_results   <- sig_results[keep_le]
  }

  if (length(leading_edges) == 0L)
    stop("No pathways remain after leading edge size filtering.")

  ## --- Intersect leading edges with universe ---------------------------------
  leading_edges <- lapply(leading_edges, intersect, universe)
  le_sizes      <- lengths(leading_edges)

  still_small <- le_sizes < min_leading_edge
  if (any(still_small)) {
    message(sprintf(
      "%d pathway(s) dropped after universe intersection (leading edge too small).",
      sum(still_small)
    ))
    leading_edges <- leading_edges[!still_small]
    sig_results   <- sig_results[!still_small]
  }

  if (length(leading_edges) == 0L)
    stop("No pathways remain after universe intersection.")

  ## --- Construct gene list from leading edge union ---------------------------
  union_genes    <- unique(unlist(leading_edges, use.names = FALSE))
  union_size     <- length(union_genes)
  union_fraction <- union_size / N

  message(sprintf(
    "Leading-edge union: %d genes (%.1f%% of universe).",
    union_size, union_fraction * 100
  ))

  if (union_fraction > max_union_fraction)
    warning(sprintf(
      "Leading-edge union covers %.1f%% of the universe (threshold: %.0f%%). ",
      union_fraction * 100, max_union_fraction * 100
    ), "MI discrimination may be reduced. Consider stricter `alpha` or ",
    "increasing `max_union_fraction` if this is expected for your dataset.",
    call. = FALSE
    )

  ## --- IT enrichment on leading edges ----------------------------------------
  message("Computing IT enrichment scores on leading edges...")
  it_res <- run_info_enrichment(
    gene_list = union_genes,
    gene_sets = leading_edges,
    universe  = universe,
    n_perm    = n_perm_it,
    min_size  = min_leading_edge,
    seed      = seed
  )

  it_scores <- setNames(it_res$info_bits, it_res$set)

  ## --- Redundancy selection --------------------------------------------------
  message("Running CMI redundancy selection on leading edges...")
  sel <- run_redundancy_selection(
    gene_list      = union_genes,
    gene_sets      = leading_edges,
    universe       = universe,
    initial_scores = it_scores,
    min_gain       = min_gain,
    fraction       = fraction,
    n_perm_gain    = n_perm_gain,
    max_pathways   = max_pathways,
    seed           = seed
  )

  if (nrow(sel) == 0L) {
    message("No pathways selected. Returning empty result.")
    return(sel)
  }

  ## --- Attach GSEA diagnostics -----------------------------------------------
  sel$p_value           <- sapply(sel$pathway,
                                  function(nm) sig_results[[nm]]$p_value)
  sel$enrichment_score  <- sapply(sel$pathway,
                                  function(nm) sig_results[[nm]]$enrichment_score)
  sel$leading_edge_size <- sapply(sel$pathway,
                                  function(nm) length(leading_edges[[nm]]))
  sel$full_set_size     <- sapply(sel$pathway,
                                  function(nm) sig_results[[nm]]$input_set_size)

  ## --- Attach metadata as attributes -----------------------------------------
  attr(sel, "min_gain")          <- attr(sel, "min_gain")
  attr(sel, "union_size")        <- union_size
  attr(sel, "universe_size")     <- N
  attr(sel, "union_fraction")    <- union_fraction
  attr(sel, "n_sig")             <- n_sig
  attr(sel, "n_dropped_small_le") <- n_dropped
  attr(sel, "it_scores") <- it_scores
  sel
}
