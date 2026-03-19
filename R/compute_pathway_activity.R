# =============================================================================
# compute_pathway_activity — per-sample pathway activity scores
# =============================================================================

#' Compute per-sample pathway activity scores
#'
#' Summarises gene expression into a single activity score per pathway per
#' sample, using the same scoring methods as \code{\link{run_info_assoc}}.
#' Returns a \strong{samples x pathways} matrix suitable for all downstream
#' biological interpretation: mediation analysis, heatmaps, group comparisons,
#' external regression, or any model that needs per-sample pathway readouts.
#'
#' @section Why this exists:
#' \code{\link{run_info_assoc}} computes pathway activity scores internally
#' but does not expose them — it only returns MI statistics. After running
#' the full Lane A2 pipeline (association + redundancy selection), users
#' typically want to visualise the selected pathways or take them into a
#' mediation model. This function provides those scores using identical
#' logic to the internal computation.
#'
#' @section Scoring methods:
#' \describe{
#'   \item{\code{"mean_z"}}{Mean of z-scored expression across pathway
#'     members. Each gene is standardised across samples before averaging,
#'     so large-variance genes do not dominate. Fast and robust.}
#'   \item{\code{"pc1"}}{First principal component of z-scored expression.
#'     Captures the dominant co-expression axis within the pathway but is
#'     slower, and PC1 sign is arbitrary (irrelevant for MI but matters for
#'     signed downstream uses such as correlation or regression
#'     coefficients).}
#' }
#'
#' @param expr Numeric matrix of expression values, \strong{samples x genes}.
#'   Column names must be gene symbols.
#' @param gene_sets Named list of character vectors defining gene set
#'   membership. Names become column names of the output matrix.
#' @param score Character. Scoring method: \code{"mean_z"} (default) or
#'   \code{"pc1"}. Must match the method used in
#'   \code{\link{run_info_assoc}} for results to be comparable.
#' @param universe Optional character vector. If supplied, both \code{expr}
#'   and gene sets are restricted to these genes before scoring.
#'   Defaults to \code{colnames(expr)}.
#' @param min_genes Integer. Minimum number of a gene set's members that
#'   must be present in \code{expr} (after universe intersection) for that
#'   pathway to be scored. Pathways below this threshold are silently
#'   dropped. Default \code{3L}.
#'
#' @return A numeric matrix with \code{nrow(expr)} rows (samples) and one
#'   column per pathway that passed the \code{min_genes} filter. Column
#'   names are pathway names from \code{gene_sets}; row names are preserved
#'   from \code{expr}. The attribute \code{"dropped"} records pathway names
#'   that were excluded by the \code{min_genes} filter, and \code{"score"}
#'   records the method used.
#'
#' @seealso \code{\link{run_info_assoc}} which computes identical scores
#'   internally during association testing,
#'   \code{\link{run_assoc_redundancy_selection}} for the upstream
#'   selection step.
#'
#' @examples
#' set.seed(1)
#' genes <- paste0("G", 1:200)
#' expr  <- matrix(
#'   rnorm(80 * 200), nrow = 80,
#'   dimnames = list(paste0("S", 1:80), genes)
#' )
#' gene_sets <- list(
#'   pathway_A = genes[1:20],
#'   pathway_B = genes[50:80]
#' )
#'
#' act <- compute_pathway_activity(expr, gene_sets)
#' dim(act)   # 80 x 2
#'
#' # Typical downstream uses:
#' # cor(act, phenotype)
#' # heatmap(act[order(phenotype), ])
#' # lm(phenotype ~ act[, "pathway_A"] + age + sex)
#'
#' @export
compute_pathway_activity <- function(
    expr,
    gene_sets,
    score     = c("mean_z", "pc1"),
    universe  = NULL,
    min_genes = 3L
) {
  ## --- Validation ------------------------------------------------------------
  if (!is.matrix(expr) || !is.numeric(expr))
    stop("`expr` must be a numeric matrix (samples x genes).")
  if (is.null(colnames(expr)))
    stop("`expr` must have gene symbols as column names.")
  if (!is.list(gene_sets) || is.null(names(gene_sets)))
    stop("`gene_sets` must be a named list.")
  if (!is.numeric(min_genes) || min_genes < 1L)
    stop("`min_genes` must be a positive integer.")

  .check_matrix_orientation(expr)
  .check_na_expr(expr)
  expr <- .remove_zero_variance(expr)

  score <- match.arg(score)

  ## --- Universe --------------------------------------------------------------
  if (!is.null(universe)) {
    universe <- intersect(universe, colnames(expr))
    if (length(universe) == 0L)
      stop("`universe` has no overlap with `colnames(expr)`.")
    expr <- expr[, universe, drop = FALSE]
  }

  ## --- Filter gene sets by coverage ------------------------------------------
  gene_sets_filt <- lapply(gene_sets, intersect, colnames(expr))
  n_genes        <- lengths(gene_sets_filt)
  keep           <- n_genes >= min_genes
  dropped        <- names(gene_sets_filt)[!keep]
  gene_sets_filt <- gene_sets_filt[keep]

  if (length(dropped) > 0L) {
    examples <- paste(head(dropped, 3L), collapse = ", ")
    message(
      length(dropped), " pathway(s) dropped: fewer than ", min_genes,
      " genes present in `expr` (e.g. ", examples,
      if (length(dropped) > 3L) ", ..." else "", ")."
    )
  }

  if (length(gene_sets_filt) == 0L)
    stop("No pathways remain after min_genes filtering. ",
         "Check that gene symbols in gene_sets match colnames(expr).")

  ## --- Scale expression (once) -----------------------------------------------
  Z <- scale(expr)

  ## --- Score function --------------------------------------------------------
  score_fun <- switch(
    score,
    mean_z = function(genes) rowMeans(Z[, genes, drop = FALSE]),
    pc1    = function(genes) {
      prcomp(Z[, genes, drop = FALSE], center = FALSE, scale. = FALSE)$x[, 1L]
    }
  )

  ## --- Compute activity matrix -----------------------------------------------
  act_list <- lapply(gene_sets_filt, score_fun)
  act_mat  <- do.call(cbind, act_list)

  rownames(act_mat) <- rownames(expr)
  colnames(act_mat) <- names(gene_sets_filt)

  attr(act_mat, "dropped") <- dropped
  attr(act_mat, "score")   <- score
  act_mat
}
