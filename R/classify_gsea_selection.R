#' Classify selected pathways after GSEA-based redundancy selection
#'
#' Adds interpretive classification columns to the output of
#' \code{\link{run_gsea_redundancy_selection}}, giving the same
#' \code{marginal_strength}, \code{signal_class}, and \code{action}
#' structure as \code{\link{classify_assoc_selection}} for Lane A3.
#'
#' @section Signal strength in GSEA context:
#' Lane A2 uses the MI z-score from \code{\link{run_info_assoc}} as the
#' signal strength axis. Lane A3 does not have a z-score, but GSEA
#' p-values serve the same role: they reflect how confidently the
#' enrichment score exceeds the permutation null. Thresholds are
#' expressed directly on p-value rather than on a transformed scale,
#' making them easier to calibrate without reference to dataset-specific
#' z-score distributions.
#'
#' The signed enrichment score (\code{enrichment_score}) is preserved in
#' the output so the user can distinguish up- from down-enriched pathways,
#' but it does not affect classification — classification is based on
#' independence (redundancy ratio) and statistical confidence (p-value).
#'
#' @section Relationship to classify_assoc_selection():
#' The classification logic is identical to
#' \code{\link{classify_assoc_selection}}. The only differences are:
#' \itemize{
#'   \item Signal strength uses \code{p_value} instead of z-score.
#'     \code{p_strong} corresponds to \code{z_strong},
#'     \code{p_moderate} to \code{z_moderate}.
#'   \item There is no \code{"suppressor"} class. Suppressors arise from
#'     residual orthogonalisation in sample space
#'     (\code{\link{run_assoc_redundancy_selection}}). GSEA redundancy
#'     selection operates in gene space (leading edges), where this
#'     artifact does not occur.
#' }
#' The output columns have the same names so downstream code (tables,
#' figures, vignettes) can treat both lanes identically.
#'
#' @section Threshold guidance:
#' Defaults correspond to conventional significance levels:
#' \itemize{
#'   \item \code{p_strong = 0.001}: result is robust, the permutation
#'     null is clearly exceeded.
#'   \item \code{p_moderate = 0.01}: reliable signal, sensitive to
#'     permutation count and dataset size.
#'   \item \code{p_value > p_moderate}: weak evidence, report with caution.
#' }
#' Adjust redundancy thresholds the same way as in
#' \code{\link{classify_assoc_selection}}: raise \code{ratio_independent}
#' to 0.6 for GO terms (many overlapping annotations), raise
#' \code{ratio_primary} to 0.8 for conservative reporting.
#'
#' @param sel_result Data frame returned by
#'   \code{\link{run_gsea_redundancy_selection}}.
#' @param p_strong Numeric. p-value ceiling for \code{"strong"} signal.
#'   Default \code{0.001}.
#' @param p_moderate Numeric. p-value ceiling for \code{"moderate"}
#'   signal. Default \code{0.01}.
#' @param ratio_primary Numeric in (0, 1). Minimum \code{redundancy_ratio}
#'   for \code{"primary"} class (requires \code{p_value <= p_strong}).
#'   Default \code{0.70}.
#' @param ratio_independent Numeric in (0, 1). Minimum
#'   \code{redundancy_ratio} for \code{"independent"} class. Default
#'   \code{0.50}.
#' @param ratio_partial Numeric in (0, 1). Minimum \code{redundancy_ratio}
#'   for \code{"partial"} class (below = \code{"redundant"}). Default
#'   \code{0.35}.
#'
#' @return The input data frame with three additional columns:
#'   \describe{
#'     \item{\code{marginal_strength}}{Ordered factor: \code{"strong"},
#'       \code{"moderate"}, \code{"weak"}. Based on \code{p_value}.}
#'     \item{\code{signal_class}}{Ordered factor: \code{"primary"},
#'       \code{"independent"}, \code{"partial"}, \code{"redundant"}.
#'       Note: no \code{"suppressor"} class (see Details).}
#'     \item{\code{action}}{Character. Plain-language reporting
#'       recommendation matching the vocabulary of
#'       \code{\link{classify_assoc_selection}}.}
#'   }
#'
#' @seealso \code{\link{run_gsea_redundancy_selection}} for the upstream
#'   selection step, \code{\link{classify_assoc_selection}} for the
#'   equivalent Lane A2 function, \code{\link{plot_gains}} to visualise
#'   the selection profile.
#'
#' @examples
#' set.seed(1)
#' universe   <- paste0("G", 1:1000)
#' gene_stats <- setNames(rnorm(1000), universe)
#' gene_stats[universe[1:70]]    <- gene_stats[universe[1:70]]    + 4
#' gene_stats[universe[500:530]] <- gene_stats[universe[500:530]] + 3
#'
#' gene_sets <- list(
#'   pathway_A1 = universe[1:60],
#'   pathway_A2 = c(universe[1:55], universe[61:65]),
#'   pathway_B1 = universe[500:530],
#'   pathway_N1 = sample(universe[200:400], 40)
#' )
#'
#' batch <- run_batch_gsea(
#'   gene_sets  = gene_sets,
#'   gene_stats = gene_stats,
#'   n_perm     = 500L,
#'   seed       = 1L
#' )
#'
#' sel <- run_gsea_redundancy_selection(
#'   gsea_results = batch$results,
#'   gene_stats   = gene_stats,
#'   alpha        = 0.05,
#'   seed         = 1L
#' )
#'
#' classified <- classify_gsea_selection(sel)
#' print(classified[, c("pathway", "marginal_strength",
#'                       "signal_class", "action")])
#'
#' @export
classify_gsea_selection <- function(
    sel_result,
    p_strong          = 0.001,
    p_moderate        = 0.01,
    ratio_primary     = 0.70,
    ratio_independent = 0.50,
    ratio_partial     = 0.35
) {
  ## --- Validation ------------------------------------------------------------
  if (!is.data.frame(sel_result))
    stop("`sel_result` must be a data frame from run_gsea_redundancy_selection().")

  required <- c("pathway", "redundancy_ratio", "p_value", "marginal_mi")
  missing  <- setdiff(required, colnames(sel_result))
  if (length(missing) > 0L)
    stop("Missing columns in `sel_result`: ", paste(missing, collapse = ", "),
         ". Pass the direct output of run_gsea_redundancy_selection().")

  if (p_strong  >= p_moderate || p_moderate >= 1 || p_strong <= 0)
    stop("Thresholds must satisfy: 0 < p_strong < p_moderate < 1.")
  if (ratio_primary     <= ratio_independent ||
      ratio_independent <= ratio_partial     ||
      ratio_partial     <= 0)
    stop("Thresholds must satisfy: ratio_primary > ratio_independent > ratio_partial > 0.")

  p     <- sel_result$p_value
  ratio <- sel_result$redundancy_ratio

  ## --- Marginal strength (p-value based) ------------------------------------
  marginal_strength <- ifelse(
    p <= p_strong,
    "strong",
    ifelse(p <= p_moderate, "moderate", "weak")
  )
  marginal_strength <- factor(
    marginal_strength,
    levels  = c("strong", "moderate", "weak"),
    ordered = TRUE
  )

  ## --- Signal class ----------------------------------------------------------
  # No "suppressor" class: leading-edge CMI selection in gene space does not
  # produce the residual-conditioning artifact that creates suppressors in
  # run_assoc_redundancy_selection's sample-space orthogonalisation.
  signal_class <- ifelse(
    ratio >= ratio_primary & p <= p_strong,
    "primary",
    ifelse(
      ratio >= ratio_independent,
      "independent",
      ifelse(ratio >= ratio_partial, "partial", "redundant")
    )
  )
  signal_class <- factor(
    signal_class,
    levels  = c("primary", "independent", "partial", "redundant"),
    ordered = TRUE
  )

  ## --- Plain-language action -------------------------------------------------
  action <- ifelse(
    signal_class == "redundant",
    "dismiss",
    ifelse(
      signal_class == "primary",
      "report",
      ifelse(
        signal_class == "independent" & marginal_strength %in% c("strong", "moderate"),
        "report as secondary",
        ifelse(
          signal_class == "independent" & marginal_strength == "weak",
          "low confidence, treat cautiously",
          ifelse(
            signal_class == "partial" & marginal_strength == "strong",
            "mention as correlated",
            "dismiss"
          )
        )
      )
    )
  )

  ## --- Append and return -----------------------------------------------------
  sel_result$marginal_strength <- marginal_strength
  sel_result$signal_class      <- signal_class
  sel_result$action            <- action

  sel_result
}
