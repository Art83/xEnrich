#' Classify selected pathways by signal strength and independence
#'
#' Adds interpretive classification columns to the output of
#' \code{\link{run_assoc_redundancy_selection}}, helping users distinguish
#' primary drivers from secondary signals, redundant shadows, and potential
#' artifacts.
#'
#' @section Classification logic:
#' Two orthogonal axes are assessed independently and then combined:
#'
#' \strong{Signal strength} (\code{marginal_strength}) reflects how robustly
#' a pathway associates with the phenotype, based on its marginal z-score
#' from \code{\link{run_info_assoc}}:
#' \itemize{
#'   \item \code{"strong"}: z-score >= \code{z_strong}. Robust signal,
#'     replicable across permutations.
#'   \item \code{"moderate"}: z-score >= \code{z_moderate}. Reliable but
#'     sensitive to sample size.
#'   \item \code{"weak"}: z-score < \code{z_moderate}. Treat cautiously
#'     regardless of independence.
#' }
#'
#' \strong{Independence class} (\code{signal_class}) reflects how much
#' independent information the pathway contributes beyond earlier-selected
#' pathways, based on \code{redundancy_ratio} and the \code{is_suppressor}
#' flag:
#' \itemize{
#'   \item \code{"primary"}: ratio >= \code{ratio_primary} AND z-score >=
#'     \code{z_strong}. Strong, independent signal. Lead finding.
#'   \item \code{"independent"}: ratio >= \code{ratio_independent}, not
#'     primary. Genuine secondary biology worth reporting.
#'   \item \code{"partial"}: ratio >= \code{ratio_partial}, not above.
#'     Contributing signal, not fully independent. Mention as correlated
#'     with earlier selections.
#'   \item \code{"redundant"}: ratio < \code{ratio_partial}. Largely
#'     captured by earlier pathways. Dismiss from primary reporting.
#'   \item \code{"suppressor"}: \code{is_suppressor == TRUE}. Conditional
#'     MI exceeded marginal MI after orthogonalisation -- selection was
#'     driven by noise interaction with already-selected pathways.
#'     Inspect before reporting.
#' }
#'
#' Note that \code{"suppressor"} takes precedence over all other classes
#' since it signals a potential artifact regardless of ratio magnitude.
#'
#' @section Threshold guidance:
#' Defaults are calibrated for typical omics datasets (n = 100-500 samples,
#' Reactome pathway gene sets). Adjust for your context:
#' \itemize{
#'   \item Noisier data (n < 100): lower \code{z_strong} to 7,
#'     \code{z_moderate} to 3.
#'   \item Many correlated pathways (e.g. GO terms): raise
#'     \code{ratio_independent} to 0.6.
#'   \item Conservative reporting: raise \code{ratio_primary} to 0.8.
#' }
#'
#' @param sel_result Data frame returned by
#'   \code{\link{run_assoc_redundancy_selection}}.
#' @param z_strong Numeric. z-score threshold for \code{"strong"} signal.
#'   Default \code{10}.
#' @param z_moderate Numeric. z-score threshold for \code{"moderate"} signal.
#'   Default \code{5}.
#' @param ratio_primary Numeric in (0, 1). Minimum \code{redundancy_ratio}
#'   for \code{"primary"} class. Default \code{0.70}.
#' @param ratio_independent Numeric in (0, 1). Minimum \code{redundancy_ratio}
#'   for \code{"independent"} class. Default \code{0.50}.
#' @param ratio_partial Numeric in (0, 1). Minimum \code{redundancy_ratio}
#'   for \code{"partial"} class (below = \code{"redundant"}). Default
#'   \code{0.35}.
#'
#' @return The input data frame with three additional columns appended:
#'   \describe{
#'     \item{\code{marginal_strength}}{Ordered factor: \code{"strong"},
#'       \code{"moderate"}, \code{"weak"}.}
#'     \item{\code{signal_class}}{Ordered factor: \code{"primary"},
#'       \code{"independent"}, \code{"partial"}, \code{"redundant"},
#'       \code{"suppressor"}.}
#'     \item{\code{action}}{Character. Plain-language reporting
#'       recommendation.}
#'   }
#'
#' @seealso \code{\link{run_assoc_redundancy_selection}} for the upstream
#'   selection step, \code{\link{plot_gains}} to visualise the selection
#'   profile.
#'
#' @examples
#' set.seed(42)
#' n_samples <- 100
#' genes     <- paste0("G", 1:300)
#' factor_A  <- rnorm(n_samples)
#' factor_B  <- rnorm(n_samples)
#' expr      <- matrix(
#'   rnorm(n_samples * length(genes)),
#'   nrow     = n_samples,
#'   dimnames = list(NULL, genes)
#' )
#' expr[, 1:30]    <- expr[, 1:30]    + factor_A
#' expr[, 200:230] <- expr[, 200:230] + factor_B
#' phenotype <- 0.7 * factor_A + 0.5 * factor_B + rnorm(n_samples, sd = 0.4)
#' gene_sets <- list(
#'   pathway_A  = genes[1:30],
#'   pathway_A2 = genes[5:35],
#'   pathway_B  = genes[200:230],
#'   pathway_C  = genes[100:130]
#' )
#' assoc <- run_info_assoc(
#'   expr = expr, phenotype = phenotype,
#'   gene_sets = gene_sets, n_perm = 500, seed = 42
#' )
#' sel <- run_assoc_redundancy_selection(
#'   assoc_results = assoc, expr = expr, phenotype = phenotype,
#'   gene_sets = gene_sets, seed = 42
#' )
#' classified <- classify_assoc_selection(sel)
#' print(classified[, c("pathway", "marginal_strength", "signal_class", "action")])
#'
#' @export
classify_assoc_selection <- function(
    sel_result,
    z_strong          = 10,
    z_moderate        = 5,
    ratio_primary     = 0.70,
    ratio_independent = 0.50,
    ratio_partial     = 0.35
) {
  ## --- Validation ------------------------------------------------------------
  if (!is.data.frame(sel_result))
    stop("`sel_result` must be a data frame from run_assoc_redundancy_selection().")

  required <- c("pathway", "redundancy_ratio", "z_score", "marginal_mi")
  missing  <- setdiff(required, colnames(sel_result))
  if (length(missing) > 0L)
    stop("Missing columns in `sel_result`: ", paste(missing, collapse = ", "),
         ". Pass the direct output of run_assoc_redundancy_selection().")

  if (!("is_suppressor" %in% colnames(sel_result))) {
    # Double residualisation (v0.0.0.9001+) does not produce the suppressor
    # artifact. Ratios > 1 are legitimate synergy: after removing pathway A's
    # contribution from both candidate and phenotype, an independent pathway B
    # may become MORE informative about the residual than it was about the
    # original phenotype (the first selection reduced noise, concentrating
    # the remaining signal). This is expected and biologically meaningful.
    sel_result$is_suppressor <- FALSE
  }

  if (ratio_primary     <= ratio_independent ||
      ratio_independent <= ratio_partial     ||
      ratio_partial     <= 0)
    stop("Thresholds must satisfy: ratio_primary > ratio_independent > ratio_partial > 0.")
  if (z_strong <= z_moderate || z_moderate <= 0)
    stop("Thresholds must satisfy: z_strong > z_moderate > 0.")

  z       <- sel_result$z_score
  ratio   <- sel_result$redundancy_ratio
  is_supp <- sel_result$is_suppressor

  # For classification, cap ratio at 1.0. Values > 1 indicate synergy
  # (conditional MI exceeds marginal MI after double residualisation) and
  # are strictly more independent than ratio = 1. The raw ratio is
  # preserved in the output for transparency; only the classification
  # thresholds operate on the capped value.
  ratio_capped <- pmin(ratio, 1.0)

  ## --- Marginal strength -----------------------------------------------------
  marginal_strength <- ifelse(
    z >= z_strong,
    "strong",
    ifelse(z >= z_moderate, "moderate", "weak")
  )
  marginal_strength <- factor(
    marginal_strength,
    levels  = c("strong", "moderate", "weak"),
    ordered = TRUE
  )

  ## --- Signal class ----------------------------------------------------------
  signal_class <- ifelse(
    is_supp,
    "suppressor",
    ifelse(
      ratio_capped >= ratio_primary & z >= z_strong,
      "primary",
      ifelse(
        ratio_capped >= ratio_independent,
        "independent",
        ifelse(ratio_capped >= ratio_partial, "partial", "redundant")
      )
    )
  )
  signal_class <- factor(
    signal_class,
    levels  = c("primary", "independent", "partial", "redundant", "suppressor"),
    ordered = TRUE
  )

  ## --- Plain-language action -------------------------------------------------
  action <- ifelse(
    signal_class == "suppressor",
    "inspect, possible artifact",
    ifelse(
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
  )

  ## --- Append and return -----------------------------------------------------
  sel_result$marginal_strength <- marginal_strength
  sel_result$signal_class      <- signal_class
  sel_result$action            <- action

  sel_result
}
