#' Non-redundant pathway selection after association analysis via
#' sample-space conditional mutual information
#'
#' Addresses the redundancy problem in pathway-phenotype association:
#' when many pathways are nominally associated with a phenotype, most
#' apparent signal is redundant because co-expressed pathways share
#' variance with the phenotype. Standard approaches report all significant
#' associations without distinguishing independent signals from correlated
#' shadows of the same underlying biology.
#'
#' This function extends \code{\link{run_redundancy_selection}} to
#' \strong{sample space}: rather than conditioning on gene membership
#' (as in the gene-list pipeline), it conditions on pathway activity
#' scores across samples. At each step, the pathway whose activity
#' contributes the most \emph{new} information about the phenotype
#' given the activity profiles of all previously selected pathways
#' is chosen.
#'
#' @section Algorithm:
#' \enumerate{
#'   \item Filter \code{assoc_results} to pathways with \code{padj <=
#'     alpha}.
#'   \item Re-compute (or accept pre-computed) pathway activity scores
#'     for each significant pathway using the same method as
#'     \code{\link{run_info_assoc}}.
#'   \item Discretise all activity scores using equal-frequency binning.
#'   \item Select the pathway with the highest
#'     I(activity ; phenotype) as the first pathway.
#'   \item For each subsequent step, compute conditional MI via
#'     \strong{double residualisation}: regress both the candidate
#'     pathway's activity and the phenotype on all previously selected
#'     pathway activities. Discretise both residual vectors and compute
#'     MI between them. A pathway redundant with already-selected ones
#'     will have its activity largely explained by the regression; its
#'     residual MI with the phenotype residual collapses toward zero.
#'   \item Select the pathway with the highest residual MI exceeding
#'     \code{min_gain}. Update the selected set. Stop when no pathway
#'     clears the threshold or \code{max_pathways} is reached.
#' }
#'
#' @section Double residualisation:
#' The conditional MI at each step is estimated by orthogonalising both
#' the candidate activity \emph{and} the phenotype against all selected
#' pathway activities, then measuring MI between the two residual
#' vectors. This is the non-parametric analog of partial correlation:
#' it asks whether the candidate explains phenotypic variance that no
#' previously selected pathway accounts for. Crucially, by
#' residualising \emph{both} sides, the method avoids the suppressor
#' artifact that arises when only the candidate is orthogonalised (in
#' which case its residual can pick up noise patterns correlated with
#' the original phenotype, inflating conditional MI above marginal MI).
#' With double residualisation, conditional MI is approximately bounded
#' by marginal MI, so \code{redundancy_ratio} stays close to \eqn{[0, 1]}.
#'
#' @section Redundancy ratio interpretation:
#' As in \code{\link{run_redundancy_selection}}, the redundancy ratio
#' \code{conditional_mi / marginal_mi} is the primary interpretive output.
#' A ratio near 1 means the pathway's phenotypic association is
#' independent of all previously selected pathways; it explains variance
#' in the phenotype that no prior selection accounts for. A ratio near 0
#' means the pathway's association is entirely mediated by previously
#' selected pathways; it is a correlated shadow, not an independent
#' signal. With double residualisation, the ratio is approximately
#' bounded by 1; values slightly exceeding 1 can occur due to
#' finite-sample discretisation noise but are not systematically
#' inflated.
#'
#' @section Connection to run_info_assoc:
#' Pass the output of \code{\link{run_info_assoc}} directly as
#' \code{assoc_results}. The \code{MI_bits} column is used as the
#' marginal MI for step 1 scoring and \code{min_gain} auto-selection,
#' avoiding redundant computation of marginal MI.
#'
#' @param assoc_results Data frame returned by \code{\link{run_info_assoc}}.
#' @param expr Numeric matrix of expression values, \strong{samples x genes}.
#'   Must be the same matrix used in \code{\link{run_info_assoc}}.
#' @param phenotype Vector of per-sample phenotype labels or values.
#'   Must be the same vector used in \code{\link{run_info_assoc}}.
#' @param gene_sets Named list of character vectors defining gene set
#'   membership. Must be the same list used in \code{\link{run_info_assoc}}.
#' @param alpha Numeric. Adjusted p-value threshold for filtering
#'   \code{assoc_results}. Applied to the \code{padj} column. Default
#'   \code{0.05}.
#' @param score Character. Activity scoring method: \code{"mean_z"}
#'   (default) or \code{"pc1"}. Should match the method used in
#'   \code{\link{run_info_assoc}}.
#' @param nbins \code{NULL} (default, auto-select) or positive integer.
#'   Bin count for equal-frequency discretisation of activity scores and
#'   phenotype. Should match the value used in \code{\link{run_info_assoc}}
#'   for consistency.
#' @param min_gain \code{NULL} (default, auto-select) or positive numeric.
#'   Minimum conditional MI in bits to select a pathway. When \code{NULL},
#'   auto-selected as the maximum of a permutation noise floor (95th
#'   percentile of max CMI under phenotype permutation) and a relevance
#'   threshold (\code{fraction * max(MI_bits)}).
#' @param fraction Numeric in (0, 1). Relevance fraction for auto
#'   \code{min_gain}. Default \code{0.15}.
#' @param n_perm_gain Integer. Permutations for auto \code{min_gain} noise
#'   floor. Default \code{200L}.
#' @param max_pathways Integer. Hard cap on selected pathways. Default
#'   \code{20L}.
#' @param seed Optional integer passed to \code{\link{set.seed}}.
#'
#' @return A data frame with one row per selected pathway in selection
#'   order. Columns:
#'   \describe{
#'     \item{\code{step}}{Selection order.}
#'     \item{\code{pathway}}{Pathway name.}
#'     \item{\code{set_size}}{Gene set size after universe intersection.}
#'     \item{\code{marginal_mi}}{I(activity ; phenotype) in bits.
#'       Taken from \code{assoc_results$MI_bits}.}
#'     \item{\code{conditional_mi}}{I(activity ; phenotype | selected
#'       activities) in bits. Equals \code{marginal_mi} at step 1.}
#'     \item{\code{redundancy_ratio}}{\code{conditional_mi /
#'       marginal_mi}. Near 1 = independent; near 0 = redundant.}
#'     \item{\code{cumulative_bits}}{Running total of
#'       \code{conditional_mi}. Interpretable as total unique
#'       information explained by the selected set.}
#'     \item{\code{p_value}}{Original \code{p_value} from
#'       \code{assoc_results}.}
#'     \item{\code{padj}}{Original \code{padj} from
#'       \code{assoc_results}.}
#'     \item{\code{z_score}}{Original \code{z_score} from
#'       \code{assoc_results}.}
#'   }
#'   Attributes: \code{min_gain}, \code{n_sig}, \code{alpha},
#'   \code{score}, \code{nbins_used}.
#'
#' @seealso
#'   \code{\link{run_info_assoc}} for the upstream association step,
#'   \code{\link{run_redundancy_selection}} for the gene-space version,
#'   \code{\link{run_gsea_redundancy_selection}} for the GSEA version,
#'   \code{\link{plot_gains}} to visualise the selection profile.
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
#'
#' # Phenotype driven by two independent latent factors
#' phenotype <- 0.7 * factor_A + 0.5 * factor_B + rnorm(n_samples, sd = 0.4)
#'
#' gene_sets <- list(
#'   pathway_A  = genes[1:30],
#'   pathway_A2 = genes[5:35],    # redundant with A
#'   pathway_B  = genes[200:230], # independent signal
#'   pathway_C  = genes[100:130]  # noise
#' )
#'
#' # Step 1: find associated pathways
#' assoc <- run_info_assoc(
#'   expr      = expr,
#'   phenotype = phenotype,
#'   gene_sets = gene_sets,
#'   n_perm    = 1000,
#'   seed      = 42
#' )
#'
#' # Step 2: select non-redundant associated pathways
#' sel <- run_assoc_redundancy_selection(
#'   assoc_results = assoc,
#'   expr          = expr,
#'   phenotype     = phenotype,
#'   gene_sets     = gene_sets,
#'   alpha         = 0.05,
#'   seed          = 42
#' )
#' print(sel)
#'
#' @export
run_assoc_redundancy_selection <- function(
    assoc_results,
    expr,
    phenotype,
    gene_sets,
    alpha        = 0.05,
    score        = c("mean_z", "pc1"),
    nbins        = NULL,
    min_gain     = NULL,
    fraction     = 0.15,
    n_perm_gain  = 200L,
    max_pathways = 20L,
    seed         = NULL
) {
  ## --- Validation ------------------------------------------------------------
  if (!is.data.frame(assoc_results))
    stop("`assoc_results` must be a data frame from run_info_assoc().")
  required_cols <- c("set", "MI_bits", "p_value", "padj", "z_score")
  missing_cols  <- setdiff(required_cols, colnames(assoc_results))
  if (length(missing_cols) > 0)
    stop("Missing columns in `assoc_results`: ",
         paste(missing_cols, collapse = ", "),
         ". Pass the direct output of run_info_assoc().")
  if (!is.matrix(expr) || !is.numeric(expr))
    stop("`expr` must be a numeric matrix (samples x genes).")
  if (length(phenotype) != nrow(expr))
    stop("`phenotype` length must equal nrow(expr).")
  if (!is.list(gene_sets) || is.null(names(gene_sets)))
    stop("`gene_sets` must be a named list.")

  .check_matrix_orientation(expr)
  .check_na_expr(expr)
  .check_sample_size(expr, phenotype)
  expr <- .remove_zero_variance(expr)

  score <- match.arg(score)
  if (!is.null(seed)) set.seed(seed)

  ## --- Filter to significant associations ------------------------------------
  is_sig <- !is.na(assoc_results$padj) & assoc_results$padj <= alpha
  n_sig  <- sum(is_sig)

  if (n_sig == 0L)
    stop(
      "No significant pathways at padj <= ", alpha, ". ",
      "Consider relaxing alpha or checking run_info_assoc() results."
    )

  message(sprintf("%d / %d pathways significant at padj <= %s.",
                  n_sig, nrow(assoc_results), alpha))

  sig_sets  <- assoc_results$set[is_sig]
  sig_assoc <- assoc_results[is_sig, ]

  ## Keep only gene sets that are both significant and present in gene_sets
  available <- intersect(sig_sets, names(gene_sets))
  if (length(available) < length(sig_sets))
    message(sprintf(
      "%d significant pathway(s) not found in `gene_sets`- dropped.",
      length(sig_sets) - length(available)
    ))
  if (length(available) == 0L)
    stop("No significant pathways found in `gene_sets`.")

  sig_sets  <- available
  sig_assoc <- sig_assoc[sig_assoc$set %in% available, ]
  gene_sets <- gene_sets[sig_sets]

  ## --- Universe and expression -----------------------------------------------
  universe <- colnames(expr)
  gene_sets <- lapply(gene_sets, intersect, universe)

  ## --- Binning ---------------------------------------------------------------
  n <- nrow(expr)
  is_continuous_pheno <- is.numeric(phenotype) && length(unique(phenotype)) > 2L

  if (is.null(nbins)) {
    nbins_y_hint <- if (is_continuous_pheno) NULL else
      length(unique(phenotype))
    nbins <- .auto_nbins(n, nbins_y = nbins_y_hint)
    message("Auto-selected nbins = ", nbins, ".")
  } else {
    nbins <- as.integer(nbins)
    if (nbins < 2L) stop("`nbins` must be >= 2.")
  }

  ## --- Pathway activity scores -----------------------------------------------
  Z <- scale(expr)

  score_fun <- switch(
    score,
    mean_z = function(G) rowMeans(Z[, G, drop = FALSE]),
    pc1    = function(G) {
      prcomp(Z[, G, drop = FALSE], center = FALSE, scale. = FALSE)$x[, 1L]
    }
  )

  activity <- lapply(gene_sets, score_fun)

  ## --- Discretise activities and phenotype -----------------------------------
  activity_disc <- lapply(activity,
                          function(a) .discretize_equalfreq(a, nbins))
  y_disc <- if (is.numeric(phenotype) && length(unique(phenotype)) > 2L) {
    .discretize_equalfreq(phenotype, nbins)
  } else {
    as.integer(as.factor(phenotype))
  }

  ## --- Marginal MI from assoc_results (avoid recomputation) -----------------
  mi_marginal        <- setNames(sig_assoc$MI_bits, sig_assoc$set)
  mi_marginal        <- mi_marginal[sig_sets]   # ensure consistent order

  ## --- min_gain auto-selection -----------------------------------------------
  if (is.null(min_gain)) {
    min_gain <- .auto_min_gain_assoc(
      activity_disc = activity_disc,
      y_disc        = y_disc,
      mi_marginal   = mi_marginal,
      n_perm        = n_perm_gain,
      fraction      = fraction
    )
    message(sprintf(
      "Auto-selected min_gain = %.5f bits (noise floor = %.5f, relevance = %.5f [%.0f%% of leading MI]).",
      as.numeric(min_gain),
      attr(min_gain, "noise_floor"),
      attr(min_gain, "relevance_threshold"),
      fraction * 100
    ))
  }

  ## --- Empty result template -------------------------------------------------
  empty_result <- data.frame(
    step             = integer(0),
    pathway          = character(0),
    set_size         = integer(0),
    marginal_mi      = numeric(0),
    conditional_mi   = numeric(0),
    redundancy_ratio = numeric(0),
    cumulative_bits  = numeric(0),
    p_value          = numeric(0),
    padj             = numeric(0),
    z_score          = numeric(0),
    stringsAsFactors = FALSE
  )

  ## --- First selection -------------------------------------------------------
  best_idx  <- which.max(mi_marginal)
  best_name <- names(mi_marginal)[best_idx]
  best_gain <- mi_marginal[best_idx]

  if (best_gain < min_gain) {
    message("No pathway exceeds min_gain = ", round(min_gain, 5),
            " bits. Returning empty selection.")
    return(empty_result)
  }

  selected   <- best_name
  remaining  <- setdiff(sig_sets, best_name)
  phenotype_num <- as.numeric(phenotype)

  selected_mat <- matrix(activity[[best_name]], ncol = 1L,
                         dimnames = list(NULL, best_name))
  cum_bits <- best_gain

  assoc_lookup <- setNames(
    split(sig_assoc, seq_len(nrow(sig_assoc))),
    sig_assoc$set
  )

  log_rows       <- vector("list", min(max_pathways, length(sig_sets)))
  log_rows[[1L]] <- data.frame(
    step             = 1L,
    pathway          = best_name,
    set_size         = length(gene_sets[[best_name]]),
    marginal_mi      = best_gain,
    conditional_mi   = best_gain,
    redundancy_ratio = 1,
    cumulative_bits  = cum_bits,
    p_value          = assoc_lookup[[best_name]]$p_value,
    padj             = assoc_lookup[[best_name]]$padj,
    z_score          = assoc_lookup[[best_name]]$z_score,
    stringsAsFactors = FALSE
  )

  ## --- Greedy forward steps --------------------------------------------------
  # Double residualisation: orthogonalise BOTH candidate activity and phenotype
  # against all selected pathway activities, then measure MI between residuals.
  # This is the non-parametric analog of partial correlation and avoids the
  # suppressor artifact that arises from single-sided residualisation.
  #
  # IMPORTANT: after regression, both residuals are continuous even if the
  # original phenotype was binary. The joint table for residual MI is
  # nbins_resid x nbins_resid (symmetric), not nbins x nbins_y. We must
  # use a bin count appropriate for the continuous-continuous case to avoid
  # sparse-cell MI inflation (which causes all ratios >> 1).

  nbins_resid <- .auto_nbins(n, nbins_y = NULL)

  n_steps <- min(max_pathways, length(sig_sets))
  if (n_steps < 2L) {
    res <- do.call(rbind, Filter(Negate(is.null), log_rows))
    attr(res, "min_gain")   <- as.numeric(min_gain)
    attr(res, "n_sig")      <- n_sig
    attr(res, "alpha")      <- alpha
    attr(res, "score")      <- score
    attr(res, "nbins_used") <- nbins
    return(res)
  }

  for (step in seq(2L, n_steps)) {
    if (length(remaining) == 0L) break

    # Phenotype residuals (shared across all candidates at this step)
    y_resid <- residuals(lm(phenotype_num ~ selected_mat))
    y_resid_disc <- .discretize_equalfreq(y_resid, nbins_resid)

    cmi_vals <- vapply(
      remaining,
      function(nm) {
        a_resid <- residuals(lm(activity[[nm]] ~ selected_mat))
        a_disc  <- .discretize_equalfreq(a_resid, nbins_resid)
        .mutual_information(a_disc, y_resid_disc)
      },
      numeric(1L)
    )

    best_idx  <- which.max(cmi_vals)
    best_name <- remaining[best_idx]
    best_cmi  <- cmi_vals[best_idx]

    if (best_cmi < min_gain) {
      message(sprintf(
        "Stopping at step %d: best CMI = %.5f bits < min_gain = %.5f.",
        step, best_cmi, as.numeric(min_gain)
      ))
      break
    }

    selected  <- c(selected, best_name)
    remaining <- setdiff(remaining, best_name)
    cum_bits  <- cum_bits + best_cmi

    selected_mat <- cbind(selected_mat, activity[[best_name]])
    colnames(selected_mat)[ncol(selected_mat)] <- best_name

    log_rows[[step]] <- data.frame(
      step             = as.integer(step),
      pathway          = best_name,
      set_size         = length(gene_sets[[best_name]]),
      marginal_mi      = mi_marginal[[best_name]],
      conditional_mi   = best_cmi,
      redundancy_ratio = best_cmi / mi_marginal[[best_name]],
      cumulative_bits  = cum_bits,
      p_value          = assoc_lookup[[best_name]]$p_value,
      padj             = assoc_lookup[[best_name]]$padj,
      z_score          = assoc_lookup[[best_name]]$z_score,
      stringsAsFactors = FALSE
    )
  }

  res <- do.call(rbind, Filter(Negate(is.null), log_rows))

  attr(res, "min_gain")    <- as.numeric(min_gain)
  attr(res, "n_sig")       <- n_sig
  attr(res, "alpha")       <- alpha
  attr(res, "score")       <- score
  attr(res, "nbins_used")  <- nbins

  res
}


## ---- Internal helpers -------------------------------------------------------

#' Sample-space conditional mutual information
#'
#' Computes I(x ; y | z) where all three vectors are discrete, by
#' stratifying on z and computing z-weighted MI within each stratum.
#' This is the sample-space analog of .cmi_binary() in the gene-space
#' pipeline.
#'
#' @param x Integer vector. Discretised pathway activity.
#' @param y Integer vector. Discretised phenotype.
#' @param z Integer vector. Discretised conditioning variable
#'   (mean activity of selected pathways).
#'
#' @return Scalar CMI in bits (>= 0).
#' @keywords internal
#' @noRd
.cmi_sample_space <- function(x, y, z) {
  n      <- length(x)
  z_vals <- unique(z)
  cmi    <- 0
  for (zv in z_vals) {
    idx <- which(z == zv)
    if (length(idx) < 2L) next
    pz  <- length(idx) / n
    cmi <- cmi + pz * .mutual_information(x[idx], y[idx])
  }
  cmi
}


#' Auto-select min_gain threshold for sample-space CMI selection
#'
#' Sample-space version of .auto_min_gain(). Estimates the noise floor
#' by permuting the phenotype labels and computing the maximum marginal
#' MI across all candidate pathways, taking the 95th percentile across
#' permutations.
#'
#' @param activity_disc Named list of discretised activity vectors.
#' @param y_disc Integer vector of discretised phenotype.
#' @param mi_marginal Named numeric vector of marginal MI values.
#' @param n_perm Integer. Permutations for noise floor.
#' @param fraction Numeric. Relevance fraction.
#'
#' @return Scalar threshold with attributes noise_floor and
#'   relevance_threshold.
#' @keywords internal
#' @noRd
.auto_min_gain_assoc <- function(activity_disc, y_disc,
                                 mi_marginal,
                                 n_perm   = 200L,
                                 fraction = 0.15) {
  null_mi <- vapply(seq_len(n_perm), function(b) {
    y_perm <- sample(y_disc)
    max(vapply(
      activity_disc,
      function(ad) .mutual_information(ad, y_perm),
      numeric(1L)
    ))
  }, numeric(1L))

  noise_floor         <- as.numeric(quantile(null_mi, 0.95))
  relevance_threshold <- fraction * max(mi_marginal)
  threshold           <- max(noise_floor, relevance_threshold)

  attr(threshold, "noise_floor")         <- noise_floor
  attr(threshold, "relevance_threshold") <- relevance_threshold
  threshold
}
