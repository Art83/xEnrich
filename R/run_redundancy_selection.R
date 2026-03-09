#' Greedy non-redundant pathway selection via conditional mutual information
#'
#' Addresses a fundamental limitation of enrichment analysis: when many
#' pathways are nominally significant, most apparent signal is redundant-
#' the same genes appearing in multiple overlapping sets. Standard approaches
#' (Jaccard clustering, representative selection, ridge penalisation) either
#' discard pathways arbitrarily or produce coefficients that are difficult to
#' interpret biologically.
#'
#' This function builds a minimal non-redundant set by greedy forward
#' selection on conditional mutual information (CMI). At each step, the
#' pathway contributing the most \emph{new} information about the gene list,
#' given all previously selected pathways,is chosen. Selection stops when
#' no remaining pathway exceeds \code{min_gain} bits.
#'
#' The \code{redundancy_ratio} column (\code{conditional_mi / marginal_mi})
#' is the primary interpretive output: a value near 1 means the pathway's
#' signal is independent of all previously selected sets; a value near 0
#' means its signal is entirely explained by earlier selections.
#'
#' @section Algorithm:
#' \enumerate{
#'   \item Compute marginal MI I(gene_list ; pathway_j) for all pathways.
#'   \item Select the pathway with the highest marginal MI as S1.
#'   \item For all remaining pathways, compute
#'     I(gene_list ; pathway_j | union(S1, ..., St)),
#'     where the conditioning set is binary (in/out of the gene union of
#'     all selected pathways so far).
#'   \item Select the pathway with the highest CMI exceeding \code{min_gain}.
#'   \item Repeat until the gain criterion fails or \code{max_pathways} is
#'     reached.
#' }
#'
#' Collapsing selected pathways into a binary union variable (rather than
#' conditioning on each separately) keeps CMI estimation tractable and avoids
#' sparse-cell problems as the selected set grows.
#'
#' @section Conditioning artifact, synergy:
#' When a large pathway is selected first, the z=0 stratum (genes outside
#' it) may be enriched for signal from small independent pathways, causing
#' their CMI to exceed their marginal MI (\code{redundancy_ratio > 1}).
#' This is mathematically valid (synergy) but inflated relative to the true
#' independent contribution. Conservative \code{min_size} filtering and a
#' well-defined \code{min_gain} threshold mitigate this in practice.
#'
#' @section min_gain selection:
#' When \code{min_gain = NULL} (default), the threshold is auto-selected as
#' the maximum of:
#' \itemize{
#'   \item A \strong{noise floor}: 95th percentile of max marginal MI under
#'     gene-list permutation, ensures selection is above chance.
#'   \item A \strong{relevance threshold}: \code{fraction * max(marginal MI)}
#'     ensures each selected pathway contributes a meaningful fraction of
#'     the dominant signal.
#' }
#' The relevance threshold is the binding constraint in most real datasets.
#' Adjust \code{fraction} (default \code{0.15}) to control stringency.
#'
#' @section Integration with run_info_enrichment():
#' Pass the \code{info_bits} column from [run_info_enrichment()] as
#' \code{initial_scores} to skip redundant marginal MI computation.
#'
#' @param gene_list Character vector of query gene symbols.
#' @param gene_sets Named list of character vectors defining gene set
#'   membership.
#' @param universe Character vector of background genes.
#' @param initial_scores Optional named numeric vector of pre-computed
#'   marginal MI scores (e.g. \code{info_bits} from
#'   [run_info_enrichment()]). Names must match \code{names(gene_sets)}.
#'   Missing entries are computed internally.
#' @param min_gain \code{NULL} (default, auto-select) or a positive numeric.
#'   Minimum conditional MI in bits required to select a pathway. See the
#'   min_gain selection section.
#' @param fraction Numeric in (0, 1). Controls the relevance component of
#'   auto \code{min_gain}: a pathway must contribute at least this fraction
#'   of the leading signal. Default \code{0.15}. Ignored when \code{min_gain}
#'   is supplied explicitly.
#' @param n_perm_gain Integer. Permutations for auto \code{min_gain} noise
#'   floor. Default \code{200L}. Ignored when \code{min_gain} is supplied.
#' @param max_pathways Integer. Hard cap on selected pathways. Default
#'   \code{20L}.
#' @param min_size Integer. Minimum gene set size after universe intersection.
#'   Default \code{5L}.
#' @param max_size Integer. Maximum gene set size. Default \code{Inf}.
#' @param seed Optional integer passed to [set.seed()].
#'
#' @return A data frame with one row per selected pathway in selection order,
#'   or a zero-row data frame if nothing exceeds \code{min_gain}. Columns:
#'   \describe{
#'     \item{\code{step}}{Selection order (1 = first selected).}
#'     \item{\code{pathway}}{Pathway name.}
#'     \item{\code{set_size}}{Members after universe intersection.}
#'     \item{\code{marginal_mi}}{Unconditional I(gene_list ; pathway) in
#'       bits.}
#'     \item{\code{conditional_mi}}{I(gene_list ; pathway | union of
#'       previously selected) in bits. Equals \code{marginal_mi} at step 1.}
#'     \item{\code{redundancy_ratio}}{\code{conditional_mi / marginal_mi}.
#'       Near 1: pathway is independent of all prior selections. Near 0:
#'       pathway is redundant. Values > 1 indicate synergy (see Details).}
#'     \item{\code{cumulative_bits}}{Running total of \code{conditional_mi}.
#'       Interpretable as total unique information explained by the selected
#'       set.}
#'   }
#'
#' @seealso [run_info_enrichment()] for marginal MI enrichment scores,
#'   [run_info_assoc()] for sample-level phenotype association,
#'   [plot_gains()] to visualise the selection profile.
#'
#' @references
#' Paninski, L. (2003). Estimation of entropy and mutual information.
#' \emph{Neural Computation}, 15(6), 1191–1253.
#'
#' @examples
#' set.seed(1)
#' universe <- paste0("G", 1:1000)
#'
#' # A and B are near-redundant (55/60 gene overlap)
#' # C is an independent signal; D is noise
#' gene_sets <- list(
#'   pathway_A = universe[1:60],
#'   pathway_B = c(universe[1:55], sample(universe[61:1000], 5)),
#'   pathway_C = universe[500:560],
#'   pathway_D = universe[200:230]
#' )
#'
#' # gene_list is enriched in A/B region and C, not elsewhere
#' gene_list <- unique(c(
#'   sample(universe[1:60],    45),
#'   sample(universe[500:560], 15),
#'   sample(universe[61:499],   3)
#' ))
#'
#' # Step 1: marginal enrichment: A, B, and C all look significant
#' mi_res <- run_info_enrichment(
#'   gene_list = gene_list,
#'   gene_sets = gene_sets,
#'   universe  = universe,
#'   n_perm    = 500,
#'   seed      = 1
#' )
#' print(mi_res[, c("set", "set_size", "info_bits", "emp_p", "padj_emp")])
#'
#' # Step 2: greedy selection: B is revealed as redundant with A
#' sel <- run_redundancy_selection(
#'   gene_list      = gene_list,
#'   gene_sets      = gene_sets,
#'   universe       = universe,
#'   initial_scores = setNames(mi_res$info_bits, mi_res$set),
#'   seed           = 1
#' )
#' print(sel)
#'
#' # Fraction of total marginal signal captured by selected pathways
#' cat(sprintf(
#'   "%d pathways capture %.0f%% of total marginal MI\n",
#'   nrow(sel),
#'   100 * max(sel$cumulative_bits) / sum(mi_res$info_bits)
#' ))
#'
#' plot_gains(sel, min_gain = attr(sel, "min_gain"))
#'
#' @export
run_redundancy_selection <- function(
    gene_list,
    gene_sets,
    universe,
    initial_scores = NULL,
    min_gain       = NULL,
    fraction       = 0.15,
    n_perm_gain    = 200L,
    max_pathways   = 20L,
    min_size       = 5L,
    max_size       = Inf,
    seed           = NULL
) {
  # --- Validation -------------------------------------------------------------
  if (!is.character(gene_list))
    stop("`gene_list` must be a character vector.")
  if (!is.list(gene_sets) || is.null(names(gene_sets)))
    stop("`gene_sets` must be a named list.")
  if (!is.null(min_gain) && min_gain <= 0)
    stop("`min_gain` must be a positive number.")
  if (fraction <= 0 || fraction >= 1)
    stop("`fraction` must be in (0, 1).")

  if (!is.null(seed)) set.seed(seed)

  # --- Intersect and filter ---------------------------------------------------
  gene_list <- intersect(gene_list, universe)
  gene_sets <- lapply(gene_sets, intersect, universe)

  set_sizes <- lengths(gene_sets)
  keep      <- set_sizes >= min_size & set_sizes <= max_size
  gene_sets <- gene_sets[keep]

  n_dropped <- sum(!keep)
  if (n_dropped > 0L)
    message(n_dropped, " gene set(s) dropped by size filter.")
  if (length(gene_sets) == 0L)
    stop("No gene sets remain after filtering.")

  # --- Binary vectors and membership matrix -----------------------------------
  x_list <- as.integer(universe %in% gene_list)
  M      <- .make_membership_matrix(gene_sets, universe)

  # --- Marginal MI ------------------------------------------------------------
  mi_marginal <- vapply(names(gene_sets), function(nm) {
    if (!is.null(initial_scores) && !is.na(initial_scores[nm]))
      return(as.numeric(initial_scores[nm]))
    .mutual_information(x_list, as.integer(M[, nm]))
  }, numeric(1L))

  # --- min_gain: auto-select or validate -------------------------------------
  if (is.null(min_gain)) {
    min_gain <- .auto_min_gain(
      x_list      = x_list,
      M           = M,
      mi_marginal = mi_marginal,
      n_perm      = n_perm_gain,
      fraction    = fraction
    )
    message(
      "Auto-selected min_gain = ", round(min_gain, 5), " bits",
      " (noise floor = ", round(attr(min_gain, "noise_floor"), 5),
      ", relevance threshold = ",
      round(attr(min_gain, "relevance_threshold"), 5),
      " [", fraction * 100, "% of leading MI])."
    )
  }

  # --- Empty-result template --------------------------------------------------
  empty_result <- data.frame(
    step             = integer(0),
    pathway          = character(0),
    set_size         = integer(0),
    marginal_mi      = numeric(0),
    conditional_mi   = numeric(0),
    redundancy_ratio = numeric(0),
    cumulative_bits  = numeric(0),
    stringsAsFactors = FALSE
  )

  # --- First selection --------------------------------------------------------
  best_idx  <- which.max(mi_marginal)
  best_name <- names(mi_marginal)[best_idx]
  best_gain <- mi_marginal[best_idx]

  if (best_gain < min_gain) {
    message("No pathway exceeds min_gain = ", round(min_gain, 5),
            " bits. Returning empty selection.")
    return(empty_result)
  }

  selected  <- best_name
  remaining <- setdiff(names(gene_sets), best_name)
  z_union   <- as.integer(M[, best_name])
  cum_bits  <- best_gain

  log_rows        <- vector("list", min(max_pathways, length(gene_sets)))
  log_rows[[1L]]  <- data.frame(
    step             = 1L,
    pathway          = best_name,
    set_size         = as.integer(Matrix::colSums(M)[best_name]),
    marginal_mi      = best_gain,
    conditional_mi   = best_gain,
    redundancy_ratio = 1,
    cumulative_bits  = cum_bits,
    stringsAsFactors = FALSE
  )

  # --- Greedy forward steps ---------------------------------------------------
  for (step in seq(2L, min(max_pathways, length(gene_sets)))) {
    if (length(remaining) == 0L) break

    cmi_vals  <- vapply(
      remaining,
      function(nm) .cmi_binary(x_list, as.integer(M[, nm]), z_union),
      numeric(1L)
    )

    best_idx  <- which.max(cmi_vals)
    best_name <- remaining[best_idx]
    best_cmi  <- cmi_vals[best_idx]

    if (best_cmi < min_gain) {
      message(
        "Stopping at step ", step,
        ": best conditional gain = ", round(best_cmi, 5),
        " bits < min_gain = ", round(min_gain, 5), "."
      )
      break
    }

    selected  <- c(selected, best_name)
    remaining <- setdiff(remaining, best_name)
    z_union   <- pmax(z_union, as.integer(M[, best_name]))
    cum_bits  <- cum_bits + best_cmi

    log_rows[[step]] <- data.frame(
      step             = as.integer(step),
      pathway          = best_name,
      set_size         = as.integer(Matrix::colSums(M)[best_name]),
      marginal_mi      = mi_marginal[best_name],
      conditional_mi   = best_cmi,
      redundancy_ratio = best_cmi / mi_marginal[best_name],
      cumulative_bits  = cum_bits,
      stringsAsFactors = FALSE
    )
  }

  res <- do.call(rbind, Filter(Negate(is.null), log_rows))

  # Attach min_gain as attribute so plot_gains() can use it without
  # requiring the user to pass it separately
  attr(res, "min_gain") <- as.numeric(min_gain)
  res
}


#' Plot conditional MI gains from run_redundancy_selection
#'
#' Visualises how much new information each selected pathway contributes,
#' with the \code{min_gain} threshold marked. Intended as a diagnostic for
#' threshold calibration and as a publication figure showing the contrast
#' between marginal and conditional signal.
#'
#' Bars are coloured by \code{redundancy_ratio}: blue indicates an
#' independent pathway (ratio near 1); grey indicates partial redundancy
#' (ratio < 1). The dashed line shows where selection stopped.
#'
#' @param sel Data frame returned by [run_redundancy_selection()].
#' @param min_gain Numeric. The threshold used in the selection run.
#'   If \code{NULL}, read from \code{attr(sel, "min_gain")} when available.
#'
#' @return Invisibly returns \code{sel}.
#' @export
plot_gains <- function(sel, min_gain = NULL) {
  if (nrow(sel) == 0L) {
    message("Nothing selected-nothing to plot.")
    return(invisible(sel))
  }

  if (is.null(min_gain)) min_gain <- attr(sel, "min_gain")

  # Colour by redundancy_ratio: independent = blue, redundant = grey
  pal   <- colorRampPalette(c("#cccccc", "#2c7bb6"))(100)
  r     <- pmin(pmax(sel$redundancy_ratio, 0), 1)
  cols  <- pal[pmax(1L, as.integer(r * 99) + 1L)]

  op <- par(mar = c(8, 4.5, 3, 1))
  on.exit(par(op))

  bp <- barplot(
    sel$conditional_mi,
    names.arg = sel$pathway,
    col       = cols,
    las       = 2,
    border    = NA,
    ylab      = "Conditional MI (bits)",
    main      = "Information gain per selection step",
    ylim      = c(0, max(sel$conditional_mi) * 1.15)
  )

  labels <- ifelse(sel$redundancy_ratio > 1,
                   sprintf("r=%.2f*", sel$redundancy_ratio),
                   sprintf("r=%.2f",  sel$redundancy_ratio))

  text(
    x      = bp,
    y      = sel$conditional_mi + max(sel$conditional_mi) * 0.03,
    labels = labels,
    cex    = 0.75,
    col    = "grey30"
  )

  if (!is.null(min_gain))
    abline(h = min_gain, lty = 2, col = "firebrick", lwd = 1.5)

  legend(
    "topright",
    legend = c("Independent (ratio \u2248 1)",
               "Redundant (ratio \u2248 0)",
               "min_gain threshold",
               "* synergy: CMI > marginal MI"),
    fill   = c("#2c7bb6", "#cccccc", NA, NA),
    border = c(NA, NA, NA, NA),
    lty    = c(NA, NA, 2, NA),
    col    = c(NA, NA, "firebrick", NA),
    bty    = "n",
    cex    = 0.85
  )

  invisible(sel)
}

#' Compute conditional mutual information I(x ; y | z)
#'
#' Stratifies observations by the binary conditioning variable \code{z}
#' (in/out of the union of already-selected pathway genes) and returns the
#' z-weighted average of MI within each stratum. Uses the
#' 0 * log2(0) = 0 convention via [.mutual_information()].
#'
#' Strata with fewer than 2 observations are skipped, they cannot
#' contribute a stable MI estimate.
#'
#' @param x Integer vector (binary: gene in query list).
#' @param y Integer vector (binary: gene in candidate pathway).
#' @param z Integer vector (binary: gene in union of selected pathways).
#'
#' @return Scalar conditional MI in bits (>= 0).
#' @keywords internal
#' @noRd
.cmi_binary <- function(x, y, z) {
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


#' Automatically select the min_gain threshold for pathway selection
#'
#' Combines two criteria and takes the stricter (larger) of the two:
#' \enumerate{
#'   \item \strong{Noise floor}: 95th percentile of the maximum marginal MI
#'     across pathways under gene-list permutation. Any threshold below this
#'     value is indistinguishable from random.
#'   \item \strong{Relevance threshold}: \code{fraction * max(mi_marginal)}.
#'     A pathway must contribute at least this fraction of the leading
#'     signal to be worth selecting.
#' }
#'
#' The relevance threshold prevents selection of pathways that, while above
#' noise, are negligible relative to the dominant biology. This is the
#' criterion that correctly suppresses the near-redundant pathway in typical
#' multi-pathway enrichment scenarios.
#'
#' @param x_list Integer vector (binary gene-list membership over universe).
#' @param M Sparse membership matrix (genes x pathways).
#' @param mi_marginal Named numeric vector of pre-computed marginal MI values.
#' @param n_perm Integer. Permutations for the noise floor estimate.
#' @param fraction Numeric in (0, 1). Fraction of leading MI defining the
#'   relevance threshold. Default \code{0.15}.
#'
#' @return Scalar threshold with attributes \code{noise_floor} and
#'   \code{relevance_threshold} recording the two components.
#' @keywords internal
#' @noRd
.auto_min_gain <- function(x_list, M, mi_marginal,
                           n_perm   = 200L,
                           fraction = 0.15) {
  # Noise floor: max null MI across pathways, 95th percentile
  null_mi <- vapply(seq_len(n_perm), function(b) {
    x_perm <- sample(x_list)
    max(vapply(
      colnames(M),
      function(nm) .mutual_information(x_perm, as.integer(M[, nm])),
      numeric(1L)
    ))
  }, numeric(1L))

  noise_floor          <- as.numeric(quantile(null_mi, 0.95))
  relevance_threshold  <- fraction * max(mi_marginal)
  threshold            <- max(noise_floor, relevance_threshold)

  attr(threshold, "noise_floor")         <- noise_floor
  attr(threshold, "relevance_threshold") <- relevance_threshold
  threshold
}
