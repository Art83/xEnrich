#' Pathway-phenotype association via mutual information
#'
#' Rather than testing whether a gene list overlaps a pathway (as in
#' [run_info_enrichment()]), this function asks: \emph{how much does a
#' pathway's aggregate activity, measured across samples, explain a
#' phenotypic outcome?}
#'
#' Each pathway is summarised into a single per-sample activity score
#' (either the mean z-score across member genes, or the first principal
#' component). Mutual information (MI) between that score and the phenotype
#' vector is computed and compared against a label-permutation null,
#' yielding an empirical p-value that is free of distributional assumptions
#' about either the score or the phenotype.
#'
#' @section Scoring methods:
#' \describe{
#'   \item{\code{"mean_z"}}{Mean of z-scored expression across pathway
#'     members. Fast and robust; the default.}
#'   \item{\code{"pc1"}}{First principal component of z-scored expression.
#'     Captures co-expression structure within the pathway but is slower.
#'     PC1 sign is arbitrary, which does not affect MI (MI is
#'     sign-invariant) but would affect any signed downstream use.}
#' }
#'
#' @section Binning and MI bias:
#' MI is estimated from equal-frequency binned pathway scores. All
#' fixed-bin estimators carry a positive bias that shrinks with sample
#' size. Because both observed and null MI values are computed under
#' identical binning, this bias cancels in relative rankings and empirical
#' p-values — but absolute \code{MI_bits} values should not be compared
#' across datasets of different sizes.
#'
#' When \code{nbins = NULL} (default), bins are chosen automatically via a
#' composite rule: Rice rule (\eqn{2n^{1/3}}) capped to ensure an
#' expected frequency of at least 5 per \strong{joint} cell (not per
#' marginal bin). For continuous phenotype the joint table has
#' \eqn{nbins^2} cells, giving \eqn{nbins \leq \lfloor\sqrt{n/5}\rfloor}.
#' For binary phenotype the joint table has only \eqn{2 \times nbins}
#' cells, allowing larger \eqn{nbins}. A hard cap of 20 applies in
#' either case. Pass an explicit integer to override.
#'
#' @section No external dependencies:
#' Discretisation and MI estimation are handled by internal helpers, so
#' \pkg{infotheo} is not required.
#'
#' @param expr Numeric matrix of expression values, \strong{samples x genes}.
#'   Column names must be gene symbols.
#' @param phenotype Vector (numeric, factor, or character) of per-sample
#'   phenotype labels or values. Length must equal \code{nrow(expr)}.
#' @param gene_sets Named list of character vectors defining gene set
#'   membership (e.g. GO terms, pathways).
#' @param universe Character vector of background genes. Defaults to
#'   \code{colnames(expr)}.
#' @param score Character string: \code{"mean_z"} (default) or \code{"pc1"}.
#'   See the Scoring methods section.
#' @param nbins \code{NULL} (default, auto-select) or a positive integer
#'   number of equal-frequency bins for score discretisation. See the
#'   Binning section.
#' @param n_perm Positive integer. Phenotype-label permutations for the null
#'   distribution. Default \code{1000}.
#' @param min_size Integer. Minimum gene set size after universe intersection.
#'   Default \code{10}.
#' @param max_size Integer. Maximum gene set size after universe intersection.
#'   Default \code{Inf}. MI naturally down-weights large sets (high pathway
#'   coverage compresses marginal entropy), so this is a computational
#'   convenience rather than a statistical necessity.
#' @param nthreads Positive integer. Number of threads for the permutation
#'   loop via \code{parallel::mclapply}. Default \code{1L} (no parallelism).
#'   Values > 1 are silently reduced to 1 on Windows where forking is
#'   unavailable. On Linux/macOS the expression matrix is shared
#'   copy-on-write between workers, so memory overhead per additional thread
#'   is negligible.
#' @param seed Optional integer passed to [set.seed()]. \code{NULL} (default)
#'   leaves the RNG state untouched.
#'
#' @return A data frame with one row per gene set (sorted by \code{padj}),
#'   containing:
#'   \describe{
#'     \item{\code{set}}{Gene set name.}
#'     \item{\code{set_size}}{Members retained after universe intersection.}
#'     \item{\code{nbins_used}}{Actual bin count used for this run.}
#'     \item{\code{MI_bits}}{Observed mutual information in bits.}
#'     \item{\code{z_score}}{MI standardised against the permutation null.
#'       \code{NA} when null SD is zero.}
#'     \item{\code{p_value}}{Empirical p-value with continuity correction
#'       \eqn{(1 + B_{\geq obs}) / (n\_perm + 1)}.}
#'     \item{\code{padj}}{BH-adjusted \code{p_value}.}
#'     \item{\code{null_mean}, \code{null_sd}}{Mean and SD of the permutation
#'       null MI distribution.}
#'   }
#'
#' @seealso [run_info_enrichment()] for the gene-list (ORA) flavour,
#'   [run_enrichment()] for classical ORA/GSEA,
#'   [compute_pathway_activity()] to extract per-sample activity scores
#'   from selected pathways for downstream analysis.
#'
#' @references
#' Paninski, L. (2003). Estimation of entropy and mutual information.
#' \emph{Neural Computation}, 15(6), 1191-1253.
#'
#' @examples
#' set.seed(1)
#' n_samples <- 60
#' genes     <- paste0("G", 1:200)
#' expr      <- matrix(
#'   rnorm(n_samples * length(genes)),
#'   nrow     = n_samples,
#'   dimnames = list(NULL, genes)
#' )
#' phenotype <- sample(c("case", "ctrl"), n_samples, replace = TRUE)
#' gene_sets <- list(
#'   pathway_A = genes[1:30],
#'   pathway_B = genes[50:120],
#'   pathway_C = genes[150:200]
#' )
#'
#' res <- run_info_assoc(
#'   expr      = expr,
#'   phenotype = phenotype,
#'   gene_sets = gene_sets,
#'   n_perm    = 200,
#'   seed      = 42
#' )
#' head(res)
#'
#' @export
run_info_assoc <- function(
    expr,
    phenotype,
    gene_sets,
    universe  = NULL,
    score     = c("mean_z", "pc1"),
    nbins     = NULL,
    n_perm    = 1000L,
    min_size  = 10L,
    max_size  = Inf,
    nthreads  = 1L,
    seed      = NULL
) {
  # --- Validation -------------------------------------------------------------
  if (!is.matrix(expr) || !is.numeric(expr))
    stop("`expr` must be a numeric matrix (samples x genes).")
  if (is.null(colnames(expr)))
    stop("`expr` must have gene symbols as column names.")
  if (length(phenotype) != nrow(expr))
    stop("`phenotype` length must equal `nrow(expr)`.")
  if (!is.list(gene_sets) || is.null(names(gene_sets)))
    stop("`gene_sets` must be a named list.")
  if (!is.numeric(n_perm) || n_perm < 1L)
    stop("`n_perm` must be a positive integer.")

  .check_matrix_orientation(expr)
  .check_na_expr(expr)
  .check_sample_size(expr, phenotype)
  .check_pathway_names(names(gene_sets))
  expr <- .remove_zero_variance(expr)

  score <- match.arg(score)
  if (!is.null(seed)) set.seed(seed)
  if (is.null(seed))
    warning("`seed` is NULL. Set `seed` for reproducible results.")

  # Windows does not support forking
  nthreads <- if (.Platform$OS.type == "windows") 1L else as.integer(nthreads)

  # --- Binning ----------------------------------------------------------------
  n <- nrow(expr)

  # Determine nbins_y first so auto_nbins can account for joint table size.
  # Binary/categorical phenotype has far fewer y-bins than continuous.
  is_continuous_pheno <- is.numeric(phenotype) && length(unique(phenotype)) > 2L

  if (is.null(nbins)) {
    nbins_y_hint <- if (is_continuous_pheno) NULL else
      length(unique(phenotype))
    nbins <- .auto_nbins(n, nbins_y = nbins_y_hint)
    message("Auto-selected nbins = ", nbins,
            " (Rice rule, capped for joint-table density and 20).")
  } else {
    nbins <- as.integer(nbins)
    if (nbins < 2L) stop("`nbins` must be >= 2.")
    eff_nbins_y <- if (is_continuous_pheno) nbins else
      length(unique(phenotype))
    exp_freq <- n / (nbins * eff_nbins_y)
    if (exp_freq < 5) {
      warning(
        "`nbins` = ", nbins, " with n = ", n, " samples gives an expected ",
        "joint-cell frequency of ", round(exp_freq, 1), " (< 5). ",
        "MI estimates may be unreliable. ",
        "Consider reducing `nbins` or using `nbins = NULL` for auto-selection."
      )
    }
  }

  # --- Universe ---------------------------------------------------------------
  if (is.null(universe)) universe <- colnames(expr)
  universe <- intersect(universe, colnames(expr))
  expr     <- expr[, universe, drop = FALSE]

  # --- Gene set filtering -----------------------------------------------------
  gene_sets <- lapply(gene_sets, intersect, universe)
  set_sizes <- lengths(gene_sets)
  keep      <- set_sizes >= min_size & set_sizes <= max_size
  gene_sets <- gene_sets[keep]

  n_dropped <- sum(!keep)
  if (n_dropped > 0L)
    message(n_dropped, " gene set(s) dropped by size filter.")
  if (length(gene_sets) == 0L)
    stop(
      "No gene sets remain after size filtering ",
      "(min_size = ", min_size, ", max_size = ", max_size, ")."
    )

  # --- Pathway scores (computed once) -----------------------------------------
  Z <- scale(expr)

  score_fun <- switch(
    score,
    mean_z = function(G) rowMeans(Z[, G, drop = FALSE]),
    pc1    = function(G) {
      # MI is sign-invariant; PC1 sign flip does not affect results
      prcomp(Z[, G, drop = FALSE], center = FALSE, scale. = FALSE)$x[, 1L]
    }
  )

  pathway_scores <- lapply(gene_sets, score_fun)

  # --- Discretise scores (fixed — done once before permutations) --------------
  scores_disc <- lapply(pathway_scores,
                        function(s) .discretize_equalfreq(s, nbins))

  # Convert to matrix for vectorised MI: shape n_samples x n_pathways.
  # This is the key layout for .mi_matrix() and is shared copy-on-write
  # across mclapply workers on Linux/macOS (no per-worker memory copy).
  scores_disc_mat <- do.call(cbind, scores_disc)

  y_disc <- if (is.numeric(phenotype) && length(unique(phenotype)) > 2L) {
    .discretize_equalfreq(phenotype, nbins)
  } else {
    as.integer(as.factor(phenotype))
  }

  # nbins_y adapts automatically: 2 for binary phenotype, nbins for
  # discretised continuous, n_levels for categorical. This keeps the
  # joint frequency table as small as it needs to be, which matters
  # for binary phenotype (2 vs nbins^2 effective cells).
  nbins_y <- max(y_disc)

  # --- Observed MI (vectorised across all pathways) --------------------------
  n_sets <- ncol(scores_disc_mat)
  MI_obs <- .mi_matrix(scores_disc_mat, y_disc, nbins, nbins_y)

  # --- Permutation null -------------------------------------------------------
  # Each permutation shuffles y_disc and recomputes MI for all pathways in
  # one vectorised call. Using mclapply: workers fork and share
  # scores_disc_mat copy-on-write, so memory overhead is negligible.

  perm_fun <- function(b) .mi_matrix(scores_disc_mat, sample(y_disc),
                                     nbins, nbins_y)

  null_list <- if (nthreads > 1L) {
    parallel::mclapply(seq_len(n_perm), perm_fun, mc.cores = nthreads)
  } else {
    lapply(seq_len(n_perm), perm_fun)
  }

  # [n_sets x n_perm] matrix
  MI_null <- matrix(unlist(null_list, use.names = FALSE),
                    nrow = n_sets, ncol = n_perm)

  # --- Summary statistics -----------------------------------------------------
  null_mean <- rowMeans(MI_null)
  # Population SD -- negligible difference from sample SD at n_perm >= 200
  null_sd   <- sqrt(rowMeans((MI_null - null_mean)^2))

  z_score <- ifelse(null_sd > 0, (MI_obs - null_mean) / null_sd, NA_real_)
  p_value <- (1L + rowSums(MI_null >= MI_obs)) / (n_perm + 1L)
  padj    <- p.adjust(p_value, method = "BH")

  # --- Output -----------------------------------------------------------------
  res <- data.frame(
    set        = names(gene_sets),
    set_size   = lengths(gene_sets),
    nbins_used = nbins,
    MI_bits    = MI_obs,
    z_score    = z_score,
    p_value    = p_value,
    padj       = padj,
    null_mean  = null_mean,
    null_sd    = null_sd,
    row.names  = NULL,
    stringsAsFactors = FALSE
  )

  res[order(res$padj), ]
}


# =============================================================================
# Internal helpers
# =============================================================================

#' Automatically select number of bins for MI estimation
#'
#' Uses the Rice rule as a base, capped by a minimum-expected-cell-frequency
#' criterion and a hard ceiling of 20 to prevent sparse-cell MI inflation.
#'
#' The frequency cap is \strong{joint-table aware}: MI is estimated from a
#' \code{nbins_x * nbins_y} contingency table, so the constraint is
#' \code{n / (nbins_x * nbins_y) >= 5} (expected count >= 5 per joint cell).
#' For continuous phenotype (\code{nbins_y = nbins_x}), this gives
#' \code{nbins <= floor(sqrt(n / 5))}. For binary phenotype
#' (\code{nbins_y = 2}), the constraint is much less restrictive:
#' \code{nbins <= floor(n / 10)}.
#'
#' @param n Integer. Number of observations (samples).
#' @param nbins_y Integer or \code{NULL}. Number of bins on the phenotype
#'   axis. \code{NULL} (default) assumes symmetric binning (continuous
#'   phenotype: \code{nbins_y = nbins_x}).
#'
#' @return A single integer bin count.
#' @keywords internal
#' @noRd
.auto_nbins <- function(n, nbins_y = NULL) {
  rice <- ceiling(2 * n^(1/3))   # Rice rule

  if (is.null(nbins_y)) {
    # Continuous phenotype: both axes use nbins -> joint table is nbins^2
    # Require n / nbins^2 >= 5  =>  nbins <= sqrt(n / 5)
    freq_cap <- floor(sqrt(n / 5))
  } else {
    # Known nbins_y (e.g. 2 for binary): joint table is nbins * nbins_y
    # Require n / (nbins * nbins_y) >= 5  =>  nbins <= n / (5 * nbins_y)
    freq_cap <- floor(n / (5L * nbins_y))
  }

  as.integer(max(2L, min(rice, freq_cap, 20L)))
}


#' Equal-frequency discretisation
#'
#' Bins a numeric vector into \code{nbins} equal-frequency (quantile) bins.
#' Ties in the quantile breaks are absorbed so the number of unique bins may
#' be less than \code{nbins} when the data are heavily discrete.
#'
#' @param x Numeric vector.
#' @param nbins Positive integer. Target number of bins.
#'
#' @return Integer vector of bin indices (same length as \code{x}).
#' @keywords internal
#' @noRd
.discretize_equalfreq <- function(x, nbins) {
  breaks <- unique(
    quantile(x, probs = seq(0, 1, length.out = nbins + 1L), na.rm = TRUE)
  )
  as.integer(cut(x, breaks = breaks, include.lowest = TRUE))
}


#' Vectorised mutual information across all pathways for one permutation
#'
#' Computes I(pathway_score ; phenotype) for every pathway simultaneously,
#' using integer joint-bin indices and \code{tabulate()} rather than
#' \code{table()}. This is the performance-critical inner loop for
#' \code{run_info_assoc()} permutation testing.
#'
#' \strong{Why this is faster than the scalar approach:}
#' The joint bin index computation is a single vectorised matrix operation
#' over all n_samples x n_pathways entries at once. \code{tabulate()} on
#' pre-computed integer indices is O(n) with no S3 dispatch or sorting
#' overhead, compared to \code{table()} which is O(n log n) plus
#' dimnames allocation. At n = 8 000 samples and 1 500 pathways the
#' per-permutation cost drops from ~O(n_pathways * n * log n) to
#' ~O(n_pathways * nbins_x * nbins_y + n_samples).
#'
#' \strong{Why nbins_y is separate from nbins_x:}
#' For binary phenotype, \code{as.integer(as.factor(phenotype))} gives
#' values 1..2, so \code{nbins_y = 2}. Using \code{nbins_x} (e.g. 20)
#' for both dimensions would allocate a 20x20 joint table where 18 columns
#' are always zero -- wasted work. Keeping them separate makes binary
#' phenotype ~10x cheaper per call than the symmetric case.
#'
#' @param scores_disc_mat Integer matrix (n_samples x n_pathways), values
#'   in 1..\code{nbins_x}. Each column is one pathway's discretised
#'   activity scores. Produced once before the permutation loop and shared
#'   copy-on-write across \code{mclapply} workers.
#' @param y_disc Integer vector (length n_samples), values in
#'   1..\code{nbins_y}. Discretised phenotype (permuted each call).
#' @param nbins_x Positive integer. Number of bins for pathway scores
#'   (the \code{nbins} argument of \code{run_info_assoc}).
#' @param nbins_y Positive integer. Number of distinct values in
#'   \code{y_disc}. Equals 2 for binary phenotype, \code{nbins_x} for
#'   discretised continuous, or \code{nlevels} for categorical.
#'
#' @return Numeric vector of length \code{ncol(scores_disc_mat)}:
#'   mutual information in bits, one value per pathway.
#' @keywords internal
#' @noRd
.mi_matrix <- function(scores_disc_mat, y_disc, nbins_x, nbins_y) {
  n          <- nrow(scores_disc_mat)
  n_pathways <- ncol(scores_disc_mat)
  nb_xy      <- nbins_x * nbins_y

  # Row-major joint bin index: 1 .. nbins_x * nbins_y
  joint_mat <- (scores_disc_mat - 1L) * nbins_y + (y_disc - 1L) + 1L

  # Marginal y bins (shared across all pathways)
  py    <- tabulate(y_disc, nbins = nbins_y) / n
  k_y   <- sum(py > 0)

  vapply(seq_len(n_pathways), function(p) {
    counts   <- tabulate(joint_mat[, p], nbins = nb_xy)
    pxy      <- matrix(counts / n, nrow = nbins_x, ncol = nbins_y, byrow = TRUE)
    px       <- rowSums(pxy)
    outer_pp <- outer(px, py)
    terms    <- ifelse(pxy <= 0, 0, pxy * log2(pxy / outer_pp))
    mi_raw   <- sum(terms)

    # Miller-Madow bias correction
    k_x        <- sum(px > 0)
    k_xy       <- sum(pxy > 0)
    correction <- (k_x + k_y - k_xy - 1) / (2 * n * log(2))

    max(0, mi_raw + correction)
  }, numeric(1L))
}


#' Empirical mutual information from two discrete vectors
#'
#' Computes I(x ; y) in bits from the empirical joint frequency table.
#' Uses the 0 * log2(0) = 0 convention for empty cells.
#' Applies the Miller-Madow bias correction:
#' \eqn{MI_{MM} = MI_{plugin} + (k_x + k_y - k_{xy} - 1) / (2n \ln 2)}
#' where \eqn{k_x}, \eqn{k_y}, \eqn{k_{xy}} are the number of non-zero
#' marginal and joint bins. The correction is always negative (subtracts
#' the positive finite-sample bias), and is clamped to zero.
#'
#' Used by \code{run_redundancy_selection} and
#' \code{run_assoc_redundancy_selection} for per-pathway scalar MI.
#' For batch permutation testing, see the faster \code{.mi_matrix()}.
#'
#' @param x Integer vector (discrete).
#' @param y Integer vector (discrete), same length as \code{x}.
#'
#' @return Scalar mutual information in bits (bias-corrected, >= 0).
#' @keywords internal
#' @noRd
.mutual_information <- function(x, y) {
  n        <- length(x)
  tbl      <- table(x, y)
  pxy      <- tbl / n
  px       <- rowSums(pxy)
  py       <- colSums(pxy)
  outer_pp <- outer(px, py)
  terms    <- ifelse(pxy <= 0, 0, pxy * log2(pxy / outer_pp))
  mi_raw   <- sum(terms)

  # Miller-Madow bias correction
  k_x  <- sum(px > 0)
  k_y  <- sum(py > 0)
  k_xy <- sum(pxy > 0)
  correction <- (k_x + k_y - k_xy - 1) / (2 * n * log(2))

  max(0, mi_raw + correction)
}
