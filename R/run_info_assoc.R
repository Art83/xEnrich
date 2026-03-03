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
#' composite rule: Rice rule (\eqn{2n^{1/3}}) capped at
#' \eqn{\lfloor n/5 \rfloor} (expected cell frequency \eqn{\geq 5}) and
#' hard-capped at 20 to limit sparse-cell inflation. Pass an explicit
#' integer to override.
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
#'   unavailable.
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
#'   [run_conditional_enrichment()] for ridge-regularised redundancy.
#'
#' @references
#' Paninski, L. (2003). Estimation of entropy and mutual information.
#' \emph{Neural Computation}, 15(6), 1191–1253.
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

  score <- match.arg(score)
  if (!is.null(seed)) set.seed(seed)

  # Windows does not support forking
  nthreads <- if (.Platform$OS.type == "windows") 1L else as.integer(nthreads)

  # --- Binning ----------------------------------------------------------------
  n <- nrow(expr)
  if (is.null(nbins)) {
    nbins <- .auto_nbins(n)
    message("Auto-selected nbins = ", nbins,
            " (Rice rule, capped at floor(n/5) and 20).")
  } else {
    nbins <- as.integer(nbins)
    if (nbins < 2L) stop("`nbins` must be >= 2.")
    exp_freq <- n / nbins
    if (exp_freq < 5) {
      warning(
        "`nbins` = ", nbins, " with n = ", n, " samples gives an expected ",
        "cell frequency of ", round(exp_freq, 1), " (< 5). ",
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
  y_disc <- as.integer(as.factor(phenotype))

  # --- Observed MI ------------------------------------------------------------
  MI_obs <- vapply(
    scores_disc,
    function(sd) .mutual_information(sd, y_disc),
    numeric(1L)
  )

  # --- Permutation null -------------------------------------------------------
  n_sets <- length(scores_disc)

  perm_one <- function(b) {
    y_perm <- sample(y_disc)
    vapply(scores_disc,
           function(sd) .mutual_information(sd, y_perm),
           numeric(1L))
  }

  null_list <- if (nthreads > 1L) {
    parallel::mclapply(seq_len(n_perm), perm_one, mc.cores = nthreads)
  } else {
    lapply(seq_len(n_perm), perm_one)
  }

  # [n_sets x n_perm] matrix
  MI_null <- matrix(unlist(null_list, use.names = FALSE),
                    nrow = n_sets, ncol = n_perm)

  # --- Summary statistics -----------------------------------------------------
  null_mean <- rowMeans(MI_null)
  # Population SD — negligible difference from sample SD at n_perm >= 200
  null_sd   <- sqrt(rowMeans((MI_null - null_mean)^2))

  z_score <- ifelse(null_sd > 0, (MI_obs - null_mean) / null_sd, NA_real_)
  p_value <- (1L + rowSums(MI_null >= MI_obs)) / (n_perm + 1L)
  padj    <- p.adjust(p_value, method = "BH")

  # --- Output -----------------------------------------------------------------
  res <- data.frame(
    set       = names(gene_sets),
    set_size  = lengths(gene_sets),
    nbins_used = nbins,
    MI_bits   = MI_obs,
    z_score   = z_score,
    p_value   = p_value,
    padj      = padj,
    null_mean = null_mean,
    null_sd   = null_sd,
    row.names = NULL,
    stringsAsFactors = FALSE
  )

  res[order(res$padj), ]
}


#' Automatically select number of bins for MI estimation
#'
#' Uses the Rice rule as a base, capped by the minimum-expected-cell-frequency
#' criterion (expected count >= 5) and a hard ceiling of 20 to prevent
#' sparse-cell MI inflation. This is called when \code{nbins = NULL} in
#' [run_info_assoc()].
#'
#' @param n Integer. Number of observations (samples).
#'
#' @return A single integer bin count.
#' @keywords internal
#' @noRd
.auto_nbins <- function(n) {
  rice     <- ceiling(2 * n^(1/3))   # Rice rule
  freq_cap <- floor(n / 5L)          # expected cell frequency >= 5
  as.integer(min(rice, freq_cap, 20L))
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


#' Empirical mutual information from two discrete vectors
#'
#' Computes I(x ; y) in bits from the empirical joint frequency table.
#' Uses the 0 * log2(0) = 0 convention for empty cells.
#'
#' Note: all fixed-bin MI estimates carry a positive bias that decreases
#' with n. Because the same binning is applied to both observed and
#' permutation-null data, this bias cancels in relative rankings and
#' empirical p-values, but absolute \code{MI_bits} values should not be
#' compared across datasets with different sample sizes.
#'
#' @param x Integer vector (discrete).
#' @param y Integer vector (discrete), same length as \code{x}.
#'
#' @return Scalar mutual information in bits.
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
  sum(terms)
}
