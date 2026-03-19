#' Run enrichment analysis using information theory
#'
#' Scores gene sets by mutual information (MI) between gene list membership
#' and pathway membership, rather than by raw overlap count alone. This
#' addresses a known limitation of classical over-representation analysis
#' (ORA): large gene sets can achieve highly significant hypergeometric
#' p-values through sheer size even when the overlap is biologically
#' unspecific. MI naturally penalises uninformative sets because when a
#' pathway covers a large fraction of the universe, membership conveys little
#' information and \code{info_bits} approaches zero.
#'
#' Two complementary p-values are returned:
#' \describe{
#'   \item{\code{p_enrich}}{Exact one-sided hypergeometric p-value on the
#'     overlap count \emph{k}. Comparable to standard ORA tools.}
#'   \item{\code{emp_p}}{Empirical p-value derived from a permutation null on
#'     \emph{information bits}. Ranks pathways by how surprising their MI is,
#'     not just their overlap size.}
#' }
#' When the two measures disagree, the disagreement is itself informative:
#' a set with low \code{p_enrich} but high \code{emp_p} is large and overlaps
#' by expectation; a set with high \code{p_enrich} but low \code{emp_p} is
#' small yet unusually precise. Users may sort by either depending on their
#' question.
#'
#' @param gene_list Character vector of query gene symbols.
#' @param gene_sets Named list of character vectors defining gene set
#'   membership (e.g. GO terms, pathways).
#' @param universe Character vector of background genes. If \code{NULL}
#'   (default), the union of all genes across \code{gene_sets} is used and a
#'   message is emitted. Supplying an explicit universe is strongly recommended
#'   for reproducibility.
#' @param n_perm Integer. Number of hypergeometric draws used to build the
#'   permutation null for \code{emp_p} and \code{z_score}. Default \code{1000}.
#' @param min_size Integer. Gene sets with fewer than \code{min_size} members
#'   after intersection with \code{universe} are silently dropped. Default
#'   \code{5}.
#' @param max_size Integer. Gene sets with more than \code{max_size} members
#'   after intersection with \code{universe} are silently dropped. Defaults to
#'   \code{Inf} (no upper limit). Note that MI naturally down-weights large,
#'   uninformative sets, so filtering is a computational convenience rather
#'   than a mathematical necessity.
#' @param seed Optional integer passed to [set.seed()] for reproducible
#'   permutations. \code{NULL} (default) leaves the RNG state untouched.
#'
#' @return A data frame with one row per gene set (after size filtering),
#'   containing:
#'   \describe{
#'     \item{\code{set}}{Gene set name.}
#'     \item{\code{overlap}}{Number of query genes in the set (\emph{k}).}
#'     \item{\code{set_size}}{Gene set size after universe intersection (\emph{m}).}
#'     \item{\code{list_size}}{Query list size after universe intersection (\emph{n}).}
#'     \item{\code{universe_size}}{Background size (\emph{N}).}
#'     \item{\code{expected}}{Expected overlap under the null (\code{n * m / N}).}
#'     \item{\code{ratio}}{Observed / expected overlap ratio.}
#'     \item{\code{info_bits}}{Mutual information I(list ; set) in bits.}
#'     \item{\code{signed_bits}}{Signed MI: positive when \code{overlap > expected}
#'       (enrichment), negative when below (depletion). The sign is derived from
#'       the direction of overlap, not from MI itself, which is always
#'       non-negative.}
#'     \item{\code{z_score}}{Standard scores of \code{info_bits} relative to the
#'       permutation null. \code{NA} when null standard deviation is zero.}
#'     \item{\code{emp_p}}{Empirical p-value: proportion of null bits >=
#'       observed, with a +1 / (n_perm + 1) continuity correction.}
#'     \item{\code{p_enrich}}{Exact hypergeometric enrichment p-value.}
#'     \item{\code{p_deplete}}{Exact hypergeometric depletion p-value.}
#'     \item{\code{null_mean}, \code{null_sd}}{Mean and SD of the permutation
#'       null bits distribution.}
#'     \item{\code{padj_enrich}, \code{padj_deplete}, \code{padj_emp}}{BH-adjusted
#'       versions of the corresponding p-values.}
#'   }
#'
#' @seealso [run_enrichment()] for classical ORA/GSEA.
#'
#' @examples
#' universe <- paste0("G", 1:500)
#' gene_list <- sample(universe, 50)
#' gene_sets <- list(
#'   pathway_A = sample(universe, 30),
#'   pathway_B = sample(universe, 100),
#'   pathway_C = sample(universe, 10)
#' )
#'
#' res <- run_info_enrichment(
#'   gene_list = gene_list,
#'   gene_sets = gene_sets,
#'   universe  = universe,
#'   n_perm    = 500,
#'   seed      = 42
#' )
#' head(res[order(res$emp_p), ])
#'
#' @export
run_info_enrichment <- function(
    gene_list,
    gene_sets,
    universe = NULL,
    n_perm   = 1000,
    min_size = 5,
    max_size = Inf,
    seed     = NULL
) {
  # --- Input validation -------------------------------------------------------
  if (!is.character(gene_list))
    stop("`gene_list` must be a character vector of gene symbols.")
  if (!is.list(gene_sets) || is.null(names(gene_sets)))
    stop("`gene_sets` must be a named list.")
  if (!is.numeric(n_perm) || n_perm < 1)
    stop("`n_perm` must be a positive integer.")

  .check_duplicates(gene_list, "gene_list")
  .check_pathway_names(names(gene_sets))

  if (!is.null(seed)) set.seed(seed)

  # --- Universe ---------------------------------------------------------------
  if (is.null(universe)) {
    message(
      "`universe` not supplied defaulting to the union of all gene sets ",
      "(", length(unique(unlist(gene_sets, use.names = FALSE))), " genes). ",
      "Supply an explicit universe for reproducible results."
    )
    universe <- unique(unlist(gene_sets, use.names = FALSE))
  }

  .check_duplicates(universe, "universe")
  .check_case_mismatch(gene_list, universe)

  gene_list <- intersect(gene_list, universe)
  gene_sets <- lapply(gene_sets, intersect, universe)

  if (length(gene_list) == 0)
    stop("No genes in `gene_list` are present in `universe`.")

  # --- Size filtering ---------------------------------------------------------
  set_sizes <- lengths(gene_sets)
  keep      <- set_sizes >= min_size & set_sizes <= max_size
  gene_sets <- gene_sets[keep]

  if (length(gene_sets) == 0)
    stop(
      "No gene sets remain after size filtering ",
      "(min_size = ", min_size, ", max_size = ", max_size, "). ",
      "Consider relaxing these thresholds."
    )

  n_dropped <- sum(!keep)
  if (n_dropped > 0)
    message(n_dropped, " gene set(s) dropped by size filter.")

  # --- Core loop --------------------------------------------------------------
  N <- length(universe)
  n <- length(gene_list)

  results <- lapply(names(gene_sets), function(set_name) {
    G <- gene_sets[[set_name]]
    m <- length(G)
    k <- length(intersect(gene_list, G))

    obs_bits <- .ite_overlap_bits(N, n, m, k)
    expected <- n * m / N
    ratio    <- if (expected > 0) k / expected else NA_real_

    p_enrich  <- phyper(k - 1, m, N - m, n, lower.tail = FALSE)
    p_deplete <- phyper(k,     m, N - m, n, lower.tail = TRUE)

    # Vectorised null: one rhyper call, one .ite_overlap_bits call
    k_null    <- rhyper(n_perm, m, N - m, n)
    bits_null <- .ite_overlap_bits(N, n, m, k_null)

    null_mean <- mean(bits_null)
    null_sd   <- sd(bits_null)
    z_score   <- if (null_sd > 0) (obs_bits - null_mean) / null_sd else NA_real_
    emp_p     <- (sum(bits_null >= obs_bits) + 1L) / (n_perm + 1L)

    data.frame(
      set           = set_name,
      overlap       = k,
      set_size      = m,
      list_size     = n,
      universe_size = N,
      expected      = expected,
      ratio         = ratio,
      info_bits     = obs_bits,
      signed_bits   = sign(k - expected) * obs_bits,
      z_score       = z_score,
      emp_p         = emp_p,
      p_enrich      = p_enrich,
      p_deplete     = p_deplete,
      null_mean     = null_mean,
      null_sd       = null_sd,
      stringsAsFactors = FALSE
    )
  })

  # --- Assemble and adjust ----------------------------------------------------
  res              <- do.call(rbind, results)
  res$padj_enrich  <- p.adjust(res$p_enrich,  method = "BH")
  res$padj_deplete <- p.adjust(res$p_deplete, method = "BH")
  res$padj_emp     <- p.adjust(res$emp_p,     method = "BH")

  res[order(res$emp_p), ]
}


#' Compute mutual information (bits) for a 2x2 gene-set contingency table
#'
#' Internal workhorse for [run_info_enrichment()]. Computes I(gene_list ; gene_set)
#' from the four joint probabilities implied by overlap count \code{k}.
#' Vectorised over \code{k} - pass a scalar for a single observed value or a
#' vector for the permutation null in one call.
#'
#' Applies the Miller-Madow bias correction. For a 2x2 table the correction
#' is at most 1/(2N ln 2) bits — negligible at typical universe sizes but
#' included for consistency with the sample-space MI functions.
#'
#' @param N Integer. Universe size.
#' @param n Integer. Gene list size.
#' @param m Integer. Gene set size.
#' @param k Integer (scalar or vector). Overlap count(s).
#'
#' @return Numeric vector of mutual information values in bits (>= 0).
#' @keywords internal
#' @noRd
.ite_overlap_bits <- function(N, n, m, k) {
  p11 <- k / N
  p10 <- (n - k) / N
  p01 <- (m - k) / N
  p00 <- (N - n - m + k) / N
  p1. <- n / N
  p0. <- 1 - p1.
  p.1 <- m / N
  p.0 <- 1 - p.1

  # Returns 0 for empty cells consistent with the 0*log(0) = 0 convention
  term <- function(p, pa, pb) {
    ifelse(p <= 0, 0, p * log2(p / (pa * pb)))
  }

  mi_raw <- term(p11, p1., p.1) +
    term(p10, p1., p.0) +
    term(p01, p0., p.1) +
    term(p00, p0., p.0)

  # Miller-Madow: k_x = k_y = 2 (always, given n>0, m>0, n<N, m<N).
  # k_xy = number of non-zero joint cells.
  k_xy <- (p11 > 0) + (p10 > 0) + (p01 > 0) + (p00 > 0)
  correction <- (3 - k_xy) / (2 * N * log(2))

  pmax(0, mi_raw + correction)
}




