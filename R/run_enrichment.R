# =============================================================================
# Core enrichment engine — hypergeometric and GSEA-like scoring
# =============================================================================


# -----------------------------------------------------------------------------
# Main function
# -----------------------------------------------------------------------------

#' Run enrichment analysis: hypergeometric test or GSEA-like scoring
#'
#' A unified interface for two enrichment approaches sharing a common output
#' structure. The active method is controlled by \code{method}; certain
#' arguments are only meaningful for one method (see Input contracts below).
#'
#' @section Input contracts:
#' \describe{
#'   \item{Hypergeometric (\code{method = "hypergeometric"})}{
#'     Requires \code{gene_set} (the pathway/biological set) and
#'     \code{gene_list} (the query). \code{universe} defaults to the union
#'     of both. \code{gene_stats}, \code{gsea_weight}, \code{n_perm},
#'     \code{adaptive}, \code{adaptive_mode}, \code{eps}, and
#'     \code{keep_trace} are ignored with a warning.}
#'   \item{GSEA (\code{method = "gsea"})}{
#'     Requires \code{gene_stats} (named numeric vector of per-gene scores —
#'     the ranked list) and \code{gene_list} (the gene set being tested).
#'     \code{universe} defaults to \code{names(gene_stats)}.
#'     \code{gene_set} is ignored.}
#' }
#'
#' @section Terminology note:
#' In this function, \code{gene_list} is the \strong{gene set} (pathway)
#' being tested for enrichment, and \code{gene_stats} supplies the
#' \strong{ranked list} of all genes. This is the reverse of the naming
#' convention in some tools (e.g. fgsea). See the \code{gene_stats} argument
#' description.
#'
#' @section GSEA weighting:
#' Gene weights are computed as \code{|gene_stats|^gsea_weight}. Taking the
#' absolute value means genes with large negative scores (e.g. strongly
#' down-regulated) receive the same weight magnitude as genes with large
#' positive scores. This is intentional and consistent with the original
#' weighted GSEA formulation (Subramanian et al. 2005). Set
#' \code{gsea_weight = 0} for the unweighted (classic KS) statistic.
#'
#' @section Adaptive permutation modes:
#' \describe{
#'   \item{\code{"pvalue"}}{Draws batches until RSE (relative standard error
#'     of p-hat) < \code{eps}. Best when an accurate p-value is needed.}
#'   \item{\code{"decision"}}{Draws batches until the Wilson CI for p-hat
#'     lies entirely above or below \code{alpha}. Best when only the
#'     sig/nonsig decision matters.}
#'   \item{\code{"refine"}}{Runs a fixed \code{n_perm} first; if the result
#'     is borderline, switches to decision-mode refinement. Best for
#'     computational efficiency across many gene sets.}
#' }
#'
#' @param gene_list Character vector of gene symbols defining the gene set
#'   under test. For hypergeometric: the query set. For GSEA: the pathway
#'   whose enrichment is being assessed in the \code{gene_stats} ranking.
#' @param universe Character vector of background genes. For hypergeometric,
#'   defaults to \code{union(gene_list, gene_set)}. For GSEA, defaults to
#'   \code{names(gene_stats)}.
#' @param method Character. \code{"hypergeometric"} (default) or
#'   \code{"gsea"}.
#' @param gene_set Character vector. The pathway / biological set for the
#'   hypergeometric test. Ignored for GSEA.
#' @param gene_stats Named numeric vector. Per-gene scores defining the
#'   ranked list for GSEA (e.g. log fold-change, t-statistic). Names are
#'   gene symbols. Ignored for hypergeometric.
#' @param gsea_weight Numeric. Exponent for weighting hits by
#'   \code{|gene_stats|^gsea_weight}. Default \code{1}. Use \code{0} for
#'   the classic unweighted KS statistic.
#' @param n_perm Integer. Permutations for the GSEA null. Default \code{1000}.
#'   Ignored for hypergeometric.
#' @param alternative Character. \code{"greater"} (default), \code{"less"},
#'   or \code{"two.sided"}. For hypergeometric with \code{"two.sided"},
#'   delegates to [stats::fisher.test()].
#' @param adaptive Logical. Use adaptive permutation stopping? Default
#'   \code{FALSE}. GSEA only; ignored with a warning for hypergeometric.
#' @param adaptive_mode Character. \code{"pvalue"} (default),
#'   \code{"decision"}, or \code{"refine"}. See Adaptive permutation modes.
#'   Ignored when \code{adaptive = FALSE}.
#' @param eps Numeric. RSE threshold for \code{adaptive_mode = "pvalue"}.
#'   Default \code{0.01}.
#' @param alpha Numeric. Significance level for CI-based decisions and
#'   \code{stability_flag}. Default \code{0.05}.
#' @param seed Optional integer. Passed to [set.seed()] before any
#'   stochastic computation.
#' @param keep_trace Logical. If \code{TRUE}, attach per-batch diagnostic
#'   traces to the \code{inference} slot. Default \code{FALSE}. Ignored for
#'   hypergeometric.
#'
#' @return A named list with:
#'   \describe{
#'     \item{\code{method}}{Method used.}
#'     \item{\code{p_value}}{Final p-value.}
#'     \item{\code{overlap}}{Number of query genes found in the set /
#'       ranking.}
#'     \item{\code{input_set_size}}{Size of \code{gene_list} after universe
#'       intersection.}
#'     \item{\code{universe_size}}{Size of \code{universe}.}
#'     \item{\code{fold_change}}{(Hypergeometric only.) Observed / expected
#'       overlap ratio.}
#'     \item{\code{group_set}}{(Hypergeometric only.) Size of
#'       \code{gene_set} in universe.}
#'     \item{\code{enrichment_score}}{(GSEA only.) Signed ES.}
#'     \item{\code{universe_size_used}}{(GSEA only.) Genes retained after
#'       NA removal and universe intersection.}
#'     \item{\code{leading_edge}}{(GSEA only.) Genes driving the ES peak.}
#'     \item{\code{inference}}{(GSEA only.) List of inference diagnostics.
#'       \code{NULL} for hypergeometric.}
#'   }
#'
#' @seealso [run_info_enrichment()] for information-theoretic scoring,
#'   [run_redundancy_selection()] for post-enrichment redundancy reduction.
#'
#' @references
#' Subramanian, A. et al. (2005). Gene set enrichment analysis.
#' \emph{PNAS}, 102(43), 15545–15550.
#'
#' @examples
#' # --- Hypergeometric ---
#' universe  <- paste0("G", 1:500)
#' gene_list <- sample(universe, 50)
#' gene_set  <- sample(universe, 40)
#'
#' run_enrichment(
#'   gene_list = gene_list,
#'   universe  = universe,
#'   method    = "hypergeometric",
#'   gene_set  = gene_set
#' )
#'
#' # --- GSEA, fixed permutations ---
#' gene_stats        <- setNames(rnorm(500), universe)
#' gene_list_gsea    <- sample(universe, 40)
#'
#' run_enrichment(
#'   gene_list  = gene_list_gsea,
#'   method     = "gsea",
#'   gene_stats = gene_stats,
#'   n_perm     = 500,
#'   seed       = 1
#' )
#'
#' # --- GSEA, adaptive (refine mode) ---
#' run_enrichment(
#'   gene_list     = gene_list_gsea,
#'   method        = "gsea",
#'   gene_stats    = gene_stats,
#'   n_perm        = 500,
#'   adaptive      = TRUE,
#'   adaptive_mode = "refine",
#'   seed          = 1
#' )
#'
#' @export
run_enrichment <- function(
    gene_list,
    universe      = NULL,
    method        = c("hypergeometric", "gsea"),
    gene_set      = NULL,
    gene_stats    = NULL,
    gsea_weight   = 1,
    n_perm        = 1000L,
    alternative   = c("greater", "less", "two.sided"),
    adaptive      = FALSE,
    adaptive_mode = c("pvalue", "decision", "refine"),
    eps           = 0.01,
    alpha         = 0.05,
    seed          = NULL,
    keep_trace    = FALSE
) {
  # --- Argument matching ------------------------------------------------------
  method        <- match.arg(method)
  alternative   <- match.arg(alternative)
  adaptive_mode <- match.arg(adaptive_mode)

  # --- Seed (always first stochastic action) ----------------------------------
  if (!is.null(seed)) set.seed(seed)

  # --- Shared input validation ------------------------------------------------
  if (is.null(gene_list) || length(gene_list) == 0L)
    stop("`gene_list` must not be empty.")

  .check_duplicates(gene_list, "gene_list")

  # ============================================================================
  # Hypergeometric path
  # ============================================================================
  if (method == "hypergeometric") {

    if (is.null(gene_set))
      stop("Provide `gene_set` (the pathway) for method = 'hypergeometric'.")

    # Warn about GSEA-only arguments being silently ignored
    if (!is.null(gene_stats))
      warning("`gene_stats` is ignored for method = 'hypergeometric'.")
    if (adaptive)
      warning("`adaptive` is ignored for method = 'hypergeometric'.")
    if (!identical(adaptive_mode, "pvalue"))
      warning("`adaptive_mode` is ignored for method = 'hypergeometric'.")
    if (keep_trace)
      warning("`keep_trace` is ignored for method = 'hypergeometric'.")
    if (n_perm != 1000L)
      warning("`n_perm` is ignored for method = 'hypergeometric'.")

    # Universe default
    if (is.null(universe))
      universe <- union(gene_list, gene_set)

    .check_duplicates(universe, "universe")
    .check_case_mismatch(gene_list, universe)

    gene_list <- intersect(gene_list, universe)
    gene_set  <- intersect(gene_set,  universe)

    if (length(gene_list) == 0L)
      stop("After intersecting with `universe`, `gene_list` is empty.")
    if (length(gene_set) == 0L)
      stop("After intersecting with `universe`, `gene_set` is empty.")

    total_pop      <- length(universe)
    pop_success    <- length(gene_set)
    pop_fail       <- total_pop - pop_success
    sample_size    <- length(gene_list)
    sample_success <- sum(gene_list %in% gene_set)

    if (alternative == "greater") {
      pval <- phyper(sample_success - 1, pop_success, pop_fail,
                     sample_size, lower.tail = FALSE)
    } else if (alternative == "less") {
      pval <- phyper(sample_success, pop_success, pop_fail,
                     sample_size, lower.tail = TRUE)
    } else {
      mat  <- matrix(
        c(sample_success,
          sample_size    - sample_success,
          pop_success    - sample_success,
          pop_fail       - (sample_size - sample_success)),
        nrow = 2, byrow = TRUE
      )
      pval <- fisher.test(mat, alternative = "two.sided")$p.value
    }

    fc <- if (pop_success == 0 || sample_size == 0) NA_real_ else
      (sample_success / sample_size) / (pop_success / total_pop)

    return(list(
      method         = "hypergeometric",
      p_value        = pval,
      overlap        = sample_success,
      input_set_size = sample_size,
      group_set      = pop_success,
      universe_size  = total_pop,
      fold_change    = fc,
      inference      = NULL          # consistent slot with GSEA output
    ))
  }

  # ============================================================================
  # GSEA path
  # ============================================================================

  if (is.null(gene_stats))
    stop("Provide `gene_stats` (named numeric ranking) for method = 'gsea'.")
  if (n_perm <= 0L)
    stop("`n_perm` must be > 0 for GSEA.")

  # Warn about hyper-only argument
  if (!is.null(gene_set))
    warning("`gene_set` is ignored for method = 'gsea'.")
  if (!adaptive && !identical(adaptive_mode, "pvalue"))
    warning("`adaptive_mode` is ignored when `adaptive = FALSE`.")

  # Clean gene_stats
  gene_stats <- gene_stats[!is.na(names(gene_stats))]
  gene_stats <- gene_stats[!is.na(gene_stats)]

  # Universe default for GSEA
  if (is.null(universe)) universe <- names(gene_stats)

  .check_duplicates(universe, "universe")
  .check_case_mismatch(gene_list, universe)

  gene_stats <- gene_stats[names(gene_stats) %in% universe]
  gene_list  <- intersect(gene_list, universe)

  if (length(gene_list) == 0L)
    stop("After intersecting with `universe`, `gene_list` is empty.")
  if (length(gene_stats) == 0L)
    stop("After restricting to `universe` and removing NAs, `gene_stats` is empty.")
  if (anyDuplicated(names(gene_stats)))
    stop("`gene_stats` must have unique names.")

  ranked_genes <- sort(gene_stats, decreasing = TRUE)
  hits         <- names(ranked_genes) %in% gene_list
  Nh           <- sum(hits)
  N_used       <- length(ranked_genes)

  # No overlap — return early with consistent structure
  if (Nh == 0L) {
    warning("No overlap between `gene_list` and ranked `gene_stats`.")
    return(list(
      method             = "gsea",
      enrichment_score   = 0,
      p_value            = NA_real_,
      overlap            = 0L,
      input_set_size     = length(gene_list),
      universe_size      = length(universe),
      universe_size_used = N_used,
      leading_edge       = character(0),
      inference          = NULL
    ))
  }

  # Weights: |stat|^p — absolute value is intentional (see ?run_enrichment)
  w_base <- abs(ranked_genes)^gsea_weight

  # ES and running score in one call — avoids duplicating the cumsum
  es_trace     <- .compute_es(hits, w_base, N_used, return_trace = TRUE)
  ES           <- es_trace$ES
  running_score <- es_trace$running_score

  if (Nh == N_used) {
    leading_edge <- character(0)
  } else {
    # Leading edge
    peak_index   <- if (ES >= 0) which.max(running_score) else which.min(running_score)
    leading_edge <- if (ES >= 0) {
      names(ranked_genes)[hits & seq_along(ranked_genes) <= peak_index]
    } else {
      names(ranked_genes)[hits & seq_along(ranked_genes) >= peak_index]
    }
  }


  # Permutation function (closure over w_base, N_used, Nh)
  perm_fun <- .make_perm_fun(w_base, N_used, Nh)

  # --- Inference --------------------------------------------------------------
  inference <- NULL

  if (!adaptive) {
    perm_ES <- replicate(n_perm, perm_fun())
    k0      <- .count_extreme(perm_ES, ES, alternative)
    inf     <- .permutation_inference(k0, n_perm, alpha)
    pval    <- inf$p_hat
    inference <- list(method = "fixed", fixed_k = k0, fixed_B = n_perm,
                      res = inf)

  } else if (adaptive_mode == "pvalue") {
    res  <- .adaptive_pvalue(ES, perm_fun, alternative,
                             eps = eps, keep_trace = keep_trace)
    pval <- res$p
    inference <- list(method = "adaptive_pvalue", res = res)

  } else if (adaptive_mode == "decision") {
    res  <- .adaptive_decision(ES, perm_fun, alternative,
                               alpha = alpha, keep_trace = keep_trace)
    pval <- res$p
    inference <- list(method = "adaptive_decision", res = res)

  } else {
    # refine: fixed run first, adaptive only if borderline
    perm_ES <- replicate(n_perm, perm_fun())
    k0      <- .count_extreme(perm_ES, ES, alternative)
    inf     <- .permutation_inference(k0, n_perm, alpha)

    if (inf$decision_alpha != "borderline") {
      pval <- inf$p_hat
      inference <- list(method = "fixed_not_refined",
                        fixed_k = k0, fixed_B = n_perm, res = inf)
    } else {
      res  <- .adaptive_refine_decision(ES, perm_fun, k0, n_perm,
                                        alternative, alpha,
                                        keep_trace = keep_trace)
      pval <- res$p
      inference <- list(method = "fixed_refined",
                        fixed_k = k0, fixed_B = n_perm,
                        initial = inf, refined = res)
    }
  }

  list(
    method             = "gsea",
    enrichment_score   = ES,
    p_value            = pval,
    overlap            = Nh,
    input_set_size     = length(gene_list),
    universe_size      = length(universe),
    universe_size_used = N_used,
    leading_edge       = leading_edge,
    inference          = inference
  )
}



# -----------------------------------------------------------------------------
# Internal helpers
# -----------------------------------------------------------------------------

#' Compute enrichment score (ES) from a hit mask and weights
#'
#' Implements a weighted Kolmogorov-Smirnov-like running sum. The score is
#' the maximum absolute deviation of the running sum, signed by the direction
#' of the peak. Optionally returns the full running score trace for leading
#' edge computation and diagnostic plots.
#'
#' @param hit_mask Logical vector. TRUE where ranked genes are in the gene set.
#' @param w_base Numeric vector. Per-gene weights (typically |stat|^p).
#' @param N_used Integer. Total number of ranked genes.
#' @param return_trace Logical. If TRUE, return a list with \code{ES} and
#'   \code{running_score}; otherwise return scalar ES. Default FALSE.
#'
#' @return Scalar ES, or a list with \code{ES} and \code{running_score} if
#'   \code{return_trace = TRUE}.
#' @keywords internal
#' @noRd
.compute_es <- function(hit_mask, w_base, N_used, return_trace = FALSE) {
  denom  <- sum(w_base[hit_mask])
  n_miss <- N_used - sum(hit_mask)

  if (denom <= 0 || n_miss <= 0) {
    if (return_trace)
      return(list(ES = 0, running_score = rep(0, N_used)))
    return(0)
  }

  step_hit  <- numeric(N_used)
  step_miss <- numeric(N_used)

  step_hit[hit_mask]   <- w_base[hit_mask] / denom
  step_miss[!hit_mask] <- 1 / n_miss

  rs     <- cumsum(step_hit - step_miss)
  ES_pos <- max(rs)
  ES_neg <- min(rs)
  ES     <- if (abs(ES_pos) >= abs(ES_neg)) ES_pos else ES_neg

  if (return_trace) list(ES = ES, running_score = rs) else ES
}


#' Build a permutation function for GSEA null distribution
#'
#' Returns a zero-argument function that, when called, randomly permutes
#' gene set membership and returns the permuted ES. \code{force()} is used
#' to capture the enclosing values eagerly and avoid lazy-evaluation bugs
#' when the caller's environment changes between calls.
#'
#' @param w_base Numeric vector of per-gene weights.
#' @param N_used Integer. Total ranked genes.
#' @param Nh Integer. Gene set size (number of hits).
#'
#' @return A zero-argument function returning a scalar ES.
#' @keywords internal
#' @noRd
.make_perm_fun <- function(w_base, N_used, Nh) {
  force(w_base); force(N_used); force(Nh)
  function() {
    idx      <- sample.int(N_used, Nh, replace = FALSE)
    hit_mask <- logical(N_used)
    hit_mask[idx] <- TRUE
    .compute_es(hit_mask, w_base, N_used)
  }
}


#' Count permutation statistics at least as extreme as the observed statistic
#'
#' @param Tb Numeric vector of permuted statistics.
#' @param obs Scalar observed statistic.
#' @param alternative Character. One of "greater", "less", "two.sided".
#'
#' @return Integer count.
#' @keywords internal
#' @noRd
.count_extreme <- function(Tb, obs,
                           alternative = c("greater", "less", "two.sided")) {
  alternative <- match.arg(alternative)
  if      (alternative == "greater")   sum(Tb >= obs)
  else if (alternative == "less")      sum(Tb <= obs)
  else                                 sum(abs(Tb) >= abs(obs))
}


#' Wilson score confidence interval for a proportion
#'
#' Preferred over the Wald interval for small counts and extreme proportions,
#' both of which arise routinely in permutation p-value estimation.
#'
#' @param x Integer. Number of successes.
#' @param n Integer. Number of trials.
#' @param conf.level Numeric. Confidence level. Default 0.95.
#'
#' @return Named list with \code{lower} and \code{upper}.
#' @keywords internal
#' @noRd
.wilson_ci <- function(x, n, conf.level = 0.95) {
  if (n <= 0) return(list(lower = NA_real_, upper = NA_real_))
  z      <- stats::qnorm(1 - (1 - conf.level) / 2)
  phat   <- x / n
  denom  <- 1 + z^2 / n
  centre <- (phat + z^2 / (2 * n)) / denom
  half   <- (z * sqrt((phat * (1 - phat) + z^2 / (4 * n)) / n)) / denom
  list(lower = max(0, centre - half),
       upper = min(1, centre + half))
}


#' Classify a p-value CI relative to alpha
#'
#' @param ci_low,ci_high Numeric. Wilson CI bounds.
#' @param alpha Numeric. Significance level.
#' @param tol Numeric. Tolerance around alpha to define the borderline zone.
#'
#' @return Character: "sig", "nonsig", or "borderline".
#' @keywords internal
#' @noRd
.decision_from_ci <- function(ci_low, ci_high, alpha = 0.05, tol = 1e-4) {
  if (ci_high <= alpha - tol) return("sig")
  if (ci_low  >= alpha + tol) return("nonsig")
  "borderline"
}


#' Summarise a fixed permutation run into inference quantities
#'
#' @param k_extreme Integer. Extreme permutation count.
#' @param n_perm Integer. Total permutations.
#' @param alpha Numeric. Significance level.
#'
#' @return Named list with p_hat, CI, MC SE, decision, and stability flag.
#' @keywords internal
#' @noRd
.permutation_inference <- function(k_extreme, n_perm, alpha = 0.05) {
  p_hat <- (k_extreme + 1) / (n_perm + 1)
  mc_se <- sqrt(p_hat * (1 - p_hat) / n_perm)
  ci    <- .wilson_ci(x = k_extreme + 1, n = n_perm + 1)

  decision_alpha <- .decision_from_ci(
    ci_low  = ci$lower,
    ci_high = ci$upper,
    alpha   = alpha
  )

  list(
    k_extreme      = k_extreme,
    n_perm         = n_perm,
    p_hat          = p_hat,
    p_ci_low       = ci$lower,
    p_ci_high      = ci$upper,
    mc_se          = mc_se,
    decision_alpha = decision_alpha,
    stability_flag = decision_alpha == "borderline"
  )
}


# -----------------------------------------------------------------------------
# Adaptive inference engines
# -----------------------------------------------------------------------------

#' Adaptive permutation p-value via RSE stopping criterion
#'
#' Draws permutations in batches until the relative standard error (RSE) of
#' the running p-value estimate falls below \code{eps}, or \code{max_perm}
#' is reached. RSE = SE / p_hat, where SE is the Monte Carlo standard error.
#' Low RSE means the p-value estimate has stabilised.
#'
#' @param obs_stat Scalar observed test statistic.
#' @param perm_fun Zero-argument function returning a single permuted stat.
#' @param alternative Character. One of "greater", "less", "two.sided".
#' @param eps Numeric. RSE threshold for early stopping. Default 0.1.
#' @param min_perm Integer. Minimum permutations before stopping. Default 200.
#' @param max_perm Integer. Hard cap on permutations. Default 1e5.
#' @param batch Integer. Permutations per batch. Default 200.
#' @param verbose Logical. Print running status. Default FALSE.
#' @param keep_trace Logical. Store running B, p_hat, RSE. Default FALSE.
#'
#' @return Named list: \code{p}, \code{k}, \code{B}, \code{RSE},
#'   \code{converged}, and optionally \code{trace}.
#' @keywords internal
#' @noRd
.adaptive_pvalue <- function(obs_stat,
                             perm_fun,
                             alternative = c("greater", "less", "two.sided"),
                             eps        = 0.1,
                             min_perm   = 200L,
                             max_perm   = 1e5,
                             batch      = 200L,
                             verbose    = FALSE,
                             keep_trace = FALSE) {
  alternative <- match.arg(alternative)
  k <- 0L
  B <- 0L
  if (keep_trace) trace <- list(B = integer(), p_hat = numeric(),
                                RSE = numeric())

  repeat {
    Tb <- replicate(batch, perm_fun())
    B  <- B + length(Tb)
    k  <- k + .count_extreme(Tb, obs_stat, alternative)

    p_hat <- (k + 1) / (B + 1)
    SE    <- sqrt(p_hat * (1 - p_hat) / B)
    RSE   <- if (p_hat > 0) SE / p_hat else Inf

    if (verbose)
      message(sprintf("B = %d | p_hat = %.4g | RSE = %.3f", B, p_hat, RSE))

    if (keep_trace) {
      trace$B     <- c(trace$B,     B)
      trace$p_hat <- c(trace$p_hat, p_hat)
      trace$RSE   <- c(trace$RSE,   RSE)
    }

    if (B >= min_perm && RSE < eps) break
    if (B >= max_perm)              break
  }

  result <- list(p = p_hat, k = k, B = B, RSE = RSE,
                 converged = RSE < eps)
  if (keep_trace) result$trace <- trace
  result
}


#' Adaptive permutation p-value via Wilson CI decision boundary
#'
#' Draws permutations until the Wilson confidence interval for the p-value
#' lies entirely above or below \code{alpha} (i.e. decision is resolved),
#' or \code{B_max} is reached.
#'
#' @param obs_stat Scalar observed test statistic.
#' @param perm_fun Zero-argument permutation function.
#' @param alternative Character. One of "greater", "less", "two.sided".
#' @param alpha Numeric. Significance level. Default 0.05.
#' @param B_min Integer. Minimum permutations before checking. Default 100.
#' @param B_max Integer. Hard cap. Default 20000.
#' @param batch_size Integer. Permutations per batch. Default 50.
#' @param keep_trace Logical. Store per-batch diagnostics. Default FALSE.
#'
#' @return Named list: \code{p}, \code{B}, \code{k}, \code{mc_se},
#'   \code{p_ci_low}, \code{p_ci_high}, \code{decision_alpha},
#'   \code{converged}, and optionally \code{trace}.
#' @keywords internal
#' @noRd
.adaptive_decision <- function(obs_stat,
                               perm_fun,
                               alternative = c("greater", "less", "two.sided"),
                               alpha      = 0.05,
                               B_min      = 100L,
                               B_max      = 20000L,
                               batch_size = 50L,
                               keep_trace = FALSE) {
  alternative <- match.arg(alternative)
  B        <- 0L
  k        <- 0L
  ci       <- NULL
  decision <- "borderline"
  if (keep_trace) trace <- list()

  repeat {
    perm_stats <- replicate(batch_size, perm_fun())
    B <- B + batch_size
    k <- k + .count_extreme(perm_stats, obs_stat, alternative)

    if (B >= B_min) {
      ci       <- .wilson_ci(x = k + 1, n = B + 1)
      decision <- .decision_from_ci(ci$lower, ci$upper, alpha)

      if (keep_trace) {
        trace[[length(trace) + 1L]] <- list(
          B        = B,
          k        = k,
          p_hat    = (k + 1) / (B + 1),
          ci_low   = ci$lower,
          ci_high  = ci$upper,
          decision = decision
        )
      }
      if (decision != "borderline") break
    }
    if (B >= B_max) break
  }

  if (is.null(ci)) ci <- .wilson_ci(x = k + 1, n = B + 1)
  p_hat <- (k + 1) / (B + 1)
  mc_se <- sqrt(p_hat * (1 - p_hat) / B)

  result <- list(
    p              = p_hat,
    B              = B,
    k              = k,
    mc_se          = mc_se,
    p_ci_low       = ci$lower,
    p_ci_high      = ci$upper,
    decision_alpha = decision,
    converged      = decision != "borderline"
  )
  if (keep_trace) result$trace <- trace
  result
}


#' Refine a borderline fixed-permutation result adaptively
#'
#' Called by [run_enrichment()] when \code{adaptive_mode = "refine"} and an
#' initial fixed run lands in the borderline zone. Continues drawing
#' permutations, accumulating on top of the fixed run, until the Wilson CI
#' resolves or \code{B_max} is reached.
#'
#' @param obs_stat Scalar observed test statistic.
#' @param perm_fun Zero-argument permutation function.
#' @param k0,B0 Integer. Extreme count and total from the initial fixed run.
#' @param alternative Character. One of "greater", "less", "two.sided".
#' @param alpha Numeric. Significance level.
#' @param B_max Integer. Hard cap on total permutations (initial + refinement).
#' @param batch_size Integer. Permutations per batch in refinement.
#' @param keep_trace Logical. Default FALSE.
#'
#' @return Named list with same structure as [.adaptive_decision()].
#' @keywords internal
#' @noRd
.adaptive_refine_decision <- function(obs_stat,
                                      perm_fun,
                                      k0,
                                      B0,
                                      alternative = c("greater", "less",
                                                      "two.sided"),
                                      alpha      = 0.05,
                                      B_max      = 20000L,
                                      batch_size = 50L,
                                      keep_trace = FALSE) {
  alternative <- match.arg(alternative)
  B        <- as.integer(B0)
  k        <- as.integer(k0)
  decision <- "borderline"
  if (keep_trace) trace <- list()

  ci <- .wilson_ci(x = k + 1, n = B + 1)

  repeat {
    if (B >= B_max) break

    perm_stats <- replicate(batch_size, perm_fun())
    B <- B + batch_size
    k <- k + .count_extreme(perm_stats, obs_stat, alternative)

    ci       <- .wilson_ci(x = k + 1, n = B + 1)
    decision <- .decision_from_ci(ci$lower, ci$upper, alpha)

    if (keep_trace) {
      trace[[length(trace) + 1L]] <- list(
        B        = B,
        k        = k,
        p_hat    = (k + 1) / (B + 1),
        ci_low   = ci$lower,
        ci_high  = ci$upper,
        decision = decision
      )
    }
    if (decision != "borderline") break
  }

  p_hat  <- (k + 1) / (B + 1)
  result <- list(
    p              = p_hat,
    B              = B,
    k              = k,
    p_ci_low       = ci$lower,
    p_ci_high      = ci$upper,
    decision_alpha = decision,
    converged      = decision != "borderline"
  )
  if (keep_trace) result$trace <- trace
  result
}
