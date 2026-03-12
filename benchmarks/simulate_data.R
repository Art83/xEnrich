# =============================================================================
# Simulation framework for xEnrich benchmarking
#
# Generates synthetic gene expression + pathway data with known ground truth.
# Three scenarios of increasing complexity, each parametrised by:
#   n_samples        number of observations
#   n_genes          size of the gene universe
#   signal_strength  effect size of latent factors on phenotype
#   seed             RNG seed for reproducibility
#
# Ground truth: which pathways are "true signals" (driven by an independent
# latent factor) is known exactly. Benchmarks measure how well each method
# recovers the true signal set.
# =============================================================================


#' Generate a synthetic expression + pathway dataset with known ground truth
#'
#' @param scenario    Integer 1, 2, or 3. See details below.
#' @param n_samples   Number of samples (rows of expr).
#' @param n_genes     Total gene universe size.
#' @param signal_strength Numeric in (0, Inf). SD of the latent factor
#'   contribution to gene expression. Higher = stronger signal.
#'   Typical values: 0.5 (weak), 1.0 (moderate), 2.0 (strong).
#' @param noise_sd    SD of gene-level noise. Default 1.
#' @param seed        RNG seed.
#'
#' @return A list with:
#'   $expr          numeric matrix (n_samples x n_genes)
#'   $phenotype     numeric vector (n_samples)
#'   $gene_sets     named list of character vectors
#'   $gene_stats    named numeric vector (for GSEA lane)
#'   $true_signals  character vector of pathway names that are true signals
#'   $scenario      integer
#'   $params        list of all parameters used
#'
#' Scenario descriptions:
#'   1 - Clean separation: 2 independent signals, 4 correlated shadows, 2 noise
#'   2 - Weak signal:      2 signals (one weak), 3 shadows, 2 noise
#'   3 - High redundancy:  1 dominant signal, 8 near-identical shadows,
#'                         2 independent moderate signals
simulate_data <- function(
    scenario       = 1L,
    n_samples      = 100L,
    n_genes        = 2000L,
    signal_strength = 1.0,
    noise_sd       = 1.0,
    seed           = 42L,
    compute_gene_stats=TRUE
) {
  set.seed(seed)
  stopifnot(scenario %in% 1:4)
  stopifnot(n_genes >= 1000L)   # scenario 3 uses indices up to 990

  genes <- paste0("G", seq_len(n_genes))

  # ---------------------------------------------------------------------------
  # Build gene sets per scenario
  # Each scenario defines:
  #   factors       list of latent factors (each a numeric vector, n_samples)
  #   pathway_defs  list of lists with $genes, $factor_idx, $is_signal
  # ---------------------------------------------------------------------------

  if (scenario == 1L) {
    # Two independent latent factors
    fA <- rnorm(n_samples)
    fB <- rnorm(n_samples)

    pathway_defs <- list(
      # True signals
      list(name = "signal_A",   genes = genes[1:60],    factor = fA, is_signal = TRUE),
      list(name = "signal_B",   genes = genes[301:360],  factor = fB, is_signal = TRUE),
      # Correlated shadows of A (high gene overlap)
      list(name = "shadow_A1",  genes = genes[1:55],    factor = fA, is_signal = FALSE),
      list(name = "shadow_A2",  genes = c(genes[5:60], genes[61:65]),  factor = fA, is_signal = FALSE),
      list(name = "shadow_A3",  genes = genes[1:50],    factor = fA, is_signal = FALSE),
      list(name = "shadow_A4",  genes = c(genes[10:60], genes[66:70]), factor = fA, is_signal = FALSE),
      # Noise (no latent factor)
      list(name = "noise_1",    genes = genes[700:750],  factor = NULL, is_signal = FALSE),
      list(name = "noise_2",    genes = genes[800:840],  factor = NULL, is_signal = FALSE)
    )

    phenotype_factors <- list(list(f = fA, w = 0.7), list(f = fB, w = 0.7))

  } else if (scenario == 2L) {
    # Two signals: one moderate, one weak
    fA <- rnorm(n_samples)
    fB <- rnorm(n_samples)

    pathway_defs <- list(
      list(name = "signal_A",   genes = genes[1:60],    factor = fA, is_signal = TRUE),
      list(name = "signal_B",   genes = genes[301:360],  factor = fB, is_signal = TRUE),  # weak
      list(name = "shadow_A1",  genes = genes[1:55],    factor = fA, is_signal = FALSE),
      list(name = "shadow_A2",  genes = c(genes[5:60], genes[61:65]), factor = fA, is_signal = FALSE),
      list(name = "shadow_A3",  genes = genes[10:60],   factor = fA, is_signal = FALSE),
      list(name = "noise_1",    genes = genes[700:750],  factor = NULL, is_signal = FALSE),
      list(name = "noise_2",    genes = genes[800:840],  factor = NULL, is_signal = FALSE)
    )

    # fB has lower phenotype weight -> weaker but still detectable signal
    phenotype_factors <- list(list(f = fA, w = 0.7), list(f = fB, w = 0.5))

  } else if (scenario == 3L) {
    # One dominant signal with 8 near-identical shadows, plus 2 independent
    fA <- rnorm(n_samples)
    fB <- rnorm(n_samples)
    fC <- rnorm(n_samples)

    # Shadows: vary start/end position slightly but >90% overlap with signal_A
    shadow_defs <- lapply(seq_len(8L), function(i) {
      offset <- (i - 1L) * 3L
      list(
        name      = paste0("shadow_A", i),
        genes     = genes[(1 + offset):(58 + offset)],
        factor    = fA,
        is_signal = FALSE
      )
    })

    pathway_defs <- c(
      list(
        list(name = "signal_A", genes = genes[1:60],    factor = fA, is_signal = TRUE),
        list(name = "signal_B", genes = genes[301:360], factor = fB, is_signal = TRUE),
        list(name = "signal_C", genes = genes[601:650], factor = fC, is_signal = TRUE)
      ),
      shadow_defs,
      list(
        list(name = "noise_1", genes = genes[900:940],  factor = NULL, is_signal = FALSE),
        list(name = "noise_2", genes = genes[950:990],  factor = NULL, is_signal = FALSE)
      )
    )

    phenotype_factors <- list(
      list(f = fA, w = 0.7),
      list(f = fB, w = 0.5),
      list(f = fC, w = 0.4)
    )

  } else {
    # Null scenario: no latent factors, all pathways are noise
    # Used to measure false positive rate — any selection is an error
    pathway_defs <- list(
      list(name = "null_1", genes = genes[1:60],    factor = NULL, is_signal = FALSE),
      list(name = "null_2", genes = genes[61:120],  factor = NULL, is_signal = FALSE),
      list(name = "null_3", genes = genes[121:180], factor = NULL, is_signal = FALSE),
      list(name = "null_4", genes = genes[181:240], factor = NULL, is_signal = FALSE),
      list(name = "null_5", genes = genes[301:360], factor = NULL, is_signal = FALSE),
      list(name = "null_6", genes = genes[361:420], factor = NULL, is_signal = FALSE),
      list(name = "null_7", genes = genes[421:480], factor = NULL, is_signal = FALSE),
      list(name = "null_8", genes = genes[481:540], factor = NULL, is_signal = FALSE)
    )
    phenotype_factors <- list()   # phenotype is pure noise
  }

  # ---------------------------------------------------------------------------
  # Build expression matrix
  # Each gene in a pathway is driven by the pathway's latent factor +
  # gene-level noise. Genes not in any signal pathway are pure noise.
  # ---------------------------------------------------------------------------
  expr <- matrix(
    rnorm(n_samples * n_genes, sd = noise_sd),
    nrow     = n_samples,
    dimnames = list(NULL, genes)
  )

  for (pd in pathway_defs) {
    if (!is.null(pd$factor)) {
      g <- intersect(pd$genes, genes)
      expr[, g] <- expr[, g] + signal_strength * pd$factor
    }
  }

  # ---------------------------------------------------------------------------
  # Phenotype: weighted sum of latent factors + noise
  # ---------------------------------------------------------------------------
  phenotype <- if (length(phenotype_factors) > 0L) {
    Reduce("+", lapply(phenotype_factors, function(x) x$w * x$f)) +
      rnorm(n_samples, sd = 0.5)
  } else {
    rnorm(n_samples, sd = 0.5)   # scenario 4: pure noise phenotype
  }

  # ---------------------------------------------------------------------------
  # gene_sets: named list
  # true_signals: names of pathways that are true independent signals
  # ---------------------------------------------------------------------------
  gene_sets    <- setNames(
    lapply(pathway_defs, `[[`, "genes"),
    sapply(pathway_defs, `[[`, "name")
  )
  true_signals <- sapply(pathway_defs, function(pd) {
    if (pd$is_signal) pd$name else NULL
  })
  true_signals <- unlist(true_signals[!sapply(true_signals, is.null)])

  # ---------------------------------------------------------------------------
  # gene_stats for GSEA lane: correlation of each gene with phenotype
  # (mimics a differential expression t-statistic)
  # Only computed when needed — skip with compute_gene_stats = FALSE for
  # Lane C benchmarks to avoid 2000-gene correlation at n=500.
  # ---------------------------------------------------------------------------
  gene_stats <- if (compute_gene_stats) {
    setNames(apply(expr, 2L, function(g) cor(g, phenotype, method="pearson")), genes)
  } else NULL
  # (callers that don't need gene_stats can ignore this field)

  list(
    expr         = expr,
    phenotype    = phenotype,
    gene_sets    = gene_sets,
    gene_stats   = gene_stats,
    true_signals = true_signals,
    scenario     = scenario,
    params       = list(
      n_samples       = n_samples,
      n_genes         = n_genes,
      signal_strength = signal_strength,
      noise_sd        = noise_sd,
      seed            = seed
    )
  )
}


#' Compute precision, recall and F1 from selected vs true pathway names
#'
#' @param selected   Character vector of selected pathway names.
#' @param truth      Character vector of true signal pathway names.
#'
#' @return Named numeric vector with precision, recall, f1.
score_selection <- function(selected, truth) {
  if (length(selected) == 0L) {
    return(c(precision = NA_real_, recall = 0, f1 = NA_real_,
             n_selected = 0L, n_true = length(truth),
             tp = 0L, fp = 0L, fn = length(truth)))
  }
  tp <- length(intersect(selected, truth))
  fp <- length(setdiff(selected, truth))
  fn <- length(setdiff(truth, selected))

  precision <- if ((tp + fp) > 0) tp / (tp + fp) else NA_real_
  recall    <- if ((tp + fn) > 0) tp / (tp + fn) else NA_real_
  pr_sum    <- sum(c(precision, recall), na.rm = FALSE)
  f1        <- if (!is.na(pr_sum) && pr_sum > 0) {
    2 * precision * recall / pr_sum
  } else NA_real_

  c(precision  = precision,
    recall     = recall,
    f1         = f1,
    n_selected = length(selected),
    n_true     = length(truth),
    tp         = tp,
    fp         = fp,
    fn         = fn)
}
