#' Run enrichment analysis using hypergeometric test or GSEA-like score
#'
#' @param gene_list Character vector of user-supplied gene symbols.
#' @param universe Character vector of all genes in the background set.
#' @param method Enrichment method to use: "hypergeometric" or "gsea".
#' @param enriched_genes Required if method = "hypergeometric".
#' @param gene_stats Required if method = "gsea". Named numeric vector of gene scores (e.g., mean pTPM).
#' @param gsea_weight Exponent for GSEA running sum (default = 1).
#' @param n_perm Number of permutations for empirical p-value (GSEA only).
#' @param seed Optional random seed for reproducibility.
#'
#' @return A named list with enrichment results.
#' @export
run_enrichment <- function(
    gene_list,
    universe,
    method = c("hypergeometric", "gsea"),
    enriched_genes = NULL,
    gene_stats = NULL,
    gsea_weight = 1,
    n_perm = 1000,
    alternative = c("greater", "less", "two.sided"),
    adaptive = FALSE,
    adaptive_mode = c("pvalue", "decision", "refine"),
    eps = 0.01,
    seed = NULL
) {
  if (length(gene_list) == 0 || is.null(gene_list)) {
    stop("gene_list must not be empty.")
  }
  if(method=="gsea" && n_perm <= 0){
    stop("Number of permutations for gsea is not specified.")
  }

  method <- match.arg(method)
  alternative <- match.arg(alternative)

  gene_list <- intersect(gene_list, universe)

  ## ---------------- Hypergeometric ----------------
  if (method == "hypergeometric") {
    if (is.null(enriched_genes)) {
      stop("You must provide 'enriched_genes' for hypergeometric test.")
    }
    if (method == "hypergeometric" && adaptive) {
      warning("adaptive inference is ignored for hypergeometric test")
    }

    enriched_genes <- intersect(enriched_genes, universe)

    N <- length(universe)
    K <- length(enriched_genes)
    n <- length(gene_list)
    k <- sum(gene_list %in% enriched_genes)

    pval <- phyper(k - 1, K, N - K, n, lower.tail = FALSE)

    return(list(
      method = "hypergeometric",
      p_value = pval,
      overlap = k,
      input_set_size = n,
      group_set = K,
      universe_size = N,
      fold_change = (k / n) / (K / N)
    ))
  }

  ## ---------------- GSEA ----------------
  if (is.null(gene_stats)) {
    stop("You must provide 'gene_stats' for GSEA.")
  }
  if (adaptive && n_perm <= 0) {
    stop("adaptive = TRUE requires n_perm > 0")
  }

  # Restrict stats to universe
  gene_stats <- gene_stats[names(gene_stats) %in% universe]
  ranked_genes <- sort(gene_stats, decreasing = TRUE)

  hits <- names(ranked_genes) %in% gene_list
  Nh <- sum(hits)
  N <- length(ranked_genes)

  if (Nh == 0) {
    warning("No overlap between gene_list and ranked gene_stats.")
    return(list(
      method = "gsea",
      enrichment_score = 0,
      p_value = NA_real_,
      overlap = 0,
      input_set_size = length(gene_list),
      universe_size = N,
      leading_edge = I(list(character(0)))
    ))
  }

  # Weighting
  weights <- abs(ranked_genes)^gsea_weight
  weights <- weights / sum(weights[hits])

  step_hit <- numeric(N)
  step_miss <- numeric(N)
  step_hit[hits] <- weights[hits]
  step_miss[!hits] <- 1 / (N - Nh)

  running_score <- cumsum(step_hit - step_miss)

  ES_pos <- max(running_score)
  ES_neg <- min(running_score)
  ES <- if (abs(ES_pos) > abs(ES_neg)) ES_pos else ES_neg

  # Permutation-based p-value
  perm_fun <- function() {
    random_genes <- sample(names(ranked_genes), Nh)
    random_hits  <- names(ranked_genes) %in% random_genes

    weights_r <- abs(ranked_genes)^gsea_weight
    denom <- sum(weights_r[random_hits])
    if (denom == 0) return(0)
    weights_r <- weights_r / denom

    step_hit_r  <- numeric(N)
    step_miss_r <- numeric(N)

    step_hit_r[random_hits]  <- weights_r[random_hits]
    step_miss_r[!random_hits] <- 1 / (N - Nh)

    rs <- cumsum(step_hit_r - step_miss_r)

    ES_pos_r <- max(rs)
    ES_neg_r <- min(rs)

    if (alternative == "two.sided") {
      max(abs(ES_pos_r), abs(ES_neg_r))
    } else {
      if (abs(ES_pos_r) > abs(ES_neg_r)) ES_pos_r else ES_neg_r
    }
  }

  pval <- NA_real_
  if (!is.null(seed)) set.seed(seed)

  perm_ES <- replicate(n_perm, perm_fun())

  k0 <- if (alternative == "greater") {
    sum(perm_ES >= ES)
  } else if (alternative == "less") {
    sum(perm_ES <= ES)
  } else {
    sum(perm_ES >= abs(ES))
  }
  B0 <- n_perm

  inf_fixed <- .permutation_inference(
    k_extreme = k0,
    n_perm    = B0
  )

  if (!adaptive) {
    pval <- inf_fixed$p_hat
    inference <- c(
      list(method = "fixed"),
      inf_fixed
    )

  } else {
    adaptive_mode <- match.arg(adaptive_mode)
    if(adaptive_mode == "pvalue") {
      res <- .adaptive_pvalue(
        obs_stat   = ES,
        perm_fun   = perm_fun,
        alternative = alternative,
        eps        = eps,
        seed       = seed
      )
      pval <- res$p
      inference <- list(
        method = "adaptive",
        res
      )
    } else if(adaptive_mode == "decision"){
      res <- .adaptive_decision(
        obs_stat    = ES,
        perm_fun    = perm_fun,
        alternative = alternative,
        seed        = seed
      )
      pval <- res$p
      inference <- list(
        method = "adaptive_decision",
        res
      )
    } else if (adaptive_mode == "refine") {
      if(inf_fixed$decision_alpha != "borderline"){
        # nothing to refine
        pval <- inf_fixed$p_hat
        inference <- c(
          list(method = "fixed_not_refined"),
          old_pval = inf_fixed$p_hat,
          old_p_ci_low = inf_fixed$p_ci_low,
          old_p_ci_high = inf_fixed$p_ci_high,
          old_decision_alpha = inf_fixed$decision_alpha,
          inf_fixed
        )
      } else {
        res <- .adaptive_refine_decision(
          obs_stat = ES,
          perm_fun = perm_fun,
          k0 = k0,
          B0 = B0,
          alternative = alternative,
          alpha = alpha,
          seed = seed
        )

        pval <- res$p
        inference <- c(
          list(method = "fixed_refine"),
          old_pval = inf_fixed$p_hat,
          old_p_ci_low = inf_fixed$p_ci_low,
          old_p_ci_high = inf_fixed$p_ci_high,
          old_decision_alpha = inf_fixed$decision_alpha,
          new_pval = res$p_hat,
          new_p_ci_low = res$p_ci_low,
          new_p_ci_high = res$p_ci_high,
          new_decision_alpha = res$decision_alpha,
          converged = res$converged,
          res,
          list(
            fixed_k = k0,
            fixed_B = B0,
            resolved = res$decision_alpha != "borderline"
          )
        )
      }
    }
  }

  # Leading edge (appropriate peak depending on ES sign)
  peak_index <- if (ES >= 0) {
    which.max(running_score)
  } else {
    which.min(running_score)
  }

  if (ES >= 0) {
    leading_edge <- names(ranked_genes)[hits & seq_along(ranked_genes) <= peak_index]
  } else {
    leading_edge <- names(ranked_genes)[hits & seq_along(ranked_genes) >= peak_index]
  }

  return(list(
    method = "gsea",
    enrichment_score = ES,
    p_value = pval,
    overlap = Nh,
    input_set_size = length(gene_list),
    universe_size = length(universe),
    leading_edge = I(list(leading_edge)),
    inference = inference
  ))
}
