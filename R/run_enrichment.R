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
    alpha=0.05,
    seed = NULL
) {
  ##---------------- Preprocessing, matching arguments-----------
  method <- match.arg(method)
  alternative <- match.arg(alternative)
  adaptive_mode <- match.arg(adaptive_mode)

  if (is.null(gene_list) || length(gene_list) == 0) stop("gene_list must not be empty.")
  gene_list <- intersect(gene_list, universe)
  if (length(gene_list) == 0) stop("After intersecting with universe, gene_list is empty.")
  if (method == "gsea" && n_perm <= 0) {
    stop("n_perm must be > 0 for GSEA.")
  }
  ## ---------------- Hypergeometric ----------------
  if (method == "hypergeometric") {
    if (is.null(enriched_genes)) {
      stop("You must provide 'enriched_genes' for hypergeometric test.")
    }
    if (adaptive) warning("adaptive inference is ignored for hypergeometric test")

    enriched_genes <- intersect(enriched_genes, universe)

    N <- length(universe)
    K <- length(enriched_genes)
    n <- length(gene_list)
    k <- sum(gene_list %in% enriched_genes)

    # alternative support
    if (alternative == "greater") {
      pval <- phyper(k - 1, K, N - K, n, lower.tail = FALSE)
    } else if (alternative == "less") {
      pval <- phyper(k, K, N - K, n, lower.tail = TRUE)
    } else {
      # Two-sided: use Fisher exact test for a clean definition
      mat <- matrix(c(k, n - k, K - k, (N - K) - (n - k)), nrow = 2, byrow = TRUE)
      pval <- fisher.test(mat, alternative = "two.sided")$p.value
    }

    fc <- if (K == 0) NA_real_ else (k / n) / (K / N)

    return(list(
      method = "hypergeometric",
      p_value = pval,
      overlap = k,
      input_set_size = n,
      group_set = K,
      universe_size = N,
      fold_change = fc
    ))
  }

  ## ---------------- GSEA ----------------
  if (is.null(gene_stats)) stop("You must provide 'gene_stats' for GSEA.")
  if (!is.null(seed)) set.seed(seed)

  # Restrict stats to universe
  gene_stats <- gene_stats[!is.na(names(gene_stats))]
  gene_stats <- gene_stats[!is.na(gene_stats)]
  gene_stats <- gene_stats[names(gene_stats) %in% universe]

  ranked_genes <- sort(gene_stats, decreasing = TRUE)
  if (length(ranked_genes) == 0) stop("After restricting gene_stats to universe (and removing NA), no genes remain to rank.")

  hits <- names(ranked_genes) %in% gene_list
  Nh <- sum(hits)
  N_used <- length(ranked_genes)

  if (Nh == 0) {
    warning("No overlap between gene_list and ranked gene_stats.")
    return(list(
      method = "gsea",
      enrichment_score = 0,
      p_value = NA_real_,
      overlap = 0,
      input_set_size = length(gene_list),
      universe_size = length(universe),
      universe_size_used = N_used,
      leading_edge = I(list(character(0)))
    ))
  }

  # Weighting
  w_base <- abs(ranked_genes)^gsea_weight

  compute_es <- function(hit_mask) {
    denom <- sum(w_base[hit_mask])
    if (denom == 0) return(0)

    step_hit <- numeric(N_used)
    step_miss <- numeric(N_used)

    step_hit[hit_mask] <- w_base[hit_mask] / denom
    step_miss[!hit_mask] <- 1 / (N_used - sum(hit_mask))

    rs <- cumsum(step_hit - step_miss)
    ES_pos <- max(rs)
    ES_neg <- min(rs)
    if (abs(ES_pos) > abs(ES_neg)) ES_pos else ES_neg
  }

  ES <- compute_es(hits)

  # Permutation-based p-value
  perm_fun <- function() {
    random_genes <- sample(names(ranked_genes), Nh, replace = FALSE)
    random_hits  <- names(ranked_genes) %in% random_genes
    compute_es(random_hits)
  }

  count_extreme <- function(Tb, obs) {
    if (alternative == "greater") sum(Tb >= obs)
    else if (alternative == "less") sum(Tb <= obs)
    else sum(abs(Tb) >= abs(obs))
  }

##----------------------Inference-----------------------
  pval <- NA_real_
  inference <- NULL

  if (!adaptive) {
    perm_ES <- replicate(n_perm, perm_fun())
    k0 <- count_extreme(perm_ES, ES)
    B0 <- n_perm
    inf_fixed <- .permutation_inference(
      k_extreme = k0,
      n_perm    = B0,
      alpha=alpha
    )
    pval <- inf_fixed$p_hat
    inference <- c(list(method = "fixed", fixed_k = k0, fixed_B = B0), inf_fixed)

  } else {

    if(adaptive_mode == "pvalue") {
      res <- .adaptive_pvalue(
        obs_stat   = ES,
        perm_fun   = perm_fun,
        alternative = alternative,
        eps        = eps
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
        alpha = alpha
      )
      pval <- res$p
      inference <- list(
        method = "adaptive_decision",
        res
      )
    } else if (adaptive_mode == "refine") {
      # do an initial fixed run then refine only if borderline
      perm_ES <- replicate(n_perm, perm_fun())
      k0 <- count_extreme(perm_ES, ES)
      B0 <- n_perm

      inf_fixed <- .permutation_inference(k_extreme = k0, n_perm = B0, alpha=alpha)

      if(inf_fixed$decision_alpha != "borderline"){
        # nothing to refine
        pval <- inf_fixed$p_hat
        inference <- c(
          list(method = "fixed_not_refined", fixed_k = k0, fixed_B = B0),
          inf_fixed
        )
      } else {
        res <- .adaptive_refine_decision(
          obs_stat = ES,
          perm_fun = perm_fun,
          k0 = k0,
          B0 = B0,
          alternative = alternative,
          alpha = alpha
        )
        pval <- res$p
        inference <- c(list(
          method = "fixed_refine",fixed_k = k0, fixed_B = B0),
          old = inf_fixed,
          new = res
        )
      }
    }
  }

  # Leading edge (appropriate peak depending on ES sign)
  denom_obs <- sum(w_base[hits])
  step_hit <- numeric(N_used)
  step_miss <- numeric(N_used)
  step_hit[hits] <- w_base[hits] / denom_obs
  step_miss[!hits] <- 1 / (N_used - Nh)
  running_score <- cumsum(step_hit - step_miss)

  peak_index <- if (ES >= 0) which.max(running_score) else which.min(running_score)

  leading_edge <- names(ranked_genes)[hits & seq_along(ranked_genes) <= peak_index]

  return(list(
    method = "gsea",
    enrichment_score = ES,
    p_value = pval,
    overlap = Nh,
    input_set_size = length(gene_list),
    universe_size = length(universe),
    universe_size_used = N_used,
    leading_edge = I(list(leading_edge)),
    inference = inference
  ))
}


#' @noRd
#' @keywords internal
.adaptive_pvalue <- function(
    obs_stat,
    perm_fun,
    alternative = c("greater", "less", "two.sided"),
    eps = 0.1,
    min_perm = 200,
    max_perm = 1e5,
    batch = 200,
    verbose = FALSE
) {
  alternative <- match.arg(alternative)

  # Counters
  k <- 0L   # number of extreme permutations
  B <- 0L   # total permutations

  # Optional trace for diagnostics / plotting later
  trace <- list(B = integer(), p_hat = numeric(), RSE = numeric())

  repeat {

    # --- draw a batch of permutations ---
    Tb <- replicate(batch, perm_fun())

    B <- B + length(Tb)

    # --- update extreme count ---
    if (alternative == "greater") {
      k <- k + sum(Tb >= obs_stat)
    } else if (alternative == "less") {
      k <- k + sum(Tb <= obs_stat)
    } else {  # two.sided
      k <- k + sum(abs(Tb) >= abs(obs_stat))
    }

    # --- p-value estimate with +1 correction ---
    p_hat <- (k + 1) / (B + 1)

    # --- Monte Carlo precision ---
    SE  <- sqrt(p_hat * (1 - p_hat) / B)
    RSE <- SE / p_hat

    if (verbose) {
      message(
        sprintf(
          "B = %d | p̂ = %.4g | RSE = %.3f",
          B, p_hat, RSE
        )
      )
    }

    # Store trace (cheap, optional)
    trace$B     <- c(trace$B, B)
    trace$p_hat <- c(trace$p_hat, p_hat)
    trace$RSE   <- c(trace$RSE, RSE)

    # --- stopping rules ---
    if (B >= min_perm && RSE < eps) {
      break
    }

    if (B >= max_perm) {
      break
    }
  }

  list(
    p        = p_hat,
    k        = k,
    B        = B,
    RSE      = RSE,
    converged = (RSE < eps),
    trace    = trace
  )
}


#' @noRd
#' @keywords internal
.adaptive_decision <- function(obs_stat,
                               perm_fun,
                               alternative = c("greater", "less", "two.sided"),
                               alpha = 0.05,
                               B_min = 100,
                               B_max = 20000,
                               batch_size = 50) {

  alternative <- match.arg(alternative)

  B <- 0L
  k <- 0L
  trace <- list()

  repeat {

    # run next batch
    perm_stats <- replicate(batch_size, perm_fun())
    B <- B + batch_size

    # update extreme count
    k_batch <- if (alternative == "greater") {
      sum(perm_stats >= obs_stat)
    } else if (alternative == "less") {
      sum(perm_stats <= obs_stat)
    } else {
      sum(abs(perm_stats) >= abs(obs_stat))
    }

    k <- k + k_batch

    # don't evaluate CI too early
    if (B >= B_min) {

      ci <- binom::binom.confint(
        x = k,
        n = B,
        methods = "wilson"
      )

      decision <- .decision_from_ci(ci_low = ci$lower, ci_high = ci$upper, alpha = alpha)

      trace[[length(trace) + 1]] <- list(
        B = B,
        k = k,
        p_hat = (k + 1) / (B + 1),
        ci_low = ci$lower,
        ci_high = ci$upper,
        decision = decision
      )

      # stop if decision resolved
      if (decision != "borderline") break
    }

    if (B >= B_max) break
  }

  # final values
  p_hat <- (k + 1) / (B + 1)
  mc_se <- sqrt(p_hat * (1 - p_hat) / B)

  ci <- binom::binom.confint(
    x = k,
    n = B,
    methods = "wilson"
  )

  decision_alpha <- .decision_from_ci(ci_low = ci$lower, ci_high = ci$upper, alpha = alpha)

  list(
    p = p_hat,
    B = B,
    k = k,
    mc_se = mc_se,
    p_ci_low = ci$lower,
    p_ci_high = ci$upper,
    decision_alpha = decision_alpha,
    converged = (decision_alpha != "borderline"),
    trace = trace
  )
}



#' @noRd
#' @keywords internal
.adaptive_refine_decision <- function(obs_stat,
                                      perm_fun,
                                      k0,
                                      B0,
                                      alternative = c("greater", "less", "two.sided"),
                                      alpha = 0.05,
                                      B_max = 20000,
                                      batch_size = 50) {

  alternative <- match.arg(alternative)
  B <- as.integer(B0)
  k <- as.integer(k0)
  trace <- list()

  repeat {

    if (B >= B_max) break

    perm_stats <- replicate(batch_size, perm_fun())
    B <- B + batch_size

    k_batch <- if (alternative == "greater") {
      sum(perm_stats >= obs_stat)
    } else if (alternative == "less") {
      sum(perm_stats <= obs_stat)
    } else {
      sum(abs(perm_stats) >= abs(obs_stat))
    }

    k <- k + k_batch

    ci <- binom::binom.confint(
      x = k,
      n = B,
      methods = "wilson"
    )

    decision <- .decision_from_ci(ci_low = ci$lower, ci_high = ci$upper, alpha = alpha)

    trace[[length(trace) + 1]] <- list(
      B = B,
      k = k,
      p_hat = (k + 1) / (B + 1),
      ci_low = ci$lower,
      ci_high = ci$upper,
      decision = decision
    )

    if (decision != "borderline") break
  }

  p_hat <- (k + 1) / (B + 1)

  ci <- binom::binom.confint(
    x = k,
    n = B,
    methods = "wilson"
  )

  decision_alpha <- .decision_from_ci(ci_low = ci$lower, ci_high = ci$upper, alpha = alpha)

  list(
    p = p_hat,
    B = B,
    k = k,
    p_ci_low = ci$lower,
    p_ci_high = ci$upper,
    decision_alpha = decision_alpha,
    converged = (decision_alpha != "borderline"),
    trace = trace
  )
}

#' @noRd
#' @keywords internal
.permutation_inference <- function(k_extreme, n_perm, alpha = 0.05) {

  p_hat <- (k_extreme + 1) / (n_perm + 1)
  mc_se <- sqrt(p_hat * (1 - p_hat) / n_perm)

  ci <- binom::binom.confint(
    x = k_extreme,
    n = n_perm,
    methods = "wilson"
  )

  decision_alpha <- if (ci$upper < alpha) {
    "sig"
  } else if (ci$lower > alpha) {
    "nonsig"
  } else {
    "borderline"
  }

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

#' @noRd
#' @keywords internal
.decision_from_ci <- function(ci_low, ci_high, alpha = 0.05, tol = 1e-4) {
  if (ci_high <= (alpha - tol)) return("sig")
  if (ci_low  >= (alpha + tol)) return("nonsig")
  "borderline"
}
