#' Non-redundant pathway selection via partial correlation
#'
#' A linear-algebraic alternative to \code{\link{run_assoc_redundancy_selection}}
#' that replaces conditional mutual information with partial correlation in the
#' redundancy selection step. The initial association test is identical —
#' permutation-based MI from \code{\link{run_info_assoc}} — but redundancy is
#' assessed by asking whether a candidate pathway's activity correlates with
#' the phenotype residual after removing the variance explained by already-
#' selected pathways.
#'
#' @section Algorithm:
#' \enumerate{
#'   \item Filter \code{assoc_results} to pathways with \code{padj <= alpha}.
#'   \item Score pathway activity (mean z-score or PC1) for each significant
#'     pathway.
#'   \item Select the pathway with the highest \code{MI_bits} as the first
#'     pathway.
#'   \item Regress the phenotype on all selected pathway activities. Compute
#'     residuals.
#'   \item For each remaining pathway, compute the Pearson correlation between
#'     its activity and the current residuals. Test via t-statistic with
#'     \eqn{df = n - k - 2} where \eqn{k} is the number of selected pathways.
#'   \item Select the pathway with the smallest p-value if it passes both
#'     \code{alpha_pcor} and \code{min_pcor}. Update residuals.
#'   \item Repeat until no candidate passes or \code{max_pathways} is reached.
#' }
#'
#' @section Comparison with run_assoc_redundancy_selection:
#' Partial correlation assumes a linear relationship between pathway activity
#' and phenotype residuals. This assumption is reasonable for mean z-score
#' activity (aggregated across genes; approximately normal by CLT) but may
#' miss nonlinear associations. The advantage is statistical efficiency at
#' small \eqn{n}: no discretisation is required, the t-test uses all
#' information in the data, and no permutations are needed for the redundancy
#' step itself.
#'
#' Use \code{run_assoc_redundancy_selection} when:
#' \itemize{
#'   \item Phenotype is categorical or non-normal
#'   \item Pathway-phenotype relationships may be nonlinear
#'   \item You want a fully non-parametric pipeline
#' }
#'
#' Use \code{run_assoc_pcor_selection} when:
#' \itemize{
#'   \item \eqn{n < 200} and statistical power in the redundancy step matters
#'   \item Phenotype is continuous and approximately normal
#'   \item Runtime is a constraint (no permutations in redundancy step)
#'   \item Interpretability of partial \eqn{r} is preferred over MI in bits
#' }
#'
#' @param assoc_results Data frame returned by \code{\link{run_info_assoc}}.
#' @param expr Numeric matrix of expression values, \strong{samples x genes}.
#'   Must be the same matrix used in \code{\link{run_info_assoc}}.
#' @param phenotype Numeric vector of per-sample phenotype values. Must be
#'   continuous. Length must equal \code{nrow(expr)}.
#' @param gene_sets Named list of character vectors defining gene set
#'   membership. Must be the same list used in \code{\link{run_info_assoc}}.
#' @param alpha Numeric. Adjusted p-value threshold for filtering
#'   \code{assoc_results}. Applied to the \code{padj} column. Default
#'   \code{0.05}.
#' @param score Character. Activity scoring method: \code{"mean_z"}
#'   (default) or \code{"pc1"}.
#' @param alpha_pcor Numeric. Significance threshold for the partial
#'   correlation t-test. Default \code{0.05}.
#' @param min_pcor Numeric. Minimum absolute partial correlation to select
#'   a pathway, regardless of p-value. Acts as an effect-size floor.
#'   Default \code{0.1}.
#' @param max_pathways Integer. Hard cap on selected pathways. Default
#'   \code{20L}.
#' @param seed Optional integer passed to \code{\link{set.seed}}.
#'
#' @return A data frame with one row per selected pathway in selection order.
#'   Columns:
#'   \describe{
#'     \item{\code{step}}{Selection order.}
#'     \item{\code{pathway}}{Pathway name.}
#'     \item{\code{set_size}}{Gene set size after universe intersection.}
#'     \item{\code{marginal_mi}}{I(activity ; phenotype) in bits from
#'       \code{assoc_results}.}
#'     \item{\code{marginal_r}}{Pearson correlation of activity with
#'       phenotype (marginal, not partial).}
#'     \item{\code{partial_r}}{Partial correlation of activity with phenotype
#'       given all previously selected pathway activities. Equals
#'       \code{marginal_r} at step 1.}
#'     \item{\code{redundancy_ratio}}{\code{abs(partial_r) /
#'       abs(marginal_r)}. Near 1 = independent; near 0 = redundant.
#'       Comparable to the redundancy ratio in
#'       \code{\link{run_assoc_redundancy_selection}}.}
#'     \item{\code{pcor_pvalue}}{Two-sided p-value for the partial
#'       correlation t-test. \code{NA} at step 1.}
#'     \item{\code{p_value}}{Original \code{p_value} from
#'       \code{assoc_results}.}
#'     \item{\code{padj}}{Original \code{padj} from \code{assoc_results}.}
#'     \item{\code{z_score}}{Original \code{z_score} from
#'       \code{assoc_results}.}
#'   }
#'   Attributes: \code{alpha}, \code{alpha_pcor}, \code{min_pcor},
#'   \code{n_sig}, \code{score}.
#'
#' @seealso
#'   \code{\link{run_info_assoc}} for the upstream association step,
#'   \code{\link{run_assoc_redundancy_selection}} for the MI-based version,
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
#' phenotype <- 0.7 * factor_A + 0.5 * factor_B + rnorm(n_samples, sd = 0.4)
#'
#' gene_sets <- list(
#'   pathway_A  = genes[1:30],
#'   pathway_A2 = genes[5:35],
#'   pathway_B  = genes[200:230],
#'   pathway_C  = genes[100:130]
#' )
#'
#' assoc <- run_info_assoc(
#'   expr      = expr,
#'   phenotype = phenotype,
#'   gene_sets = gene_sets,
#'   n_perm    = 500,
#'   seed      = 42
#' )
#'
#' sel <- run_assoc_pcor_selection(
#'   assoc_results = assoc,
#'   expr          = expr,
#'   phenotype     = phenotype,
#'   gene_sets     = gene_sets,
#'   seed          = 42
#' )
#' print(sel)
#'
#' @export
run_assoc_pcor_selection <- function(
    assoc_results,
    expr,
    phenotype,
    gene_sets,
    alpha        = 0.05,
    score        = c("mean_z", "pc1"),
    alpha_pcor   = 0.05,
    min_pcor     = 0.1,
    max_pathways = 20L,
    seed         = NULL
) {
  ## --- Validation ------------------------------------------------------------
  if (!is.data.frame(assoc_results))
    stop("`assoc_results` must be a data frame from run_info_assoc().")
  required_cols <- c("set", "MI_bits", "p_value", "padj", "z_score")
  missing_cols  <- setdiff(required_cols, colnames(assoc_results))
  if (length(missing_cols) > 0L)
    stop("Missing columns in `assoc_results`: ",
         paste(missing_cols, collapse = ", "),
         ". Pass the direct output of run_info_assoc().")
  if (!is.matrix(expr) || !is.numeric(expr))
    stop("`expr` must be a numeric matrix (samples x genes).")
  if (!is.numeric(phenotype))
    stop("`phenotype` must be numeric. For categorical phenotypes use ",
         "run_assoc_redundancy_selection() instead.")
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
    stop("No significant pathways at padj <= ", alpha, ". ",
         "Consider relaxing alpha or checking run_info_assoc() results.")

  message(sprintf("%d / %d pathways significant at padj <= %s.",
                  n_sig, nrow(assoc_results), alpha))

  sig_sets  <- assoc_results$set[is_sig]
  sig_assoc <- assoc_results[is_sig, ]

  available <- intersect(sig_sets, names(gene_sets))
  if (length(available) < length(sig_sets))
    message(sprintf("%d significant pathway(s) not found in `gene_sets` - dropped.",
                    length(sig_sets) - length(available)))
  if (length(available) == 0L)
    stop("No significant pathways found in `gene_sets`.")

  sig_sets  <- available
  sig_assoc <- sig_assoc[sig_assoc$set %in% available, ]
  gene_sets <- gene_sets[sig_sets]

  ## --- Universe and activity scores ------------------------------------------
  universe  <- colnames(expr)
  gene_sets <- lapply(gene_sets, intersect, universe)
  Z         <- scale(expr)

  score_fun <- switch(
    score,
    mean_z = function(G) rowMeans(Z[, G, drop = FALSE]),
    pc1    = function(G) {
      prcomp(Z[, G, drop = FALSE], center = FALSE, scale. = FALSE)$x[, 1L]
    }
  )

  activity <- lapply(gene_sets, score_fun)

  ## --- Marginal correlations and MI -----------------------------------------
  phenotype_num <- as.numeric(phenotype)
  n             <- length(phenotype_num)

  marginal_r  <- vapply(activity,
                        function(a) cor(a, phenotype_num),
                        numeric(1L))
  mi_marginal <- setNames(sig_assoc$MI_bits, sig_assoc$set)[sig_sets]

  ## --- Empty result template -------------------------------------------------
  empty_result <- data.frame(
    step             = integer(0),
    pathway          = character(0),
    set_size         = integer(0),
    marginal_mi      = numeric(0),
    marginal_r       = numeric(0),
    partial_r        = numeric(0),
    redundancy_ratio = numeric(0),
    pcor_pvalue      = numeric(0),
    p_value          = numeric(0),
    padj             = numeric(0),
    z_score          = numeric(0),
    stringsAsFactors = FALSE
  )

  ## --- First selection (highest MI_bits) ------------------------------------
  best_idx  <- which.max(mi_marginal)
  best_name <- names(mi_marginal)[best_idx]

  selected    <- best_name
  remaining   <- setdiff(sig_sets, best_name)
  selected_df <- data.frame(a1 = activity[[best_name]])
  fit         <- lm(phenotype_num ~ ., data = selected_df)
  resid_curr  <- residuals(fit)

  assoc_lookup <- setNames(
    split(sig_assoc, seq_len(nrow(sig_assoc))),
    sig_assoc$set
  )

  log_rows       <- vector("list", min(max_pathways, length(sig_sets)))
  log_rows[[1L]] <- data.frame(
    step             = 1L,
    pathway          = best_name,
    set_size         = length(gene_sets[[best_name]]),
    marginal_mi      = mi_marginal[[best_name]],
    marginal_r       = marginal_r[[best_name]],
    partial_r        = marginal_r[[best_name]],   # partial = marginal at step 1
    redundancy_ratio = 1,
    pcor_pvalue      = NA_real_,
    p_value          = assoc_lookup[[best_name]]$p_value,
    padj             = assoc_lookup[[best_name]]$padj,
    z_score          = assoc_lookup[[best_name]]$z_score,
    stringsAsFactors = FALSE
  )

  ## --- Greedy forward steps --------------------------------------------------
  n_steps <- min(max_pathways, length(sig_sets))
  if (n_steps < 2L) {
    res <- do.call(rbind, Filter(Negate(is.null), log_rows))
    attr(res, "alpha")      <- alpha
    attr(res, "alpha_pcor") <- alpha_pcor
    attr(res, "min_pcor")   <- min_pcor
    attr(res, "n_sig")      <- n_sig
    attr(res, "score")      <- score
    return(res)
  }

  for (step in seq(2L, n_steps)) {
    if (length(remaining) == 0L) break

    k          <- length(selected)
    df_resid   <- n - k - 2L

    if (df_resid < 1L) {
      message(sprintf("Stopping at step %d: degrees of freedom exhausted.", step))
      break
    }

    # Partial correlation of each candidate with current residuals
    pcor_vals <- vapply(remaining, function(nm) {
      cor(activity[[nm]], resid_curr)
    }, numeric(1L))

    # t-statistic and two-sided p-value
    t_vals  <- pcor_vals * sqrt(df_resid / pmax(1 - pcor_vals^2, 1e-10))
    p_vals  <- 2 * pt(-abs(t_vals), df = df_resid)

    # Select best candidate that passes both thresholds
    passes   <- p_vals < alpha_pcor & abs(pcor_vals) >= min_pcor
    if (!any(passes)) {
      message(sprintf(
        "Stopping at step %d: no candidate passes alpha_pcor=%.3f and min_pcor=%.2f.",
        step, alpha_pcor, min_pcor
      ))
      break
    }

    best_idx  <- which.min(p_vals * (!passes) + p_vals * passes *
                             (p_vals == min(p_vals[passes])) +
                             (1 - as.numeric(passes)))
    # simpler: just pick lowest p among passing
    best_idx  <- which(passes)[which.min(p_vals[passes])]
    best_name <- remaining[best_idx]
    best_pcor <- pcor_vals[best_idx]
    best_pval <- p_vals[best_idx]

    selected  <- c(selected, best_name)
    remaining <- setdiff(remaining, best_name)

    selected_df[[paste0("a", length(selected))]] <- activity[[best_name]]
    fit        <- lm(phenotype_num ~ ., data = selected_df)
    resid_curr <- residuals(fit)

    redundancy_ratio <- if (abs(marginal_r[[best_name]]) > 1e-10) {
      abs(best_pcor) / abs(marginal_r[[best_name]])
    } else NA_real_

    log_rows[[step]] <- data.frame(
      step             = as.integer(step),
      pathway          = best_name,
      set_size         = length(gene_sets[[best_name]]),
      marginal_mi      = mi_marginal[[best_name]],
      marginal_r       = marginal_r[[best_name]],
      partial_r        = best_pcor,
      redundancy_ratio = redundancy_ratio,
      pcor_pvalue      = best_pval,
      p_value          = assoc_lookup[[best_name]]$p_value,
      padj             = assoc_lookup[[best_name]]$padj,
      z_score          = assoc_lookup[[best_name]]$z_score,
      stringsAsFactors = FALSE
    )
  }

  res <- do.call(rbind, Filter(Negate(is.null), log_rows))

  attr(res, "alpha")      <- alpha
  attr(res, "alpha_pcor") <- alpha_pcor
  attr(res, "min_pcor")   <- min_pcor
  attr(res, "n_sig")      <- n_sig
  attr(res, "score")      <- score

  res
}
