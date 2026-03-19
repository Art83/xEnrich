## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 8,
  fig.height = 5
)
library(xEnrich)

## ----synthetic-data-----------------------------------------------------------
set.seed(42)
n_samples <- 120
n_genes   <- 500
genes     <- paste0("G", seq_len(n_genes))

# Two independent latent factors drive two clusters of genes
factor_A <- rnorm(n_samples)
factor_B <- rnorm(n_samples)

expr <- matrix(
  rnorm(n_samples * n_genes),
  nrow     = n_samples,
  dimnames = list(paste0("S", seq_len(n_samples)), genes)
)

# Inject signal into two gene blocks
expr[, 1:40]    <- expr[, 1:40]    + 2.5 * factor_A
expr[, 300:340] <- expr[, 300:340] + 2.0 * factor_B

# Phenotype: weighted sum of both factors + noise
phenotype <- 0.7 * factor_A + 0.5 * factor_B + rnorm(n_samples, sd = 0.4)

# Pathway definitions: two true signals, one shadow, one noise
gene_sets <- list(
  immune_response   = genes[1:40],           # true signal (factor A)
  immune_regulation = genes[1:35],           # redundant shadow of immune_response
  solute_transport  = genes[300:340],        # true signal (factor B)
  cell_adhesion     = genes[200:230],        # noise
  lipid_metabolism  = genes[400:430]         # noise
)

## ----assoc--------------------------------------------------------------------
assoc <- run_info_assoc(
  expr      = expr,
  phenotype = phenotype,
  gene_sets = gene_sets,
  n_perm    = 500,
  seed      = 42
)

assoc[, c("set", "set_size", "MI_bits", "z_score", "p_value", "padj")]

## ----redundancy---------------------------------------------------------------
sel <- run_assoc_redundancy_selection(
  assoc_results = assoc,
  expr          = expr,
  phenotype     = phenotype,
  gene_sets     = gene_sets,
  alpha         = 0.05,
  seed          = 42
)

sel[, c("step", "pathway", "marginal_mi", "conditional_mi",
        "redundancy_ratio", "cumulative_bits")]

## ----classify-----------------------------------------------------------------
classified <- classify_assoc_selection(sel)

classified[, c("pathway", "marginal_strength", "signal_class", "action")]

## ----activity-----------------------------------------------------------------
act <- compute_pathway_activity(
  expr      = expr,
  gene_sets = gene_sets[sel$pathway],
  score     = "mean_z"
)

cat("Activity matrix:", nrow(act), "samples x", ncol(act), "pathways\n")
cat("Correlation with phenotype:\n")
round(cor(act, phenotype), 3)

## ----gains-plot---------------------------------------------------------------
plot_gains(sel)

## ----gene-list-lane-----------------------------------------------------------
# Simulate a DE gene list
gene_list <- unique(c(
  sample(genes[1:40],    30),   # enriched in immune_response
  sample(genes[300:340], 20),   # enriched in solute_transport
  sample(genes[100:200], 5)     # a few noise genes
))

# MI-scored enrichment
enrich <- run_info_enrichment(
  gene_list = gene_list,
  gene_sets = gene_sets,
  universe  = genes,
  n_perm    = 500,
  seed      = 42
)

enrich[, c("set", "overlap", "set_size", "info_bits", "ratio",
           "padj_emp", "padj_enrich")]

## ----gene-space-redundancy----------------------------------------------------
# Redundancy selection in gene space
sel_gl <- run_redundancy_selection(
  gene_list      = gene_list,
  gene_sets      = gene_sets,
  universe       = genes,
  initial_scores = setNames(enrich$info_bits, enrich$set),
  seed           = 42
)

sel_gl[, c("step", "pathway", "marginal_mi", "conditional_mi",
           "redundancy_ratio")]

## ----gsea-lane----------------------------------------------------------------
# Gene-level statistics: correlation with phenotype
gene_stats <- setNames(
  apply(expr, 2, function(g) cor(g, phenotype)),
  genes
)

batch <- run_batch_gsea(
  gene_sets  = gene_sets,
  gene_stats = gene_stats,
  n_perm     = 500,
  adaptive   = TRUE,
  seed       = 42
)

batch$summary[, c("pathway", "enrichment_score", "p_value", "padj",
                   "leading_edge_size")]

## ----gsea-redundancy, eval = nrow(batch$summary[!is.na(batch$summary$padj) & batch$summary$padj < 0.05, ]) >= 2----
# Redundancy on leading edges
sel_gsea <- run_gsea_redundancy_selection(
  gsea_results = batch$results,
  gene_stats   = gene_stats,
  alpha        = 0.05,
  seed         = 42
)

# Classify
classified_gsea <- classify_gsea_selection(sel_gsea)
classified_gsea[, c("pathway", "signal_class", "action")]

## ----adaptive-detail----------------------------------------------------------
res <- run_enrichment(
  gene_list     = gene_sets$immune_response,
  gene_stats    = gene_stats,
  method        = "gsea",
  n_perm        = 500,
  adaptive      = TRUE,
  adaptive_mode = "decision",
  alpha         = 0.05,
  seed          = 42,
  keep_trace    = TRUE
)

inf <- res$inference$res
cat(sprintf(
  "ES = %.3f\np = %.4f\nCI = [%.4f, %.4f]\nPermutations used: %d\nDecision: %s\n",
  res$enrichment_score, inf$p,
  inf$p_ci_low, inf$p_ci_high,
  inf$B, inf$decision_alpha
))

## ----session------------------------------------------------------------------
sessionInfo()

