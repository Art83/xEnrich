# xEnrich: Multi-Scale Enrichment for Mechanistic and Biomarker Discovery

[![R-CMD-check](https://github.com/Art83/xEnrich/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Art83/xEnrich/actions/workflows/R-CMD-check.yaml)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19108669.svg)](https://doi.org/10.5281/zenodo.19108669)
[![Codecov](https://codecov.io/gh/Art83/xEnrich/branch/master/graph/badge.svg)](https://codecov.io/gh/Art83/xEnrich)

**xEnrich** is an R package that translates differential omics results into
biological context — answering not just *which* pathways are altered but
*where* in the body and *in which cell type* those changes are most relevant.

## The problem

Standard enrichment analysis has three persistent weaknesses:

1. **P-values from fixed permutation schemes carry unquantified Monte Carlo
   uncertainty.** You run 1 000 permutations and get p = 0.048 — but how
   confident are you that the true p is below 0.05?

2. **Significant pathway lists conflate independent signals with redundant
   co-annotations of the same biology.** GO:Fatty acid oxidation,
   GO:Fatty acid beta-oxidation, and KEGG:Fatty acid metabolism are not
   three independent findings — they are one finding described three times.

3. **After you know *which* pathways are enriched, you still don't know
   *where in the body* those pathways are active**, or *which cell type*
   is the origin.

## Two tracks, one interface

### Track A — What is altered?

| Function | Purpose |
|----------|---------|
| `run_info_enrichment()` | MI-scored gene-list enrichment (replaces standard ORA) |
| `run_info_assoc()` | Pathway activity vs phenotype association via MI |
| `run_enrichment()` | Classical hypergeometric or GSEA with adaptive permutation and Wilson CIs |
| `run_batch_gsea()` | Batch GSEA across all pathways with adaptive stopping |
| `run_redundancy_selection()` | Gene-space CMI: identifies independent signals |
| `run_assoc_redundancy_selection()` | Sample-space CMI: same, for expression-based associations |
| `run_assoc_pcor_selection()` | Partial correlation alternative for small n |
| `run_gsea_redundancy_selection()` | Leading-edge CMI for GSEA results |
| `classify_assoc_selection()` | Plain-language classification: report / mention / dismiss |
| `classify_gsea_selection()` | Same, for GSEA lane |
| `compute_pathway_activity()` | Per-sample pathway scores for downstream analysis |
| `summarize_selection()` | Variance-explained decomposition (R² per pathway) |

### Track B — Where in the body?

| Function | Purpose |
|----------|---------|
| `enrich_by_hpa()` | Tissue enrichment via HPA RNA data |
| `enrich_by_hpa_ihc()` | Tissue enrichment via HPA IHC protein data |
| `enrich_by_tabula()` | Cell-type enrichment via Tabula Sapiens scRNA |
| `build_tissue_sets()` | Construct tissue gene sets from reference data |
| `run_location_pipeline()` | All-in-one: organ → cell type |
| `summarize_expression()` | Explore reference data before choosing cutoffs |

### Convergence

The gene list bridges both tracks:

```
gene list → Track A: run_info_assoc()          → WHAT  (pathway + biology)
         → Track B: run_location_pipeline()   → WHERE (organ + cell type)
         ↓
         Together: "solute transport / proximal tubule cells / kidney"
```

## Installation

```r
# install.packages("remotes")
remotes::install_github("Art83/xEnrich")
```

Reference data (Human Protein Atlas and Tabula Sapiens) are downloaded on
first use from [Zenodo](https://doi.org/10.5281/zenodo.19106109) and cached
locally via `load_reference()`.

## Quick start

### Pathway association with redundancy selection

```r
library(xEnrich)

# expr: samples × genes matrix; phenotype: numeric or factor
assoc <- run_info_assoc(
  expr      = expr_matrix,
  phenotype = phenotype_vec,
  gene_sets = reactome_sets,
  n_perm    = 1000,
  seed      = 42
)

# Which significant pathways are independent?
sel <- run_assoc_redundancy_selection(
  assoc_results = assoc,
  expr          = expr_matrix,
  phenotype     = phenotype_vec,
  gene_sets     = reactome_sets,
  seed          = 42
)

# Plain-language classification
classified <- classify_assoc_selection(sel)

# Variance-explained decomposition
result <- summarize_selection(classified, expr_matrix, phenotype_vec, reactome_sets)
# "2 pathways explain 34% of phenotype variance"

# Per-sample activity for downstream models
act <- compute_pathway_activity(expr_matrix, reactome_sets[sel$pathway])
```

### Gene-list enrichment

```r
res <- run_info_enrichment(
  gene_list = my_de_genes,
  gene_sets = reactome_sets,
  universe  = all_measured_genes,
  n_perm    = 1000,
  seed      = 42
)

sel <- run_redundancy_selection(
  gene_list      = my_de_genes,
  gene_sets      = reactome_sets,
  universe       = all_measured_genes,
  initial_scores = setNames(res$info_bits, res$set),
  seed           = 42
)

plot_gains(sel)
```

### Adaptive GSEA with confidence intervals

```r
batch <- run_batch_gsea(
  gene_sets  = reactome_sets,
  gene_stats = log2fc_vector,
  adaptive   = TRUE,
  seed       = 42
)

# batch$summary includes padj, perms_used, and adaptive savings
```

### Location pipeline

```r
hpa_rna <- load_reference("rna_hpa", source = "hpa")[[1]]
tabula  <- load_reference("kidney",  source = "tabula")[[1]]

result <- run_location_pipeline(
  gene_list   = my_de_genes,
  hpa_data    = hpa_rna,
  tabula_data = list(kidney = tabula),
  hpa_type    = "rna",
  seed        = 42
)

print(result)     # organ-level summary
summary(result)   # flat data frame of top cell types
```

## Vignette

For a full walkthrough with synthetic data, see:

```r
vignette("xEnrich")
```

## Citation

If you use xEnrich in published work, please cite:

> Shvetcov A (2026). *xEnrich: Multi-scale enrichment with information theory
> and location context.* R package version 0.1.0.
> https://doi.org/10.5281/zenodo.19108669

```bibtex
@software{xEnrich,
  author  = {Shvetcov, Artur},
  title   = {{xEnrich}: Multi-scale enrichment with information theory
             and location context},
  year    = {2026},
  version = {0.1.0},
  doi     = {10.5281/zenodo.19108669},
  url     = {https://github.com/Art83/xEnrich}
}
```

## License

MIT © Artur Shvetcov
