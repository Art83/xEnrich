# xEnrich 0.0.0.9000

* Initial development version.

## Track A — Information-theoretic enrichment

* `run_info_enrichment()`: MI-scored gene-list enrichment.
* `run_info_assoc()`: pathway activity–phenotype association via MI.
  Joint-table-aware binning ensures stable MI estimates at small sample sizes.
* `run_enrichment()`: hypergeometric test and GSEA with adaptive permutation
  (RSE, decision, and refine modes) and Wilson confidence intervals on p-values.
* `run_batch_gsea()`: batch GSEA with adaptive stopping and BH adjustment.
* `run_redundancy_selection()`: gene-space conditional MI for non-redundant
  pathway selection.
* `run_assoc_redundancy_selection()`: sample-space conditional MI via double
  residualisation — identifies independent pathway–phenotype associations.
* `run_assoc_pcor_selection()`: partial correlation alternative for small n
  or continuous phenotype.
* `run_gsea_redundancy_selection()`: leading-edge-based CMI selection for
  GSEA results.
* `classify_assoc_selection()`, `classify_gsea_selection()`: plain-language
  classification of selected pathways (primary / independent / partial /
  redundant).
* `compute_pathway_activity()`: per-sample pathway scores for downstream
  analysis.
* `summarize_selection()`: variance-explained decomposition — translates MI
  selection into marginal, incremental, and cumulative R-squared values.

## Track B — Location-aware enrichment

* `enrich_by_hpa()`, `enrich_by_hpa_ihc()`: tissue-level enrichment via
  Human Protein Atlas RNA and IHC data.
* `enrich_by_tabula()`: cell-type enrichment via Tabula Sapiens scRNA data.
* `build_tissue_sets()`: construct tissue gene sets from reference data.
* `run_location_pipeline()`: all-in-one organ → cell-type resolution.
* `summarize_expression()`: explore reference data before choosing cutoffs.
* `load_reference()`: download and cache reference datasets from Zenodo.

## MI estimation

* All mutual information estimators (gene-space and sample-space) now apply
  the Miller-Madow bias correction, reducing the positive finite-sample bias
  inherent in plug-in MI estimation. The correction is small at typical sample
  sizes but improves accuracy of absolute MI values.
