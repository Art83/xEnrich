
# tests/testthat/test_run_enrichment.R

test_that("run_enrichment works for hypergeometric test", {
  gene_list <- c("A", "B", "C")
  universe <- LETTERS[1:20]
  enriched_genes <- c("A", "D", "E", "F")

  result <- run_enrichment(
    gene_list = gene_list,
    universe = universe,
    method = "hypergeometric",
    enriched_genes = enriched_genes
  )

  expect_type(result$p_value, "double")
  expect_true(result$p_value >= 0 && result$p_value <= 1)
  expect_equal(result$method, "hypergeometric")
})


test_that("run_enrichment works for hypergeometric test", {
  gene_list <- paste0("Gene", c(5, 10, 15, 20, 25, 30, 35, 45))
  universe <- paste0("Gene", 1:100)
  enriched_genes <- paste0("Gene", c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100))

  result <- run_enrichment(
    gene_list = gene_list,
    universe = universe,
    method = "hypergeometric",
    enriched_genes = enriched_genes
  )

  expect_type(result$p_value, "double")
  expect_true(round(result$p_value,2) == 0.03)
  expect_equal(result$method, "hypergeometric")
})


test_that("run_enrichment works for GSEA with minimal example", {
  gene_list <- c("gene1", "gene3", "gene5")
  universe <- paste0("gene", 1:10)
  gene_stats <- setNames(rev(seq_along(universe)), universe)

  result <- run_enrichment(
    gene_list = gene_list,
    universe = universe,
    method = "gsea",
    gene_stats = gene_stats,
    n_perm = 0
  )

  expect_type(result$enrichment_score, "double")
  expect_equal(result$method, "gsea")
  expect_true(result$enrichment_score >= -1 && result$enrichment_score <= 1)
})



test_that("run_enrichment handles no overlap", {
  result <- run_enrichment(
    gene_list = c("X", "Y"),
    universe = LETTERS[1:10],
    method = "hypergeometric",
    enriched_genes = c("A", "B", "C")
  )

  expect_true(result$p_value > 0.5)
})


test_that("gsea validation does not affect hypergeometric mode", {
  result <- run_enrichment(
    gene_list = c("A", "B"),
    universe = LETTERS[1:10],
    method = "hypergeometric",
    enriched_genes = c("A", "C"),
    n_perm = 0
  )

  expect_equal(result$method, "hypergeometric")
})
