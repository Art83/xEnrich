
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
    n_perm = 50
  )

  expect_type(result$enrichment_score, "double")
  expect_equal(result$method, "gsea")
  expect_true(result$enrichment_score >= -1 && result$enrichment_score <= 1)
})



test_that("run_enrichment handles no overlap", {
  expect_error(result <- run_enrichment(
    gene_list = c("X", "Y"),
    universe = LETTERS[1:10],
    method = "hypergeometric",
    enriched_genes = c("A", "B", "C")
  ))
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


test_that("run_enrichment errors for GSEA when n_perm is not positive", {
  gene_list <- c("gene1", "gene3", "gene5")
  universe <- paste0("gene", 1:10)
  gene_stats <- setNames(rev(seq_along(universe)), universe)

  expect_error(
    run_enrichment(
      gene_list = gene_list,
      universe = universe,
      method = "gsea",
      gene_stats = gene_stats,
      n_perm = 0
    ),
    "n_perm must be > 0"
  )
})

test_that("run_enrichment handles GSEA case where all ranked genes are hits", {
  universe <- paste0("gene", 1:5)
  gene_stats <- setNames(c(5, 4, 3, 2, 1), universe)

  result <- run_enrichment(
    gene_list = universe,
    universe = universe,
    method = "gsea",
    gene_stats = gene_stats,
    n_perm = 20,
    seed = 1
  )

  expect_equal(result$method, "gsea")
  expect_true(is.list(result$leading_edge))
  expect_length(result$leading_edge[[1]], 0)
  expect_true(is.finite(result$enrichment_score))
})


test_that("Hypergeometric: Basic functionality and correct calculation", {
  gene_list <- paste0("Gene", c(5, 10, 15, 20, 25, 30, 35, 45))
  universe <- paste0("Gene", 1:100)
  enriched_genes <- paste0("Gene", c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100))

  result <- run_enrichment(
    gene_list = gene_list,
    universe = universe,
    method = "hypergeometric",
    enriched_genes = enriched_genes,
    alternative = "greater"
  )

  expect_type(result$p_value, "double")
  # Based on dhyper/phyper logic, the p-value for k=3, n=8, K=10, N=100 is ~0.031
  expect_equal(round(result$p_value, 3), 0.031)
})

test_that("GSEA: Leading Edge for Negative Enrichment (Drivers at the bottom)", {
  universe <- paste0("gene", 1:10)
  # Strong negative signal: hits are at the bottom of the ranked list
  gene_stats <- setNames(10:1, universe)
  gene_list <- c("gene9", "gene10")

  result <- run_enrichment(
    gene_list = gene_list,
    universe = universe,
    method = "gsea",
    gene_stats = gene_stats,
    n_perm = 100,
    seed = 42
  )

  # For negative ES, leading edge should be genes AT or AFTER the trough
  # gene9 and gene10 should be present
  le <- result$leading_edge[[1]]
  expect_true(all(c("gene9", "gene10") %in% le))
  expect_true(result$enrichment_score < 0)
})

test_that("GSEA: Leading Edge for Positive Enrichment (Drivers at the top)", {
  universe <- paste0("gene", 1:10)
  gene_stats <- setNames(10:1, universe)
  gene_list <- c("gene1", "gene2")

  result <- run_enrichment(
    gene_list = gene_list,
    universe = universe,
    method = "gsea",
    gene_stats = gene_stats,
    n_perm = 100
  )

  le <- result$leading_edge[[1]]
  expect_true(all(c("gene1", "gene2") %in% le))
  expect_true(result$enrichment_score > 0)
})

test_that("Edge Case: All genes in universe are hits", {
  # If every gene is a hit, n_miss is 0.
  # The code should return ES=0 and an empty leading edge rather than crashing.
  universe <- paste0("gene", 1:5)
  gene_stats <- setNames(c(5, 4, 3, 2, 1), universe)

    result <- run_enrichment(
      gene_list = universe,
      universe = universe,
      method = "gsea",
      gene_stats = gene_stats,
      n_perm = 10
    )

  expect_equal(result$enrichment_score, 0)
  expect_length(result$leading_edge[[1]], 0)
})

test_that("Edge Case: gene_stats contain zeros or negatives", {
  # GSEA uses abs(stats), so zeros can be problematic if gsea_weight > 0
  universe <- c("A", "B", "C")
  gene_stats <- setNames(c(0, 0, 0), universe)
  gene_list <- c("A")

  result <- run_enrichment(
    gene_list = gene_list,
    universe = universe,
    method = "gsea",
    gene_stats = gene_stats,
    n_perm = 10
  )

  # Should handle zero weights gracefully (returns 0 score)
  expect_equal(result$enrichment_score, 0)
})

test_that("Adaptive Mode: pvalue convergence", {
  universe <- paste0("gene", 1:100)
  gene_stats <- setNames(rnorm(100), universe)
  gene_list <- sample(universe, 10)

  result <- run_enrichment(
    gene_list = gene_list,
    universe = universe,
    method = "gsea",
    gene_stats = gene_stats,
    adaptive = TRUE,
    adaptive_mode = "pvalue",
    eps = 0.2 # Loose eps for fast test
  )

  expect_true(result$inference$res$converged)
  expect_true("trace" %in% names(result$inference$res) )
})
