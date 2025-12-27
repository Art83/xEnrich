test_that("ES is identical for fixed and adaptive inference", {

  set.seed(1)

  gene_stats <- rnorm(1000)
  names(gene_stats) <- paste0("g", seq_along(gene_stats))

  gene_list <- sample(names(gene_stats), 50)
  universe  <- names(gene_stats)

  res_fixed <- run_enrichment(
    gene_list = gene_list,
    universe = universe,
    gene_stats = gene_stats,
    method = "gsea",
    n_perm = 500,
    adaptive = FALSE,
    seed = 1
  )

  res_adapt <- run_enrichment(
    gene_list = gene_list,
    universe = universe,
    gene_stats = gene_stats,
    method = "gsea",
    n_perm = 500,
    adaptive = TRUE,
    eps = 0.01,
    seed = 1
  )

  expect_equal(res_fixed$enrichment_score, res_adapt$enrichment_score)
})




test_that("Adaptive p-value matches fixed permutation p-value", {

  set.seed(2)

  gene_stats <- rnorm(1000)
  names(gene_stats) <- paste0("g", seq_along(gene_stats))

  gene_list <- sample(names(gene_stats), 40)
  universe  <- names(gene_stats)

  res_fixed <- run_enrichment(
    gene_list,
    universe,
    gene_stats = gene_stats,
    method = "gsea",
    n_perm = 5000,
    adaptive = FALSE,
    seed = 2
  )

  res_adapt <- run_enrichment(
    gene_list,
    universe,
    gene_stats = gene_stats,
    method = "gsea",
    n_perm = 5000,
    adaptive = TRUE,
    eps = 0.02,
    seed = 2
  )

  expect_true(abs(res_fixed$p_value - res_adapt$p_value) < 0.01)
})


test_that("Adaptive stops early for non-significant gene sets", {

  set.seed(3)

  gene_stats <- rnorm(2000)
  names(gene_stats) <- paste0("g", seq_along(gene_stats))

  # Random gene list → null-like
  gene_list <- sample(names(gene_stats), 20)
  universe  <- names(gene_stats)

  res <- run_enrichment(
    gene_list,
    universe,
    gene_stats = gene_stats,
    method = "gsea",
    n_perm = 10000,
    adaptive = TRUE,
    eps = 0.1,
    seed = 3
  )

  expect_true(res$inference$B < 3000)
  expect_true(res$p_value > 0.1)
})



test_that("Adaptive runs longer for significant enrichment", {

  set.seed(4)

  gene_stats <- sort(rnorm(2000), decreasing = TRUE)
  names(gene_stats) <- paste0("g", seq_along(gene_stats))

  gene_list <- names(gene_stats)[1:30]  # strong signal
  universe  <- names(gene_stats)

  res <- run_enrichment(
    gene_list,
    universe,
    gene_stats = gene_stats,
    method = "gsea",
    n_perm = 20000,
    adaptive = TRUE,
    eps = 0.1,
    seed = 4
  )

  expect_true(res$p_value < 0.05)
  expect_true(res$inference$B > 5000)
})


test_that("Adaptive inference is reproducible with a seed", {

  gene_stats <- rnorm(1500)
  names(gene_stats) <- paste0("g", seq_along(gene_stats))

  gene_list <- sample(names(gene_stats), 30)
  universe  <- names(gene_stats)

  res1 <- run_enrichment(
    gene_list,
    universe,
    gene_stats = gene_stats,
    method = "gsea",
    n_perm = 8000,
    adaptive = TRUE,
    eps = 0.05,
    seed = 123
  )

  res2 <- run_enrichment(
    gene_list,
    universe,
    gene_stats = gene_stats,
    method = "gsea",
    n_perm = 8000,
    adaptive = TRUE,
    eps = 0.05,
    seed = 123
  )

  expect_equal(res1$p_value, res2$p_value)
  expect_equal(res1$inference$B, res2$inference$B)
})


test_that("Inference structure is correct", {

  gene_stats <- rnorm(1000)
  names(gene_stats) <- paste0("g", seq_along(gene_stats))

  gene_list <- sample(names(gene_stats), 20)
  universe  <- names(gene_stats)

  res <- run_enrichment(
    gene_list,
    universe,
    gene_stats = gene_stats,
    method = "gsea",
    n_perm = 500,
    adaptive = TRUE
  )

  expect_true("inference" %in% names(res))
  expect_true(res$inference$method == "adaptive")
  expect_true(all(c("B", "RSE", "trace", "p", "k", "converged", "trace") %in% names(res$inference)))
})
