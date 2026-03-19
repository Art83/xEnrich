# tests/testthat/test-run_info_assoc.R
#
# Tests for run_info_assoc() -- expression-based pathway association via MI.
# Organised as:
#   1. Input validation
#   2. Output structure
#   3. Statistical correctness
#   4. Edge cases
#   5. Reproducibility

# ---------------------------------------------------------------------------
# Shared fixtures
# ---------------------------------------------------------------------------

make_assoc_data <- function(n = 80, seed = 1) {
  set.seed(seed)
  genes  <- paste0("G", 1:200)
  factor <- rnorm(n)
  expr   <- matrix(rnorm(n * 200), nrow = n, dimnames = list(NULL, genes))
  # Strong signal in pathway_A
  expr[, 1:20] <- expr[, 1:20] + 3 * factor
  phenotype <- 0.8 * factor + rnorm(n, sd = 0.3)
  gene_sets <- list(
    pathway_A = genes[1:20],    # true signal
    pathway_B = genes[101:120], # noise
    pathway_C = genes[151:170]  # noise
  )
  list(expr = expr, phenotype = phenotype, gene_sets = gene_sets,
       genes = genes)
}

# ---------------------------------------------------------------------------
# 1. Input validation
# ---------------------------------------------------------------------------

test_that("run_info_assoc errors when expr is not a numeric matrix", {
  d <- make_assoc_data()
  bad_expr <- as.data.frame(d$expr)
  expect_error(
    run_info_assoc(bad_expr, d$phenotype, d$gene_sets, n_perm = 50, seed = 1),
    regexp = "numeric matrix"
  )
})

test_that("run_info_assoc errors when expr has no column names", {
  d <- make_assoc_data()
  colnames(d$expr) <- NULL
  expect_error(
    run_info_assoc(d$expr, d$phenotype, d$gene_sets, n_perm = 50, seed = 1),
    regexp = "column names"
  )
})

test_that("run_info_assoc errors when phenotype length mismatches expr rows", {
  d <- make_assoc_data()
  expect_error(
    run_info_assoc(d$expr, d$phenotype[-1], d$gene_sets, n_perm = 50, seed = 1),
    regexp = "length"
  )
})

test_that("run_info_assoc errors when gene_sets is not a named list", {
  d <- make_assoc_data()
  expect_error(
    run_info_assoc(d$expr, d$phenotype, unname(d$gene_sets), n_perm = 50, seed = 1),
    regexp = "named list"
  )
})

test_that("run_info_assoc errors when n_perm < 1", {
  d <- make_assoc_data()
  expect_error(
    run_info_assoc(d$expr, d$phenotype, d$gene_sets, n_perm = 0, seed = 1),
    regexp = "n_perm"
  )
})

test_that("run_info_assoc errors when nbins < 2", {
  d <- make_assoc_data()
  expect_error(
    run_info_assoc(d$expr, d$phenotype, d$gene_sets,
                   nbins = 1, n_perm = 50, seed = 1),
    regexp = "nbins"
  )
})

test_that("run_info_assoc warns when nbins gives sparse cells", {
  d <- make_assoc_data(n = 60)   # n=60 avoids sample-size warning, still sparse at nbins=20
  expect_warning(
    suppressMessages(
      run_info_assoc(d$expr, d$phenotype, d$gene_sets,
                     nbins = 20, n_perm = 50, seed = 1)
    ),
    regexp = "expected"
  )
})

test_that("run_info_assoc errors when all gene sets are too small", {
  d    <- make_assoc_data()
  tiny <- list(p1 = d$genes[1:2], p2 = d$genes[3:4])
  expect_error(
    run_info_assoc(d$expr, d$phenotype, tiny,
                   min_size = 10, n_perm = 50, seed = 1),
    regexp = "No gene sets remain"
  )
})

# ---------------------------------------------------------------------------
# 2. Output structure
# ---------------------------------------------------------------------------

test_that("run_info_assoc returns a data frame with correct columns", {
  d   <- make_assoc_data()
  res <- suppressMessages(
    run_info_assoc(d$expr, d$phenotype, d$gene_sets, n_perm = 50, seed = 1)
  )
  expect_s3_class(res, "data.frame")
  expected_cols <- c("set", "set_size", "nbins_used", "MI_bits",
                     "z_score", "p_value", "padj", "null_mean", "null_sd")
  expect_true(all(expected_cols %in% colnames(res)))
})

test_that("run_info_assoc returns one row per gene set passing size filter", {
  d   <- make_assoc_data()
  res <- suppressMessages(
    run_info_assoc(d$expr, d$phenotype, d$gene_sets, n_perm = 50, seed = 1)
  )
  expect_equal(nrow(res), length(d$gene_sets))
})

test_that("run_info_assoc output is sorted by padj ascending", {
  d   <- make_assoc_data()
  res <- suppressMessages(
    run_info_assoc(d$expr, d$phenotype, d$gene_sets, n_perm = 50, seed = 1)
  )
  expect_true(all(diff(res$padj) >= 0))
})

test_that("run_info_assoc p_values are in [0, 1]", {
  d   <- make_assoc_data()
  res <- suppressMessages(
    run_info_assoc(d$expr, d$phenotype, d$gene_sets, n_perm = 50, seed = 1)
  )
  expect_true(all(res$p_value >= 0 & res$p_value <= 1))
  expect_true(all(res$padj    >= 0 & res$padj    <= 1))
})

test_that("run_info_assoc MI_bits are non-negative", {
  d   <- make_assoc_data()
  res <- suppressMessages(
    run_info_assoc(d$expr, d$phenotype, d$gene_sets, n_perm = 50, seed = 1)
  )
  expect_true(all(res$MI_bits >= 0))
})

test_that("run_info_assoc set_size matches actual intersection sizes", {
  d   <- make_assoc_data()
  res <- suppressMessages(
    run_info_assoc(d$expr, d$phenotype, d$gene_sets, n_perm = 50, seed = 1)
  )
  universe       <- colnames(d$expr)
  expected_sizes <- sapply(d$gene_sets, function(gs) length(intersect(gs, universe)))
  actual_sizes   <- setNames(res$set_size, res$set)
  expect_equal(actual_sizes[names(expected_sizes)], expected_sizes)
})

# ---------------------------------------------------------------------------
# 3. Statistical correctness
# ---------------------------------------------------------------------------

test_that("run_info_assoc detects a strong signal pathway", {
  d   <- make_assoc_data(n = 120, seed = 42)
  res <- suppressMessages(
    run_info_assoc(d$expr, d$phenotype, d$gene_sets,
                   n_perm = 200, seed = 42)
  )
  # pathway_A should have the lowest padj
  expect_equal(res$set[1], "pathway_A")
  expect_true(res$padj[1] < 0.05)
})

test_that("run_info_assoc noise pathways are not significant under strong signal", {
  d   <- make_assoc_data(n = 120, seed = 42)
  res <- suppressMessages(
    run_info_assoc(d$expr, d$phenotype, d$gene_sets,
                   n_perm = 200, seed = 42)
  )
  noise_padj <- res$padj[res$set %in% c("pathway_B", "pathway_C")]
  expect_true(all(noise_padj > 0.05))
})

test_that("run_info_assoc MI_bits for signal pathway exceeds noise pathways", {
  d   <- make_assoc_data(n = 120, seed = 42)
  res <- suppressMessages(
    run_info_assoc(d$expr, d$phenotype, d$gene_sets,
                   n_perm = 200, seed = 42)
  )
  mi_signal <- res$MI_bits[res$set == "pathway_A"]
  mi_noise  <- res$MI_bits[res$set %in% c("pathway_B", "pathway_C")]
  expect_true(all(mi_signal > mi_noise))
})

test_that("run_info_assoc z_score is positive for signal pathway", {
  d   <- make_assoc_data(n = 120, seed = 42)
  res <- suppressMessages(
    run_info_assoc(d$expr, d$phenotype, d$gene_sets,
                   n_perm = 200, seed = 42)
  )
  expect_true(res$z_score[res$set == "pathway_A"] > 0)
})

test_that("run_info_assoc size filtering drops pathways outside bounds", {
  d <- make_assoc_data()
  # Gene sets are size 20. min_size=18/max_size=22 keeps them.
  # Add a tiny pathway outside bounds to verify it gets dropped.
  gene_sets_extra <- c(d$gene_sets, list(tiny = d$genes[1:5]))
  res <- suppressMessages(
    run_info_assoc(d$expr, d$phenotype, gene_sets_extra,
                   min_size = 18, max_size = 22, n_perm = 50, seed = 1)
  )
  expect_true(all(res$set_size >= 18 & res$set_size <= 22))
  expect_false("tiny" %in% res$set)
})

test_that("run_info_assoc pc1 score method runs without error", {
  d   <- make_assoc_data()
  expect_no_error(
    suppressMessages(
      run_info_assoc(d$expr, d$phenotype, d$gene_sets,
                     score = "pc1", n_perm = 50, seed = 1)
    )
  )
})

test_that("run_info_assoc accepts factor phenotype", {
  d         <- make_assoc_data()
  pheno_fac <- factor(ifelse(d$phenotype > 0, "case", "ctrl"))
  expect_no_error(
    suppressMessages(
      run_info_assoc(d$expr, pheno_fac, d$gene_sets, n_perm = 50, seed = 1)
    )
  )
})

# ---------------------------------------------------------------------------
# 4. Edge cases
# ---------------------------------------------------------------------------

test_that("run_info_assoc handles universe subset smaller than expr columns", {
  d        <- make_assoc_data()
  universe <- colnames(d$expr)[1:50]
  res      <- suppressMessages(
    run_info_assoc(d$expr, d$phenotype, d$gene_sets,
                   universe = universe, n_perm = 50, seed = 1)
  )
  expect_true(all(res$set_size <= 50))
})

test_that("run_info_assoc drops gene sets with no universe intersection silently", {
  d         <- make_assoc_data()
  ghost_set <- list(pathway_A = d$gene_sets$pathway_A,
                    ghost      = paste0("FAKE", 1:20))
  res       <- suppressMessages(
    run_info_assoc(d$expr, d$phenotype, ghost_set, n_perm = 50, seed = 1)
  )
  expect_false("ghost" %in% res$set)
})

test_that("run_info_assoc with explicit nbins overrides auto-selection", {
  d   <- make_assoc_data()
  res <- suppressMessages(
    run_info_assoc(d$expr, d$phenotype, d$gene_sets,
                   nbins = 5, n_perm = 50, seed = 1)
  )
  expect_true(all(res$nbins_used == 5))
})

# ---------------------------------------------------------------------------
# 5. Reproducibility
# ---------------------------------------------------------------------------

test_that("run_info_assoc is reproducible with same seed", {
  d    <- make_assoc_data()
  res1 <- suppressMessages(
    run_info_assoc(d$expr, d$phenotype, d$gene_sets, n_perm = 50, seed = 99)
  )
  res2 <- suppressMessages(
    run_info_assoc(d$expr, d$phenotype, d$gene_sets, n_perm = 50, seed = 99)
  )
  expect_equal(res1$p_value, res2$p_value)
  expect_equal(res1$MI_bits, res2$MI_bits)
})

test_that("run_info_assoc null distribution varies with seed", {
  d    <- make_assoc_data()
  res1 <- suppressMessages(
    run_info_assoc(d$expr, d$phenotype, d$gene_sets, n_perm = 200, seed = 1)
  )
  res2 <- suppressMessages(
    run_info_assoc(d$expr, d$phenotype, d$gene_sets, n_perm = 200, seed = 2)
  )
  # MI_bits are deterministic (no randomness in scoring)
  expect_equal(res1$MI_bits, res2$MI_bits)
  # Null distribution statistics (null_mean, null_sd) should differ
  # between seeds because permutation draws differ
  expect_false(identical(res1$null_mean, res2$null_mean))
})

test_that("binary integer phenotype works (not discretized to constant)", {
  # Regression test: as.integer(status == "case") is numeric with only
  # two unique values, so it must be treated as categorical (not passed
  # through discretize_equalfreq which collapses it to a constant bin).
  # n = 60 matches bin_pheno length (30 + 30).
  d         <- make_assoc_data(n = 60, seed = 7)
  bin_pheno <- as.integer(c(rep(0L, 30L), rep(1L, 30L)))
  expect_no_error(
    res <- suppressMessages(
      run_info_assoc(d$expr, bin_pheno, d$gene_sets, n_perm = 100, seed = 1)
    )
  )
  expect_gt(max(res$MI_bits), 0)
})
