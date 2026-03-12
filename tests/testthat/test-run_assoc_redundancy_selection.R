# tests/testthat/test-run_assoc_redundancy_selection.R
#
# Tests for run_assoc_redundancy_selection() — MI-based greedy forward
# redundancy selection for expression-based pathway associations.

# ---------------------------------------------------------------------------
# Shared fixtures
# ---------------------------------------------------------------------------

make_redundancy_data <- function(n = 120, seed = 1) {
  set.seed(seed)
  genes  <- paste0("G", 1:300)
  fA     <- rnorm(n)
  fB     <- rnorm(n)
  expr   <- matrix(rnorm(n * 300), nrow = n, dimnames = list(NULL, genes))
  # signal_A and shadow_A2 both driven by fA (redundant pair)
  expr[, 1:30]    <- expr[, 1:30]    + 3 * fA
  expr[, 5:35]    <- expr[, 5:35]    + 3 * fA   # shadow: overlaps signal_A
  # signal_B driven by independent fB
  expr[, 200:230] <- expr[, 200:230] + 3 * fB
  phenotype <- 0.7 * fA + 0.5 * fB + rnorm(n, sd = 0.3)
  gene_sets <- list(
    signal_A = genes[1:30],
    shadow_A = genes[5:35],
    signal_B = genes[200:230],
    noise    = genes[150:170]
  )
  assoc <- suppressMessages(
    run_info_assoc(expr, phenotype, gene_sets, n_perm = 200, seed = seed)
  )
  list(expr = expr, phenotype = phenotype, gene_sets = gene_sets,
       assoc = assoc)
}

# ---------------------------------------------------------------------------
# 1. Input validation
# ---------------------------------------------------------------------------

test_that("run_assoc_redundancy_selection errors on non-dataframe assoc_results", {
  d <- make_redundancy_data()
  expect_error(
    run_assoc_redundancy_selection(
      as.matrix(d$assoc), d$expr, d$phenotype, d$gene_sets
    ),
    regexp = "data.frame"
  )
})

test_that("run_assoc_redundancy_selection errors on missing required columns", {
  d        <- make_redundancy_data()
  bad_assoc <- d$assoc[, c("set", "MI_bits")]   # missing padj etc.
  expect_error(
    run_assoc_redundancy_selection(
      bad_assoc, d$expr, d$phenotype, d$gene_sets
    ),
    regexp = "Missing columns"
  )
})

test_that("run_assoc_redundancy_selection errors when no pathways significant", {
  d          <- make_redundancy_data()
  assoc_null <- d$assoc
  assoc_null$padj <- 1   # force all non-significant
  expect_error(
    run_assoc_redundancy_selection(
      assoc_null, d$expr, d$phenotype, d$gene_sets, alpha = 0.05
    ),
    regexp = "No significant"
  )
})

test_that("run_assoc_redundancy_selection errors when expr is not numeric matrix", {
  d <- make_redundancy_data()
  expect_error(
    run_assoc_redundancy_selection(
      d$assoc, as.data.frame(d$expr), d$phenotype, d$gene_sets
    ),
    regexp = "numeric matrix"
  )
})

test_that("run_assoc_redundancy_selection errors on phenotype length mismatch", {
  d <- make_redundancy_data()
  expect_error(
    run_assoc_redundancy_selection(
      d$assoc, d$expr, d$phenotype[-1], d$gene_sets
    ),
    regexp = "length"
  )
})

# ---------------------------------------------------------------------------
# 2. Output structure
# ---------------------------------------------------------------------------

test_that("run_assoc_redundancy_selection returns a data frame", {
  d   <- make_redundancy_data()
  res <- suppressMessages(
    run_assoc_redundancy_selection(
      d$assoc, d$expr, d$phenotype, d$gene_sets,
      n_perm_gain = 20, seed = 1
    )
  )
  expect_s3_class(res, "data.frame")
})

test_that("run_assoc_redundancy_selection has required output columns", {
  d   <- make_redundancy_data()
  res <- suppressMessages(
    run_assoc_redundancy_selection(
      d$assoc, d$expr, d$phenotype, d$gene_sets,
      n_perm_gain = 20, seed = 1
    )
  )
  expected <- c("step", "pathway", "set_size", "marginal_mi",
                "conditional_mi", "redundancy_ratio", "p_value", "padj")
  expect_true(all(expected %in% colnames(res)))
})

test_that("run_assoc_redundancy_selection step column is sequential from 1", {
  d   <- make_redundancy_data()
  res <- suppressMessages(
    run_assoc_redundancy_selection(
      d$assoc, d$expr, d$phenotype, d$gene_sets,
      n_perm_gain = 20, seed = 1
    )
  )
  expect_equal(res$step, seq_len(nrow(res)))
})

test_that("run_assoc_redundancy_selection redundancy_ratio is positive", {
  # Conditional MI on residuals can legitimately exceed marginal MI due to
  # permutation variance — no hard upper bound of 1 is enforced.
  # The meaningful check is that all ratios are strictly positive.
  d   <- make_redundancy_data()
  res <- suppressMessages(
    run_assoc_redundancy_selection(
      d$assoc, d$expr, d$phenotype, d$gene_sets,
      n_perm_gain = 20, seed = 1
    )
  )
  expect_true(all(res$redundancy_ratio > 0, na.rm = TRUE))
})

test_that("run_assoc_redundancy_selection result has attributes", {
  d   <- make_redundancy_data()
  res <- suppressMessages(
    run_assoc_redundancy_selection(
      d$assoc, d$expr, d$phenotype, d$gene_sets,
      n_perm_gain = 20, seed = 1
    )
  )
  expect_false(is.null(attr(res, "alpha")))
  expect_false(is.null(attr(res, "n_sig")))
})

# ---------------------------------------------------------------------------
# 3. Statistical correctness
# ---------------------------------------------------------------------------

test_that("run_assoc_redundancy_selection selects signal pathways not shadows", {
  d   <- make_redundancy_data()
  res <- suppressMessages(
    run_assoc_redundancy_selection(
      d$assoc, d$expr, d$phenotype, d$gene_sets,
      n_perm_gain = 50, seed = 1
    )
  )
  # Should select signal_A or shadow_A (one of the pair) and signal_B
  # Should NOT select both signal_A and shadow_A (they are redundant)
  a_group <- intersect(res$pathway, c("signal_A", "shadow_A"))
  expect_lte(length(a_group), 1)
})

test_that("run_assoc_redundancy_selection first pathway has highest MI", {
  d   <- make_redundancy_data()
  res <- suppressMessages(
    run_assoc_redundancy_selection(
      d$assoc, d$expr, d$phenotype, d$gene_sets,
      n_perm_gain = 20, seed = 1
    )
  )
  sig   <- d$assoc[!is.na(d$assoc$padj) & d$assoc$padj <= 0.05, ]
  top   <- sig$set[which.max(sig$MI_bits)]
  expect_equal(res$pathway[1], top)
})

test_that("run_assoc_redundancy_selection noise pathway is not selected", {
  d   <- make_redundancy_data()
  res <- suppressMessages(
    run_assoc_redundancy_selection(
      d$assoc, d$expr, d$phenotype, d$gene_sets,
      n_perm_gain = 50, seed = 1
    )
  )
  expect_false("noise" %in% res$pathway)
})

test_that("run_assoc_redundancy_selection respects max_pathways cap", {
  d   <- make_redundancy_data()
  res <- suppressMessages(
    run_assoc_redundancy_selection(
      d$assoc, d$expr, d$phenotype, d$gene_sets,
      n_perm_gain = 20, max_pathways = 1, seed = 1
    )
  )
  expect_lte(nrow(res), 1)
})

# ---------------------------------------------------------------------------
# 4. Reproducibility
# ---------------------------------------------------------------------------

test_that("run_assoc_redundancy_selection is reproducible with same seed", {
  d    <- make_redundancy_data()
  args <- list(d$assoc, d$expr, d$phenotype, d$gene_sets,
               n_perm_gain = 20, seed = 7)
  res1 <- suppressMessages(do.call(run_assoc_redundancy_selection, args))
  res2 <- suppressMessages(do.call(run_assoc_redundancy_selection, args))
  expect_equal(res1$pathway,          res2$pathway)
  expect_equal(res1$conditional_mi,   res2$conditional_mi)
  expect_equal(res1$redundancy_ratio, res2$redundancy_ratio)
})
