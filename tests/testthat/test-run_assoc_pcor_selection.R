# tests/testthat/test-run_assoc_pcor_selection.R
#
# Tests for run_assoc_pcor_selection(), partial-correlation-based greedy
# forward redundancy selection. Mirrors test-run_assoc_redundancy_selection.R
# so the two functions can be compared directly.

# ---------------------------------------------------------------------------
# Shared fixture (identical to MI test file for direct comparability)
# ---------------------------------------------------------------------------

make_pcor_data <- function(n = 120, seed = 1) {
  set.seed(seed)
  genes  <- paste0("G", 1:300)
  fA     <- rnorm(n)
  fB     <- rnorm(n)
  expr   <- matrix(rnorm(n * 300), nrow = n, dimnames = list(NULL, genes))
  expr[, 1:30]    <- expr[, 1:30]    + 3 * fA
  expr[, 5:35]    <- expr[, 5:35]    + 3 * fA
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

test_that("run_assoc_pcor_selection errors on non-dataframe assoc_results", {
  d <- make_pcor_data()
  expect_error(
    run_assoc_pcor_selection(
      as.matrix(d$assoc), d$expr, d$phenotype, d$gene_sets
    ),
    regexp = "data.frame"
  )
})

test_that("run_assoc_pcor_selection errors on missing required columns", {
  d         <- make_pcor_data()
  bad_assoc <- d$assoc[, c("set", "MI_bits")]
  expect_error(
    run_assoc_pcor_selection(
      bad_assoc, d$expr, d$phenotype, d$gene_sets
    ),
    regexp = "Missing columns"
  )
})

test_that("run_assoc_pcor_selection errors when no pathways significant", {
  d             <- make_pcor_data()
  assoc_null    <- d$assoc
  assoc_null$padj <- 1
  expect_error(
    run_assoc_pcor_selection(
      assoc_null, d$expr, d$phenotype, d$gene_sets, alpha = 0.05
    ),
    regexp = "No significant"
  )
})

test_that("run_assoc_pcor_selection errors for non-numeric phenotype", {
  d <- make_pcor_data()
  expect_error(
    run_assoc_pcor_selection(
      d$assoc, d$expr,
      factor(ifelse(d$phenotype > 0, "case", "ctrl")),
      d$gene_sets
    ),
    regexp = "numeric"
  )
})

test_that("run_assoc_pcor_selection errors on phenotype length mismatch", {
  d <- make_pcor_data()
  expect_error(
    run_assoc_pcor_selection(
      d$assoc, d$expr, d$phenotype[-1], d$gene_sets
    ),
    regexp = "length"
  )
})

# ---------------------------------------------------------------------------
# 2. Output structure
# ---------------------------------------------------------------------------

test_that("run_assoc_pcor_selection returns a data frame", {
  d   <- make_pcor_data()
  res <- suppressMessages(
    run_assoc_pcor_selection(
      d$assoc, d$expr, d$phenotype, d$gene_sets, seed = 1
    )
  )
  expect_s3_class(res, "data.frame")
})

test_that("run_assoc_pcor_selection has required output columns", {
  d   <- make_pcor_data()
  res <- suppressMessages(
    run_assoc_pcor_selection(
      d$assoc, d$expr, d$phenotype, d$gene_sets, seed = 1
    )
  )
  expected <- c("step", "pathway", "set_size", "marginal_mi",
                "marginal_r", "partial_r", "redundancy_ratio",
                "pcor_pvalue", "p_value", "padj")
  expect_true(all(expected %in% colnames(res)))
})

test_that("run_assoc_pcor_selection step column is sequential from 1", {
  d   <- make_pcor_data()
  res <- suppressMessages(
    run_assoc_pcor_selection(
      d$assoc, d$expr, d$phenotype, d$gene_sets, seed = 1
    )
  )
  expect_equal(res$step, seq_len(nrow(res)))
})

test_that("run_assoc_pcor_selection redundancy_ratio is positive", {
  # partial_r / marginal_r can exceed 1 in small samples due to sampling
  # variance in the residuals. No hard upper bound of 1 is enforced.
  d   <- make_pcor_data()
  res <- suppressMessages(
    run_assoc_pcor_selection(
      d$assoc, d$expr, d$phenotype, d$gene_sets, seed = 1
    )
  )
  expect_true(all(res$redundancy_ratio > 0, na.rm = TRUE))
})

test_that("run_assoc_pcor_selection step 1 has NA pcor_pvalue", {
  d   <- make_pcor_data()
  res <- suppressMessages(
    run_assoc_pcor_selection(
      d$assoc, d$expr, d$phenotype, d$gene_sets, seed = 1
    )
  )
  expect_true(is.na(res$pcor_pvalue[1]))
})

test_that("run_assoc_pcor_selection step 1 partial_r equals marginal_r", {
  d   <- make_pcor_data()
  res <- suppressMessages(
    run_assoc_pcor_selection(
      d$assoc, d$expr, d$phenotype, d$gene_sets, seed = 1
    )
  )
  expect_equal(res$partial_r[1], res$marginal_r[1])
})

test_that("run_assoc_pcor_selection result has expected attributes", {
  d   <- make_pcor_data()
  res <- suppressMessages(
    run_assoc_pcor_selection(
      d$assoc, d$expr, d$phenotype, d$gene_sets, seed = 1
    )
  )
  expect_false(is.null(attr(res, "alpha")))
  expect_false(is.null(attr(res, "alpha_pcor")))
  expect_false(is.null(attr(res, "min_pcor")))
  expect_false(is.null(attr(res, "n_sig")))
})

# ---------------------------------------------------------------------------
# 3. Statistical correctness
# ---------------------------------------------------------------------------

test_that("run_assoc_pcor_selection selects at most one of redundant pair", {
  d   <- make_pcor_data()
  res <- suppressMessages(
    run_assoc_pcor_selection(
      d$assoc, d$expr, d$phenotype, d$gene_sets, seed = 1
    )
  )
  a_group <- intersect(res$pathway, c("signal_A", "shadow_A"))
  expect_lte(length(a_group), 1)
})

test_that("run_assoc_pcor_selection first pathway has highest MI among significant", {
  d   <- make_pcor_data()
  res <- suppressMessages(
    run_assoc_pcor_selection(
      d$assoc, d$expr, d$phenotype, d$gene_sets, seed = 1
    )
  )
  sig <- d$assoc[!is.na(d$assoc$padj) & d$assoc$padj <= 0.05, ]
  top <- sig$set[which.max(sig$MI_bits)]
  expect_equal(res$pathway[1], top)
})

test_that("run_assoc_pcor_selection noise pathway not selected", {
  d   <- make_pcor_data()
  res <- suppressMessages(
    run_assoc_pcor_selection(
      d$assoc, d$expr, d$phenotype, d$gene_sets, seed = 1
    )
  )
  expect_false("noise" %in% res$pathway)
})

test_that("run_assoc_pcor_selection respects max_pathways cap", {
  d   <- make_pcor_data()
  res <- suppressMessages(
    run_assoc_pcor_selection(
      d$assoc, d$expr, d$phenotype, d$gene_sets,
      max_pathways = 1, seed = 1
    )
  )
  expect_lte(nrow(res), 1)
})

test_that("run_assoc_pcor_selection marginal_r reflects true phenotype relationship", {
  d   <- make_pcor_data()
  res <- suppressMessages(
    run_assoc_pcor_selection(
      d$assoc, d$expr, d$phenotype, d$gene_sets, seed = 1
    )
  )
  # All selected pathways should have |marginal_r| > min_pcor (default 0.1)
  expect_true(all(abs(res$marginal_r) >= 0.1))
})

# ---------------------------------------------------------------------------
# 4. Agreement with MI version on redundancy decisions
# ---------------------------------------------------------------------------

test_that("pcor and MI versions agree on first selected pathway", {
  d      <- make_pcor_data()
  res_mi <- suppressMessages(
    run_assoc_redundancy_selection(
      d$assoc, d$expr, d$phenotype, d$gene_sets,
      n_perm_gain = 50, seed = 1
    )
  )
  res_pc <- suppressMessages(
    run_assoc_pcor_selection(
      d$assoc, d$expr, d$phenotype, d$gene_sets, seed = 1
    )
  )
  expect_equal(res_mi$pathway[1], res_pc$pathway[1])
})

# ---------------------------------------------------------------------------
# 5. Reproducibility
# ---------------------------------------------------------------------------

test_that("run_assoc_pcor_selection is reproducible with same seed", {
  d    <- make_pcor_data()
  args <- list(d$assoc, d$expr, d$phenotype, d$gene_sets, seed = 7)
  res1 <- suppressMessages(do.call(run_assoc_pcor_selection, args))
  res2 <- suppressMessages(do.call(run_assoc_pcor_selection, args))
  expect_equal(res1$pathway,          res2$pathway)
  expect_equal(res1$partial_r,        res2$partial_r)
  expect_equal(res1$redundancy_ratio, res2$redundancy_ratio)
})

test_that("run_assoc_pcor_selection is deterministic (no permutations in redundancy step)", {
  # pcor redundancy step is closed-form — results should be identical
  # regardless of seed (unlike MI version)
  d    <- make_pcor_data()
  res1 <- suppressMessages(
    run_assoc_pcor_selection(d$assoc, d$expr, d$phenotype, d$gene_sets, seed = 1)
  )
  res2 <- suppressMessages(
    run_assoc_pcor_selection(d$assoc, d$expr, d$phenotype, d$gene_sets, seed = 99)
  )
  # Redundancy decisions should be identical — seed only affects RNG state, not computation
  expect_equal(res1$pathway,    res2$pathway)
  expect_equal(res1$partial_r,  res2$partial_r)
})
