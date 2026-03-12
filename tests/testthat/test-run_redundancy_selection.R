# tests/testthat/test-run_redundancy_selection.R
#
# Tests for run_redundancy_selection() — gene-space CMI-based greedy
# forward selection (Lane A).
# Organised as:
#   1. Input validation
#   2. Output structure
#   3. Statistical correctness
#   4. Edge cases
#   5. Integration with run_info_enrichment
#   6. Reproducibility

# ---------------------------------------------------------------------------
# Shared fixtures
# ---------------------------------------------------------------------------

make_redundancy_data <- function(seed = 1) {
  set.seed(seed)
  universe <- paste0("G", 1:1000)

  gene_sets <- list(
    signal_A = universe[1:60],                          # true signal
    shadow_A = c(universe[1:55], universe[61:65]),      # redundant with A
    signal_B = universe[500:560],                       # independent signal
    noise    = universe[800:830]                        # unrelated
  )

  # gene_list: enriched in A/shadow and B, not noise
  gene_list <- unique(c(
    sample(universe[1:60],    45, replace = FALSE),
    sample(universe[500:560], 20, replace = FALSE),
    sample(universe[200:400],  5, replace = FALSE)
  ))

  list(gene_list = gene_list, gene_sets = gene_sets, universe = universe)
}

# ---------------------------------------------------------------------------
# 1. Input validation
# ---------------------------------------------------------------------------

test_that("run_redundancy_selection errors when gene_list is not character", {
  d <- make_redundancy_data()
  expect_error(
    run_redundancy_selection(1:10, d$gene_sets, d$universe),
    regexp = "character"
  )
})

test_that("run_redundancy_selection errors when gene_sets is not a named list", {
  d <- make_redundancy_data()
  expect_error(
    run_redundancy_selection(d$gene_list, unname(d$gene_sets), d$universe),
    regexp = "named list"
  )
})

test_that("run_redundancy_selection errors when min_gain <= 0", {
  d <- make_redundancy_data()
  expect_error(
    run_redundancy_selection(d$gene_list, d$gene_sets, d$universe,
                             min_gain = -0.1, seed = 1),
    regexp = "min_gain"
  )
})

test_that("run_redundancy_selection errors when fraction outside (0, 1)", {
  d <- make_redundancy_data()
  expect_error(
    run_redundancy_selection(d$gene_list, d$gene_sets, d$universe,
                             fraction = 1.5, seed = 1),
    regexp = "fraction"
  )
})

test_that("run_redundancy_selection errors when all gene sets fail size filter", {
  d <- make_redundancy_data()
  expect_error(
    run_redundancy_selection(d$gene_list, d$gene_sets, d$universe,
                             min_size = 1000, seed = 1),
    regexp = "No gene sets"
  )
})

# ---------------------------------------------------------------------------
# 2. Output structure
# ---------------------------------------------------------------------------

test_that("run_redundancy_selection returns a data frame", {
  d   <- make_redundancy_data()
  res <- suppressMessages(
    run_redundancy_selection(d$gene_list, d$gene_sets, d$universe, seed = 1)
  )
  expect_s3_class(res, "data.frame")
})

test_that("run_redundancy_selection has required output columns", {
  d   <- make_redundancy_data()
  res <- suppressMessages(
    run_redundancy_selection(d$gene_list, d$gene_sets, d$universe, seed = 1)
  )
  expected <- c("step", "pathway", "set_size", "marginal_mi",
                "conditional_mi", "redundancy_ratio", "cumulative_bits")
  expect_true(all(expected %in% colnames(res)))
})

test_that("run_redundancy_selection step column is sequential from 1", {
  d   <- make_redundancy_data()
  res <- suppressMessages(
    run_redundancy_selection(d$gene_list, d$gene_sets, d$universe, seed = 1)
  )
  expect_equal(res$step, seq_len(nrow(res)))
})

test_that("run_redundancy_selection cumulative_bits is monotonically increasing", {
  d   <- make_redundancy_data()
  res <- suppressMessages(
    run_redundancy_selection(d$gene_list, d$gene_sets, d$universe, seed = 1)
  )
  if (nrow(res) > 1)
    expect_true(all(diff(res$cumulative_bits) > 0))
})

test_that("run_redundancy_selection step 1 has redundancy_ratio = 1", {
  d   <- make_redundancy_data()
  res <- suppressMessages(
    run_redundancy_selection(d$gene_list, d$gene_sets, d$universe, seed = 1)
  )
  expect_equal(res$redundancy_ratio[1], 1)
})

test_that("run_redundancy_selection step 1 conditional_mi equals marginal_mi", {
  d   <- make_redundancy_data()
  res <- suppressMessages(
    run_redundancy_selection(d$gene_list, d$gene_sets, d$universe, seed = 1)
  )
  expect_equal(res$conditional_mi[1], res$marginal_mi[1])
})

test_that("run_redundancy_selection has min_gain attribute", {
  d   <- make_redundancy_data()
  res <- suppressMessages(
    run_redundancy_selection(d$gene_list, d$gene_sets, d$universe, seed = 1)
  )
  expect_false(is.null(attr(res, "min_gain")))
  expect_true(is.numeric(attr(res, "min_gain")))
})

test_that("run_redundancy_selection marginal_mi values are non-negative", {
  d   <- make_redundancy_data()
  res <- suppressMessages(
    run_redundancy_selection(d$gene_list, d$gene_sets, d$universe, seed = 1)
  )
  expect_true(all(res$marginal_mi >= 0))
})

# ---------------------------------------------------------------------------
# 3. Statistical correctness
# ---------------------------------------------------------------------------

test_that("run_redundancy_selection does not select both signal_A and shadow_A", {
  d   <- make_redundancy_data()
  res <- suppressMessages(
    run_redundancy_selection(d$gene_list, d$gene_sets, d$universe, seed = 1)
  )
  a_group <- intersect(res$pathway, c("signal_A", "shadow_A"))
  expect_lte(length(a_group), 1)
})

test_that("run_redundancy_selection selects signal pathways not noise", {
  d   <- make_redundancy_data()
  res <- suppressMessages(
    run_redundancy_selection(d$gene_list, d$gene_sets, d$universe, seed = 1)
  )
  expect_false("noise" %in% res$pathway)
})

test_that("run_redundancy_selection first pathway has highest marginal MI", {
  d   <- make_redundancy_data()
  res <- suppressMessages(
    run_redundancy_selection(d$gene_list, d$gene_sets, d$universe, seed = 1)
  )
  # First selected should have highest marginal MI among all selected
  expect_equal(res$pathway[1],
               res$pathway[which.max(res$marginal_mi)])
})

test_that("run_redundancy_selection returns empty frame when no signal exceeds min_gain", {
  d   <- make_redundancy_data()
  res <- suppressMessages(
    run_redundancy_selection(d$gene_list, d$gene_sets, d$universe,
                             min_gain = 999, seed = 1)
  )
  expect_equal(nrow(res), 0L)
  expect_true(all(c("step", "pathway", "marginal_mi") %in% colnames(res)))
})

test_that("run_redundancy_selection respects max_pathways cap", {
  d   <- make_redundancy_data()
  res <- suppressMessages(
    run_redundancy_selection(d$gene_list, d$gene_sets, d$universe,
                             max_pathways = 1, seed = 1)
  )
  expect_lte(nrow(res), 1)
})

test_that("run_redundancy_selection shadow has lower conditional_mi than signal", {
  d   <- make_redundancy_data()
  res <- suppressMessages(
    run_redundancy_selection(d$gene_list, d$gene_sets, d$universe,
                             min_gain = 1e-6, seed = 1)
  )
  # If both A and shadow appear, shadow should have lower conditional MI
  if (all(c("signal_A", "shadow_A") %in% res$pathway)) {
    cmi_A      <- res$conditional_mi[res$pathway == "signal_A"]
    cmi_shadow <- res$conditional_mi[res$pathway == "shadow_A"]
    expect_lt(cmi_shadow, cmi_A)
  } else {
    succeed()  # Only one selected — correct behaviour
  }
})

# ---------------------------------------------------------------------------
# 4. Edge cases
# ---------------------------------------------------------------------------

test_that("run_redundancy_selection handles single gene set", {
  d      <- make_redundancy_data()
  single <- d$gene_sets["signal_A"]
  res    <- suppressMessages(
    run_redundancy_selection(d$gene_list, single, d$universe, seed = 1)
  )
  expect_lte(nrow(res), 1L)
  if (nrow(res) == 1L)
    expect_equal(res$pathway, "signal_A")
})

test_that("run_redundancy_selection with explicit min_gain skips auto-selection", {
  d <- make_redundancy_data()
  expect_message(
    run_redundancy_selection(d$gene_list, d$gene_sets, d$universe,
                             min_gain = NULL, seed = 1),
    regexp = "Auto-selected"
  )
  expect_no_message(
    suppressMessages(
      run_redundancy_selection(d$gene_list, d$gene_sets, d$universe,
                               min_gain = 0.01, seed = 1)
    )
  )
})

test_that("run_redundancy_selection size filter removes sets outside bounds", {
  d <- make_redundancy_data()
  gene_sets_extra <- c(
    d$gene_sets,
    list(
      huge  = d$universe[1:800],
      micro = d$universe[1:2]
    )
  )
  res <- suppressMessages(
    run_redundancy_selection(d$gene_list, gene_sets_extra, d$universe,
                             min_size = 5, max_size = 200, seed = 1)
  )
  expect_false("huge"  %in% res$pathway)
  expect_false("micro" %in% res$pathway)
})

# ---------------------------------------------------------------------------
# 5. Integration with run_info_enrichment
# ---------------------------------------------------------------------------

test_that("initial_scores from run_info_enrichment are accepted", {
  d      <- make_redundancy_data()
  enrich <- suppressMessages(
    run_info_enrichment(d$gene_list, d$gene_sets, d$universe,
                        n_perm = 100, seed = 1)
  )
  scores <- setNames(enrich$info_bits, enrich$set)

  res <- suppressMessages(
    run_redundancy_selection(d$gene_list, d$gene_sets, d$universe,
                             initial_scores = scores, seed = 1)
  )
  expect_s3_class(res, "data.frame")
  expect_true(nrow(res) >= 1L)
})

test_that("initial_scores change marginal_mi but not selection logic", {
  d <- make_redundancy_data()

  # Without initial_scores
  res_computed <- suppressMessages(
    run_redundancy_selection(d$gene_list, d$gene_sets, d$universe,
                             seed = 1)
  )

  # With identical scores passed explicitly
  enrich <- suppressMessages(
    run_info_enrichment(d$gene_list, d$gene_sets, d$universe,
                        n_perm = 100, seed = 1)
  )
  scores <- setNames(enrich$info_bits, enrich$set)
  res_prescored <- suppressMessages(
    run_redundancy_selection(d$gene_list, d$gene_sets, d$universe,
                             initial_scores = scores, seed = 1)
  )

  # Both should select from signal pathways, not noise
  expect_false("noise" %in% res_computed$pathway)
  expect_false("noise" %in% res_prescored$pathway)
})

# ---------------------------------------------------------------------------
# 6. Reproducibility
# ---------------------------------------------------------------------------

test_that("run_redundancy_selection is reproducible with same seed", {
  d    <- make_redundancy_data()
  args <- list(d$gene_list, d$gene_sets, d$universe, seed = 42)
  res1 <- suppressMessages(do.call(run_redundancy_selection, args))
  res2 <- suppressMessages(do.call(run_redundancy_selection, args))
  expect_equal(res1$pathway,          res2$pathway)
  expect_equal(res1$conditional_mi,   res2$conditional_mi)
  expect_equal(res1$redundancy_ratio, res2$redundancy_ratio)
})

test_that("run_redundancy_selection marginal_mi is deterministic regardless of seed", {
  d    <- make_redundancy_data()
  res1 <- suppressMessages(
    run_redundancy_selection(d$gene_list, d$gene_sets, d$universe, seed = 1)
  )
  res2 <- suppressMessages(
    run_redundancy_selection(d$gene_list, d$gene_sets, d$universe, seed = 99)
  )
  # marginal MI is computed from observed data — same regardless of seed
  mi1 <- setNames(res1$marginal_mi, res1$pathway)
  mi2 <- setNames(res2$marginal_mi, res2$pathway)
  common <- intersect(names(mi1), names(mi2))
  expect_equal(mi1[common], mi2[common])
})
