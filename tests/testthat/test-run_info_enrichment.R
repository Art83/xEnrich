# tests/testthat/test-run_info_enrichment.R
#
# Tests for run_info_enrichment() — MI-based over-representation analysis.
# Organised as:
#   1. Input validation
#   2. Output structure
#   3. Statistical correctness
#   4. Edge cases
#   5. Reproducibility

# ---------------------------------------------------------------------------
# Shared fixtures
# ---------------------------------------------------------------------------

make_enrich_data <- function(seed = 1) {
  set.seed(seed)
  universe  <- paste0("G", 1:500)
  # gene_list contains all members of pathway_A — perfect overlap
  gene_list <- universe[1:20]
  gene_sets <- list(
    pathway_A = universe[1:20],    # perfect overlap with gene_list
    pathway_B = universe[101:150], # no overlap
    pathway_C = universe[1:10]     # partial overlap (subset of A)
  )
  list(gene_list = gene_list, gene_sets = gene_sets, universe = universe)
}

# ---------------------------------------------------------------------------
# 1. Input validation
# ---------------------------------------------------------------------------

test_that("run_info_enrichment errors when gene_list is not a character vector", {
  d <- make_enrich_data()
  expect_error(
    run_info_enrichment(1:10, d$gene_sets, d$universe, n_perm = 50),
    regexp = "character"
  )
})

test_that("run_info_enrichment errors when gene_sets is not a named list", {
  d <- make_enrich_data()
  expect_error(
    run_info_enrichment(d$gene_list, unname(d$gene_sets), d$universe, n_perm = 50),
    regexp = "named list"
  )
})

test_that("run_info_enrichment errors when n_perm < 1", {
  d <- make_enrich_data()
  expect_error(
    run_info_enrichment(d$gene_list, d$gene_sets, d$universe, n_perm = 0),
    regexp = "n_perm"
  )
})

test_that("run_info_enrichment errors when all gene sets fail size filter", {
  d    <- make_enrich_data()
  tiny <- list(p1 = d$universe[1:2], p2 = d$universe[3:4])
  expect_error(
    run_info_enrichment(d$gene_list, tiny, d$universe,
                        min_size = 10, n_perm = 50),
    regexp = "No gene sets"
  )
})

test_that("run_info_enrichment warns then errors when gene_list has no overlap with universe", {
  d <- make_enrich_data()
  # Function warns about identifier mismatch, then errors as no genes remain
  expect_error(
    suppressWarnings(
      run_info_enrichment(paste0("FAKE", 1:10), d$gene_sets, d$universe,
                          n_perm = 50, seed = 1)
    ),
    regexp = "No genes"
  )
})

# ---------------------------------------------------------------------------
# 2. Output structure
# ---------------------------------------------------------------------------

test_that("run_info_enrichment returns a data frame", {
  d   <- make_enrich_data()
  res <- suppressMessages(
    run_info_enrichment(d$gene_list, d$gene_sets, d$universe,
                        n_perm = 100, seed = 1)
  )
  expect_s3_class(res, "data.frame")
})

test_that("run_info_enrichment returns expected columns", {
  d   <- make_enrich_data()
  res <- suppressMessages(
    run_info_enrichment(d$gene_list, d$gene_sets, d$universe,
                        n_perm = 100, seed = 1)
  )
  expected <- c("set", "overlap", "set_size", "info_bits",
                "p_enrich", "emp_p", "z_score")
  expect_true(all(expected %in% colnames(res)))
})

test_that("run_info_enrichment returns one row per passing gene set", {
  d   <- make_enrich_data()
  res <- suppressMessages(
    run_info_enrichment(d$gene_list, d$gene_sets, d$universe,
                        n_perm = 100, seed = 1)
  )
  expect_equal(nrow(res), length(d$gene_sets))
})

test_that("run_info_enrichment p_enrich and emp_p are in [0, 1]", {
  d   <- make_enrich_data()
  res <- suppressMessages(
    run_info_enrichment(d$gene_list, d$gene_sets, d$universe,
                        n_perm = 100, seed = 1)
  )
  expect_true(all(res$p_enrich >= 0 & res$p_enrich <= 1))
  expect_true(all(res$emp_p    >= 0 & res$emp_p    <= 1))
})

test_that("run_info_enrichment info_bits are non-negative", {
  d   <- make_enrich_data()
  res <- suppressMessages(
    run_info_enrichment(d$gene_list, d$gene_sets, d$universe,
                        n_perm = 100, seed = 1)
  )
  expect_true(all(res$info_bits >= 0))
})

test_that("run_info_enrichment set_size reflects post-universe intersection", {
  d   <- make_enrich_data()
  res <- suppressMessages(
    run_info_enrichment(d$gene_list, d$gene_sets, d$universe,
                        n_perm = 100, seed = 1)
  )
  expected_sizes <- sapply(d$gene_sets,
                           function(gs) length(intersect(gs, d$universe)))
  actual_sizes   <- setNames(res$set_size, res$set)
  expect_equal(actual_sizes[names(expected_sizes)], expected_sizes)
})

# ---------------------------------------------------------------------------
# 3. Statistical correctness
# ---------------------------------------------------------------------------

test_that("run_info_enrichment assigns lowest p_enrich to perfectly overlapping set", {
  d   <- make_enrich_data()
  res <- suppressMessages(
    run_info_enrichment(d$gene_list, d$gene_sets, d$universe,
                        n_perm = 200, seed = 1)
  )
  expect_equal(res$set[which.min(res$p_enrich)], "pathway_A")
})

test_that("run_info_enrichment assigns highest info_bits to perfectly overlapping set", {
  d   <- make_enrich_data()
  res <- suppressMessages(
    run_info_enrichment(d$gene_list, d$gene_sets, d$universe,
                        n_perm = 200, seed = 1)
  )
  expect_equal(res$set[which.max(res$info_bits)], "pathway_A")
})

test_that("run_info_enrichment gives zero overlap to non-overlapping pathway", {
  d   <- make_enrich_data()
  res <- suppressMessages(
    run_info_enrichment(d$gene_list, d$gene_sets, d$universe,
                        n_perm = 100, seed = 1)
  )
  expect_equal(res$overlap[res$set == "pathway_B"], 0)
})

test_that("run_info_enrichment p_enrich=1 for non-overlapping pathway", {
  d   <- make_enrich_data()
  res <- suppressMessages(
    run_info_enrichment(d$gene_list, d$gene_sets, d$universe,
                        n_perm = 100, seed = 1)
  )
  # No overlap = p_enrich should be 1 (one-sided test, no enrichment)
  expect_equal(res$p_enrich[res$set == "pathway_B"], 1)
})

test_that("run_info_enrichment size filter removes large and small sets", {
  d <- make_enrich_data()
  gene_sets_extra <- c(
    d$gene_sets,
    list(
      huge  = d$universe[200:400],  # 201 genes — too large
      micro = d$universe[1:2]       # 2 genes — too small
    )
  )
  res <- suppressMessages(
    run_info_enrichment(d$gene_list, gene_sets_extra, d$universe,
                        min_size = 5, max_size = 100, n_perm = 50, seed = 1)
  )
  expect_false("huge"  %in% res$set)
  expect_false("micro" %in% res$set)
  expect_true(all(res$set_size >= 5 & res$set_size <= 100))
})

test_that("run_info_enrichment uses union of gene_sets as universe when NULL", {
  d <- make_enrich_data()
  expect_message(
    run_info_enrichment(d$gene_list, d$gene_sets, universe = NULL,
                        n_perm = 50, seed = 1),
    regexp = "universe"
  )
})

# ---------------------------------------------------------------------------
# 4. Edge cases
# ---------------------------------------------------------------------------

test_that("run_info_enrichment handles single gene set", {
  d      <- make_enrich_data()
  single <- d$gene_sets["pathway_A"]
  res    <- suppressMessages(
    run_info_enrichment(d$gene_list, single, d$universe,
                        n_perm = 50, seed = 1)
  )
  expect_equal(nrow(res), 1L)
  expect_equal(res$set, "pathway_A")
})

test_that("run_info_enrichment handles gene_list larger than any single set", {
  universe  <- paste0("G", 1:500)
  gene_list <- universe[1:200]
  gene_sets <- list(small_A = universe[1:20], small_B = universe[100:120])
  res <- suppressMessages(
    run_info_enrichment(gene_list, gene_sets, universe,
                        n_perm = 50, seed = 1)
  )
  expect_equal(nrow(res), 2L)
})

# ---------------------------------------------------------------------------
# 5. Reproducibility
# ---------------------------------------------------------------------------

test_that("run_info_enrichment is reproducible with same seed", {
  d    <- make_enrich_data()
  res1 <- suppressMessages(
    run_info_enrichment(d$gene_list, d$gene_sets, d$universe,
                        n_perm = 100, seed = 42)
  )
  res2 <- suppressMessages(
    run_info_enrichment(d$gene_list, d$gene_sets, d$universe,
                        n_perm = 100, seed = 42)
  )
  expect_equal(res1$emp_p,     res2$emp_p)
  expect_equal(res1$info_bits, res2$info_bits)
  expect_equal(res1$z_score,   res2$z_score)
})

test_that("run_info_enrichment info_bits are deterministic regardless of seed", {
  d    <- make_enrich_data()
  res1 <- suppressMessages(
    run_info_enrichment(d$gene_list, d$gene_sets, d$universe,
                        n_perm = 100, seed = 1)
  )
  res2 <- suppressMessages(
    run_info_enrichment(d$gene_list, d$gene_sets, d$universe,
                        n_perm = 100, seed = 99)
  )
  # info_bits computed from observed data — deterministic
  expect_equal(res1$info_bits, res2$info_bits)
  # emp_p derived from permutations — varies with seed
  expect_false(identical(res1$emp_p, res2$emp_p))
})
