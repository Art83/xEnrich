# tests/testthat/test-fixes-and-edge-cases.R
#
# Tests covering:
#   1. .auto_nbins joint-table density fix
#   2. Double residualisation in run_assoc_redundancy_selection
#   3. classify_assoc_selection
#   4. classify_gsea_selection (basic)
#   5. compute_pathway_activity
#   6. Defensive checks (utils_checks.R)
#   7. Edge cases from problems-to-address

# ===========================================================================
# 1. .auto_nbins — joint-table density awareness
# ===========================================================================

test_that(".auto_nbins with continuous phenotype caps for joint table", {
  # n=100, continuous: old code gave 10 (100 cells, 1/cell).
  # New code: sqrt(100/5) = 4.47 -> 4 (16 cells, 6.25/cell).
  nbins <- xEnrich:::.auto_nbins(100, nbins_y = NULL)
  expect_lte(nbins, 5L)
  expect_gte(nbins, 2L)
  # Verify: n / nbins^2 >= 5
  expect_gte(100 / nbins^2, 5)
})

test_that(".auto_nbins with binary phenotype allows more bins", {
  # n=100, binary (nbins_y=2): floor(100/10) = 10.
  # Joint table is 10*2=20 cells, 5/cell. OK.
  nbins <- xEnrich:::.auto_nbins(100, nbins_y = 2L)
  expect_gte(nbins, 4L)
  # Should allow more bins than the continuous case
  nbins_cont <- xEnrich:::.auto_nbins(100, nbins_y = NULL)
  expect_gte(nbins, nbins_cont)
})

test_that(".auto_nbins at n=500 gives reasonable bin count", {
  nbins <- xEnrich:::.auto_nbins(500, nbins_y = NULL)
  # sqrt(500/5) = 10. Rice = ceiling(2*500^(1/3)) = 16. min(10,16,20) = 10.
  expect_equal(nbins, 10L)
  expect_gte(500 / nbins^2, 5)
})

test_that(".auto_nbins never returns less than 2", {
  # Very small n
  expect_gte(xEnrich:::.auto_nbins(10, nbins_y = NULL), 2L)
  expect_gte(xEnrich:::.auto_nbins(5,  nbins_y = NULL), 2L)
})

test_that(".auto_nbins caps at 20", {
  expect_lte(xEnrich:::.auto_nbins(100000, nbins_y = NULL), 20L)
  expect_lte(xEnrich:::.auto_nbins(100000, nbins_y = 2L),   20L)
})

test_that("run_info_assoc uses joint-aware binning for binary phenotype", {
  set.seed(1)
  n     <- 80
  genes <- paste0("G", 1:100)
  expr  <- matrix(rnorm(n * 100), nrow = n, dimnames = list(NULL, genes))
  pheno <- rep(c(0L, 1L), each = n / 2)
  gs    <- list(pw = genes[1:20])

  res <- suppressMessages(suppressWarnings(
    run_info_assoc(expr, pheno, gs, n_perm = 50, seed = 1)
  ))
  # With binary phenotype (nbins_y=2), nbins should be larger than
  # the continuous-phenotype case
  nbins_used <- res$nbins_used[1]
  nbins_cont <- xEnrich:::.auto_nbins(n, nbins_y = NULL)
  expect_gte(nbins_used, nbins_cont)
})


# ===========================================================================
# 2. Double residualisation — suppressor artifacts
# ===========================================================================

test_that("double residualisation selects independent signals correctly", {
  set.seed(42)
  n     <- 150
  genes <- paste0("G", 1:300)
  fA    <- rnorm(n)
  fB    <- rnorm(n)
  expr  <- matrix(rnorm(n * 300), nrow = n, dimnames = list(NULL, genes))
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
    run_info_assoc(expr, phenotype, gene_sets, n_perm = 300, seed = 42)
  )

  sel <- suppressMessages(
    run_assoc_redundancy_selection(
      assoc, expr, phenotype, gene_sets,
      n_perm_gain = 50, seed = 42
    )
  )

  # Should not select both signal_A and shadow_A (redundant pair)
  a_group <- intersect(sel$pathway, c("signal_A", "shadow_A"))
  expect_lte(length(a_group), 1)

  # Should select signal_B (independent from A)
  expect_true("signal_B" %in% sel$pathway)

  # Noise should not be selected
  expect_false("noise" %in% sel$pathway)

  # All ratios positive, cumulative_bits increasing
  expect_true(all(sel$redundancy_ratio > 0))
  if (nrow(sel) >= 2) {
    expect_true(all(diff(sel$cumulative_bits) > 0))
  }
})

test_that("cumulative_bits is monotonically increasing", {
  set.seed(1)
  n     <- 120
  genes <- paste0("G", 1:300)
  fA    <- rnorm(n)
  fB    <- rnorm(n)
  expr  <- matrix(rnorm(n * 300), nrow = n, dimnames = list(NULL, genes))
  expr[, 1:30]    <- expr[, 1:30]    + 3 * fA
  expr[, 200:230] <- expr[, 200:230] + 3 * fB
  phenotype <- 0.7 * fA + 0.5 * fB + rnorm(n, sd = 0.3)

  gene_sets <- list(
    signal_A = genes[1:30],
    shadow_A = genes[5:35],
    signal_B = genes[200:230],
    noise    = genes[150:170]
  )

  assoc <- suppressMessages(
    run_info_assoc(expr, phenotype, gene_sets, n_perm = 200, seed = 1)
  )
  sel <- suppressMessages(
    run_assoc_redundancy_selection(
      assoc, expr, phenotype, gene_sets,
      n_perm_gain = 50, seed = 1
    )
  )

  if (nrow(sel) >= 2) {
    expect_true(all(diff(sel$cumulative_bits) > 0))
  }
})


# ===========================================================================
# 3. classify_assoc_selection
# ===========================================================================

test_that("classify_assoc_selection adds expected columns", {
  set.seed(42)
  n     <- 120
  genes <- paste0("G", 1:300)
  fA    <- rnorm(n)
  fB    <- rnorm(n)
  expr  <- matrix(rnorm(n * 300), nrow = n, dimnames = list(NULL, genes))
  expr[, 1:30]    <- expr[, 1:30]    + 3 * fA
  expr[, 200:230] <- expr[, 200:230] + 3 * fB
  phenotype <- 0.7 * fA + 0.5 * fB + rnorm(n, sd = 0.3)

  gene_sets <- list(
    signal_A = genes[1:30],
    shadow_A = genes[5:35],
    signal_B = genes[200:230]
  )

  assoc <- suppressMessages(
    run_info_assoc(expr, phenotype, gene_sets, n_perm = 200, seed = 42)
  )
  sel <- suppressMessages(
    run_assoc_redundancy_selection(
      assoc, expr, phenotype, gene_sets,
      n_perm_gain = 50, seed = 42
    )
  )

  classified <- classify_assoc_selection(sel)

  expect_true("marginal_strength" %in% colnames(classified))
  expect_true("signal_class"      %in% colnames(classified))
  expect_true("action"            %in% colnames(classified))
  expect_s3_class(classified$marginal_strength, "factor")
  expect_s3_class(classified$signal_class,      "factor")
})

test_that("classify_assoc_selection errors on bad threshold ordering", {
  fake <- data.frame(
    pathway = "A", redundancy_ratio = 1, z_score = 10,
    marginal_mi = 1, is_suppressor = FALSE
  )
  # primary <= independent is invalid
  expect_error(
    classify_assoc_selection(fake, ratio_primary = 0.3, ratio_independent = 0.5),
    regexp = "ratio_primary"
  )
})

test_that("classify_assoc_selection first step is always primary or independent", {
  fake <- data.frame(
    step = 1L, pathway = "A", redundancy_ratio = 1, z_score = 15,
    marginal_mi = 0.5
  )
  classified <- classify_assoc_selection(fake)
  # Ratio=1 and z=15 -> should be primary
  expect_equal(as.character(classified$signal_class), "primary")
  expect_equal(classified$action, "report")
})


# ===========================================================================
# 4. compute_pathway_activity
# ===========================================================================

test_that("compute_pathway_activity returns correct dimensions", {
  set.seed(1)
  expr <- matrix(rnorm(80 * 100), nrow = 80,
                 dimnames = list(paste0("S", 1:80), paste0("G", 1:100)))
  gs   <- list(pw_A = paste0("G", 1:20), pw_B = paste0("G", 50:70))
  act  <- suppressMessages(compute_pathway_activity(expr, gs))
  expect_equal(nrow(act), 80)
  expect_equal(ncol(act), 2)
  expect_equal(colnames(act), c("pw_A", "pw_B"))
})

test_that("compute_pathway_activity drops pathways below min_genes", {
  set.seed(1)
  expr <- matrix(rnorm(80 * 100), nrow = 80,
                 dimnames = list(NULL, paste0("G", 1:100)))
  gs   <- list(big = paste0("G", 1:20), tiny = paste0("G", 1:2))
  act  <- suppressMessages(compute_pathway_activity(expr, gs, min_genes = 5))
  expect_equal(ncol(act), 1)
  expect_equal(colnames(act), "big")
  expect_true("tiny" %in% attr(act, "dropped"))
})

test_that("compute_pathway_activity mean_z scores have zero mean per pathway", {
  set.seed(1)
  expr <- matrix(rnorm(200 * 50), nrow = 200,
                 dimnames = list(NULL, paste0("G", 1:50)))
  gs   <- list(pw = paste0("G", 1:20))
  act  <- suppressMessages(compute_pathway_activity(expr, gs, score = "mean_z"))
  # Mean of z-scored means should be ~0
  expect_lt(abs(mean(act[, 1])), 0.3)
})

test_that("compute_pathway_activity errors on non-matrix input", {
  df <- data.frame(G1 = rnorm(10), G2 = rnorm(10))
  gs <- list(pw = c("G1", "G2"))
  expect_error(compute_pathway_activity(df, gs), regexp = "numeric matrix")
})


# ===========================================================================
# 5. Defensive checks (utils_checks.R)
# ===========================================================================

test_that(".check_duplicates warns on duplicate genes", {
  expect_warning(
    xEnrich:::.check_duplicates(c("TP53", "BRCA1", "TP53"), "gene_list"),
    regexp = "duplicate"
  )
})

test_that(".check_duplicates is silent on unique input", {
  expect_silent(
    xEnrich:::.check_duplicates(c("TP53", "BRCA1"), "gene_list")
  )
})

test_that(".check_case_mismatch warns on low overlap", {
  expect_warning(
    xEnrich:::.check_case_mismatch(
      paste0("gene", 1:20),     # lowercase
      paste0("GENE", 1:1000)    # uppercase
    ),
    regexp = "case"
  )
})

test_that(".check_case_mismatch is silent on good overlap", {
  genes <- paste0("G", 1:100)
  expect_silent(
    xEnrich:::.check_case_mismatch(genes[1:20], genes)
  )
})

test_that(".remove_zero_variance removes constant columns", {
  mat <- matrix(c(1, 1, 1, 1, 2, 3, 4, 5), nrow = 4,
                dimnames = list(NULL, c("const", "vary")))
  expect_warning(
    result <- xEnrich:::.remove_zero_variance(mat),
    regexp = "zero variance"
  )
  expect_equal(ncol(result), 1L)
  expect_equal(colnames(result), "vary")
})

test_that(".check_matrix_orientation warns when nrow >> ncol", {
  # 1000 rows x 10 cols looks like genes-as-rows (transposed)
  mat <- matrix(rnorm(10000), nrow = 1000, ncol = 10)
  expect_warning(
    xEnrich:::.check_matrix_orientation(mat),
    regexp = "transpose"
  )
})

test_that(".check_matrix_orientation is silent for typical orientation", {
  # 100 rows x 500 cols is normal (samples x genes)
  mat <- matrix(rnorm(50000), nrow = 100, ncol = 500)
  expect_silent(
    xEnrich:::.check_matrix_orientation(mat)
  )
})

test_that(".check_sample_size warns on small n", {
  mat <- matrix(rnorm(30 * 100), nrow = 30, ncol = 100)
  pheno <- rep(c("A", "B"), each = 15)
  expect_warning(
    xEnrich:::.check_sample_size(mat, pheno),
    regexp = "fewer than"
  )
})

test_that(".check_na_expr warns on NAs", {
  mat <- matrix(rnorm(100), nrow = 10, ncol = 10)
  mat[1, 1] <- NA
  expect_warning(
    xEnrich:::.check_na_expr(mat),
    regexp = "NA"
  )
})


# ===========================================================================
# 6. Edge cases from problems-to-address
# ===========================================================================

test_that("run_info_enrichment handles gene_list = entire universe", {
  universe  <- paste0("G", 1:200)
  gene_list <- universe   # degenerate: everything
  gene_sets <- list(pw = universe[1:30])
  # Should not crash; MI should be ~0 since knowing list membership
  # is trivial when list = universe
  res <- suppressMessages(
    run_info_enrichment(gene_list, gene_sets, universe,
                        n_perm = 50, seed = 1)
  )
  expect_equal(nrow(res), 1L)
  # MI should be 0 or very near 0 (every gene is in the list)
  expect_lt(res$info_bits, 0.01)
})

test_that("run_info_enrichment handles gene_list > 50% of universe", {
  universe  <- paste0("G", 1:200)
  gene_list <- universe[1:150]  # 75% of universe
  gene_sets <- list(pw = universe[1:30])
  # Should work, possibly with a warning in future
  res <- suppressMessages(
    run_info_enrichment(gene_list, gene_sets, universe,
                        n_perm = 50, seed = 1)
  )
  expect_equal(nrow(res), 1L)
})

test_that("run_redundancy_selection handles single significant pathway", {
  set.seed(1)
  universe  <- paste0("G", 1:500)
  gene_list <- c(universe[1:30], sample(universe[60:500], 5))
  gene_sets <- list(pathway_A = universe[1:50])

  res <- suppressMessages(
    run_redundancy_selection(
      gene_list, gene_sets, universe, seed = 1, min_gain = 1e-10
    )
  )
  # Should return 1 row, not error
  expect_equal(nrow(res), 1L)
  expect_equal(res$pathway, "pathway_A")
  expect_equal(res$step, 1L)
})

test_that("run_info_assoc gives MI=0 for constant phenotype", {
  set.seed(1)
  expr  <- matrix(rnorm(80 * 50), nrow = 80,
                  dimnames = list(NULL, paste0("G", 1:50)))
  pheno <- rep(1, 80)  # all same value — MI is zero by definition
  gs    <- list(pw = paste0("G", 1:20))
  res <- suppressMessages(suppressWarnings(
    run_info_assoc(expr, pheno, gs, n_perm = 50, seed = 1)
  ))
  # MI should be 0 or near-0; p_value should be non-significant
  expect_lt(res$MI_bits[1], 0.01)
})

test_that("run_enrichment hypergeometric handles zero overlap", {
  universe <- paste0("G", 1:500)
  gl       <- universe[1:20]
  gs       <- universe[100:120]    # no overlap with gl
  res      <- run_enrichment(gl, universe, "hypergeometric", gene_set = gs)
  expect_equal(res$overlap, 0L)
  expect_equal(res$p_value, 1)
})

test_that("run_enrichment GSEA handles minimal overlap", {
  set.seed(1)
  universe   <- paste0("G", 1:500)
  gene_stats <- setNames(rnorm(500), universe)
  gl         <- universe[1:20]
  res <- suppressWarnings(
    run_enrichment(gl, method = "gsea", gene_stats = gene_stats,
                   n_perm = 50, seed = 1)
  )
  expect_true(is.numeric(res$enrichment_score))
  expect_true(is.numeric(res$p_value))
})

test_that("run_info_assoc handles pathway with genes not in expr", {
  set.seed(1)
  expr <- matrix(rnorm(80 * 50), nrow = 80,
                 dimnames = list(NULL, paste0("G", 1:50)))
  pheno <- rnorm(80)
  # pathway_A has 15 genes in expr + 15 not in expr
  gs <- list(pathway_A = c(paste0("G", 1:15), paste0("FAKE", 1:15)))
  res <- suppressMessages(suppressWarnings(
    run_info_assoc(expr, pheno, gs, n_perm = 50, seed = 1, min_size = 5)
  ))
  # set_size should be 15 (only genes in expr), not 30
  expect_equal(res$set_size, 15L)
})

test_that("plot_gains works on empty selection", {
  empty <- data.frame(
    step = integer(0), pathway = character(0), set_size = integer(0),
    marginal_mi = numeric(0), conditional_mi = numeric(0),
    redundancy_ratio = numeric(0), cumulative_bits = numeric(0)
  )
  expect_message(plot_gains(empty), regexp = "Nothing")
})


# ===========================================================================
# 7. summarize_selection — variance explained decomposition
# ===========================================================================

test_that("summarize_selection adds R2 columns", {
  set.seed(42)
  n <- 120; genes <- paste0("G", 1:300)
  fA <- rnorm(n); fB <- rnorm(n)
  expr <- matrix(rnorm(n * 300), nrow = n, dimnames = list(NULL, genes))
  expr[, 1:30]    <- expr[, 1:30]    + 3 * fA
  expr[, 200:230] <- expr[, 200:230] + 3 * fB
  phenotype <- 0.7 * fA + 0.5 * fB + rnorm(n, sd = 0.3)
  gene_sets <- list(
    signal_A = genes[1:30], signal_B = genes[200:230],
    noise    = genes[100:130]
  )

  assoc <- suppressMessages(
    run_info_assoc(expr, phenotype, gene_sets, n_perm = 200, seed = 42)
  )
  sel <- suppressMessages(
    run_assoc_redundancy_selection(
      assoc, expr, phenotype, gene_sets, seed = 42
    )
  )

  result <- summarize_selection(sel, expr, phenotype, gene_sets)

  expect_true("marginal_r2"    %in% colnames(result))
  expect_true("incremental_r2" %in% colnames(result))
  expect_true("cumulative_r2"  %in% colnames(result))
  expect_false(is.null(attr(result, "total_r2")))
  expect_false(is.null(attr(result, "adj_r2")))
  expect_false(is.null(attr(result, "residual_pct")))
})

test_that("summarize_selection R2 values are in [0, 1]", {
  set.seed(1)
  n <- 120; genes <- paste0("G", 1:200)
  fA <- rnorm(n)
  expr <- matrix(rnorm(n * 200), nrow = n, dimnames = list(NULL, genes))
  expr[, 1:30] <- expr[, 1:30] + 3 * fA
  phenotype <- 0.8 * fA + rnorm(n, sd = 0.3)
  gene_sets <- list(signal_A = genes[1:30], noise = genes[100:130])

  assoc <- suppressMessages(
    run_info_assoc(expr, phenotype, gene_sets, n_perm = 200, seed = 1)
  )
  sel <- suppressMessages(
    run_assoc_redundancy_selection(
      assoc, expr, phenotype, gene_sets, seed = 1
    )
  )

  result <- summarize_selection(sel, expr, phenotype, gene_sets)

  expect_true(all(result$marginal_r2    >= 0 & result$marginal_r2    <= 1))
  expect_true(all(result$incremental_r2 >= 0))
  expect_true(all(result$cumulative_r2  >= 0 & result$cumulative_r2  <= 1))
})

test_that("summarize_selection cumulative_r2 is non-decreasing", {
  set.seed(42)
  n <- 120; genes <- paste0("G", 1:300)
  fA <- rnorm(n); fB <- rnorm(n)
  expr <- matrix(rnorm(n * 300), nrow = n, dimnames = list(NULL, genes))
  expr[, 1:30]    <- expr[, 1:30]    + 3 * fA
  expr[, 200:230] <- expr[, 200:230] + 3 * fB
  phenotype <- 0.7 * fA + 0.5 * fB + rnorm(n, sd = 0.3)
  gene_sets <- list(
    signal_A = genes[1:30], signal_B = genes[200:230],
    noise    = genes[100:130]
  )

  assoc <- suppressMessages(
    run_info_assoc(expr, phenotype, gene_sets, n_perm = 200, seed = 42)
  )
  sel <- suppressMessages(
    run_assoc_redundancy_selection(
      assoc, expr, phenotype, gene_sets, seed = 42
    )
  )

  result <- summarize_selection(sel, expr, phenotype, gene_sets)

  if (nrow(result) >= 2) {
    expect_true(all(diff(result$cumulative_r2) >= -1e-10))
  }
})

test_that("summarize_selection signal pathway has higher marginal_r2 than noise", {
  set.seed(1)
  n <- 120; genes <- paste0("G", 1:200)
  fA <- rnorm(n)
  expr <- matrix(rnorm(n * 200), nrow = n, dimnames = list(NULL, genes))
  expr[, 1:30] <- expr[, 1:30] + 3 * fA
  phenotype <- 0.8 * fA + rnorm(n, sd = 0.3)
  gene_sets <- list(signal_A = genes[1:30], noise = genes[100:130])

  assoc <- suppressMessages(
    run_info_assoc(expr, phenotype, gene_sets, n_perm = 200, seed = 1)
  )
  sel <- suppressMessages(
    run_assoc_redundancy_selection(
      assoc, expr, phenotype, gene_sets, seed = 1
    )
  )

  result <- summarize_selection(sel, expr, phenotype, gene_sets)

  # Signal pathway should be step 1 (highest MI) and have decent R2
  expect_equal(result$pathway[1], "signal_A")
  expect_gt(result$marginal_r2[1], 0.1)
})

test_that("summarize_selection errors on missing pathways", {
  fake <- data.frame(pathway = "nonexistent", step = 1L,
                     marginal_mi = 1, conditional_mi = 1, redundancy_ratio = 1)
  expr <- matrix(rnorm(40), nrow = 10, dimnames = list(NULL, paste0("G", 1:4)))
  expect_error(
    summarize_selection(fake, expr, rnorm(10), list(other = "G1")),
    regexp = "not found"
  )
})


# ===========================================================================
# 8. Miller-Madow bias correction
# ===========================================================================

test_that("Miller-Madow corrected MI is less than or equal to plugin MI", {
  # For independent variables, plugin MI has positive bias.
  # Corrected MI should be smaller (closer to zero).
  set.seed(1)
  x <- sample(1:5, 200, replace = TRUE)
  y <- sample(1:5, 200, replace = TRUE)  # independent of x

  mi_corrected <- xEnrich:::.mutual_information(x, y)

  # Manual plugin (no correction)
  tbl      <- table(x, y)
  pxy      <- tbl / length(x)
  px       <- rowSums(pxy)
  py       <- colSums(pxy)
  outer_pp <- outer(px, py)
  terms    <- ifelse(pxy <= 0, 0, pxy * log2(pxy / outer_pp))
  mi_plugin <- sum(terms)

  expect_lte(mi_corrected, mi_plugin)
})

test_that("Miller-Madow corrected MI is non-negative", {
  set.seed(42)
  x <- sample(1:3, 50, replace = TRUE)
  y <- sample(1:3, 50, replace = TRUE)
  mi <- xEnrich:::.mutual_information(x, y)
  expect_gte(mi, 0)
})

test_that("Miller-Madow correction is negligible for large N in gene space", {
  # For a 2x2 table at N=1000, correction < 0.001 bits
  N <- 1000; n <- 100; m <- 80; k <- 15
  mi <- xEnrich:::.ite_overlap_bits(N, n, m, k)

  # Same computation without correction
  p11 <- k/N; p10 <- (n-k)/N; p01 <- (m-k)/N; p00 <- (N-n-m+k)/N
  p1. <- n/N; p0. <- 1-p1.; p.1 <- m/N; p.0 <- 1-p.1
  term <- function(p, pa, pb) ifelse(p <= 0, 0, p * log2(p / (pa * pb)))
  mi_raw <- term(p11, p1., p.1) + term(p10, p1., p.0) +
    term(p01, p0., p.1) + term(p00, p0., p.0)

  # Difference should be tiny
  expect_lt(abs(mi - mi_raw), 0.002)
})

