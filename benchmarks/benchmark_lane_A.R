# =============================================================================
# Benchmark runner — Lane A (gene list enrichment)
#
# Compares:
#   xEnrich      run_info_enrichment -> run_redundancy_selection
#   baseline_ORA all significant hypergeometric (no redundancy selection)
#   jaccard      Jaccard clustering at 0.7 (clusterProfiler::simplify surrogate)
#
# Grid parameters that actually vary the difficulty of redundancy suppression:
#   signal_fraction  fraction of signal pathway genes in the gene list
#   n_shadows        number of redundant shadow pathways
#   shadow_overlap   Jaccard similarity between each shadow and its signal
#
# No expression data needed — Lane A is purely about gene list vs gene sets.
# All functions are self-contained inside run_one_A so worker serialisation
# is clean. Universe and pathway constants are exported as global objects.
#
# Run from package root:
#   devtools::load_all()
#   source("benchmarks/benchmark_lane_A.R")
# =============================================================================

library(xEnrich)
library(parallel)

N_CORES <- max(1L, detectCores() - 1L)
cat(sprintf("Using %d cores\n", N_CORES))


# -----------------------------------------------------------------------------
# Fixed universe and pathway regions — exported to workers as globals
# -----------------------------------------------------------------------------

UNIVERSE     <- paste0("G", 1:2000)
SIGNAL_A     <- UNIVERSE[1:100]           # true signal 1
SIGNAL_B     <- UNIVERSE[501:600]         # true signal 2 (well separated)
SHADOW_POOL  <- UNIVERSE[101:500]         # 400 genes for building shadows
NOISE_POOL   <- UNIVERSE[801:1000]        # noise pathways drawn from here
TRUE_SIGNALS <- c("signal_A", "signal_B")


# -----------------------------------------------------------------------------
# Build shadow gene sets with a specified Jaccard overlap to SIGNAL_A
# Each shadow is 100 genes: `overlap_k` from SIGNAL_A + unique genes from pool
# -----------------------------------------------------------------------------

make_shadow_sets <- function(n_shadows, shadow_overlap) {
  # Jaccard J with |A| = |shadow| = 100:  J = k / (200 - k)  =>  k = 200J/(1+J)
  overlap_k  <- round(200 * shadow_overlap / (1 + shadow_overlap))
  unique_per <- 100L - as.integer(overlap_k)

  if (n_shadows * unique_per > length(SHADOW_POOL))
    stop(sprintf(
      "Shadow pool too small for %d shadows at overlap %.2f (need %d genes, have %d)",
      n_shadows, shadow_overlap, n_shadows * unique_per, length(SHADOW_POOL)
    ))

  sets <- lapply(seq_len(n_shadows), function(i) {
    pool_idx     <- ((i - 1L) * unique_per + 1L):(i * unique_per)
    unique_genes <- SHADOW_POOL[pool_idx]
    c(SIGNAL_A[seq_len(overlap_k)], unique_genes)
  })
  setNames(sets, paste0("shadow_A", seq_len(n_shadows)))
}


# -----------------------------------------------------------------------------
# Jaccard-based clustering (clusterProfiler::simplify surrogate)
# Fixed: early exit is a plain if/else, not return() inside system.time
# -----------------------------------------------------------------------------

jaccard_select <- function(gene_sets, pvalues, threshold = 0.7) {
  ord       <- order(pvalues)
  nms       <- names(gene_sets)[ord]
  remaining <- nms
  selected  <- character(0)

  while (length(remaining) > 0L) {
    best     <- remaining[1L]
    selected <- c(selected, best)
    gs_best  <- gene_sets[[best]]

    if (length(remaining) > 1L) {
      jacc <- sapply(remaining[-1L], function(nm) {
        gs  <- gene_sets[[nm]]
        int <- length(intersect(gs_best, gs))
        uni <- length(union(gs_best, gs))
        if (uni == 0L) 0 else int / uni
      })
      to_drop   <- names(jacc)[jacc >= threshold]
      remaining <- setdiff(remaining[-1L], to_drop)
    } else {
      remaining <- character(0)
    }
  }
  selected
}


# -----------------------------------------------------------------------------
# Scoring (redefined here so workers don't need simulate_data.R)
# -----------------------------------------------------------------------------

score_selection <- function(selected, truth) {
  if (length(selected) == 0L) {
    return(c(precision = NA_real_, recall = 0, f1 = NA_real_,
             n_selected = 0L, n_true = length(truth),
             tp = 0L, fp = 0L, fn = length(truth)))
  }
  tp <- length(intersect(selected, truth))
  fp <- length(setdiff(selected, truth))
  fn <- length(setdiff(truth, selected))
  precision <- tp / (tp + fp)
  recall    <- tp / (tp + fn)
  f1 <- if ((precision + recall) > 0) 2 * precision * recall / (precision + recall) else 0
  c(precision = precision, recall = recall, f1 = f1,
    n_selected = length(selected), n_true = length(truth),
    tp = tp, fp = fp, fn = fn)
}


# -----------------------------------------------------------------------------
# Single replicate — no function arguments, uses exported globals directly
# -----------------------------------------------------------------------------

run_one_A <- function(params) {
  signal_fraction <- params$signal_fraction
  n_shadows       <- as.integer(params$n_shadows)
  shadow_overlap  <- params$shadow_overlap
  seed            <- params$seed

  set.seed(seed)

  # Build gene sets
  shadow_sets <- make_shadow_sets(n_shadows, shadow_overlap)
  gene_sets <- c(
    list(signal_A = SIGNAL_A, signal_B = SIGNAL_B),
    shadow_sets,
    list(noise_1 = NOISE_POOL[1:50], noise_2 = NOISE_POOL[51:100])
  )

  # Gene list: signal_fraction of each true signal + 20 random noise genes
  n_sig <- round(100 * signal_fraction)
  gene_list <- unique(c(
    sample(SIGNAL_A, n_sig),
    sample(SIGNAL_B, n_sig),
    sample(UNIVERSE[401:500], 20L)   # neutral zone, not in any pathway
  ))

  truth    <- TRUE_SIGNALS
  universe <- UNIVERSE
  N_u      <- length(universe)
  n_q      <- length(intersect(gene_list, universe))

  # ---- xEnrich --------------------------------------------------------------
  t_xe <- system.time({
    mi_res <- tryCatch(
      suppressMessages(suppressWarnings(
        xEnrich::run_info_enrichment(
          gene_list = gene_list,
          gene_sets = gene_sets,
          universe  = universe,
          n_perm    = 500L,
          seed      = seed
        )
      )),
      error = function(e) NULL
    )

    xe_selected <- character(0)
    if (!is.null(mi_res) && any(!is.na(mi_res$padj_emp) & mi_res$padj_emp <= 0.05)) {
      sel <- tryCatch(
        suppressMessages(suppressWarnings(
          xEnrich::run_redundancy_selection(
            gene_list      = gene_list,
            gene_sets      = gene_sets,
            universe       = universe,
            initial_scores = setNames(mi_res$info_bits, mi_res$set),
            seed           = seed
          )
        )),
        error = function(e) NULL
      )
      if (!is.null(sel) && nrow(sel) > 0L)
        xe_selected <- sel$pathway
    }
  })["elapsed"]

  # ---- Baseline ORA ---------------------------------------------------------
  ora_pvals <- sapply(gene_sets, function(gs) {
    m <- length(intersect(gs, universe))
    k <- length(intersect(gene_list, gs))
    if (m == 0L || n_q == 0L) return(1)
    stats::phyper(k - 1L, m, N_u - m, n_q, lower.tail = FALSE)
  })
  ora_padj          <- stats::p.adjust(ora_pvals, method = "BH")
  baseline_selected <- names(ora_padj[!is.na(ora_padj) & ora_padj <= 0.05])

  # ---- Jaccard clustering ---------------------------------------------------
  t_jc <- system.time({
    jc_selected <- if (length(baseline_selected) == 0L) {
      character(0)
    } else {
      tryCatch(
        jaccard_select(gene_sets[baseline_selected],
                       ora_pvals[baseline_selected],
                       threshold = 0.7),
        error = function(e) character(0)
      )
    }
  })["elapsed"]

  # ---- Score ----------------------------------------------------------------
  data.frame(
    signal_fraction = signal_fraction,
    n_shadows       = n_shadows,
    shadow_overlap  = shadow_overlap,
    seed            = seed,
    method          = c("xEnrich", "baseline_ORA", "jaccard"),
    rbind(
      score_selection(xe_selected,       truth),
      score_selection(baseline_selected, truth),
      score_selection(jc_selected,       truth)
    ),
    runtime_s       = c(t_xe, NA_real_, t_jc),
    stringsAsFactors = FALSE
  )
}


# -----------------------------------------------------------------------------
# Parameter grid — parameters that actually vary the redundancy problem
# -----------------------------------------------------------------------------

grid <- expand.grid(
  signal_fraction = c(0.5, 0.7, 0.9),
  n_shadows       = c(2L, 4L, 8L),
  shadow_overlap = c(0.40, 0.55, 0.65, 0.75, 0.85),
  seed            = seq_len(50L),
  stringsAsFactors = FALSE
)

grid_list <- lapply(seq_len(nrow(grid)), function(i) as.list(grid[i, ]))
cat(sprintf("Grid: %d conditions\n", nrow(grid)))


# -----------------------------------------------------------------------------
# Parallel execution
# -----------------------------------------------------------------------------

cl <- makeCluster(N_CORES)

pkg_path <- "D:/AZ/xEnrich"
clusterExport(cl, varlist = c(
  "UNIVERSE", "SIGNAL_A", "SIGNAL_B", "SHADOW_POOL", "NOISE_POOL",
  "TRUE_SIGNALS", "make_shadow_sets", "jaccard_select",
  "score_selection", "run_one_A", "pkg_path"
))
clusterEvalQ(cl, pkgload::load_all(pkg_path, quiet = TRUE))

cat("Running...\n")
t_total <- system.time({
  results_list <- parLapply(cl, grid_list, function(params) {
    tryCatch(
      run_one_A(params),
      error = function(e) NULL
    )
  })
})["elapsed"]

stopCluster(cl)
cat(sprintf("Done in %.1f min\n", t_total / 60))


# -----------------------------------------------------------------------------
# Assemble and save
# -----------------------------------------------------------------------------

results <- do.call(rbind, Filter(Negate(is.null), results_list))
saveRDS(results, "benchmarks/results/lane_A_results.rds")

summary_cols <- c("precision", "recall", "f1", "n_selected", "runtime_s")

agg <- do.call(rbind, lapply(
  split(results, paste(results$method, results$signal_fraction,
                       results$n_shadows, results$shadow_overlap, sep = "_")),
  function(df) {
    row <- df[1L, c("method", "signal_fraction", "n_shadows", "shadow_overlap")]
    for (col in summary_cols) {
      row[[paste0(col, "_mean")]] <- mean(df[[col]], na.rm = TRUE)
      row[[paste0(col, "_sd")]]   <- sd(df[[col]],   na.rm = TRUE)
    }
    row$n_replicates <- nrow(df)
    row
  }
))

rownames(agg) <- NULL
write.csv(agg, "benchmarks/results/lane_A_summary.csv", row.names = FALSE)

cat(sprintf("Raw:     lane_A_results.rds  (%d rows)\n", nrow(results)))
cat(sprintf("Summary: lane_A_summary.csv  (%d rows)\n", nrow(agg)))
