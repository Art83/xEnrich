# =============================================================================
# Benchmark runner — Lane C (expression + phenotype)
#
# Compares:
#   xEnrich     run_info_assoc -> run_assoc_redundancy_selection
#   baseline    no selection (all significant at alpha=0.05)
#   cor_cluster mean_z activity -> pairwise correlation -> greedy cluster r>0.7
#
# Parallelised over grid rows via parallel::parLapply (Windows-compatible).
# Set N_CORES below. A safe default is detectCores() - 1.
#
# Run from package root:
#   devtools::load_all()
#   source("benchmarks/benchmark_lane_C.R")
# =============================================================================

library(xEnrich)
library(parallel)
source("benchmarks/simulate_data.R")

N_CORES <- max(1L, detectCores() - 1L)
cat(sprintf("Using %d cores\n", N_CORES))


# -----------------------------------------------------------------------------
# Competitor: correlation-based clustering
# -----------------------------------------------------------------------------

cluster_select <- function(Z, phenotype, gene_sets, cor_threshold = 0.7) {
  genes <- colnames(Z)

  activity <- sapply(gene_sets, function(gs) {
    g <- intersect(gs, genes)
    if (length(g) == 0L) return(rep(NA_real_, nrow(Z)))
    rowMeans(Z[, g, drop = FALSE])
  })

  cor_pheno <- apply(activity, 2L, function(a) {
    if (all(is.na(a))) return(NA_real_)
    cor(a, phenotype, use = "complete.obs")
  })

  act_cor   <- cor(activity, use = "pairwise.complete.obs")
  remaining <- names(sort(abs(cor_pheno), decreasing = TRUE))
  remaining <- remaining[!is.na(cor_pheno[remaining])]
  selected  <- character(0)

  while (length(remaining) > 0L) {
    best      <- remaining[1L]
    selected  <- c(selected, best)
    rest      <- setdiff(remaining, best)   # never check best against itself
    if (length(rest) == 0L) {
      remaining <- character(0)
    } else {
      correlated <- names(which(abs(act_cor[best, rest]) >= cor_threshold))
      remaining  <- setdiff(rest, correlated)
    }
  }

  if (length(selected) > 0L) {
    act_cors <- sapply(selected, function(nm) {
      g <- intersect(gene_sets[[nm]], genes)
      if (length(g) == 0L) return(0)
      cor(rowMeans(Z[, g, drop = FALSE]), phenotype)
    })
    selected <- selected[abs(act_cors) > 0.1]
  }

  selected
}


# -----------------------------------------------------------------------------
# Single replicate — must be self-contained for export to workers
# -----------------------------------------------------------------------------

run_one <- function(params, simulate_data_fn, score_selection_fn,
                    cluster_select_fn) {
  scenario        <- params$scenario
  n_samples       <- params$n_samples
  signal_strength <- params$signal_strength
  seed            <- params$seed

  dat <- simulate_data_fn(
    scenario        = scenario,
    n_samples       = n_samples,
    n_genes         = 1000L,
    signal_strength = signal_strength,
    compute_gene_stats = FALSE,
    seed            = seed
  )

  truth     <- dat$true_signals
  gene_sets <- dat$gene_sets
  expr      <- dat$expr
  phenotype <- dat$phenotype
  Z         <- scale(expr)

  # ---- xEnrich --------------------------------------------------------------
  t_xe <- system.time({
    assoc <- tryCatch(
      suppressMessages(suppressWarnings(
        xEnrich::run_info_assoc(
          expr      = expr,
          phenotype = phenotype,
          gene_sets = gene_sets,
          n_perm    = 50L,
          seed      = seed
        )
      )),
      error = function(e) NULL
    )

    xe_selected <- character(0)
    if (!is.null(assoc) && any(!is.na(assoc$padj) & assoc$padj <= 0.05)) {
      sel <- tryCatch(
        suppressMessages(suppressWarnings(
          xEnrich::run_assoc_redundancy_selection(
            assoc_results = assoc,
            expr          = expr,
            phenotype     = phenotype,
            gene_sets     = gene_sets,
            alpha         = 0.05,
            n_perm_gain   = 20L,
            seed          = seed
          )
        )),
        error = function(e) NULL
      )
      if (!is.null(sel) && nrow(sel) > 0L)
        xe_selected <- sel$pathway
    }
  })["elapsed"]

  # ---- Baseline -------------------------------------------------------------
  baseline_selected <- if (!is.null(assoc)) {
    assoc$set[!is.na(assoc$padj) & assoc$padj <= 0.05]
  } else character(0)

  # ---- Correlation clustering -----------------------------------------------
  t_cc <- system.time({
    cc_selected <- tryCatch(
      cluster_select_fn(Z, phenotype, gene_sets, cor_threshold = 0.7),
      error = function(e) character(0)
    )
  })["elapsed"]

  # ---- xEnrich pcor ---------------------------------------------------------
  t_xp <- system.time({
    xp_selected <- character(0)
    if (!is.null(assoc) && any(!is.na(assoc$padj) & assoc$padj <= 0.05)) {
      sel <- tryCatch(
        suppressMessages(suppressWarnings(
          xEnrich::run_assoc_pcor_selection(
            assoc_results = assoc,
            expr          = expr,
            phenotype     = phenotype,
            gene_sets     = gene_sets,
            alpha         = 0.05,
            seed          = seed
          )
        )),
        error = function(e) NULL
      )
      if (!is.null(sel) && nrow(sel) > 0L)
        xp_selected <- sel$pathway
    }
  })["elapsed"]

  # ---- Score ----------------------------------------------------------------
  data.frame(
    scenario        = scenario,
    n_samples       = n_samples,
    signal_strength = signal_strength,
    seed            = seed,
    method          = c("xEnrich", "xEnrich_pcor", "baseline", "cor_cluster"),
    rbind(
      score_selection_fn(xe_selected,       truth),
      score_selection_fn(xp_selected,       truth),
      score_selection_fn(baseline_selected, truth),
      score_selection_fn(cc_selected,       truth)
    ),
    runtime_s       = c(t_xe, t_xp, NA_real_, t_cc),
    stringsAsFactors = FALSE
  )
}


# -----------------------------------------------------------------------------
# Parameter grid
# -----------------------------------------------------------------------------

grid <- expand.grid(
  scenario        = 1:4,
  n_samples       = c(100L, 200L, 500L),
  signal_strength = 3.0,   # fixed — varying this adds nothing above detection floor
  seed            = seq_len(20L),
  stringsAsFactors = FALSE
)

grid_list <- lapply(seq_len(nrow(grid)), function(i) as.list(grid[i, ]))

cat(sprintf("Grid: %d conditions\n", nrow(grid)))


# -----------------------------------------------------------------------------
# To time a single run before committing, call:
#   system.time(run_one(grid_list[[1]], simulate_data, score_selection, cluster_select))
# Then call run_lane_C() to execute the full benchmark.
# -----------------------------------------------------------------------------

run_lane_C <- function(n_cores = N_CORES, pkg = "D:/AZ/xEnrich") {

  cl <- makeCluster(n_cores)
  on.exit(stopCluster(cl))

  clusterExport(cl, varlist = c("simulate_data", "score_selection",
                                "cluster_select", "run_one", "pkg"),
                envir = environment())
  clusterEvalQ(cl, {
    pkgload::load_all(pkg, quiet = TRUE)
    source(file.path(pkg, "benchmarks/simulate_data.R"))
  })

  cat("Running...\n")
  t_total <- system.time({
    results_list <- parLapply(cl, grid_list, function(params) {
      tryCatch(
        run_one(params, simulate_data, score_selection, cluster_select),
        error = function(e) NULL
      )
    })
  })["elapsed"]

  cat(sprintf("Done in %.1f min\n", t_total / 60))

  results <- do.call(rbind, Filter(Negate(is.null), results_list))
  saveRDS(results, "benchmarks/results/lane_C_results.rds")

  summary_cols <- c("precision", "recall", "f1", "n_selected", "runtime_s")

  agg <- do.call(rbind, lapply(
    split(results, paste(results$method, results$scenario,
                         results$n_samples, results$signal_strength, sep = "_")),
    function(df) {
      row <- df[1L, c("method", "scenario", "n_samples", "signal_strength")]
      for (col in summary_cols) {
        row[[paste0(col, "_mean")]] <- mean(df[[col]], na.rm = TRUE)
        row[[paste0(col, "_sd")]]   <- sd(df[[col]],   na.rm = TRUE)
      }
      row$n_replicates <- nrow(df)
      row
    }
  ))

  rownames(agg) <- NULL
  write.csv(agg, "benchmarks/results/lane_C_summary.csv", row.names = FALSE)

  cat(sprintf("Raw:     lane_C_results.rds  (%d rows)\n", nrow(results)))
  cat(sprintf("Summary: lane_C_summary.csv  (%d rows)\n", nrow(agg)))
  invisible(agg)
}

cat("Functions loaded. Time one run first:\n")
cat("  system.time(run_one(grid_list[[1]], simulate_data, score_selection, cluster_select))\n")
cat("Then run the full benchmark:\n")
cat("  run_lane_C()\n")
