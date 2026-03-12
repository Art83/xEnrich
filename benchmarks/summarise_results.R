# =============================================================================
# Benchmark results visualisation — xEnrich
#
# Produces publication-quality figures for Lane A and Lane C separately.
# Lane A and C have different grid structures and cannot be combined.
#
# Figures saved to benchmarks/results/figures/
#
# Run from package root after both benchmark scripts have completed:
#   source("benchmarks/summarise_results.R")
# =============================================================================

library(ggplot2)

dir.create("benchmarks/results/figures", showWarnings = FALSE, recursive = TRUE)

METHOD_COLOURS <- c(
  "xEnrich (MI)"   = "#2166AC",
  "xEnrich (pcor)" = "#67A9CF",
  "Baseline"       = "#B2182B",
  "Jaccard"        = "#4DAC26",
  "Cor. cluster"   = "#E08214"
)

theme_bench <- function() {
  theme_bw(base_size = 11) +
    theme(
      panel.grid.minor  = element_blank(),
      strip.background  = element_rect(fill = "grey92", colour = NA),
      legend.position   = "bottom",
      legend.title      = element_blank()
    )
}

save_fig <- function(p, filename, width = 7, height = 5) {
  path <- file.path("benchmarks/results/figures", filename)
  png(path, width = width, height = height, units = "in", res = 300)
  print(p)
  dev.off()
  cat(sprintf("Saved: %s\n", path))
}


# =============================================================================
# LANE A
# =============================================================================

lane_A <- read.csv("benchmarks/results/lane_A_summary.csv",
                   stringsAsFactors = FALSE)

lane_A$method_label <- factor(
  lane_A$method,
  levels = c("xEnrich", "baseline_ORA", "jaccard"),
  labels = c("xEnrich", "Baseline", "Jaccard")
)

lane_A$shadow_label <- paste0("overlap = ", lane_A$shadow_overlap)
lane_A$shadows_label <- paste0(lane_A$n_shadows, " shadows")

cat("\n=== LANE A: F1 by method ===\n")
print(aggregate(f1_mean ~ method_label, data = lane_A,
                FUN = function(x) round(mean(x, na.rm = TRUE), 3)))

# ---- Figure A1: F1 vs shadow_overlap, faceted by n_shadows -----------------
p_A1 <- ggplot(lane_A,
               aes(x = shadow_overlap, y = f1_mean,
                   colour = method_label, group = method_label)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  geom_ribbon(aes(ymin = f1_mean - f1_sd,
                  ymax = f1_mean + f1_sd,
                  fill = method_label),
              alpha = 0.15, colour = NA) +
  facet_wrap(~ shadows_label) +
  scale_colour_manual(values = METHOD_COLOURS) +
  scale_fill_manual(values   = METHOD_COLOURS) +
  scale_x_continuous(breaks = c(0.40, 0.55, 0.65, 0.75, 0.85)) +
  labs(
    title   = "Lane A: redundancy suppression in gene lists",
    x       = "Shadow-signal Jaccard overlap",
    y       = "F1 score (mean \u00b1 SD, 50 replicates)",
    caption = paste("Vertical dashed line marks Jaccard clustering threshold (0.70).",
                    "signal_fraction averaged across 0.5, 0.7, 0.9.")
  ) +
  geom_vline(xintercept = 0.70, linetype = "dashed",
             colour = "grey40", linewidth = 0.5) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_bench()

save_fig(p_A1, "lane_A_f1_vs_overlap.png", width = 8, height = 4)

# ---- Figure A2: F1 vs signal_fraction, faceted by shadow_overlap -----------
p_A2 <- ggplot(lane_A,
               aes(x = signal_fraction, y = f1_mean,
                   colour = method_label, group = method_label)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  facet_wrap(~ shadow_label) +
  scale_colour_manual(values = METHOD_COLOURS) +
  scale_x_continuous(breaks = c(0.5, 0.7, 0.9)) +
  labs(
    title = "Lane A: effect of gene list composition",
    x     = "Fraction of signal genes in gene list",
    y     = "F1 score (mean, 50 replicates)"
  ) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_bench()

save_fig(p_A2, "lane_A_f1_vs_signal_fraction.png", width = 8, height = 4)


# =============================================================================
# LANE C
# =============================================================================

lane_C <- read.csv("benchmarks/results/lane_C_summary.csv",
                   stringsAsFactors = FALSE)

# Separate null scenario
lane_C_null   <- lane_C[lane_C$scenario == 4, ]
lane_C_signal <- lane_C[lane_C$scenario != 4, ]

lane_C_signal$method_label <- factor(
  lane_C_signal$method,
  levels = c("xEnrich", "xEnrich_pcor", "baseline", "cor_cluster"),
  labels = c("xEnrich (MI)", "xEnrich (pcor)", "Baseline", "Cor. cluster")
)

lane_C_signal$scenario_label <- factor(
  lane_C_signal$scenario,
  levels = 1:3,
  labels = c("1: Clean separation\n(2 signals, 4 shadows)",
             "2: Asymmetric signals\n(strong + weak)",
             "3: High redundancy\n(1 signal, 8 shadows)")
)

lane_C_null$method_label <- factor(
  lane_C_null$method,
  levels = c("xEnrich", "xEnrich_pcor", "baseline", "cor_cluster"),
  labels = c("xEnrich (MI)", "xEnrich (pcor)", "Baseline", "Cor. cluster")
)

cat("\n=== LANE C: F1 by method and scenario ===\n")
print(aggregate(f1_mean ~ method_label + scenario_label,
                data = lane_C_signal,
                FUN = function(x) round(mean(x, na.rm = TRUE), 3)))

cat("\n=== LANE C: False positives in null scenario (n_selected) ===\n")
print(lane_C_null[, c("method", "n_samples", "n_selected_mean", "n_selected_sd")])

# ---- Figure C1: F1 vs n_samples, faceted by scenario -----------------------
p_C1 <- ggplot(lane_C_signal,
               aes(x = n_samples, y = f1_mean,
                   colour = method_label, group = method_label)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  geom_ribbon(aes(ymin = f1_mean - f1_sd,
                  ymax = f1_mean + f1_sd,
                  fill = method_label),
              alpha = 0.15, colour = NA) +
  facet_wrap(~ scenario_label) +
  scale_colour_manual(values = METHOD_COLOURS) +
  scale_fill_manual(values   = METHOD_COLOURS) +
  scale_x_continuous(breaks = c(100, 200, 500)) +
  labs(
    title   = "Lane C: expression-based pathway association",
    x       = "Sample size",
    y       = "F1 score (mean \u00b1 SD, 20 replicates)",
    caption = "signal_strength = 3.0 (fixed)"
  ) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_bench()

save_fig(p_C1, "lane_C_f1_vs_nsamples.png", width = 9, height = 4)

# ---- Figure C2: Precision-recall by method and scenario --------------------
p_C2 <- ggplot(lane_C_signal,
               aes(x = recall_mean, y = precision_mean,
                   colour = method_label, shape = factor(n_samples))) +
  geom_point(size = 3, alpha = 0.8) +
  facet_wrap(~ scenario_label) +
  scale_colour_manual(values = METHOD_COLOURS) +
  scale_shape_manual(values = c(16, 17, 15),
                     name = "n_samples") +
  labs(
    title = "Lane C: precision vs recall tradeoff",
    x     = "Recall",
    y     = "Precision"
  ) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", colour = "grey60") +
  theme_bench()

save_fig(p_C2, "lane_C_precision_recall.png", width = 9, height = 4)

# ---- Figure C3: False positive rate (null scenario) ------------------------
p_C3 <- ggplot(lane_C_null,
               aes(x = n_samples, y = n_selected_mean,
                   colour = method_label, group = method_label)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  geom_ribbon(aes(ymin = pmax(0, n_selected_mean - n_selected_sd),
                  ymax = n_selected_mean + n_selected_sd,
                  fill = method_label),
              alpha = 0.15, colour = NA) +
  scale_colour_manual(values = METHOD_COLOURS) +
  scale_fill_manual(values   = METHOD_COLOURS) +
  scale_x_continuous(breaks = c(100, 200, 500)) +
  labs(
    title   = "Lane C: false positive rate (null scenario)",
    x       = "Sample size",
    y       = "Mean pathways selected (truth = 0)",
    caption = "Any selection under the null is a false positive."
  ) +
  coord_cartesian(ylim = c(0, NA)) +
  theme_bench()

save_fig(p_C3, "lane_C_null_fp_rate.png", width = 5, height = 4)

# ---- Figure C4: Runtime MI vs pcor -----------------------------------------
runtime_C <- lane_C_signal[
  lane_C_signal$method %in% c("xEnrich", "xEnrich_pcor"), ]
runtime_C$method_label <- factor(
  runtime_C$method,
  levels = c("xEnrich", "xEnrich_pcor"),
  labels = c("xEnrich (MI)", "xEnrich (pcor)")
)

p_C4 <- ggplot(runtime_C,
               aes(x = n_samples, y = runtime_s_mean,
                   colour = method_label, group = method_label)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  geom_ribbon(aes(ymin = pmax(0, runtime_s_mean - runtime_s_sd),
                  ymax = runtime_s_mean + runtime_s_sd,
                  fill = method_label),
              alpha = 0.15, colour = NA) +
  facet_wrap(~ scenario_label) +
  scale_colour_manual(values = METHOD_COLOURS) +
  scale_fill_manual(values   = METHOD_COLOURS) +
  scale_x_continuous(breaks = c(100, 200, 500)) +
  labs(
    title   = "Lane C: runtime — MI vs partial correlation",
    x       = "Sample size",
    y       = "Runtime (seconds per dataset)",
    caption = "pcor redundancy step requires no permutations."
  ) +
  theme_bench()

save_fig(p_C4, "lane_C_runtime_mi_vs_pcor.png", width = 9, height = 4)

cat("\nDone. Figures in benchmarks/results/figures/\n")
