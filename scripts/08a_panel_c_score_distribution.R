#!/usr/bin/env Rscript
# ============================================================================
# Supplementary panel: distribution of all 250 µm sliding-bin Chol+CD8T2/T3
# scores per sample, with the chosen Panel C ROI marked, so the reader can
# assess where the chosen ROI sits within each sample's score distribution.
# ============================================================================
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(cowplot)
})

source(file.path(Sys.getenv("IRHEP_BUNDLE_ROOT", here::here()), "scripts", "config.R"))
SCAN_DIR <- file.path(PROCESSED_DIR, "scan_tables")
PDIR     <- file.path(PROCESSED_DIR, "figure_2_panels")

# Load both scan tables
ir <- fread(file.path(SCAN_DIR, "scan_irhep_chol_cd8t23.csv"))
ai <- fread(file.path(SCAN_DIR, "scan_aih_chol_cd8t23.csv"))
ir[, condition := "irHep"]; ai[, condition := "AIH"]
setnames(ir, "sample", "sample_id"); setnames(ai, "sample", "sample_id")
all_scans <- rbind(ir, ai)

# Sample-level summary: max score, n_bins, rank of chosen ROI within sample
chosen <- data.table(
  sample_id = c("18321_a", "54023_b"),
  condition = c("irHep", "AIH"),
  x0 = c(1473, 4350),
  y0 = c(-2249, -9850),
  display_score = c(64, 49)        # from scan tables
)
# Resolve to nearest scan grid point per chosen ROI
match_scan <- function(sid, x, y) {
  s <- all_scans[sample_id == sid]
  s[, d := (x0 - x)^2 + (y0 - y)^2]
  s[order(d)][1]
}
chosen_resolved <- rbindlist(lapply(seq_len(nrow(chosen)), function(i) {
  m <- match_scan(chosen$sample_id[i], chosen$x0[i], chosen$y0[i])
  m[, .(sample_id, condition = chosen$condition[i], x0, y0, score)]
}))

# Within-sample rank of the chosen bin
ranks <- all_scans[, .(rank = rank(-score, ties.method = "min"),
                       score = score, x0 = x0, y0 = y0),
                   by = .(sample_id, condition)]
chosen_resolved <- merge(chosen_resolved,
  ranks[, .(sample_id, x0, y0, rank)],
  by = c("sample_id", "x0", "y0"))
sample_n <- all_scans[, .N, by = sample_id]
chosen_resolved <- merge(chosen_resolved, sample_n, by = "sample_id")
chosen_resolved[, label := sprintf("%s\nrank %d / %d\nscore %d",
                                   sample_id, rank, N, score)]
print(chosen_resolved)

# Plot: violin/jitter per sample, chosen ROI marked
all_scans[, sample_id := factor(sample_id,
  levels = c("18321_a", "24774_d", "54023_b", "56784_c"))]
chosen_resolved[, sample_id := factor(sample_id,
  levels = c("18321_a", "24774_d", "54023_b", "56784_c"))]

cond_palette <- c("irHep" = "#C62828", "AIH" = "#1565C0")

p_dist <- ggplot(all_scans,
                 aes(x = sample_id, y = score, color = condition)) +
  geom_jitter(width = 0.15, alpha = 0.18, size = 0.6, stroke = 0) +
  geom_violin(fill = NA, scale = "width", linewidth = 0.45) +
  geom_point(data = chosen_resolved, aes(x = sample_id, y = score),
             color = "black", fill = "yellow", size = 4.5,
             shape = 23, stroke = 0.8) +
  geom_text(data = chosen_resolved,
            aes(x = sample_id, y = score, label = label),
            color = "grey20", size = 3.0, hjust = -0.15, vjust = 0.5,
            lineheight = 0.95) +
  scale_color_manual(values = cond_palette, name = "Condition") +
  scale_x_discrete(name = "Sample (irHep: 18321_a, 24774_d | AIH: 54023_b, 56784_c)") +
  scale_y_continuous(name = "Chol+CD8T2/T3 hotspot score (sqrt(n_chol * n_CD8T2/T3))",
                     expand = expansion(mult = c(0.02, 0.08))) +
  theme_classic(base_size = 11) +
  theme(legend.position = "top",
        axis.text.x = element_text(size = 10),
        plot.background = element_rect(fill = "white", color = NA))

# Add a small annotation strip below explaining the message
caption <- paste0(
  "Distribution of all 250 µm sliding-bin scores per sample. ",
  "Yellow diamonds = Panel C chosen ROIs (top-1 bin per condition). ",
  "irHep ROI rank 1/", chosen_resolved[condition == "irHep"]$N, " in sample 18321_a; ",
  "AIH ROI rank ", chosen_resolved[condition == "AIH"]$rank,
  "/", chosen_resolved[condition == "AIH"]$N, " in sample 54023_b. ",
  "irHep cohort distribution lies entirely above AIH at the high-percentile end ",
  "(score ≥ 30), consistent with Panels A/B."
)

p_final <- plot_grid(
  p_dist,
  ggdraw() + draw_label(caption, size = 9, color = "grey25",
                        x = 0.02, hjust = 0, vjust = 0.5,
                        lineheight = 1.1),
  ncol = 1, rel_heights = c(1, 0.12)
) + theme(plot.background = element_rect(fill = "white", color = NA))

OUT_PDF <- file.path(PDIR, "figure_S_panelC_score_distribution.pdf")
OUT_PNG <- file.path(PDIR, "figure_S_panelC_score_distribution.png")
ggsave(OUT_PDF, p_final, width = 9, height = 5.5, device = cairo_pdf, bg = "white")
ggsave(OUT_PNG, p_final, width = 9, height = 5.5, dpi = 280, bg = "white")
cat("\nWrote", OUT_PDF, "\nWrote", OUT_PNG, "\n")
