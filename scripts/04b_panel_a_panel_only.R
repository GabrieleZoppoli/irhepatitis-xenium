#!/usr/bin/env Rscript
# ============================================================================
# Figure 2 — Panel A: differential cell–cell proximity heatmap.
#
# Reads `giotto_proximity_per_sample.csv` (per-sample Delaunay-graph proximity
# enrichment from `04_panel_a_proximity_compute.R`), aggregates per condition
# under the coherent-significance rule (both samples in a condition must agree
# in direction at FDR < 0.10), and writes:
#   fig2_a_heatmap_only_panel.{pdf,png} — heatmap with the colorbar inset
#                                         top-right, axis-label-size title.
#   fig2_a_legend_block.{pdf,png}       — horizontal legend strip placed by
#                                         the composite assembler.
# Cluster naming, ordering and palette come from 02_cluster_config.R.
# ============================================================================
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(scales); library(cowplot); library(patchwork)
})

source(file.path(Sys.getenv("IRHEP_BUNDLE_ROOT", here::here()), "scripts", "config.R"))
source(file.path(SCRIPTS_DIR, "02_cluster_config.R"))

GIOTTO <- file.path(PROCESSED_DIR, "spatial_graph_analysis_giotto")
PDIR   <- file.path(PROCESSED_DIR, "figure_2_panels")
dir.create(PDIR, recursive = TRUE, showWarnings = FALSE)

# ----------------------------------------------------------------------------
# 1. Build per-condition aggregate from per-sample proximity table
# ----------------------------------------------------------------------------
ps <- fread(file.path(GIOTTO, "giotto_proximity_per_sample.csv"))
ps[, sig_attract_strict := !is.na(p.adj_higher) & p.adj_higher < 0.10]
ps[, sig_avoid_strict   := !is.na(p.adj_lower)  & p.adj_lower  < 0.10]

sig_tab <- ps[, .(
  mean_enrich          = mean(enrichm, na.rm = TRUE),
  n_sig_attract_strict = sum(sig_attract_strict, na.rm = TRUE),
  n_sig_avoid_strict   = sum(sig_avoid_strict,   na.rm = TRUE)
), by = .(condition, from_label, to_label)]
sig_tab[, coherent_attract := n_sig_attract_strict >= 1 & n_sig_avoid_strict == 0]
sig_tab[, coherent_avoid   := n_sig_avoid_strict   >= 1 & n_sig_attract_strict == 0]
sig_tab[, sig_strict       := coherent_attract | coherent_avoid]

sig_strict_wide <- dcast(sig_tab[, .(condition, from_label, to_label, sig_strict)],
  from_label + to_label ~ condition, value.var = "sig_strict", fill = FALSE)
setnames(sig_strict_wide, c("irHepatitis", "AIH"), c("sig_irh_strict", "sig_aih_strict"))

per_cond_wide <- dcast(sig_tab[, .(condition, from_label, to_label, mean_enrich)],
  from_label + to_label ~ condition, value.var = "mean_enrich")
setnames(per_cond_wide, c("irHepatitis", "AIH"), c("enrich_irh", "enrich_aih"))
per_cond_wide[, delta_enrichm := enrich_irh - enrich_aih]
per_cond_wide[, abs_delta     := abs(delta_enrichm)]

plot_dt <- merge(per_cond_wide, sig_strict_wide, by = c("from_label","to_label"))
for (c in c("sig_irh_strict", "sig_aih_strict")) plot_dt[is.na(get(c)), (c) := FALSE]
plot_dt[, robust_count_strict := as.integer(sig_irh_strict) + as.integer(sig_aih_strict)]

# Lower-triangle layout: arrange (from, to) so col index <= row index
plot_dt[, idx_a := match(from_label, cluster_name_order)]
plot_dt[, idx_b := match(to_label, cluster_name_order)]
plot_dt[, col_label := ifelse(idx_a <= idx_b, from_label, to_label)]
plot_dt[, row_label := ifelse(idx_a <= idx_b, to_label,   from_label)]

# Short labels (matching Panel B nomenclature)
SHORT_LABELS <- c(
  "Hepatocyte 1" = "Hep1", "Hepatocyte 2" = "Hep2",
  "Hepatocyte 3" = "Hep3", "Hepatocyte 4" = "Hep4",
  "Ig+ hepatocyte" = "Ig+Hep",
  "CD8+ T cell 1" = "CD8T1", "CD8+ T cell 2" = "CD8T2",
  "CD8+ T cell 3" = "CD8T3", "CD4+ T cell" = "CD4T",
  "IgG+ plasmablast / plasma cell" = "IgG+ PB",
  "IgG+ B cell" = "IgG+B", "IgM+ B cell" = "IgM+B",
  "Macrophage 1" = "Mac1", "Macrophage 2" = "Mac2",
  "Stellate cell" = "Stellate",
  "Cholangiocyte" = "Chol",
  "Sinusoidal endothelial cell" = "SinEC",
  "Vascular endothelial cell" = "VascEC"
)
short_order <- unname(SHORT_LABELS[cluster_name_order])
plot_dt[, col_label := factor(SHORT_LABELS[col_label], levels = short_order)]
plot_dt[, row_label := factor(SHORT_LABELS[row_label], levels = short_order)]

# Top-10 differential pairs by |Δ enrichment| among coherent-significant tiles
N_TOP <- 10
plot_dt[, is_top := FALSE]; plot_dt[, top_dir := NA_character_]
if (any(plot_dt$robust_count_strict >= 1)) {
  plot_dt[robust_count_strict >= 1,
          rank_in := frank(-abs_delta, ties.method = "first")]
  plot_dt[!is.na(rank_in) & rank_in <= N_TOP, is_top := TRUE]
  plot_dt[is_top == TRUE,
          top_dir := ifelse(delta_enrichm > 0, "irHep", "AIH")]
}

plot_dt[, robust_factor := factor(robust_count_strict, levels = c(0, 1, 2))]
delta_lim <- c(-max(abs(plot_dt$delta_enrichm), na.rm = TRUE),
                max(abs(plot_dt$delta_enrichm), na.rm = TRUE))

# ----------------------------------------------------------------------------
# 2. Build heatmap. Colorbar legend lives inline (top-right empty triangle);
#    shape / box legends go in a separate bottom strip below.
# ----------------------------------------------------------------------------
p_base <- ggplot(plot_dt,
       aes(x = col_label, y = row_label, fill = delta_enrichm)) +
  geom_tile(color = "grey92", linewidth = 0.12) +
  geom_tile(data = plot_dt[is_top == TRUE],
            aes(x = col_label, y = row_label, color = top_dir),
            fill = NA, linewidth = 0.95, inherit.aes = FALSE,
            show.legend = FALSE) +
  scale_color_manual(values = c("irHep" = "#67001F", "AIH" = "#053061"),
                     na.value = NA, guide = "none") +
  geom_point(data = plot_dt[robust_count_strict >= 1],
             aes(shape = robust_factor),
             size = 1.3, color = "black", stroke = 0.4, fill = "black") +
  scale_shape_manual(
    values = c("0" = NA, "1" = 1, "2" = 16),
    name = "Per-condition\nsig (coherent\nat FDR<0.10)",
    labels = c("1 condition", "Both conditions"),
    breaks = c("1", "2"), drop = FALSE, na.translate = FALSE) +
  scale_y_discrete(limits = rev(levels(plot_dt$row_label)), drop = FALSE) +
  scale_x_discrete(drop = FALSE) +
  coord_equal() +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        axis.text.y = element_text(size = 9),
        panel.grid = element_blank(),
        axis.title = element_blank(),
        legend.background = element_rect(fill = "white", color = NA),
        plot.margin = margin(4, 4, 4, 4))

p_inline <- p_base +
  guides(shape = "none",
         fill = guide_colorbar(direction = "horizontal",
                               title.position = "top",
                               title.hjust = 0.5,
                               barwidth  = unit(2.4, "cm"),
                               barheight = unit(0.32, "cm"),
                               ticks.colour = "grey25",
                               frame.colour = "grey25",
                               frame.linewidth = 0.3)) +
  scale_fill_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
    limits = delta_lim,
    name = "Δ log2 enrichment (SR-irHep − AIH)",
    breaks = pretty(delta_lim, n = 4)) +
  theme(legend.position = c(0.98, 0.96),
        legend.justification = c(1, 1),
        legend.direction = "horizontal",
        legend.title = element_text(size = 8.5, face = "bold"),
        legend.text  = element_text(size = 8))

ggsave(file.path(PDIR, "fig2_a_heatmap_only_panel.pdf"),
       p_inline, width = 6, height = 5.5, device = cairo_pdf, bg = "white")
ggsave(file.path(PDIR, "fig2_a_heatmap_only_panel.png"),
       p_inline, width = 6, height = 5.5, dpi = 320, bg = "white")

# ----------------------------------------------------------------------------
# 3. Build the bottom legend strip (shape + top-10 box legends only; the
#    colorbar moved inline into the heatmap above).
# ----------------------------------------------------------------------------
p_shape <- ggplot() +
  annotate("text", x = 0.5, y = 1.65,
           label = "Per-condition sig",
           size = 4.5, fontface = "bold", hjust = 0.5) +
  annotate("point", x = 0.18, y = 1.05, shape = 1, size = 4.5, stroke = 0.7) +
  annotate("text",  x = 0.28, y = 1.05,
           label = "1 condition", hjust = 0, size = 4.0) +
  annotate("point", x = 0.18, y = 0.55, shape = 16, size = 4.5) +
  annotate("text",  x = 0.28, y = 0.55,
           label = "Both conditions", hjust = 0, size = 4.0) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0.0, 2), expand = FALSE) +
  theme_void() +
  theme(plot.margin = margin(8, 6, 8, 6))

p_box <- ggplot() +
  annotate("text", x = 0.5, y = 1.65,
           label = "Top 10 differential",
           size = 4.5, fontface = "bold", hjust = 0.5) +
  annotate("tile", x = 0.20, y = 1.05, width = 0.09, height = 0.32,
           fill = NA, color = "#67001F", linewidth = 1.4) +
  annotate("text", x = 0.28, y = 1.05,
           label = "SR-irHep-enriched", hjust = 0, size = 4.0) +
  annotate("tile", x = 0.20, y = 0.55, width = 0.09, height = 0.32,
           fill = NA, color = "#053061", linewidth = 1.4) +
  annotate("text", x = 0.28, y = 0.55,
           label = "AIH-enriched", hjust = 0, size = 4.0) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0.0, 2), expand = FALSE) +
  theme_void() +
  theme(plot.margin = margin(8, 6, 8, 6))

# Pad NULL columns flank the mini-plots so the inner grid is centered. The
# trailing NULL row reserves whitespace so cowplot's aspect-fit in the
# composite cannot bleed text descenders into the row below.
legend_block_inner <- cowplot::plot_grid(NULL, p_shape, p_box, NULL,
                                         ncol = 4,
                                         rel_widths = c(0.4, 1, 1, 0.4))
legend_block <- cowplot::plot_grid(legend_block_inner, NULL,
                                   ncol = 1, rel_heights = c(1, 0.18))

ggsave(file.path(PDIR, "fig2_a_legend_block.pdf"),
       legend_block, width = 6.2, height = 2.20, device = cairo_pdf, bg = "white")
ggsave(file.path(PDIR, "fig2_a_legend_block.png"),
       legend_block, width = 6.2, height = 2.20, dpi = 320, bg = "white")

cat("Wrote fig2_a_heatmap_only_panel.{pdf,png}\n")
cat("Wrote fig2_a_legend_block.{pdf,png}\n")
