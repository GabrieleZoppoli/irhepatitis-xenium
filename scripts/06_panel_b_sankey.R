#!/usr/bin/env Rscript
# ============================================================================
# Figure 2 — Panel B: differential cell–cell-communication Sankey.
#
# Reads `cellchat_coherent_interactions.csv` (top coherent CellChat flows
# per condition produced by `05_panel_b_cellchat.R`), filters out flows
# involving clusters that fail the minimum-cell-count rule (n < 10 in any
# sample), takes the top-K coherent flows per direction by |log2FC|, draws
# the alluvial diagram, and writes:
#   fig2_b_sankey.{pdf,png}        — the alluvium (no inline legend).
#   fig2_b_legend_block.{pdf,png}  — a single "Enriched in" mini-plot for
#                                    the composite assembler.
# ============================================================================
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(ggalluvial); library(scales)
  library(cowplot)
})

source(file.path(Sys.getenv("IRHEP_BUNDLE_ROOT", here::here()), "scripts", "config.R"))
source(file.path(SCRIPTS_DIR, "02_cluster_config.R"))

PDIR <- file.path(PROCESSED_DIR, "figure_2_panels")
dir.create(PDIR, recursive = TRUE, showWarnings = FALSE)

CC_CSV <- if (file.exists(file.path(PDIR, "cellchat_coherent_interactions.csv"))) {
  file.path(PDIR, "cellchat_coherent_interactions.csv")
} else {
  file.path(PROCESSED_DIR, "cellchat", "cellchat_coherent_interactions.csv")
}
OUT_PDF     <- file.path(PDIR, "fig2_b_sankey.pdf")
OUT_PNG     <- file.path(PDIR, "fig2_b_sankey.png")
OUT_LEG_PDF <- file.path(PDIR, "fig2_b_legend_block.pdf")
OUT_LEG_PNG <- file.path(PDIR, "fig2_b_legend_block.png")

cc <- fread(CC_CSV)
cc[, abs_lfc := abs(log2FC)]

# Minimum-cell-count filter for CCC reliability. Per-sample mean communication
# probability is unstable for clusters with n < 10 in any sample because near-
# zero denominators inflate log2FC. Ig+ hepatocyte (cluster 20) fails this
# rule in sample 24774_d (12 / 2 / 243 / 478 cells); no other active cluster
# fails the filter.
low_n_clusters <- c("Ig+ hepatocyte")
cc_clean <- cc[!source %in% low_n_clusters & !target %in% low_n_clusters]

# Top K coherent flows per condition by |log2FC|.
K <- 10
top_ir <- cc_clean[direction == "irHep"][order(-abs_lfc)][1:min(K, .N)]
top_ai <- cc_clean[direction == "AIH"  ][order(-abs_lfc)][1:min(K, .N)]
m <- unique(rbind(top_ir, top_ai))
m <- m[!is.na(source) & !is.na(target) & !is.na(pathway_label)]

cat(sprintf("Displayed: top %d per direction by |log2FC| (irHep=%d, AIH=%d, total=%d)\n",
            K, nrow(top_ir), nrow(top_ai), nrow(m)))
cat("\n--- Displayed interactions ---\n")
print(m[, .(source, target, pathway_label, direction, log2FC)])

# MHC-I direction normalisation: antigen presentation goes presenter -> CD8 T.
# CellChat sometimes labels CD8 T as source; swap only in that case so MHC-I
# arrows are biologically oriented.
cd8_cells <- c("CD8+ T cell 1", "CD8+ T cell 2", "CD8+ T cell 3")
needs_swap <- m$pathway_label %in% c("MHC-I_CD8A","MHC-I_CD8B") &
              m$source %in% cd8_cells
m[needs_swap, c("source","target") := .(target, source)]

short_cell <- c(
  "Cholangiocyte" = "Chol", "CD8+ T cell 1" = "CD8T1",
  "CD8+ T cell 2" = "CD8T2", "CD8+ T cell 3" = "CD8T3",
  "CD4+ T cell" = "CD4", "Vascular endothelial cell" = "VascEC",
  "Sinusoidal endothelial cell" = "SinEC", "Hepatocyte 1" = "Hep1",
  "Hepatocyte 2" = "Hep2", "Hepatocyte 3" = "Hep3",
  "IgM+ B cell" = "IgM+B", "IgG+ B cell" = "IgG+B",
  "IgG+ plasmablast / plasma cell" = "PlasmaB"
)
# Receptor-only pathway labels keep stratum text legible. Ligand family is
# implicit (CXCL9/10/11 -> CXCR3, CCL3/4/5 -> CCR1, etc.) and is reported
# fully in Methods / Supplementary Table S1.
short_pw <- c(
  "MHC-I_CD8A" = "MHC-I", "MHC-I_CD8B" = "MHC-I",
  "VEGF" = "VEGF",
  "CXCL_CXCR3" = "CXCR3", "CXCL_CXCR4" = "CXCR4",
  "CXCL_CXCR2" = "CXCR2", "CXCL_CXCR6" = "CXCR6",
  "CCL_CCR1" = "CCR1", "CCL_CCR4" = "CCR4",
  "FN1" = "FN1", "CD86" = "CD86", "ICOS" = "ICOS",
  "CSF" = "CSF", "CD40" = "CD40",
  "FASLG" = "FASLG", "CEACAM" = "CEACAM"
)

cell_levels <- c("Chol", "VascEC", "Hep1", "Hep2", "Hep3",
                 "SinEC", "IgM+B", "PlasmaB",
                 "CD8T1", "CD8T2", "CD8T3", "CD4")
m[, Source  := factor(short_cell[source],
  levels = intersect(cell_levels, unique(short_cell[source])))]
m[, Pathway := factor(short_pw[pathway_label],
  levels = intersect(c("MHC-I","VEGF","CXCR3","CXCR4","CXCR2","CXCR6",
                       "CCR1","CCR4","CSF","FN1",
                       "CD86","ICOS","CD40","FASLG","CEACAM"),
                     unique(short_pw[pathway_label])))]
m[, Target  := factor(short_cell[target],
  levels = intersect(cell_levels, unique(short_cell[target])))]
stopifnot(
  "Source NA"  = all(!is.na(m$Source)),
  "Target NA"  = all(!is.na(m$Target)),
  "Pathway NA" = all(!is.na(m$Pathway))
)

m[, flow_color := ifelse(direction == "irHep", "#C62828", "#1565C0")]
m[, flow_w     := scales::rescale(abs_lfc, to = c(1.5, 5.5))]
m[, flow_id    := seq_len(.N)]

wide <- as.data.frame(
  m[, .(flow_id, Source, Pathway, Target, direction, flow_color, flow_w)])

p <- ggplot(wide,
            aes(y = flow_w, axis1 = Source, axis2 = Pathway, axis3 = Target,
                fill = flow_color)) +
  geom_alluvium(width = 1/4, alpha = 0.72, color = NA,
                knot.pos = 0.3, curve_type = "quintic") +
  geom_stratum(width = 1/4, fill = "grey95", color = "grey30", linewidth = 0.35) +
  geom_text(stat = "stratum",
            aes(label = after_stat(stratum)),
            size = 2.9, fontface = "bold") +
  scale_fill_identity(guide = "none") +
  scale_x_discrete(limits = c("Source", "Pathway", "Target"),
                   expand = c(0.08, 0.08), position = "top") +
  theme_void(base_size = 11) +
  theme(axis.text.x = element_text(size = 10.5, face = "bold",
                                   margin = margin(b = 4)),
        plot.margin = margin(6, 6, 6, 6),
        plot.background = element_rect(fill = "white", color = NA))

ggsave(OUT_PDF, p, width = 7.8, height = 5.4, device = cairo_pdf, bg = "white")
ggsave(OUT_PNG, p, width = 7.8, height = 5.4, dpi = 300, bg = "white")
cat("\nWrote", OUT_PDF, "\nWrote", OUT_PNG, "\n")

# ----------------------------------------------------------------------------
# Legend strip: "Enriched in" colour-key only. The MHC-I source/target
# convention is documented in Supplementary Note 2; raw and redrawn
# directions are reported in Supplementary Table S5.
# ----------------------------------------------------------------------------
p_b_leg_enriched <- ggplot() +
  annotate("text", x = 0.18, y = 1.65,
           label = "Enriched in",
           size = 4.5, fontface = "bold", hjust = 0) +
  annotate("tile", x = 0.18, y = 1.05, width = 0.09, height = 0.32,
           fill = "#1565C0", color = "grey20", linewidth = 0.3) +
  annotate("text", x = 0.28, y = 1.05,
           label = "AIH", hjust = 0, size = 4.0) +
  annotate("tile", x = 0.18, y = 0.55, width = 0.09, height = 0.32,
           fill = "#C62828", color = "grey20", linewidth = 0.3) +
  annotate("text", x = 0.28, y = 0.55,
           label = "SR-irHep", hjust = 0, size = 4.0) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0.0, 2), expand = FALSE) +
  theme_void() +
  theme(plot.margin = margin(8, 6, 8, 6))

# Centre the mini-plot under the Sankey via NULL pad columns; bottom NULL row
# reserves whitespace so cowplot's aspect-fit cannot bleed descenders.
b_legend_inner <- cowplot::plot_grid(NULL, p_b_leg_enriched, NULL,
                                     ncol = 3,
                                     rel_widths = c(1, 1, 1))
b_legend_block <- cowplot::plot_grid(b_legend_inner, NULL,
                                     ncol = 1, rel_heights = c(1, 0.18))

ggsave(OUT_LEG_PDF, b_legend_block, width = 7.8, height = 2.20,
       device = cairo_pdf, bg = "white")
ggsave(OUT_LEG_PNG, b_legend_block, width = 7.8, height = 2.20,
       dpi = 320, bg = "white")
cat("Wrote", OUT_LEG_PDF, "\n")
cat("Wrote", OUT_LEG_PNG, "\n")

cat("\nDone.\n")
