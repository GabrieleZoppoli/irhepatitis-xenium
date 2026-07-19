#!/usr/bin/env Rscript
# ============================================================================
# Figure 2 — final composite (7.25 × 9.35 in, journal-page format).
#
# Reads the panel PNGs produced by 04a / 06 / 08 (heatmap + legend; Sankey +
# legend; ROI panel C), places them via ggdraw + draw_plot with explicit
# absolute positioning, and writes:
#   figures/main/Figure_2.{pdf,png}
#
# Layout: top row carries Panel A (heatmap + bottom legend) and Panel B
# (Sankey + bottom legend) side-by-side; bottom row carries Panel C
# (full width, enlarged so its right edge aligns with Panel B's).
# Panel labels (A, B, C) at 12 pt bold to match Figure 1.
# ============================================================================
suppressPackageStartupMessages({
  library(cowplot); library(ggplot2); library(png); library(grid)
})

source(file.path(Sys.getenv("IRHEP_BUNDLE_ROOT", here::here()), "scripts", "config.R"))

PDIR     <- file.path(PROCESSED_DIR, "figure_2_panels")
FIG_MAIN <- file.path(FIGURES_DIR, "main")
dir.create(FIG_MAIN, recursive = TRUE, showWarnings = FALSE)

PG <- function(p) ggdraw() + draw_image(readPNG(p))

pA_hm   <- PG(file.path(PDIR, "fig2_a_heatmap_only_panel.png"))
pA_leg  <- PG(file.path(PDIR, "fig2_a_legend_block.png"))
pB_sk   <- PG(file.path(PDIR, "fig2_b_sankey.png"))
pB_leg  <- PG(file.path(PDIR, "fig2_b_legend_block.png"))
pC      <- PG(file.path(PDIR, "fig2_c_rois.png"))

# Canvas dimensions (inches). Panel C slot is sized to match its source
# aspect (14.5/9.4 = 1.543) so the panel fills its slot without horizontal
# whitespace, and its right edge aligns with Panel B's.
H     <- 9.35
W     <- 7.25
BOT_H <- 4.70
GAP   <- 0.45
TOP_H <- H - BOT_H - GAP

LEG_H <- 1.20
A_W <- 0.86 / 2.00 * W
B_W <- W - A_W

y_panelC_top <- BOT_H / H
y_legend_bot <- (BOT_H + GAP) / H
y_legend_top <- (BOT_H + GAP + LEG_H) / H

fig2 <- ggdraw() +
  draw_plot(pC,     x = 0,         y = 0,
            width  = 1,            height = BOT_H / H) +
  draw_plot(pA_hm,  x = 0,         y = y_legend_top,
            width  = A_W / W,      height = (TOP_H - LEG_H) / H) +
  draw_plot(pB_sk,  x = A_W / W,   y = y_legend_top,
            width  = B_W / W,      height = (TOP_H - LEG_H) / H) +
  # Legend strips drawn AFTER bottom row so they sit on top of z-order
  # (prevents any white-background bleed from Panel C above its content).
  draw_plot(pA_leg, x = 0,         y = y_legend_bot,
            width  = A_W / W,      height = LEG_H / H) +
  draw_plot(pB_leg, x = A_W / W,   y = y_legend_bot,
            width  = B_W / W,      height = LEG_H / H) +
  # Panel labels at 12 pt bold (matches Figure 1 label_size).
  draw_label("A", x = 0.005, y = 0.995, fontface = "bold",
             size = 12, hjust = 0, vjust = 1) +
  draw_label("B", x = A_W / W + 0.005, y = 0.995, fontface = "bold",
             size = 12, hjust = 0, vjust = 1) +
  draw_label("C", x = 0.005, y = y_panelC_top - 0.005, fontface = "bold",
             size = 12, hjust = 0, vjust = 1) +
  theme(plot.background = element_rect(fill = "white", color = NA))

ggsave(file.path(FIG_MAIN, "Figure_2.pdf"),
       fig2, width = W, height = H, device = cairo_pdf, bg = "white")
ggsave(file.path(FIG_MAIN, "Figure_2.png"),
       fig2, width = W, height = H, dpi = 320, bg = "white")
cat("Wrote", file.path(FIG_MAIN, "Figure_2.pdf"), "\n")
cat("Wrote", file.path(FIG_MAIN, "Figure_2.png"), "\n")
