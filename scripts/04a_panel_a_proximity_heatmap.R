#!/usr/bin/env Rscript
# ============================================================================
# Figure 2 Panel A: differential proximity heatmap + top-15 lollipop.
# Renders both strict and relaxed coherent-rule variants.
# ============================================================================
#
# INPUT:  giotto_proximity_per_sample.csv (from 04_panel_a_proximity_compute.R)
#
# OUTPUTS (in data/processed/figure_2_panels/):
#   For both strict and relaxed coherent rules:
#   - fig2_a_proximity_heatmap_{strict|relaxed}.{pdf|png}      composite
#   - fig2_a_proximity_heatmap_{strict|relaxed}_only.{pdf|png} heatmap alone
#   - fig2_a_proximity_lollipop_{strict|relaxed}_only.{pdf|png} lollipop alone
#   And one shared:
#   - fig2_a_proximity_differential_table.csv (sig categories from both rules)
#
# COHERENT RULES:
#   strict : both samples of a condition agree on direction at FDR < 0.05.
#   relaxed: at least one sample of a condition is significant in the direction.
#
# DESIGN CHOICES:
#   - Lower triangle only (matrix is symmetric, half is enough)
#   - No lineage separator gridlines (ordering preserved but not visually grouped)
#   - Cell-type ordering follows cluster_name_order from 02_cluster_config.R
#
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

have_patchwork <- requireNamespace("patchwork", quietly = TRUE)
have_cowplot   <- requireNamespace("cowplot",   quietly = TRUE)
if (have_patchwork) {
  suppressPackageStartupMessages(library(patchwork))
} else if (have_cowplot) {
  suppressPackageStartupMessages(library(cowplot))
} else {
  stop("Neither patchwork nor cowplot is installed.")
}

# --- Paths ---
source(file.path(Sys.getenv("IRHEP_BUNDLE_ROOT", here::here()), "scripts", "config.R"))
source(file.path(SCRIPTS_DIR, "02_cluster_config.R"))

giotto_dir <- file.path(PROCESSED_DIR, "spatial_graph_analysis_giotto")
out_dir    <- file.path(PROCESSED_DIR, "figure_2_panels")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

SIG_THRESHOLD_STRICT  <- 0.10  # for the "coherent" rule (see below)
SIG_THRESHOLD_RELAXED <- 0.05  # for the "any sample sig" rule
N_TOP_HIGHLIGHT <- 10  # number of strongest pairs to bold-outline on the heatmap
                       # (filtered to coherent-supported pairs only — see make_heatmap)
N_TOP_LOLLIPOP  <- 15  # standalone lollipop file


# ============================================================================
# 1. Load + validate Giotto per-sample data
# ============================================================================
cat("Loading Giotto per-sample output...\n")
ps <- fread(file.path(giotto_dir, "giotto_proximity_per_sample.csv"))

stopifnot(
  "from_label"   %in% names(ps),
  "to_label"     %in% names(ps),
  "enrichm"      %in% names(ps),
  "p.adj_higher" %in% names(ps),
  "p.adj_lower"  %in% names(ps),
  "condition"    %in% names(ps),
  "sample_id"    %in% names(ps)
)
cat(sprintf("  %d rows; samples: %s\n", nrow(ps),
            paste(sort(unique(ps$sample_id)), collapse = ", ")))

# ============================================================================
# 2. Per-sample significance flags (TWO thresholds)
# ============================================================================
# Strict-rule per-sample sig: uses SIG_THRESHOLD_STRICT (0.10)
ps[, sig_attract_strict := !is.na(p.adj_higher) & p.adj_higher < SIG_THRESHOLD_STRICT]
ps[, sig_avoid_strict   := !is.na(p.adj_lower)  & p.adj_lower  < SIG_THRESHOLD_STRICT]

# Relaxed-rule per-sample sig: uses SIG_THRESHOLD_RELAXED (0.05)
ps[, sig_attract_relaxed := !is.na(p.adj_higher) & p.adj_higher < SIG_THRESHOLD_RELAXED]
ps[, sig_avoid_relaxed   := !is.na(p.adj_lower)  & p.adj_lower  < SIG_THRESHOLD_RELAXED]
ps[, sig_any_relaxed     := sig_attract_relaxed | sig_avoid_relaxed]

# ============================================================================
# 3. Compute BOTH consensus rules
# ============================================================================
sig_tab <- ps[, .(
  n_samples              = .N,
  n_sig_attract_strict   = sum(sig_attract_strict,  na.rm = TRUE),
  n_sig_avoid_strict     = sum(sig_avoid_strict,    na.rm = TRUE),
  n_sig_any_relaxed      = sum(sig_any_relaxed,     na.rm = TRUE)
), by = .(condition, from_label, to_label)]

sig_tab[, coherent_attract := (n_sig_attract_strict >= 1) & (n_sig_avoid_strict == 0)]
sig_tab[, coherent_avoid   := (n_sig_avoid_strict   >= 1) & (n_sig_attract_strict == 0)]
sig_tab[, sig_strict       := coherent_attract | coherent_avoid]

# Relaxed (p<0.05): at least one sample has any direction sig
sig_tab[, sig_relaxed := n_sig_any_relaxed >= 1]

# Wide format: one column per condition for each rule
sig_strict_wide <- dcast(
  sig_tab[, .(condition, from_label, to_label, sig_strict)],
  from_label + to_label ~ condition,
  value.var = "sig_strict", fill = FALSE
)
setnames(sig_strict_wide, c("irHepatitis", "AIH"),
         c("sig_irh_strict", "sig_aih_strict"))

sig_relaxed_wide <- dcast(
  sig_tab[, .(condition, from_label, to_label, sig_relaxed)],
  from_label + to_label ~ condition,
  value.var = "sig_relaxed", fill = FALSE
)
setnames(sig_relaxed_wide, c("irHepatitis", "AIH"),
         c("sig_irh_relaxed", "sig_aih_relaxed"))

# ============================================================================
# 4. Per-condition mean enrichment + differential
# ============================================================================
per_cond <- ps[, .(mean_enrichm = mean(enrichm, na.rm = TRUE)),
               by = .(condition, from_label, to_label)]

per_cond_wide <- dcast(per_cond,
                       from_label + to_label ~ condition,
                       value.var = "mean_enrichm", fill = 0)
setnames(per_cond_wide, c("irHepatitis", "AIH"), c("irHep_enrich", "AIH_enrich"))
per_cond_wide[, delta_enrichm := irHep_enrich - AIH_enrich]

# Merge BOTH significance versions
plot_dt <- merge(per_cond_wide, sig_strict_wide,
                 by = c("from_label", "to_label"), all.x = TRUE)
plot_dt <- merge(plot_dt, sig_relaxed_wide,
                 by = c("from_label", "to_label"), all.x = TRUE)
for (c in c("sig_irh_strict", "sig_aih_strict",
            "sig_irh_relaxed", "sig_aih_relaxed")) {
  plot_dt[is.na(get(c)), (c) := FALSE]
}

# Categorize for both rules (kept for table output and lollipop labels)
plot_dt[, sig_category_strict := fcase(
  sig_irh_strict &  sig_aih_strict, "both",
  sig_irh_strict & !sig_aih_strict, "irHep_only",
 !sig_irh_strict &  sig_aih_strict, "AIH_only",
  default                          = "neither"
)]
plot_dt[, sig_category_relaxed := fcase(
  sig_irh_relaxed &  sig_aih_relaxed, "both",
  sig_irh_relaxed & !sig_aih_relaxed, "irHep_only",
 !sig_irh_relaxed &  sig_aih_relaxed, "AIH_only",
  default                            = "neither"
)]

# Robustness count: number of conditions (0, 1, or 2) where the strict/relaxed
# consensus is met. This is the marker-encoding variable for the heatmap.
# It deliberately ignores DIRECTION — direction is encoded by cell color
# (Δ enrichment, red = irHep-enriched, blue = AIH-enriched). A pair where
# AIH robustly avoids contact (sig in AIH, avoidance) and irHep doesn't
# avoid will have a positive Δ (red cell) AND a "1 condition robust" marker
# from the AIH side. This avoids the misleading "AIH_only" label that
# implied AIH-attraction.
plot_dt[, robust_count_strict  := as.integer(sig_irh_strict)  + as.integer(sig_aih_strict)]
plot_dt[, robust_count_relaxed := as.integer(sig_irh_relaxed) + as.integer(sig_aih_relaxed)]

cat("\nSignificance breakdown — STRICT consensus (sig_category):\n")
print(table(plot_dt$sig_category_strict))
cat("\nRobustness counts — STRICT (0 / 1 / 2 conditions robust):\n")
print(table(plot_dt$robust_count_strict))

cat("\nSignificance breakdown — RELAXED consensus (sig_category):\n")
print(table(plot_dt$sig_category_relaxed))
cat("\nRobustness counts — RELAXED (0 / 1 / 2 conditions robust):\n")
print(table(plot_dt$robust_count_relaxed))

# Save full table with both versions
fwrite(plot_dt[order(-delta_enrichm)],
       file.path(out_dir, "fig2_a_proximity_differential_table.csv"))
cat("\nWrote: fig2_a_proximity_differential_table.csv\n")

# ============================================================================
# 5. Build lower-triangle plot data
# ============================================================================
# Order cell types by cluster_name_order from cluster_config.
# For each unordered pair (A,B), assign:
#   col_label = the one with lower index in cluster_name_order
#   row_label = the one with higher index
# This places everything in the LOWER triangle (col_idx <= row_idx).
plot_dt[, idx_a := match(from_label, cluster_name_order)]
plot_dt[, idx_b := match(to_label,   cluster_name_order)]
plot_dt[, col_label := ifelse(idx_a <= idx_b, from_label, to_label)]
plot_dt[, row_label := ifelse(idx_a <= idx_b, to_label,   from_label)]
plot_dt[, col_label := factor(col_label, levels = cluster_name_order)]
plot_dt[, row_label := factor(row_label, levels = cluster_name_order)]

# Pre-compute |Δ| for downstream ranking inside make_heatmap.
# is_top is now computed PER VERSION (strict vs relaxed) inside make_heatmap,
# because the filter "robust_count >= 1" depends on the version.
plot_dt[, abs_delta := abs(delta_enrichm)]
setorder(plot_dt, from_label, to_label)

delta_max <- max(abs(plot_dt$delta_enrichm), na.rm = TRUE)
delta_lim <- c(-delta_max, delta_max)
cat(sprintf("\nColor scale limits: (%.2f, %.2f)\n", delta_lim[1], delta_lim[2]))

# ============================================================================
# 6. Plotting helpers
# ============================================================================
make_heatmap <- function(dt, robust_col, title_suffix) {
  dt2 <- copy(dt)
  setnames(dt2, robust_col, "robust_count")
  dt2[, robust_factor := factor(robust_count, levels = c(0, 1, 2))]

  # Top-N highlight: filter to coherent-supported pairs (robust_count >= 1),
  # then take top N by |Δ enrichment|. Direction is encoded by border color.
  dt2[, is_top := FALSE]
  dt2[, top_direction := NA_character_]
  if (any(dt2$robust_count >= 1)) {
    dt2[robust_count >= 1,
        rank_in_filtered := frank(-abs_delta, ties.method = "first")]
    dt2[!is.na(rank_in_filtered) & rank_in_filtered <= N_TOP_HIGHLIGHT,
        is_top := TRUE]
    dt2[is_top == TRUE,
        top_direction := ifelse(delta_enrichm > 0, "irHep", "AIH")]
    dt2[, rank_in_filtered := NULL]
  }
  cat(sprintf("  Highlighted %d coherent top pairs (of %d coherent total)\n",
              sum(dt2$is_top), sum(dt2$robust_count >= 1)))

  ggplot(dt2,
         aes(x = col_label, y = row_label, fill = delta_enrichm)) +
    geom_tile(color = "grey92", linewidth = 0.12) +
    # Bold-outline highlight on the top-N coherent strongest pairs.
    # Border color = direction of differential (red = irHep>AIH, blue = AIH>irHep).
    geom_tile(
      data = dt2[is_top == TRUE],
      aes(x = col_label, y = row_label, color = top_direction),
      fill = NA, linewidth = 0.95,
      inherit.aes = FALSE,
      show.legend = FALSE
    ) +
    scale_color_manual(
      values = c("irHep" = "#67001F", "AIH" = "#053061"),
      na.value = NA,
      guide = "none"
    ) +
    scale_fill_gradient2(
      low      = "#2166AC",
      mid      = "white",
      high     = "#B2182B",
      midpoint = 0,
      limits   = delta_lim,
      name     = "Δ log2\nenrichment\n(irHep − AIH)",
      breaks   = pretty(delta_lim, n = 5)
    ) +
    # Markers: open circle = 1 condition robust, filled circle = both robust.
    # Direction is encoded by COLOR (Δ), not by marker.
    geom_point(
      data   = dt2[robust_count >= 1],
      aes(shape = robust_factor),
      size   = 1.3,
      color  = "black",
      stroke = 0.4,
      fill   = "black"
    ) +
    scale_shape_manual(
      values = c("0" = NA, "1" = 1, "2" = 16),
      name   = "Per-condition\nsig (coherent\nat FDR<0.10)",
      labels = c("1 condition", "Both conditions"),
      breaks = c("1", "2"),
      drop   = FALSE,
      na.translate = FALSE
    ) +
    # Reverse y so the matrix reads top-down (Hep1 at top)
    scale_y_discrete(limits = rev(levels(dt2$row_label)),
                     drop   = FALSE) +
    scale_x_discrete(drop = FALSE) +
    coord_equal() +
    theme_minimal(base_size = 8) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
      axis.text.y = element_text(size = 6),
      panel.grid  = element_blank(),
      legend.position   = "right",
      legend.key.height = unit(0.7, "lines"),
      legend.text       = element_text(size = 6),
      legend.title      = element_text(size = 7),
      plot.title    = element_text(size = 9, face = "bold"),
      plot.subtitle = element_text(size = 6.5, color = "grey40", lineheight = 1.1)
    ) +
    labs(
      title    = if (nzchar(title_suffix))
                   sprintf("Differential cell–cell proximity (irHep − AIH) — %s",
                           title_suffix) else NULL,
      subtitle = if (nzchar(title_suffix)) paste0(
        "Giotto cellProximityEnrichment, log2((obs+1)/(exp+1)); 1000 sims/sample; n=2 per condition. ",
        "Color = Δ enrichment (red = irHep>AIH).\n",
        "Markers = per-condition coherent sig (≥1 sample FDR<0.10, no opposite-sig). ",
        sprintf("Bold colored borders = top %d coherent |Δ| (red = irHep-enriched, blue = AIH-enriched).",
                N_TOP_HIGHLIGHT)
      ) else NULL,
      x = NULL, y = NULL
    )
}

make_lollipop <- function(dt, robust_col, title_suffix) {
  dt2 <- copy(dt)
  setnames(dt2, robust_col, "robust_count")
  dt2[, robust_factor := factor(robust_count, levels = c(0, 1, 2))]

  top_dt <- dt2[robust_count >= 1][order(-abs(delta_enrichm))][
    seq_len(min(N_TOP_LOLLIPOP, .N))
  ]
  if (nrow(top_dt) == 0) {
    return(ggplot() + labs(title = "No robust pairs",
                           subtitle = title_suffix) + theme_void())
  }
  top_dt[, pair_label := paste(from_label, "  ⟷  ", to_label)]
  top_dt <- top_dt[order(delta_enrichm)]
  top_dt[, pair_label := factor(pair_label, levels = pair_label)]
  top_dt[, direction := ifelse(delta_enrichm > 0, "irHep", "AIH")]

  ggplot(top_dt, aes(x = delta_enrichm, y = pair_label)) +
    geom_segment(
      aes(x = 0, xend = delta_enrichm,
          y = pair_label, yend = pair_label,
          color = direction),
      linewidth = 0.7
    ) +
    geom_point(aes(color = direction, shape = robust_factor),
               size = 3, fill = "white", stroke = 1) +
    scale_color_manual(values = c("irHep" = "#B2182B", "AIH" = "#2166AC"),
                       guide  = "none") +
    scale_shape_manual(
      values = c("0" = NA, "1" = 1, "2" = 16),
      name   = "Per-condition sig",
      labels = c("1 condition", "Both"),
      breaks = c("1", "2"),
      drop   = FALSE,
      na.translate = FALSE
    ) +
    geom_vline(xintercept = 0, linewidth = 0.3, color = "grey40") +
    theme_minimal(base_size = 9) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.text.y        = element_text(size = 7),
      plot.title         = element_text(size = 10, face = "bold"),
      legend.position    = "bottom"
    ) +
    labs(
      title = sprintf("Top %d differential pairs (%s, robust in ≥1 cond)",
                      nrow(top_dt), title_suffix),
      x = "Δ log2 enrichment (irHep − AIH)",
      y = NULL
    )
}

save_panel <- function(p_hm, p_ll, suffix) {
  # Primary output: heatmap alone (with bold-outline highlights for top hits).
  # No composite with the lollipop — the heatmap conveys the same information
  # via the highlighted tiles, and a smaller single panel matches Figure 1 sizing.
  ggsave(file.path(out_dir, sprintf("fig2_a_proximity_heatmap_%s.pdf", suffix)),
         p_hm, width = 7, height = 6.5, device = cairo_pdf)
  ggsave(file.path(out_dir, sprintf("fig2_a_proximity_heatmap_%s.png", suffix)),
         p_hm, width = 7, height = 6.5, dpi = 300)

  # Lollipop kept as a separate optional file (supplement / backup view)
  ggsave(file.path(out_dir,
                   sprintf("fig2_a_proximity_lollipop_%s_only.pdf", suffix)),
         p_ll, width = 5, height = 6, device = cairo_pdf)
  ggsave(file.path(out_dir,
                   sprintf("fig2_a_proximity_lollipop_%s_only.png", suffix)),
         p_ll, width = 5, height = 6, dpi = 300)
}

# ============================================================================
# 7. Build + save BOTH versions
# ============================================================================
cat("\n=== COHERENT rule (≥1 sample p<0.10, no opposite sig) ===\n")
hm_strict <- make_heatmap(plot_dt, "robust_count_strict", "coherent rule")
ll_strict <- make_lollipop(plot_dt, "robust_count_strict", "coherent")
save_panel(hm_strict, ll_strict, "strict")
cat("  Saved coherent-rule PDFs/PNGs (file suffix: _strict)\n")

# Panel-only version (no title/subtitle) for main figure assembly.
# Use short cluster labels (Hep1, CD8T2, etc.) that match Panel B short_names
# so the heatmap reads cleanly at composite display size.
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
plot_dt_short <- copy(plot_dt)
plot_dt_short[, col_label := factor(SHORT_LABELS[as.character(col_label)],
                                    levels = unname(SHORT_LABELS[cluster_name_order]))]
plot_dt_short[, row_label := factor(SHORT_LABELS[as.character(row_label)],
                                    levels = unname(SHORT_LABELS[cluster_name_order]))]

hm_panel <- make_heatmap(plot_dt_short, "robust_count_strict", "")
hm_panel <- hm_panel +
  theme(axis.text.x       = element_text(angle = 45, hjust = 1, size = 9),
        axis.text.y       = element_text(size = 9),
        legend.text       = element_text(size = 9),
        legend.title      = element_text(size = 9.5),
        legend.key.height = unit(0.9, "lines"),
        legend.key.width  = unit(0.45, "lines"),
        plot.margin       = margin(4, 4, 4, 4))
ggsave(file.path(out_dir, "fig2_a_proximity_heatmap_panel.pdf"),
       hm_panel, width = 6.2, height = 5.4, device = cairo_pdf)
ggsave(file.path(out_dir, "fig2_a_proximity_heatmap_panel.png"),
       hm_panel, width = 6.2, height = 5.4, dpi = 320)
cat("  Saved panel-only (short labels): fig2_a_proximity_heatmap_panel.{pdf,png}\n")

cat("\n=== RELAXED rule (any sample p<0.05) ===\n")
hm_relaxed <- make_heatmap(plot_dt, "robust_count_relaxed", "relaxed rule")
ll_relaxed <- make_lollipop(plot_dt, "robust_count_relaxed", "relaxed")
save_panel(hm_relaxed, ll_relaxed, "relaxed")
cat("  Saved relaxed-rule PDFs/PNGs (file suffix: _relaxed)\n")

# ============================================================================
# 8. Console report
# ============================================================================
print_top <- function(label, sig_col) {
  cat(sprintf("\n--- Top 15 sig pairs (%s) ---\n", label))
  show <- plot_dt[get(sig_col) != "neither"][order(-abs(delta_enrichm))][
    seq_len(min(N_TOP_LOLLIPOP, .N))
  ]
  for (i in seq_len(nrow(show))) {
    r <- show[i]
    cat(sprintf("  %-33s <-> %-33s  Δ=%+6.3f  irHep=%+6.2f  AIH=%+6.2f  [%s]\n",
                substr(r$from_label, 1, 33),
                substr(r$to_label, 1, 33),
                r$delta_enrichm, r$irHep_enrich, r$AIH_enrich,
                r[[sig_col]]))
  }
}
print_top("STRICT",  "sig_category_strict")
print_top("RELAXED", "sig_category_relaxed")

cat("\n================================================================\n")
cat("Panel A differential proximity heatmap COMPLETE\n")
cat("================================================================\n")
