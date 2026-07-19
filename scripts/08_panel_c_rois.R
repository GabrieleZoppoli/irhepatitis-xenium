#!/usr/bin/env Rscript
# ============================================================================
# Figure 2 — Panel C: 120 µm ROIs of bile-duct and parenchymal niches in
# each condition, rendered as paired H&E + segmentation overlay.
#
# For every ROI: cells passing Xenium QC (nucleus_area ≥ 8 µm², total counts
# ≥ 10) are recoloured by lineage; partner-type cells whose centroid lies
# within INTERACTION_DIST µm of any anatomy-partner cell are flagged "in
# contact" with a magenta border (#AD1457) and an in-pane jittered arrow
# (yellow on H&E, white on overlay; deterministic per ROI via set.seed).
# Transcripts of the displayed marker genes are plotted at single-molecule
# resolution within the ROI box without filtering by cell membership (gene
# specificity is then visually apparent rather than asserted by filter).
#
# Outputs: fig2_c_rois.{pdf,png} (rendered at 2× composite size so the
# downstream assembler shrinks back to ~50%, matching Panel A and B scaling).
# ============================================================================
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork); library(cowplot)
  library(arrow); library(png); library(grid); library(FNN); library(magick)
})

source(file.path(Sys.getenv("IRHEP_BUNDLE_ROOT", here::here()), "scripts", "config.R"))
source(file.path(SCRIPTS_DIR, "02_cluster_config.R"))

PDIR <- file.path(PROCESSED_DIR, "figure_2_panels")
dir.create(PDIR, recursive = TRUE, showWarnings = FALSE)

# ---- Cluster constants for the four lineages displayed in Panel C ----------
HEP   <- c(1L, 6L, 10L, 18L); CHOL <- 8L
CD8T1 <- 2L; CD8T2 <- 17L; CD8T3 <- 19L
CD4   <- 7L

# ---- Display colours (matching the manuscript palette) ---------------------
COL_CHOL           <- "#FFB300"
COL_HEP            <- "#E64B35"
COL_CD8T1          <- "#4DBBD5"
COL_CD8T2          <- "#1976D2"
COL_CD8T3          <- "#0D47A1"
COL_CD4            <- "#7BC8A4"
COL_OTHER          <- "grey93"
COL_CONTACT_BORDER <- "#AD1457"
LW_CONTACT_BORDER  <- 0.7
COL_NEUTRAL_SWATCH <- "grey85"

# ---- Analysis constants ----------------------------------------------------
INTERACTION_DIST <- 10    # µm: centroid-to-boundary threshold for "in contact"
MIN_QV           <- 20    # Xenium transcript-quality cutoff
ROI_SIZE_UM      <- 120
SCALE_BAR_UM     <- 25

# ---- Displayed transcripts (cell-specific markers) -------------------------
# SDC1 was evaluated as a second hepatocyte marker and dropped because its
# in-tissue distribution is broader than ARG1's; ARG1 alone gives a cleaner
# hepatocyte signal in the transcripts-everywhere view.
GENE_COL <- c(
  "EPCAM" = "#00695C", "CXCL6" = "#26A69A",
  "CD8A"  = "#E65100", "KLRK1" = "#FB8C00",
  "ARG1"  = "#4527A0",
  "IL7R"  = "#880E4F", "CD2"   = "#EC407A"
)
GENES_BY_TYPE <- list(
  chol = c("EPCAM",  "CXCL6"),
  cd8  = c("CD8A",   "KLRK1"),
  hep  = c("ARG1"),
  cd4  = c("IL7R",   "CD2")
)
GENES_CHOL_ROW <- c(GENES_BY_TYPE$chol, GENES_BY_TYPE$cd8)
GENES_HEP_ROW  <- c(GENES_BY_TYPE$hep,  GENES_BY_TYPE$cd4)

# ---- Selected representative ROIs (id, anatomy, condition, sample, centre) -
rois <- data.table(
  id        = c(306, 410, 601, 701),
  row       = c("Bile duct", "Bile duct", "Parenchyma", "Parenchyma"),
  col       = c("irHep", "AIH", "AIH", "irHep"),
  partner   = c("CD8", "CD8", "CD4", "CD4"),
  sample_id = c("18321_a", "54023_b", "54023_b", "18321_a"),
  cx        = c(1572, 3195, 3207, 8228),
  cy        = c(-2146, -9694, -9960, -1175),
  he_dir    = c("supp_panelC_v6_candidates_norm",
                "supp_panelC_v6_candidates_norm",
                "supp_panelC_v8_hep120_norm",
                "supp_panelC_v8_hep120_norm")
)
rois[, roi_x := cx - ROI_SIZE_UM/2]
rois[, roi_y := cy - ROI_SIZE_UM/2]

ok_cell <- function(leiden, ca, na, tc) {
  if (is.na(ca)) return(FALSE)
  if (is.na(na) || na < 8) return(FALSE)
  return(!is.na(tc) && tc >= 10)
}

# Per-sample x/y shifts used in the original multi-sample composite plotting.
samples_info <- as.data.table(SAMPLE_SHIFTS)[, .(sample_id, x_shift, y_shift)]

# ---- Load Giotto session and join cell areas / QC --------------------------
cat("Loading session...\n")
load(SESSION_SNAPSHOT)
suppressPackageStartupMessages({ library(Giotto); library(GiottoClass) })
m <- pDataDT(gobject); s <- getSpatialLocations(gobject, output = "data.table")
coords <- merge(s[, .(cell_ID, sdimx, sdimy)],
                m[!leiden %in% excluded_clusters,
                  .(cell_ID, sample_id, leiden)],
                by = "cell_ID")
coords[, leiden := as.integer(leiden)]
coords[, cell_id_raw := sub("^[0-9]+_[a-d]_", "", cell_ID)]
rm(gobject); gc(verbose = FALSE)

areas <- rbindlist(lapply(unique(rois$sample_id), function(sid) {
  a <- as.data.table(read_parquet(file.path(sample_raw_dir(sid), "cells.parquet")))
  a[, sample_id := sid]
  a[, .(cell_id_raw = cell_id, sample_id, cell_area, nucleus_area, total_counts)]
}))
coords <- merge(coords, areas, by = c("cell_id_raw", "sample_id"), all.x = TRUE)
coords[, ok := mapply(ok_cell, leiden, cell_area, nucleus_area, total_counts)]

# ---- Cache per-sample boundaries and transcripts ---------------------------
bd_cache <- list(); nc_cache <- list(); tx_cache <- list()
for (sid in unique(rois$sample_id)) {
  xo <- samples_info[sample_id == sid]$x_shift
  yo <- samples_info[sample_id == sid]$y_shift
  cb <- as.data.table(read_parquet(file.path(sample_raw_dir(sid), "cell_boundaries.parquet")))
  cb[, x_plot := vertex_x + xo]; cb[, y_plot := -vertex_y + yo]
  bd_cache[[sid]] <- cb
  nb <- as.data.table(read_parquet(file.path(sample_raw_dir(sid), "nucleus_boundaries.parquet")))
  nb[, x_plot := vertex_x + xo]; nb[, y_plot := -vertex_y + yo]
  nc_cache[[sid]] <- nb

  cat(sprintf("Loading transcripts for %s...\n", sid))
  tx <- as.data.table(read_parquet(
    file.path(sample_raw_dir(sid), "transcripts.parquet"),
    col_select = c("cell_id", "feature_name", "x_location", "y_location", "qv", "is_gene")))
  tx <- tx[is_gene == TRUE & qv >= MIN_QV &
           feature_name %in% names(GENE_COL)]
  tx[, x_plot := x_location + xo]
  tx[, y_plot := -y_location + yo]
  tx[, sample_id := sid]
  tx_cache[[sid]] <- tx
  cat(sprintf("  %s: %d transcripts of interest\n", sid, nrow(tx)))
}

# ---- Build a single ROI pair (H&E + overlay) -------------------------------
build_roi <- function(i) {
  r <- rois[i]
  rx <- r$roi_x; ry <- r$roi_y; rs <- ROI_SIZE_UM
  sid <- r$sample_id
  xlims <- c(rx, rx + rs); ylims <- c(ry, ry + rs); margin <- 5

  he_file <- file.path(PDIR, r$he_dir, sprintf("cluster_%d_he_norm.png", r$id))
  he_img <- if (file.exists(he_file)) readPNG(he_file) else NULL

  cb <- bd_cache[[sid]]; nb <- nc_cache[[sid]]; tx <- tx_cache[[sid]]
  cb_roi <- cb[x_plot >= (xlims[1] - margin) & x_plot <= (xlims[2] + margin) &
               y_plot >= (ylims[1] - margin) & y_plot <= (ylims[2] + margin)]
  nb_roi <- nb[x_plot >= (xlims[1] - margin) & x_plot <= (xlims[2] + margin) &
               y_plot >= (ylims[1] - margin) & y_plot <= (ylims[2] + margin)]
  cb_cent <- cb_roi[, .(cx = mean(x_plot), cy = mean(y_plot)), by = cell_id]
  roi_cells <- coords[sample_id == sid &
                      sdimx >= (xlims[1] - margin) & sdimx <= (xlims[2] + margin) &
                      sdimy >= (ylims[1] - margin) & sdimy <= (ylims[2] + margin)]
  cb_cent[, leiden := {
    d <- (roi_cells$sdimx - cx)^2 + (roi_cells$sdimy - cy)^2
    if (length(d) == 0) NA_integer_ else roi_cells$leiden[which.min(d)]
  }, by = cell_id]
  cb_cent[, ok := {
    d <- (roi_cells$sdimx - cx)^2 + (roi_cells$sdimy - cy)^2
    if (length(d) == 0) FALSE else roi_cells$ok[which.min(d)]
  }, by = cell_id]
  cb_cent_ok <- cb_cent[ok == TRUE]
  cb_j <- merge(cb_roi, cb_cent_ok[, .(cell_id, leiden)], by = "cell_id")

  if (r$partner == "CD8") {
    anatomy_clusters <- CHOL
    partner_clusters <- c(CD8T2, CD8T3)
    genes_show <- GENES_CHOL_ROW
  } else {
    anatomy_clusters <- HEP
    partner_clusters <- CD4
    genes_show <- GENES_HEP_ROW
  }

  # In-contact = partner cell whose centroid sits within INTERACTION_DIST µm
  # of any anatomy-partner cell-boundary vertex.
  partner_cents <- cb_cent_ok[leiden %in% partner_clusters]
  anatomy_cell_ids <- cb_cent_ok[leiden %in% anatomy_clusters]$cell_id
  interact_ids <- character(0)
  if (nrow(partner_cents) > 0 && length(anatomy_cell_ids) > 0) {
    av <- cb_roi[cell_id %in% anatomy_cell_ids, .(x_plot, y_plot)]
    if (nrow(av) > 0) {
      nn <- FNN::get.knnx(as.matrix(av), as.matrix(partner_cents[, .(cx, cy)]), k = 1)
      partner_cents[, dist := nn$nn.dist[, 1]]
      interact_ids <- partner_cents[dist <= INTERACTION_DIST]$cell_id
    }
  }
  cb_j[, cluster_role := fcase(
    leiden == CHOL,    "chol",
    leiden %in% HEP,   "hep",
    leiden == CD8T1,   "cd8t1",
    leiden == CD8T2,   "cd8t2",
    leiden == CD8T3,   "cd8t3",
    leiden == CD4,     "cd4",
    default = "other")]
  cb_j[, in_contact := cell_id %in% interact_ids &
                       leiden %in% partner_clusters]

  he_int <- cb_j[cluster_role %in% c("chol","hep","cd8t1","cd8t2","cd8t3","cd4")]
  draw_he_layer <- function(p, role_name, fill_color, lw = 0.5) {
    sub <- he_int[cluster_role == role_name]
    if (nrow(sub) == 0) return(p)
    p + geom_polygon(data = sub, aes(x = x_plot, y = y_plot, group = cell_id),
                     fill = NA, color = fill_color, linetype = "dashed",
                     linewidth = lw)
  }
  p_he <- ggplot() +
    (if (!is.null(he_img))
       annotation_raster(he_img, xmin = xlims[1], xmax = xlims[2],
                         ymin = ylims[1], ymax = ylims[2])
     else annotate("rect", xmin = xlims[1], xmax = xlims[2],
                   ymin = ylims[1], ymax = ylims[2], fill = "grey90"))
  p_he <- draw_he_layer(p_he, "cd8t1", COL_CD8T1, 0.5)
  p_he <- draw_he_layer(p_he, "cd8t2", COL_CD8T2, 0.55)
  p_he <- draw_he_layer(p_he, "cd8t3", COL_CD8T3, 0.55)
  p_he <- draw_he_layer(p_he, "cd4",   COL_CD4,   0.55)
  p_he <- draw_he_layer(p_he, "chol",  COL_CHOL,  0.65)
  p_he <- draw_he_layer(p_he, "hep",   COL_HEP,   0.6)
  ic <- cb_j[in_contact == TRUE]
  if (nrow(ic) > 0)
    p_he <- p_he + geom_polygon(data = ic, aes(x = x_plot, y = y_plot, group = cell_id),
                                fill = NA, color = COL_CONTACT_BORDER,
                                linewidth = 0.7)

  # Jittered in-pane arrows pointing to each in-contact cell. Tip lands 4 µm
  # outside the centroid (just past typical 5 µm cell radius) so the cell
  # itself is unobscured; tail at 10 µm gives a 6 µm visible shaft. Jitter
  # angle is seeded per-ROI so the same arrow pattern reproduces across runs.
  ic_centers <- if (nrow(ic) > 0)
    ic[, .(cx = mean(x_plot), cy = mean(y_plot)), by = cell_id] else
    data.table(cell_id = character(), cx = numeric(), cy = numeric())
  if (nrow(ic_centers) > 0) {
    set.seed(42L + as.integer(r$id))
    # Arrow orientations are NOT random: each arrow points away from the
    # inverse-distance-weighted vector to its k nearest in-contact neighbours
    # (k = min(3, n-1)). This naturally spreads arrows outward so they don't
    # overlap with each other when in-contact cells are clustered. A small
    # angular jitter (±0.10 rad) breaks ties for symmetric configurations
    # and preserves reproducibility via the per-ROI seed.
    coords_mat <- as.matrix(ic_centers[, .(cx, cy)])
    n_ic <- nrow(coords_mat)
    theta_vec <- numeric(n_ic)
    if (n_ic >= 2) {
      k_nn <- min(3L, n_ic - 1L)
      jitter_eps <- 0.10
      for (i in seq_len(n_ic)) {
        dx <- coords_mat[i, 1] - coords_mat[-i, 1, drop = FALSE]
        dy <- coords_mat[i, 2] - coords_mat[-i, 2, drop = FALSE]
        d  <- sqrt(dx^2 + dy^2)
        ord <- order(d)[seq_len(k_nn)]
        wts <- 1 / (d[ord] + 1e-6)
        vx <- sum(dx[ord] * wts) / sum(wts)
        vy <- sum(dy[ord] * wts) / sum(wts)
        theta_vec[i] <- atan2(vy, vx) + runif(1, -jitter_eps, jitter_eps)
      }
    } else {
      theta_vec <- runif(n_ic, 0, 2 * pi)
    }
    ic_centers[, theta := theta_vec]
    ic_centers[, tip_x  := cx + 4 * cos(theta)]
    ic_centers[, tip_y  := cy + 4 * sin(theta)]
    ic_centers[, tail_x := cx + 10 * cos(theta)]
    ic_centers[, tail_y := cy + 10 * sin(theta)]
    arrow_g <- arrow(length = unit(0.12, "inches"), type = "closed", angle = 26)
    # H&E view: yellow core + black halo (high contrast on pink/purple H&E).
    p_he <- p_he + geom_segment(data = ic_centers,
                                aes(x = tail_x, y = tail_y, xend = tip_x, yend = tip_y),
                                arrow = arrow_g, color = "black",
                                linewidth = 2.2, lineend = "round") +
                   geom_segment(data = ic_centers,
                                aes(x = tail_x, y = tail_y, xend = tip_x, yend = tip_y),
                                arrow = arrow_g, color = "#FFD600",
                                linewidth = 1.1, lineend = "round")
  }

  bar <- SCALE_BAR_UM
  add_scale_bar <- function(p) {
    p +
      annotate("segment", x = xlims[2] - 5 - bar, xend = xlims[2] - 5,
               y = ylims[1] + 5, yend = ylims[1] + 5,
               color = "black", linewidth = 1.6) +
      annotate("text", x = xlims[2] - 5 - bar/2, y = ylims[1] + 12,
               label = sprintf("%d µm", bar),
               color = "black", size = 5.0, fontface = "bold")
  }
  # Per-ROI sample-identifier text (top-left of each ROI).
  # White on H&E, black on segmentation overlay — palette matches Figure 1
  # sample-ID titles (plain bold, no halo).
  add_sample_label_he <- function(p) {
    p +
      annotate("text", x = xlims[1] + 5, y = ylims[2] - 5,
               label = sid, color = "white", size = 4.0,
               fontface = "bold", hjust = 0, vjust = 1)
  }
  add_sample_label_ov <- function(p) {
    p +
      annotate("text", x = xlims[1] + 5, y = ylims[2] - 5,
               label = sid, color = "black", size = 4.0,
               fontface = "bold", hjust = 0, vjust = 1)
  }
  p_he <- p_he +
    coord_fixed(xlim = xlims, ylim = ylims, expand = FALSE) +
    theme_void(base_size = 8) +
    theme(plot.margin = margin(1, 1, 1, 1),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6))
  p_he <- add_scale_bar(p_he)
  p_he <- add_sample_label_he(p_he)

  p_ov <- ggplot()
  for (rv in c("other","cd8t1","cd8t2","cd8t3","cd4","chol","hep")) {
    sub <- cb_j[cluster_role == rv & in_contact == FALSE]
    if (nrow(sub) == 0) next
    fill_col <- if (rv=="other") COL_OTHER else
                if (rv=="cd8t1") COL_CD8T1 else
                if (rv=="cd8t2") COL_CD8T2 else
                if (rv=="cd8t3") COL_CD8T3 else
                if (rv=="cd4")   COL_CD4 else
                if (rv=="chol")  COL_CHOL else COL_HEP
    border_col <- if (rv=="other") "grey70" else "grey25"
    border_w   <- if (rv=="other") 0.1 else 0.2
    alpha_v    <- if (rv=="other") 1 else 0.78
    p_ov <- p_ov + geom_polygon(data = sub, aes(x = x_plot, y = y_plot, group = cell_id),
                                fill = fill_col, color = border_col,
                                linewidth = border_w, alpha = alpha_v)
  }
  for (rv in c("cd8t2","cd8t3","cd4")) {
    sub <- cb_j[cluster_role == rv & in_contact == TRUE]
    if (nrow(sub) == 0) next
    fill_col <- if (rv=="cd8t2") COL_CD8T2 else
                if (rv=="cd8t3") COL_CD8T3 else COL_CD4
    p_ov <- p_ov + geom_polygon(data = sub, aes(x = x_plot, y = y_plot, group = cell_id),
                                fill = fill_col, color = COL_CONTACT_BORDER,
                                linewidth = LW_CONTACT_BORDER, alpha = 0.85)
  }
  if (nrow(nb_roi) > 0)
    p_ov <- p_ov + geom_polygon(data = nb_roi, aes(x = x_plot, y = y_plot, group = cell_id),
                                fill = NA, color = "grey35", linewidth = 0.08)

  # Transcripts are drawn for every cell in the ROI box (no cell-membership
  # filter). If the displayed genes are cell-specific, the natural spatial
  # distribution will show the dots clustering inside the relevant lineage.
  # Layer order: transcript dots BEFORE arrows so arrows sit on top of dots.
  tx_roi <- tx[x_plot >= xlims[1] & x_plot <= xlims[2] &
               y_plot >= ylims[1] & y_plot <= ylims[2] &
               feature_name %in% genes_show]
  if (nrow(tx_roi) > 0) {
    p_ov <- p_ov + geom_point(data = tx_roi,
                              aes(x = x_plot, y = y_plot, color = feature_name),
                              size = 1.6, shape = 19, alpha = 1) +
      scale_color_manual(values = GENE_COL[genes_show], guide = "none")
  }

  # Overlay view arrows: white core + black halo — neutral against the
  # cell-colour palette (yellow would blend into Chol gold). Drawn AFTER
  # transcripts so they remain visible over dots.
  if (nrow(ic_centers) > 0) {
    p_ov <- p_ov + geom_segment(data = ic_centers,
                                aes(x = tail_x, y = tail_y, xend = tip_x, yend = tip_y),
                                arrow = arrow_g, color = "black",
                                linewidth = 2.4, lineend = "round") +
                   geom_segment(data = ic_centers,
                                aes(x = tail_x, y = tail_y, xend = tip_x, yend = tip_y),
                                arrow = arrow_g, color = "white",
                                linewidth = 1.2, lineend = "round")
  }
  p_ov <- p_ov +
    coord_fixed(xlim = xlims, ylim = ylims, expand = FALSE) +
    theme_void(base_size = 8) +
    theme(plot.margin = margin(1, 1, 1, 1),
          panel.background = element_rect(fill = "white", color = NA),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6))
  p_ov <- add_scale_bar(p_ov)
  p_ov <- add_sample_label_ov(p_ov)

  list(he = p_he, ov = p_ov)
}

cat("\n=== Building Panel C mini-panels ===\n")
panels <- lapply(seq_len(nrow(rois)), build_roi)
names(panels) <- paste0("roi_", rois$id)

make_pair <- function(p) plot_grid(p$he, p$ov, ncol = 2, rel_widths = c(1, 1))
p_306 <- make_pair(panels$roi_306)
p_410 <- make_pair(panels$roi_410)
p_601 <- make_pair(panels$roi_601)
p_701 <- make_pair(panels$roi_701)

# ---- Column / row labels ---------------------------------------------------
hdr_irhep <- ggdraw() + draw_label("SR-irHep", fontface = "bold", size = 14, color = "grey15")
hdr_aih   <- ggdraw() + draw_label("AIH",   fontface = "bold", size = 14, color = "grey15")
col_hdrs  <- plot_grid(hdr_irhep, hdr_aih, ncol = 2, rel_widths = c(1, 1))

row_lbl_bd <- ggdraw() + draw_label("Bile duct",  fontface = "bold",
                                    size = 13, angle = 90, color = "grey15")
row_lbl_pa <- ggdraw() + draw_label("Parenchyma", fontface = "bold",
                                    size = 13, angle = 90, color = "grey15")

LBL_W <- 0.045

row1_body <- plot_grid(p_306, p_410, ncol = 2, rel_widths = c(1, 1))
row1      <- plot_grid(row_lbl_bd, row1_body, ncol = 2, rel_widths = c(LBL_W, 1))
row2_body <- plot_grid(p_701, p_601, ncol = 2, rel_widths = c(1, 1))
row2      <- plot_grid(row_lbl_pa, row2_body, ncol = 2, rel_widths = c(LBL_W, 1))

col_hdrs_aligned <- plot_grid(ggdraw(), col_hdrs, ncol = 2, rel_widths = c(LBL_W, 1))

# ---- Below-row cell-type legend strips ------------------------------------
# Strip rel_height kept generous (0.22) so the dots (size 7 ≈ 0.275 in
# diameter) have full vertical clearance from the ROI row above them.
build_row_legend <- function(entries, in_contact_idx = integer(0),
                             total_slots = NULL) {
  n <- length(entries$label)
  if (is.null(total_slots)) total_slots <- n
  slot_w <- 1.85
  offset <- (total_slots - n) / 2
  dt <- data.table(label = entries$label,
                   color = entries$color,
                   x = (seq_len(n) - 1 + offset) * slot_w + 0.5,
                   y = 1.00)
  dt[, in_contact := seq_len(.N) %in% in_contact_idx]

  ggplot(dt) +
    geom_point(aes(x = x, y = y, fill = label), shape = 21,
               size = 7,
               stroke = ifelse(dt$in_contact, LW_CONTACT_BORDER * 2, 0.6),
               color = ifelse(dt$in_contact, COL_CONTACT_BORDER, "grey20")) +
    geom_text(aes(x = x + 0.20, y = y, label = label),
              hjust = 0, size = 4.0, color = "grey15") +
    scale_fill_manual(values = setNames(dt$color, dt$label), guide = "none") +
    scale_x_continuous(limits = c(0, total_slots * slot_w + 0.5),
                       expand = c(0, 0)) +
    scale_y_continuous(limits = c(0.5, 1.5), expand = c(0, 0)) +
    theme_void() +
    theme(plot.margin = margin(6, 4, 6, 4))
}

TOTAL_SLOTS <- 5

chol_legend <- build_row_legend(
  entries = list(
    label = c("Chol", "CD8 T1", "CD8 T2", "CD8 T3", "in contact"),
    color = c(COL_CHOL, COL_CD8T1, COL_CD8T2, COL_CD8T3, COL_NEUTRAL_SWATCH)),
  in_contact_idx = 5L,
  total_slots = TOTAL_SLOTS)

hep_legend <- build_row_legend(
  entries = list(
    label = c("Hep", "CD4 T", "in contact"),
    color = c(COL_HEP, COL_CD4, COL_NEUTRAL_SWATCH)),
  in_contact_idx = 3L,
  total_slots = TOTAL_SLOTS)

chol_legend_aligned <- plot_grid(ggdraw(), chol_legend, ncol = 2, rel_widths = c(LBL_W, 1))
hep_legend_aligned  <- plot_grid(ggdraw(), hep_legend,  ncol = 2, rel_widths = c(LBL_W, 1))

# ---- Right-side transcript legend ------------------------------------------
# y coordinates align with the centres of row1 and row2 inside panel_block
# (row1 centre = 3.655, row2 centre = -0.193 in tx_legend coords).
build_tx_legend <- function() {
  bd_dt <- data.table(
    label = sprintf("%s · %s", names(GENE_COL[GENES_CHOL_ROW]),
                    c("Chol", "Chol", "CD8 T", "CD8 T")),
    color = GENE_COL[GENES_CHOL_ROW],
    y = c(4.755, 4.155, 3.155, 2.555))
  pa_dt <- data.table(
    label = sprintf("%s · %s", names(GENE_COL[GENES_HEP_ROW]),
                    c("Hep", "CD4 T", "CD4 T")),
    color = GENE_COL[GENES_HEP_ROW],
    y = c(0.607, -0.393, -0.993))

  ggplot() +
    annotate("text", x = -0.05, y = 5.35,
             label = "Transcripts",
             hjust = 0, fontface = "bold", size = 4.5, color = "grey15") +
    geom_point(data = bd_dt, aes(x = 0.05, y = y), shape = 21, size = 6,
               stroke = 0.5, fill = bd_dt$color, color = "grey20") +
    geom_text(data = bd_dt, aes(x = 0.21, y = y, label = label),
              hjust = 0, size = 4.0, color = "grey15") +
    geom_point(data = pa_dt, aes(x = 0.05, y = y), shape = 21, size = 6,
               stroke = 0.5, fill = pa_dt$color, color = "grey20") +
    geom_text(data = pa_dt, aes(x = 0.21, y = y, label = label),
              hjust = 0, size = 4.0, color = "grey15") +
    scale_x_continuous(limits = c(-0.1, 1.25), expand = c(0, 0)) +
    scale_y_continuous(limits = c(-2.4, 5.4), expand = c(0, 0)) +
    theme_void() +
    theme(plot.margin = margin(0, 1, 0, 5))
}
p_tx_legend <- build_tx_legend()

panel_block <- plot_grid(col_hdrs_aligned,
                         row1, chol_legend_aligned,
                         row2, hep_legend_aligned,
                         ncol = 1,
                         rel_heights = c(0.04, 1.28, 0.22, 1.28, 0.22)) +
  theme(plot.background = element_rect(fill = "white", color = NA))

fig <- plot_grid(panel_block, p_tx_legend, ncol = 2,
                 rel_widths = c(1, 0.17)) +
  theme(plot.background = element_rect(fill = "white", color = NA))

# Render at 2× the composite slot so the assembler's downstream shrink to
# ~7.25 × 4.7 in matches Panel A's and Panel B's source-to-composite scaling.
OUT_PDF <- file.path(PDIR, "fig2_c_rois.pdf")
OUT_PNG <- file.path(PDIR, "fig2_c_rois.png")
ggsave(OUT_PDF, fig, width = 14.5, height = 9.4, device = cairo_pdf, bg = "white")
ggsave(OUT_PNG, fig, width = 14.5, height = 9.4, dpi = 320, bg = "white")
cat(sprintf("Wrote %s\nWrote %s\n", OUT_PDF, OUT_PNG))
