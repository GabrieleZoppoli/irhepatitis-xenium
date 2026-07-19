#!/usr/bin/env Rscript
# ============================================================================
# 03_figure_1_main.R — Figure 1 build: cellular landscape and cluster
# definition. Composes UMAP, per-condition composition, marker dot-plot,
# H&E lens ROIs and spatial overviews into a multi-panel PDF. Reads the
# Giotto session snapshot and cluster config; writes panels to figures/main/.
# ============================================================================

suppressPackageStartupMessages({
  library(Giotto); library(GiottoClass); library(data.table); library(ggplot2)
  library(patchwork); library(cowplot); library(scales); library(ggrepel)
  library(pheatmap); library(grid); library(arrow); library(terra)
  library(circlize); library(jsonlite); library(png)
})


source(file.path(Sys.getenv("IRHEP_BUNDLE_ROOT", here::here()), "scripts", "config.R"))
# --- Setup ---
working_dir  <- RAW_DIR
table_dir    <- file.path(PROCESSED_DIR, "tables")
main_fig_dir <- file.path(PROCESSED_DIR, "main_figures")
panel_dir    <- file.path(main_fig_dir, "panels")
rds_dir      <- file.path(main_fig_dir, "rds")
for (d in c(main_fig_dir, panel_dir, rds_dir))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

samples_info <- data.table(
  sample_id = c("18321_a", "24774_d", "54023_b", "56784_c"),
  condition = c("irHepatitis", "irHepatitis", "AIH", "AIH"),
  x_shift   = c(0, 12000, 0, 12000),
  y_shift   = c(0, 0, -8000, -8000)
)

roi_coords <- data.table(
  sample_id = c("18321_a", "24774_d", "54023_b", "56784_c"),
  roi_x = c(1649, 15300, 4700, 13450),
  roi_y = c(-2124, -1200, -9800, -9930)
)

# --- Cluster config ---
source(file.path(SCRIPTS_DIR, "02_cluster_config.R"))

# Re-assign table_dir; the cluster config above overwrites it.
table_dir    <- file.path(PROCESSED_DIR, "tables")

# Active (non-excluded) cluster labels and colors for figure legends
active_labels <- cell_type_labels[!names(cell_type_labels) %in% excluded_clusters]
active_colors <- setNames(
  cluster_colors_by_number[names(active_labels)],
  active_labels
)

sig_colors <- c("Up in irHepatitis" = "#E41A1C", "Up in AIH" = "#377EB8", "NS" = "grey70")

# ============================================================================
# LOAD DATA
# ============================================================================
cat("=== Loading data ===\n")

gobj <- readRDS(file.path(PROCESSED_DIR, "gobj_with_spatial_networks.rds"))

meta <- pDataDT(gobj)
leiden_col <- grep("leiden", names(meta), value = TRUE)[1]
if (leiden_col != "leiden") meta[, leiden := get(leiden_col)]
meta[, cluster_name := cell_type_labels[as.character(leiden)]]
meta <- meta[!leiden %in% as.integer(excluded_clusters)]
meta[, cluster_name := factor(cluster_name, levels = cluster_name_order)]
if (!"condition" %in% names(meta)) {
  meta <- merge(meta, samples_info[, c("sample_id", "condition"), with = FALSE],
                by = "sample_id")
}

coords_dt <- merge(
  getSpatialLocations(gobj, output = "data.table"),
  meta[, .(cell_ID, leiden, cluster_name, sample_id, condition)],
  by = "cell_ID"
)
coords_dt[, cluster_name := factor(cluster_name, levels = cluster_name_order)]

expr_raw <- getExpression(gobj, values = "raw", output = "matrix")
expr_norm <- tryCatch(
  getExpression(gobj, values = "normalized", output = "matrix"),
  error = function(e) NULL
)
if (is.null(expr_norm)) expr_norm <- expr_raw

# Load UMAP coordinates
umap_dt <- tryCatch({
  umap_locs <- getDimReduction(gobj, reduction = "cells", reduction_method = "umap",
                                name = "umap", output = "data.table")
  umap_locs <- as.data.table(umap_locs, keep.rownames = "cell_ID")
  umap_locs
}, error = function(e) {
  cat("  UMAP extraction via getDimReduction failed, trying alternative...\n")
  NULL
})

# Load analysis tables (use file= to handle spaces in paths)
de_thresh    <- fread(file = file.path(table_dir, "de_thresh_edgeR.csv"))
svg_thresh   <- fread(file = file.path(table_dir, "svg_thresh_per_sample.csv"))
evidence     <- fread(file = file.path(table_dir, "synthesis_gene_evidence_matrix.csv"))
cluster_markers <- fread(file = file.path(PROCESSED_DIR, "tables", "cluster_markers.csv"))
cluster_comp <- fread(file = file.path(PROCESSED_DIR, "tables", "cluster_composition.csv"))

# Proximity data (permutation-tested)
perm_overall <- fread(file = file.path(PROCESSED_DIR, "spatial_graph_analysis",
                                       "proximity_permutation_results.csv"))
perm_cond    <- fread(file = file.path(PROCESSED_DIR, "spatial_graph_analysis",
                                       "proximity_permutation_per_condition.csv"))
perm_diff    <- fread(file = file.path(PROCESSED_DIR, "spatial_graph_analysis",
                                       "proximity_permutation_differential.csv"))

# Condition comparison proximity (for immune-stroma heatmap)
prox_cond <- fread(file = file.path(PROCESSED_DIR, "spatial_graph_analysis",
                                    "proximity_condition_comparison.csv"))

# Per-cluster DE
per_cluster_de <- fread(file = file.path(PROCESSED_DIR, "pseudobulk_de",
                                         "edgeR_per_cluster_de.csv"))
per_cluster_de <- per_cluster_de[!cluster %in% as.integer(excluded_clusters)]

# Filter cluster 21 from composition and marker tables
cluster_comp <- cluster_comp[!leiden %in% as.integer(excluded_clusters)]
cluster_markers <- cluster_markers[!cluster %in% as.integer(excluded_clusters)]

cat("All data loaded.\n")


# ============================================================================
# FIGURE 1: CELLULAR LANDSCAPE & CLUSTER DEFINITION
# ============================================================================
# Layout:  A (dot plot + composition bar chart, aligned x-axes)
#          B (condition UMAPs side by side + cluster legend)
#          C (irHep H&E + lens) | D (AIH H&E + lens)
cat("\n=== FIGURE 1: Cellular Landscape & Cluster Definition ===\n")

# --- Legend labels: "1 Hepatocyte 1", "2 CD8+ T cell 1", etc. ---
# display_number maps original Leiden numbers to sequential 1-18
# Cross-reference abbreviations shown in Figure 2 and the main text
cluster_abbrev <- c(
  "Hepatocyte 1" = "Hep1", "Hepatocyte 2" = "Hep2", "Hepatocyte 3" = "Hep3",
  "Hepatocyte 4" = "Hep4", "Ig+ hepatocyte" = "Ig+Hep",
  "CD8+ T cell 1" = "CD8T1", "CD8+ T cell 2" = "CD8T2", "CD8+ T cell 3" = "CD8T3",
  "CD4+ T cell" = "CD4T", "Cholangiocyte" = "Chol",
  "IgG+ plasmablast / plasma cell" = "IgG+PB", "IgG+ B cell" = "IgG+B",
  "IgM+ B cell" = "IgM+B", "Macrophage 1" = "Mac1", "Macrophage 2" = "Mac2",
  "Stellate cell" = "Stellate", "Sinusoidal endothelial cell" = "SinEC",
  "Vascular endothelial cell" = "VascEC"
)
abbr_suffix <- ifelse(is.na(cluster_abbrev[active_labels]), "",
                      paste0(" (", cluster_abbrev[active_labels], ")"))
legend_labels <- setNames(
  paste0(display_number[names(active_labels)], " ", active_labels, abbr_suffix),
  active_labels
)
legend_colors <- setNames(active_colors[names(legend_labels)], legend_labels)
seq_order <- active_labels[as.character(sort(as.integer(names(active_labels))))]
legend_order <- legend_labels[seq_order]

# Cluster order for dot plot & composition: numerical (retained clusters only)
num_cluster_order <- active_labels[as.character(sort(as.integer(names(active_labels))))]
# Mapping cluster name -> display number (for dotplot x-axis labels)
name_to_number <- setNames(display_number[names(active_labels)], active_labels)

# --- UMAP data ---
umap_plot <- merge(umap_dt, meta[, .(cell_ID, leiden, cluster_name, sample_id, condition)],
                   by = "cell_ID")
umap_plot[, legend_label := legend_labels[as.character(cluster_name)]]
umap_plot[, legend_label := factor(legend_label, levels = legend_order)]

umap_centroids_dt <- umap_plot[, .(x = median(Dim.1), y = median(Dim.2)), by = leiden]
umap_centroids_dt[, display_num := display_number[as.character(leiden)]]

# --- Condition UMAPs (Panel B bottom) — clean, no halos ---
make_condition_umap <- function(cond) {
  bg <- umap_plot[, .(Dim.1, Dim.2)]
  fg <- umap_plot[condition == cond]

  # Display label: harmonised short form used in the manuscript figures
  display_cond <- c("irHepatitis" = "SR-irHep", "AIH" = "AIH")[cond]
  if (is.na(display_cond)) display_cond <- cond

  ggplot() +
    geom_point(data = bg, aes(x = Dim.1, y = Dim.2),
               color = "grey88", size = 0.01, alpha = 0.3) +
    geom_point(data = fg, aes(x = Dim.1, y = Dim.2, color = legend_label),
               size = 0.04, alpha = 0.5) +
    geom_text(data = umap_centroids_dt, aes(x = x, y = y, label = display_num),
              size = 1.8, fontface = "bold", color = "black") +
    scale_color_manual(values = legend_colors, name = "Cluster", drop = FALSE) +
    coord_fixed() +
    labs(x = "UMAP 1", y = "UMAP 2", title = display_cond) +
    theme_minimal(base_size = 7) +
    theme(legend.position = "none",
          plot.title = element_text(face = "bold", hjust = 0.5, size = 9),
          axis.text = element_text(size = 5),
          axis.title = element_text(size = 6))
}

p_1b_irHep <- make_condition_umap("irHepatitis")
p_1b_AIH   <- make_condition_umap("AIH")

# --- Shared legend (for Panels B, C, D — right column) ---
legend_dummy_dt <- data.table(
  x = seq_along(legend_order), y = 1,
  legend_label = factor(legend_order, levels = legend_order)
)
p_legend <- ggplot(legend_dummy_dt, aes(x = x, y = y, color = legend_label)) +
  geom_point(size = 3) +
  scale_color_manual(values = legend_colors, name = "Cluster", drop = FALSE) +
  guides(color = guide_legend(ncol = 1, byrow = TRUE,
                              override.aes = list(size = 3, alpha = 1))) +
  theme_void(base_size = 7) +
  theme(legend.position = "right",
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 8, face = "bold"),
        legend.key.height = unit(0.40, "cm"),
        legend.key.width = unit(0.35, "cm"),
        legend.spacing.y = unit(0.04, "cm"),
        legend.margin = margin(0, 0, 0, 4))
shared_legend <- cowplot::get_legend(p_legend)

# --- Panel A: Dot plot (marker gene characterization) ---
# Greedy 1-per-cluster non-redundant DE marker selection.
# Cluster display numbers above dot plot.

cat("  Panel A: Selecting 1 non-redundant marker per cluster...\n")
claimed_genes <- character(0)
marker_selection <- list()
for (cl_num in sort(as.integer(names(active_labels)))) {
  cl_de <- cluster_markers[cluster == cl_num & FDR < 0.05 & logFC > 0 &
                              !feats %in% claimed_genes]
  cl_de <- cl_de[order(-logFC)]
  if (nrow(cl_de) > 0) {
    for (j in seq_len(nrow(cl_de))) {
      if (cl_de$feats[j] %in% rownames(expr_raw)) {
        claimed_genes <- c(claimed_genes, cl_de$feats[j])
        marker_selection[[as.character(cl_num)]] <- cl_de$feats[j]
        break
      }
    }
  }
}
dotplot_genes <- unlist(marker_selection[as.character(sort(as.integer(names(marker_selection))))])
n_dotplot_genes <- length(dotplot_genes)
cat(sprintf("  Panel A: %d non-redundant markers selected\n", n_dotplot_genes))

if (n_dotplot_genes > 0) {
  dot_list <- list()
  for (g in dotplot_genes) {
    for (cl in unique(meta$cluster_name)) {
      cl_cells <- meta[cluster_name == cl, cell_ID]
      g_expr <- as.numeric(expr_raw[g, cl_cells])
      dot_list[[paste(g, cl)]] <- data.table(
        gene = g, cluster_name = cl,
        pct_expressing = 100 * mean(g_expr > 0),
        mean_expr = mean(g_expr)
      )
    }
  }
  dot_dt <- rbindlist(dot_list)
  dot_dt[, cluster_name := factor(cluster_name, levels = num_cluster_order)]
  # Gene order: cluster 1's marker at bottom, cluster 20's at top
  dot_dt[, gene := factor(gene, levels = dotplot_genes)]

  p_1d_dot <- ggplot(dot_dt, aes(x = cluster_name, y = gene,
                                   size = pct_expressing, color = mean_expr)) +
    geom_point() +
    scale_size_continuous(range = c(0.3, 2.2), name = "% Expr.") +
    scale_color_viridis_c(option = "inferno", name = "Mean\nExpr.") +
    scale_x_discrete(drop = FALSE, position = "top",
                     labels = function(x) name_to_number[x]) +
    labs(x = NULL, y = NULL, tag = "A") +
    theme_minimal(base_size = 7) +
    theme(axis.text.x.top = element_text(size = 6),
          axis.text.y = element_text(size = 5.5, face = "italic"),
          legend.text = element_text(size = 5),
          legend.title = element_text(size = 5.5),
          legend.key.size = unit(0.25, "cm"),
          plot.tag = element_text(size = 12, face = "bold", hjust = 0, vjust = 1),
          plot.tag.position = c(0.005, 0.995),
          plot.margin = margin(t = 10, b = 5, l = 5, r = 5))
} else {
  p_1d_dot <- ggplot() + annotate("text", x=0.5, y=0.5, label="No markers") + theme_void()
}

# --- Panel B (top): Composition bar chart (per-sample-averaged, irHep first) ---
comp_dt <- copy(cluster_comp)
comp_dt[, cluster_name := cell_type_labels[as.character(leiden)]]
comp_dt <- comp_dt[!is.na(cluster_name)]
comp_dt[, sample_total := sum(N), by = .(sample_id)]
comp_dt[, sample_prop := N / sample_total]

# Average within condition (each sample weighted equally)
comp_mean <- comp_dt[, .(proportion = mean(sample_prop)), by = .(leiden, cluster_name, condition)]
comp_mean[, cluster_name := factor(cluster_name, levels = num_cluster_order)]
comp_mean[, condition := factor(condition, levels = c("irHepatitis", "AIH"))]

p_1d_comp <- ggplot(comp_mean, aes(x = cluster_name, y = proportion, fill = condition)) +
  geom_col(position = "dodge", width = 0.7) +
  scale_fill_manual(values = condition_colors,
                    labels = c("irHepatitis" = "SR-irHep", "AIH" = "AIH")) +
  scale_x_discrete(drop = FALSE, labels = function(x) name_to_number[x]) +
  labs(x = "Cluster", y = "Proportion", fill = "Condition") +
  theme_minimal(base_size = 7) +
  theme(axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        legend.position = "bottom",
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 6.5),
        legend.key.size = unit(0.25, "cm"),
        plot.margin = margin(t = 5, b = 5, l = 5, r = 5))

# --- Save proportion tables ---
v3_table_dir <- file.path(PROCESSED_DIR, "tables")
dir.create(v3_table_dir, recursive = TRUE, showWarnings = FALSE)

# Table 1: Full per-sample transparency
prop_per_sample <- comp_dt[, .(leiden, cluster_name, sample_id, condition,
                                n_cells = N, sample_total, proportion = sample_prop)]
setorder(prop_per_sample, leiden, sample_id)
fwrite(prop_per_sample, file.path(v3_table_dir, "cluster_proportions_per_sample.csv"))
cat(sprintf("  Saved cluster_proportions_per_sample.csv: %d rows\n", nrow(prop_per_sample)))

# Table 2: Condition-level summary with individual sample values
prop_summary <- comp_dt[, .(leiden, cluster_name, sample_id, condition, pct = sample_prop * 100)]
prop_wide_samples <- dcast(prop_summary, leiden + cluster_name ~ sample_id,
                           value.var = "pct", fill = 0)
prop_wide_cond <- dcast(
  prop_summary[, .(mean_pct = mean(pct)), by = .(leiden, cluster_name, condition)],
  leiden + cluster_name ~ condition, value.var = "mean_pct", fill = 0
)
if ("irHepatitis" %in% names(prop_wide_cond))
  setnames(prop_wide_cond, "irHepatitis", "irHep_mean_pct")
if ("AIH" %in% names(prop_wide_cond))
  setnames(prop_wide_cond, "AIH", "AIH_mean_pct")
prop_final <- merge(prop_wide_cond, prop_wide_samples, by = c("leiden", "cluster_name"))
sample_cols <- setdiff(names(prop_final), c("leiden", "cluster_name", "irHep_mean_pct", "AIH_mean_pct"))
for (sc in sample_cols) {
  cond_label <- samples_info[sample_id == sc, condition]
  if (length(cond_label) > 0) {
    new_name <- paste0(ifelse(cond_label == "irHepatitis", "irHep_", "AIH_"), sc, "_pct")
    setnames(prop_final, sc, new_name)
  }
}
setorder(prop_final, leiden)
fwrite(prop_final, file.path(v3_table_dir, "cluster_proportions_summary.csv"))
cat(sprintf("  Saved cluster_proportions_summary.csv: %d rows\n", nrow(prop_final)))

# --- Panels C & D: H&E overview + cell/nucleus contour lens below ---
# C: 24774_d (irHep), D: 54023_b (AIH)
# Lens uses cell + nuclear boundary polygons (same style as Fig 3 ROIs)
# with cluster number labels on centroids.
# ROI coordinates are LOWER-LEFT corner of the 600x600 µm square.
he_overview_dir <- file.path(main_fig_dir, "he_overviews")
roi_size <- 600  # micrometers

# ROI coordinates optimized by pathology-specific scan (scan_rois_pathology.R):
# irHep: T-cell enrichment, ductulopenia, interface hepatitis
# AIH: TLS-like B/plasma aggregates, depleted parenchyma, emperipolesis
roi_e <- c(15200, -1400)   # 24774_d (irHep) — user-selected ROI
roi_f <- c(794, -11063)    # 54023_b (AIH) — 9.4% hep (depleted), 20% B/plasma, TLS-like

cat(sprintf("  Panel E lens ROI lower-left: (%.0f, %.0f)\n", roi_e[1], roi_e[2]))
cat(sprintf("  Panel F lens ROI lower-left: (%.0f, %.0f)\n", roi_f[1], roi_f[2]))

he_extents <- list(
  "24774_d" = c(xmin = 12643.8947, xmax = 17261.7471, ymin = -2402.7649, ymax = -630.037),
  "54023_b" = c(xmin = 43.9531, xmax = 7020.4067, ymin = -11233.6331, ymax = -8430.6443)
)

# Generate the H&E overview PNGs via the Python producer if missing.
# 03 only needs overviews for the two ROIs displayed in Fig 1 Panels C/D
# (24774_d and 54023_b). The bundle ships them tracked under
# data/processed/main_figures/he_overviews/, so this step is a no-op on a
# clean clone; it is the from-raw rebuild path.
he_pngs_missing <- !all(file.exists(file.path(he_overview_dir,
                                              paste0("he_overview_",
                                                     c("24774_d", "54023_b"),
                                                     ".png"))))
if (he_pngs_missing) {
  cat("  H&E overview PNGs missing; running 03b_sample_he_overviews.py...\n")
  dir.create(he_overview_dir, recursive = TRUE, showWarnings = FALSE)
  he_cfg <- list(
    working_dir = RAW_DIR,
    output_dir  = he_overview_dir,
    samples = list(
      list(sample_id = "24774_d",
           extent = as.list(he_extents$`24774_d`)),
      list(sample_id = "54023_b",
           extent = as.list(he_extents$`54023_b`))
    )
  )
  cfg_path <- tempfile(fileext = ".json")
  jsonlite::write_json(he_cfg, cfg_path, auto_unbox = TRUE, pretty = TRUE)
  py_script <- file.path(SCRIPTS_DIR, "03b_sample_he_overviews.py")
  system2("python3", args = c(shQuote(py_script), shQuote(cfg_path)))
  unlink(cfg_path)
}

# Load parquet boundary data for lens panels
cat("  Loading parquet boundary data for H&E lens panels...\n")
fig1_boundary <- list()
for (sid in c("24774_d", "54023_b")) {
  dir_name <- paste0("0027420_", sid)
  sample_path <- file.path(working_dir, dir_name)
  x_off <- samples_info[sample_id == sid]$x_shift
  y_off <- samples_info[sample_id == sid]$y_shift

  cb_file <- file.path(sample_path, "cell_boundaries.parquet")
  if (file.exists(cb_file)) {
    cb <- as.data.table(read_parquet(cb_file))
    cb[, `:=`(x_plot = vertex_x + x_off, y_plot = -vertex_y + y_off)]
    fig1_boundary[[sid]] <- cb
  }

  nb_file <- file.path(sample_path, "nucleus_boundaries.parquet")
  if (file.exists(nb_file)) {
    nb <- as.data.table(read_parquet(nb_file))
    nb[, `:=`(x_plot = vertex_x + x_off, y_plot = -vertex_y + y_off)]
    attr(fig1_boundary[[sid]], "nucleus") <- nb
  }
}
cat("  Parquet boundary data loaded.\n")

# Build H&E overview + cell/nucleus contour lens below
make_he_with_contour_lens <- function(roi_corner, sample_id_val, sample_title,
                                       img_extent, boundary_data) {
  # ROI: lower-left corner to corner + roi_size
  x_lo <- roi_corner[1]; x_hi <- roi_corner[1] + roi_size
  y_lo <- roi_corner[2]; y_hi <- roi_corner[2] + roi_size

  # H&E overview
  he_png <- file.path(he_overview_dir, paste0("he_overview_", sample_id_val, ".png"))
  img <- readPNG(he_png)
  ext_xmin <- img_extent["xmin"]; ext_xmax <- img_extent["xmax"]
  ext_ymin <- img_extent["ymin"]; ext_ymax <- img_extent["ymax"]

  he_bar_len <- 1000; he_bar_pad <- 150
  he_bar_x0 <- ext_xmax - he_bar_len - he_bar_pad
  he_bar_y0 <- ext_ymin + he_bar_pad

  p_he <- ggplot() +
    annotation_raster(img, xmin = ext_xmin, xmax = ext_xmax,
                      ymin = ext_ymin, ymax = ext_ymax) +
    annotate("rect", xmin = x_lo, xmax = x_hi,
             ymin = y_lo, ymax = y_hi,
             fill = NA, color = "black", linewidth = 0.5) +
    annotate("segment", x = he_bar_x0, xend = he_bar_x0 + he_bar_len,
             y = he_bar_y0, yend = he_bar_y0,
             color = "black", linewidth = 0.8) +
    annotate("text", x = he_bar_x0 + he_bar_len / 2, y = he_bar_y0 + 120,
             label = "1 mm", color = "black", size = 2.5, fontface = "bold") +
    coord_fixed(xlim = c(ext_xmin, ext_xmax),
                ylim = c(ext_ymin, ext_ymax), expand = FALSE) +
    labs(title = sample_title) +
    theme_void() +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 7),
          plot.margin = margin(1, 1, 1, 1))

  # Contour lens: cell boundaries colored by cluster + nuclear outlines
  cb <- boundary_data
  nuc <- attr(boundary_data, "nucleus")
  margin_pq <- 10

  p_lens <- ggplot()
  roi_cluster_names_local <- character(0)

  if (!is.null(cb) && nrow(cb) > 0) {
    cb_roi <- cb[x_plot >= (x_lo - margin_pq) & x_plot <= (x_hi + margin_pq) &
                   y_plot >= (y_lo - margin_pq) & y_plot <= (y_hi + margin_pq)]

    if (nrow(cb_roi) > 0) {
      # Match cell boundaries to cluster assignments via nearest centroid
      cb_centroids <- cb_roi[, .(cx = mean(x_plot), cy = mean(y_plot)), by = cell_id]
      roi_cells <- coords_dt[sample_id == sample_id_val &
                               sdimx >= x_lo & sdimx <= x_hi &
                               sdimy >= y_lo & sdimy <= y_hi]

      if (nrow(roi_cells) > 0) {
        cb_centroids[, `:=`(cluster_name = NA_character_, match_dist = NA_real_)]
        cb_centroids[, c("cluster_name", "match_dist") := {
          dists <- (roi_cells$sdimx - cx)^2 + (roi_cells$sdimy - cy)^2
          idx <- which.min(dists)
          list(as.character(roi_cells$cluster_name[idx]), sqrt(dists[idx]))
        }, by = cell_id]

        # Only assign cluster if match is close (within ~cell diameter, 15 um)
        # Distant matches are unmatched cells (excluded clusters or edge artifacts)
        cb_centroids[match_dist > 15, cluster_name := NA_character_]

        cb_roi_matched <- merge(cb_roi, cb_centroids[, .(cell_id, cluster_name)], by = "cell_id")
        cb_roi_colored <- cb_roi_matched[!is.na(cluster_name)]
        cb_roi_grey    <- cb_roi_matched[is.na(cluster_name)]

        # Draw unmatched/excluded cells as thin grey outlines first
        if (nrow(cb_roi_grey) > 0) {
          p_lens <- p_lens +
            geom_polygon(data = cb_roi_grey,
                         aes(x = x_plot, y = y_plot, group = cell_id),
                         fill = NA, color = "grey75", linewidth = 0.08)
        }

        # Draw matched cells with cluster-colored fill
        p_lens <- p_lens +
          geom_polygon(data = cb_roi_colored,
                       aes(x = x_plot, y = y_plot, group = cell_id,
                           fill = cluster_name),
                       color = "grey40", linewidth = 0.1, alpha = 0.6) +
          scale_fill_manual(values = cluster_colors, guide = "none")

        roi_cluster_names_local <- unique(as.character(cb_roi_colored$cluster_name))
      } else {
        p_lens <- p_lens +
          geom_polygon(data = cb_roi,
                       aes(x = x_plot, y = y_plot, group = cell_id),
                       fill = NA, color = "grey70", linewidth = 0.08)
      }
    }
  }

  if (!is.null(nuc) && nrow(nuc) > 0) {
    nuc_roi <- nuc[x_plot >= (x_lo - margin_pq) & x_plot <= (x_hi + margin_pq) &
                     y_plot >= (y_lo - margin_pq) & y_plot <= (y_hi + margin_pq)]
    if (nrow(nuc_roi) > 0) {
      p_lens <- p_lens +
        geom_polygon(data = nuc_roi,
                     aes(x = x_plot, y = y_plot, group = cell_id),
                     fill = NA, color = "grey50", linewidth = 0.05)
    }
  }

  p_lens <- p_lens +
    coord_fixed(xlim = c(x_lo, x_hi), ylim = c(y_lo, y_hi), expand = FALSE) +
    theme_void() +
    theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
          plot.margin = margin(2, 2, 2, 2))

  # Stack: H&E on top, contour lens below
  # Give lens more height so the square fills the same width as the wide H&E
  panel <- cowplot::plot_grid(p_he, p_lens, ncol = 1, rel_heights = c(0.45, 1))
  list(panel = panel, roi_clusters = roi_cluster_names_local)
}

res_e <- make_he_with_contour_lens(roi_e, "24774_d", "24774_d (SR-irHep)",
                                    he_extents[["24774_d"]], fig1_boundary[["24774_d"]])
res_f <- make_he_with_contour_lens(roi_f, "54023_b", "54023_b (AIH)",
                                    he_extents[["54023_b"]], fig1_boundary[["54023_b"]])
p_1e <- res_e$panel
p_1f <- res_f$panel
roi_clusters_ef <- unique(c(res_e$roi_clusters, res_f$roi_clusters))

rm(fig1_boundary); gc(verbose = FALSE)

# --- Assemble Figure 1 ---
# Layout: Row 1 = A (dotplot + bar chart, aligned x-axes)
#         Rows 2-3 = B (UMAPs) + C/D (H&E + lens) | shared legend column

# Row 1 (Panel A): Dotplot stacked with composition bar chart, aligned x-axes
# Remove x-axis text from dotplot (bar chart below will have the numbers)
p_1a_dot_clean <- p_1d_dot +
  theme(axis.text.x.top = element_blank())

panel_a_patchwork <- patchwork::patchworkGrob(
  (p_1a_dot_clean / p_1d_comp) + plot_layout(heights = c(1.2, 1))
)

# Rows 2-3 left column: B (UMAPs) stacked with C/D (H&E + lens)
umap_pair <- cowplot::plot_grid(p_1b_irHep, p_1b_AIH, ncol = 2)

he_panels <- cowplot::plot_grid(p_1e, p_1f, ncol = 2,
                                 labels = c("C", "D"),
                                 label_size = 12, label_fontface = "bold",
                                 label_x = 0.005, label_y = 0.995,
                                 hjust = 0, vjust = 1)

bcd_left <- cowplot::plot_grid(
  umap_pair, NULL, he_panels,
  ncol = 1, rel_heights = c(1, 0.03, 1.6),
  labels = c("B", "", ""), label_size = 12, label_fontface = "bold",
  label_x = 0.005, label_y = 0.995, hjust = 0, vjust = 1
)

# Right column: single shared legend spanning B, C, D
legend_col <- cowplot::plot_grid(NULL, shared_legend, NULL,
                                  ncol = 1, rel_heights = c(0.05, 3, 0.05))

bcd_row <- cowplot::plot_grid(bcd_left, legend_col,
                                ncol = 2, rel_widths = c(3.5, 1))

# Final assembly: A (dotplot + bar) / B+C+D (UMAPs + H&E + shared legend)
fig1 <- cowplot::plot_grid(panel_a_patchwork, NULL, bcd_row,
                            ncol = 1, rel_heights = c(1.8, 0.21, 2.8))

ggsave(file.path(main_fig_dir, "figure_1.9_landscape.pdf"), fig1,
       width = 210, height = 297, units = "mm", limitsize = FALSE)
saveRDS(list(p_1a_dot = p_1d_dot, p_1a_comp = p_1d_comp,
             p_1b_irHep = p_1b_irHep, p_1b_AIH = p_1b_AIH,
             shared_legend = shared_legend, p_1c = p_1e, p_1d = p_1f,
             comp_dt = comp_dt, comp_mean = comp_mean),
        file.path(rds_dir, "fig1_panels.rds"))
cat("Saved: figure_1.9_landscape.pdf\n")

# 03 produces Figure 1. Figure 2 is built by 09_figure_2_assembled.R.
cat("\n03 complete.\n")
