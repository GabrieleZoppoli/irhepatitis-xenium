#!/usr/bin/env Rscript
# ============================================================================
# Per-cluster supplementary ROI figures.
# For each cluster, finds the densest 200×200 µm region across all samples
# and writes a 4-panel figure: H&E, DAPI + scale bar, cell boundaries coloured
# by cluster (target highlighted), top-marker transcript dots + nucleus
# outlines.
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(arrow)
  library(png)
  library(cowplot)
  library(patchwork)
  library(jsonlite)
  library(grid)
  library(Giotto)
  library(GiottoClass)
})


source(file.path(Sys.getenv("IRHEP_BUNDLE_ROOT", here::here()), "scripts", "config.R"))
cat("================================================================\n")
cat(" SUPPLEMENTARY ROI FIGURES\n")
cat("================================================================\n\n")

# --- Paths ---
working_dir <- RAW_DIR
PROCESSED_DIR <- PROCESSED_DIR
PROCESSED_DIR <- PROCESSED_DIR
roi_output   <- file.path(PROCESSED_DIR, "supplementary_rois")
crop_dir     <- file.path(roi_output, "crops")
dir.create(crop_dir, recursive = TRUE, showWarnings = FALSE)

# Source cluster config
source(file.path(SCRIPTS_DIR, "02_cluster_config.R"))

# Sample metadata
samples_info <- data.table(
  sample_id = c("18321_a", "24774_d", "54023_b", "56784_c"),
  condition = c("irHepatitis", "irHepatitis", "AIH", "AIH"),
  x_shift   = c(0, 12000, 0, 12000),
  y_shift   = c(0, 0, -8000, -8000)
)


# ============================================================================
# 1. TOP MARKER GENE PER CLUSTER (for transcript panel)
# ============================================================================

roi_genes <- c(
  "1"  = "ARG1",      # Hepatocyte 1 — canonical hepatocyte marker
  "2"  = "KLRK1",     # NKG2D+ CD8+ T cell — KLRK1 = NKG2D
  "3"  = "CD163",     # Resident macrophage
  "4"  = "DCN",       # Stellate cell
  "5"  = "S100A9",    # Myeloid cell
  "6"  = "SDC1",      # Interface hepatocyte — distinguishes from Cl 1 (ARG1)
  "7"  = "IL7R",      # CD4+ T cell — CD127, canonical naive CD4 marker
  "8"  = "CXCL6",     # Cholangiocyte
  "9"  = "IGHGP",     # IgG+ plasma cell 1 — matches "IgG+"
  "10" = "APOE",      # Hepatocyte 3 — distinguishes from Cl 1, 6
  "11" = "FCGR2B",    # Sinusoidal endothelial cell
  "12" = "IGKC",      # Uncharacterized 1
  "13" = "IGHG1",     # IgG+ plasma cell 2 — distinguishes from Cl 9
  "14" = "SPARCL1",   # Vascular endothelial cell
  "15" = "MARCO",     # Uncharacterized 2
  "16" = "IGHM",      # IgM+ B cell — matches "IgM+"
  "17" = "CD8A",      # CD8+ T cell — matches "CD8+"
  "18" = "ARG1",      # Proliferating hepatocyte — shows hepatocyte lineage
  "19" = "GZMB",      # Proliferating CD8+ T cell — shows cytotoxic CD8 identity
  "20" = "IGHG1"      # IgG+ plasma cell 3 — matches "IgG+"
)


# ============================================================================
# 2. LOAD GIOTTO OBJECT -> SPATIAL COORDINATES
# ============================================================================

cat("Loading Giotto object...\n")
gobj <- readRDS(file.path(PROCESSED_DIR, "gobj_with_spatial_networks.rds"))

meta <- pDataDT(gobj)
coords_dt <- getSpatialLocations(gobj, output = "data.table")
coords_dt <- merge(coords_dt, meta[, .(cell_ID, leiden, sample_id, condition)],
                   by = "cell_ID")
coords_dt <- coords_dt[!leiden %in% as.integer(excluded_clusters)]

cat(sprintf("  %s cells across %d clusters, %d samples\n",
            format(nrow(coords_dt), big.mark = ","),
            length(unique(coords_dt$leiden)),
            length(unique(coords_dt$sample_id))))

rm(gobj)
gc(verbose = FALSE)


# ============================================================================
# 3. FIND OPTIMAL 200x200 ROI PER CLUSTER
# ============================================================================

cat("\nFinding optimal 200x200 um ROIs per cluster...\n")

ROI_SIZE <- 200
GRID_STEP <- 50

find_best_roi <- function(cl_id, coords, grid_step = GRID_STEP) {
  cl_cells <- coords[leiden == cl_id]
  best <- list(count = 0, sample = NA, x = NA, y = NA)

  for (sid in unique(cl_cells$sample_id)) {
    sid_cells <- cl_cells[sample_id == sid]
    if (nrow(sid_cells) < 3) next

    x_range <- range(sid_cells$sdimx)
    y_range <- range(sid_cells$sdimy)

    x_starts <- seq(x_range[1] - ROI_SIZE / 2, x_range[2], by = grid_step)
    y_starts <- seq(y_range[1] - ROI_SIZE / 2, y_range[2], by = grid_step)

    for (xs in x_starts) {
      for (ys in y_starts) {
        n <- sum(sid_cells$sdimx >= xs & sid_cells$sdimx < (xs + ROI_SIZE) &
                   sid_cells$sdimy >= ys & sid_cells$sdimy < (ys + ROI_SIZE))
        if (n > best$count) {
          best <- list(count = n, sample = sid, x = xs, y = ys)
        }
      }
    }
  }
  best
}

roi_specs <- data.table(
  cluster    = integer(),
  cluster_name = character(),
  sample_id  = character(),
  condition  = character(),
  roi_x      = numeric(),
  roi_y      = numeric(),
  roi_size   = numeric(),
  n_target   = integer(),
  n_total    = integer(),
  gene       = character()
)

for (cl in sort(as.integer(names(cell_type_labels)))) {
  cat(sprintf("  Cluster %2d (%s)... ", cl, cell_type_labels[as.character(cl)]))
  best <- find_best_roi(cl, coords_dt)

  if (is.na(best$sample)) {
    cat("SKIPPED (no cells found)\n")
    next
  }

  n_total <- nrow(coords_dt[
    sample_id == best$sample &
      sdimx >= best$x & sdimx < (best$x + ROI_SIZE) &
      sdimy >= best$y & sdimy < (best$y + ROI_SIZE)
  ])

  cond <- samples_info[sample_id == best$sample]$condition
  gene <- roi_genes[as.character(cl)]

  cat(sprintf("sample=%s, n_target=%d, n_total=%d, gene=%s\n",
              best$sample, best$count, n_total, gene))

  roi_specs <- rbindlist(list(roi_specs, data.table(
    cluster = cl,
    cluster_name = cell_type_labels[as.character(cl)],
    sample_id = best$sample,
    condition = cond,
    roi_x = best$x,
    roi_y = best$y,
    roi_size = ROI_SIZE,
    n_target = best$count,
    n_total = n_total,
    gene = gene
  )))
}

fwrite(roi_specs, file.path(roi_output, "roi_summary.csv"))
cat(sprintf("\n  ROI summary saved: %d clusters\n", nrow(roi_specs)))


# ============================================================================
# 4. WRITE JSON CONFIG + CALL PYTHON FOR IMAGE EXTRACTION
# ============================================================================

cat("\nPreparing image extraction...\n")

roi_list <- lapply(1:nrow(roi_specs), function(i) {
  list(
    cluster_id = roi_specs$cluster[i],
    sample_id  = roi_specs$sample_id[i],
    roi_x      = roi_specs$roi_x[i],
    roi_y      = roi_specs$roi_y[i],
    roi_size   = roi_specs$roi_size[i]
  )
})

config_json <- list(
  working_dir = working_dir,
  output_dir  = crop_dir,
  pixel_size  = 0.2125,
  rois        = roi_list
)

json_path <- file.path(crop_dir, "roi_config.json")
write_json(config_json, json_path, auto_unbox = TRUE, pretty = TRUE)
cat(sprintf("  Config written: %s\n", json_path))

# Call Python ROI-crop extractor (shared with 07_panel_c_he_extract.R).
py_script <- file.path(SCRIPTS_DIR, "07a_panel_c_he_extract.py")
cat("  Running Python image extraction...\n")
py_cmd <- sprintf("python3 '%s' '%s'", py_script, json_path)
py_result <- system(py_cmd, intern = FALSE)

if (py_result != 0) {
  cat("  WARNING: Python extraction returned non-zero exit code.\n")
  cat("  Will proceed with available crops (missing panels shown as grey).\n")
}


# ============================================================================
# 5. LOAD PARQUET DATA PER SAMPLE (cell boundaries + transcripts)
# ============================================================================

cat("\nLoading parquet data...\n")

needed_samples <- unique(roi_specs$sample_id)
boundary_cache <- list()
transcript_cache <- list()

for (sid in needed_samples) {
  dir_name <- paste0("0027420_", sid)
  sample_path <- file.path(working_dir, dir_name)
  x_off <- samples_info[sample_id == sid]$x_shift
  y_off <- samples_info[sample_id == sid]$y_shift

  cb_file <- file.path(sample_path, "cell_boundaries.parquet")
  if (file.exists(cb_file)) {
    cat(sprintf("  Loading cell boundaries: %s...\n", sid))
    cb <- as.data.table(read_parquet(cb_file))
    cb[, `:=`(x_plot = vertex_x + x_off, y_plot = -vertex_y + y_off)]
    boundary_cache[[sid]] <- cb
  }

  nb_file <- file.path(sample_path, "nucleus_boundaries.parquet")
  if (file.exists(nb_file)) {
    cat(sprintf("  Loading nucleus boundaries: %s...\n", sid))
    nb <- as.data.table(read_parquet(nb_file))
    nb[, `:=`(x_plot = vertex_x + x_off, y_plot = -vertex_y + y_off)]
    attr(boundary_cache[[sid]], "nucleus") <- nb
  }

  tx_file <- file.path(sample_path, "transcripts.parquet")
  if (file.exists(tx_file)) {
    needed_genes <- roi_specs[sample_id == sid]$gene
    cat(sprintf("  Loading transcripts: %s (genes: %s)...\n",
                sid, paste(unique(needed_genes), collapse = ", ")))
    tx <- as.data.table(read_parquet(tx_file))
    tx <- tx[feature_name %in% unique(needed_genes)]
    tx[, `:=`(x_plot = x_location + x_off, y_plot = -y_location + y_off)]
    transcript_cache[[sid]] <- tx
  }
}

cat("  Parquet data loaded.\n")


# ============================================================================
# 6. ASSEMBLE 4-PANEL FIGURES
# ============================================================================

cat("\nGenerating supplementary ROI figures...\n")

theme_roi <- theme_void() +
  theme(
    plot.title = element_text(size = 11, face = "bold", hjust = 0.5,
                              margin = margin(b = 3)),
    plot.margin = margin(2, 2, 2, 2)
  )

for (i in 1:nrow(roi_specs)) {
  cl      <- roi_specs$cluster[i]
  cl_name <- roi_specs$cluster_name[i]
  sid     <- roi_specs$sample_id[i]
  cond    <- roi_specs$condition[i]
  rx      <- roi_specs$roi_x[i]
  ry      <- roi_specs$roi_y[i]
  rs      <- roi_specs$roi_size[i]
  gene    <- roi_specs$gene[i]
  n_tgt   <- roi_specs$n_target[i]
  n_tot   <- roi_specs$n_total[i]

  x_off <- samples_info[sample_id == sid]$x_shift
  y_off <- samples_info[sample_id == sid]$y_shift

  cat(sprintf("  [%2d/20] Cluster %2d: %s (%s, %s)...\n",
              i, cl, cl_name, sid, gene))

  xlims <- c(rx, rx + rs)
  ylims <- c(ry, ry + rs)

  # Panel 1: H&E
  he_file <- file.path(crop_dir, sprintf("cluster_%02d_he.png", cl))
  if (file.exists(he_file)) {
    he_img <- readPNG(he_file)
    p_he <- ggplot() +
      annotation_raster(he_img,
                        xmin = xlims[1], xmax = xlims[2],
                        ymin = ylims[1], ymax = ylims[2]) +
      coord_fixed(xlim = xlims, ylim = ylims, expand = FALSE) +
      ggtitle("H&E") +
      theme_roi
  } else {
    p_he <- ggplot() +
      annotate("rect", xmin = xlims[1], xmax = xlims[2],
               ymin = ylims[1], ymax = ylims[2], fill = "grey90") +
      annotate("text", x = mean(xlims), y = mean(ylims),
               label = "H&E\nnot available", size = 4, color = "grey50") +
      coord_fixed(xlim = xlims, ylim = ylims, expand = FALSE) +
      ggtitle("H&E") +
      theme_roi
  }

  # Panel 2: DAPI + scale bar
  dapi_file <- file.path(crop_dir, sprintf("cluster_%02d_dapi.png", cl))
  if (file.exists(dapi_file)) {
    dapi_img <- readPNG(dapi_file)
    p_dapi <- ggplot() +
      annotation_raster(dapi_img,
                        xmin = xlims[1], xmax = xlims[2],
                        ymin = ylims[1], ymax = ylims[2]) +
      annotate("segment",
               x = xlims[2] - 60, xend = xlims[2] - 10,
               y = ylims[1] + 10, yend = ylims[1] + 10,
               color = "white", linewidth = 1.5) +
      annotate("text",
               x = xlims[2] - 35, y = ylims[1] + 20,
               label = "50 \u00b5m", color = "white",
               size = 2.8, fontface = "bold") +
      coord_fixed(xlim = xlims, ylim = ylims, expand = FALSE) +
      ggtitle("DAPI") +
      theme_roi
  } else {
    p_dapi <- ggplot() +
      annotate("rect", xmin = xlims[1], xmax = xlims[2],
               ymin = ylims[1], ymax = ylims[2], fill = "grey20") +
      annotate("text", x = mean(xlims), y = mean(ylims),
               label = "DAPI\nnot available", size = 4, color = "grey70") +
      coord_fixed(xlim = xlims, ylim = ylims, expand = FALSE) +
      ggtitle("DAPI") +
      theme_roi
  }

  # Panel 3: Clusters + Cell Borders
  p_clusters <- ggplot()

  cb <- boundary_cache[[sid]]
  nuc <- attr(boundary_cache[[sid]], "nucleus")

  if (!is.null(cb)) {
    margin <- 10
    cb_roi <- cb[x_plot >= (xlims[1] - margin) & x_plot <= (xlims[2] + margin) &
                   y_plot >= (ylims[1] - margin) & y_plot <= (ylims[2] + margin)]

    if (nrow(cb_roi) > 0) {
      cb_centroids <- cb_roi[, .(cx = mean(x_plot), cy = mean(y_plot)),
                              by = cell_id]

      roi_cells <- coords_dt[sample_id == sid &
                               sdimx >= xlims[1] & sdimx <= xlims[2] &
                               sdimy >= ylims[1] & sdimy <= ylims[2]]

      if (nrow(roi_cells) > 0) {
        cb_centroids[, leiden := {
          dists <- (roi_cells$sdimx - cx)^2 + (roi_cells$sdimy - cy)^2
          roi_cells$leiden[which.min(dists)]
        }, by = cell_id]

        cb_roi_cl <- merge(cb_roi, cb_centroids[, .(cell_id, leiden)],
                           by = "cell_id")
        cb_roi_cl[, is_target := leiden == cl]

        cb_other <- cb_roi_cl[is_target == FALSE]
        if (nrow(cb_other) > 0) {
          p_clusters <- p_clusters +
            geom_polygon(data = cb_other,
                         aes(x = x_plot, y = y_plot, group = cell_id),
                         fill = NA, color = "grey70", linewidth = 0.08)
        }

        cb_target <- cb_roi_cl[is_target == TRUE]
        cl_color <- cluster_colors_by_number[as.character(cl)]
        if (nrow(cb_target) > 0) {
          p_clusters <- p_clusters +
            geom_polygon(data = cb_target,
                         aes(x = x_plot, y = y_plot, group = cell_id),
                         fill = cl_color, color = "grey40",
                         linewidth = 0.1, alpha = 0.6)
        }
      }
    }

    if (!is.null(nuc)) {
      nuc_roi <- nuc[x_plot >= (xlims[1] - margin) &
                       x_plot <= (xlims[2] + margin) &
                       y_plot >= (ylims[1] - margin) &
                       y_plot <= (ylims[2] + margin)]
      if (nrow(nuc_roi) > 0) {
        p_clusters <- p_clusters +
          geom_polygon(data = nuc_roi,
                       aes(x = x_plot, y = y_plot, group = cell_id),
                       fill = NA, color = "grey50", linewidth = 0.05)
      }
    }
  }

  p_clusters <- p_clusters +
    coord_fixed(xlim = xlims, ylim = ylims, expand = FALSE) +
    ggtitle(paste0(cl_name, " (n=", n_tgt, ")")) +
    theme_roi

  # Panel 4: Gene Transcripts + nucleus outlines
  p_transcripts <- ggplot()

  if (!is.null(nuc)) {
    nuc_roi2 <- nuc[x_plot >= (xlims[1] - margin) &
                      x_plot <= (xlims[2] + margin) &
                      y_plot >= (ylims[1] - margin) &
                      y_plot <= (ylims[2] + margin)]
    if (nrow(nuc_roi2) > 0) {
      p_transcripts <- p_transcripts +
        geom_polygon(data = nuc_roi2,
                     aes(x = x_plot, y = y_plot, group = cell_id),
                     fill = "grey92", color = "grey70", linewidth = 0.08)
    }
  }

  if (!is.null(cb)) {
    cb_roi2 <- cb[x_plot >= (xlims[1] - margin) & x_plot <= (xlims[2] + margin) &
                    y_plot >= (ylims[1] - margin) & y_plot <= (ylims[2] + margin)]
    if (nrow(cb_roi2) > 0) {
      p_transcripts <- p_transcripts +
        geom_polygon(data = cb_roi2,
                     aes(x = x_plot, y = y_plot, group = cell_id),
                     fill = NA, color = "grey50", linewidth = 0.12)
    }
  }

  tx <- transcript_cache[[sid]]
  n_tx <- 0
  if (!is.null(tx)) {
    tx_roi <- tx[feature_name == gene &
                   x_plot >= xlims[1] & x_plot <= xlims[2] &
                   y_plot >= ylims[1] & y_plot <= ylims[2]]
    n_tx <- nrow(tx_roi)
    if (n_tx > 0) {
      p_transcripts <- p_transcripts +
        geom_point(data = tx_roi,
                   aes(x = x_plot, y = y_plot),
                   color = "red", size = 0.8, alpha = 0.8)
    }
  }

  p_transcripts <- p_transcripts +
    coord_fixed(xlim = xlims, ylim = ylims, expand = FALSE) +
    ggtitle(paste0(gene, " (n=", n_tx, ")")) +
    theme_roi

  # Assemble 2x2 grid with title
  top_row <- (p_he | p_dapi)
  bottom_row <- (p_clusters | p_transcripts)
  combined <- top_row / bottom_row

  combined <- combined + plot_annotation(
    title = sprintf("Cluster %d: %s", cl, cl_name),
    subtitle = sprintf("Sample %s (%s) | 200 \u00d7 200 \u00b5m | %d cells (%d %s)",
                       sid, cond, n_tot, n_tgt, cl_name),
    theme = theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey40")
    )
  )

  out_file <- file.path(roi_output,
                         sprintf("SuppFig_ROI_Cluster%02d_%s.pdf",
                                 cl, gsub("[^A-Za-z0-9]", "_", cl_name)))
  cairo_pdf(out_file, width = 10, height = 10)
  print(combined)
  dev.off()

  cat(sprintf("    -> %s\n", basename(out_file)))
}


# ============================================================================
# 7. SUMMARY
# ============================================================================

cat("\n================================================================\n")
cat(" SUPPLEMENTARY ROI FIGURES COMPLETE\n")
cat("================================================================\n\n")
cat(sprintf("Output: %s\n", roi_output))
cat(sprintf("  %d PDF figures generated\n", nrow(roi_specs)))
cat(sprintf("  ROI summary: %s\n", file.path(roi_output, "roi_summary.csv")))
cat(sprintf("  Image crops: %s\n\n", crop_dir))

cat("ROI Summary:\n")
cat(sprintf("  %-4s %-28s %-10s %-12s %6s %6s %-8s\n",
            "Cl", "Name", "Sample", "Condition", "Target", "Total", "Gene"))
cat(paste(rep("-", 82), collapse = ""), "\n")
for (i in 1:nrow(roi_specs)) {
  cat(sprintf("  %-4d %-28s %-10s %-12s %6d %6d %-8s\n",
              roi_specs$cluster[i], roi_specs$cluster_name[i],
              roi_specs$sample_id[i], roi_specs$condition[i],
              roi_specs$n_target[i], roi_specs$n_total[i],
              roi_specs$gene[i]))
}
cat("\n")
