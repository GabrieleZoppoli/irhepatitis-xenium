# ============================================================================
# 01_giotto_normalize_cluster.R
# Takes the per-sample Xenium Analyzer outputs (already-segmented cells with
# decoded transcripts, in `data/raw/0027420_<sample_id>/`) and produces the
# integrated Giotto session that all downstream analyses consume.
#
# Per-sample inputs (Xenium-Analyzer outputs, NOT raw imaging):
#   cell_feature_matrix.h5, cells.parquet, transcripts.parquet,
#   cell_boundaries.parquet, nucleus_boundaries.parquet
#
# Pipeline:
#   1. Read the four per-sample cell-feature matrices into Giotto objects.
#   2. QC-filter (>= 4 transcripts/cell, >= 3 cells/gene, qv >= 20).
#   3. Log-normalize (library-size, scale factor 5000).
#   4. PCA (50 PCs) -> Harmony batch correction (theta = 2, sigma = 0.1,
#      lambda = 1).
#   5. UMAP and SNN-Leiden clustering (k = 15, resolution = 0.5).
#   6. Save the integrated multi-sample object to
#      `data/processed/session_snapshot.RData`.
# ============================================================================
source(file.path(Sys.getenv("IRHEP_BUNDLE_ROOT", here::here()), "scripts", "config.R"))

# Explicit script-level seed (belt-and-braces). Giotto's stochastic functions
# (runPCA, runGiottoHarmony, runUMAP, doLeidenCluster) also default to
# `set_seed = TRUE, seed_number = 1234`, but we set the R-side seed here so
# any non-Giotto random call earlier in the script (rare) is also reproducible.
set.seed(1234)

# ==============================================================================
# ENVIRONMENT SETUP
# ==============================================================================

library(Giotto)
library(GiottoClass)
library(data.table)
library(Matrix)
library(rhdf5)
library(arrow)
library(harmony)
library(ggplot2)
library(terra)
library(jsonlite)
library(RColorBrewer)
library(scales)

# Python environment
python_path <- Sys.getenv("RETICULATE_PYTHON",
  unset = system("which python", intern = TRUE)[1])
reticulate::use_python(python_path, required = TRUE)

# Parallel processing
library(future)
plan("multicore", workers = 8)
options(future.globals.maxSize = 16000 * 1024^2)

# ==============================================================================
# CONFIGURATION
# ==============================================================================

root_dir <- RAW_DIR
out_dir  <- PROCESSED_DIR

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
for (d in c("plots", "tables", "objects")) {
  dir.create(file.path(out_dir, d), showWarnings = FALSE)
}

instrs <- createGiottoInstructions(
  python_path = python_path,
  save_dir = file.path(out_dir, "plots"),
  save_plot = TRUE,
  show_plot = FALSE
)


params <- list(
  # Filtering
  min_transcripts_per_cell = 4,
  min_cells_per_gene = 3,
  qv_threshold = 20,

  # Normalization
  scale_factor = 5000,

  # Dimensionality reduction (optimized for IO panel ~380 genes)
  n_pcs = 50,
  harmony_dims = 15,

  # Harmony batch correction
  harmony_theta = 2,        # Diversity penalty (increase if samples don't mix)
  harmony_sigma = 0.1,      # Soft k-means width
  harmony_lambda = 1,       # Ridge regression penalty

  # UMAP (tuned for visual clarity)
  umap_dims = 15,
  umap_neighbors = 15,
  umap_min_dist = 0.0005,

  # Clustering
  snn_k = 15,
  leiden_resolution = 0.5
)


cluster_palette_21 <- c(
  "1"  = "#E64B35", "2"  = "#4DBBD5", "3"  = "#00A087", "4"  = "#3C5488",
  "5"  = "#F39B7F", "6"  = "#8491B4", "7"  = "#91D1C2", "8"  = "#DC0000",
  "9"  = "#7E6148", "10" = "#B09C85", "11" = "#00A6A6", "12" = "#E9C46A",
  "13" = "#F4A261", "14" = "#E76F51", "15" = "#264653", "16" = "#2A9D8F",
  "17" = "#A8DADC", "18" = "#457B9D", "19" = "#1D3557", "20" = "#606060",
  "21" = "#9B5DE5"
)

condition_colors <- c("irHepatitis" = "#E64B35", "AIH" = "#4DBBD5")

sample_colors <- c(
  "18321_a" = "#E64B35", "24774_d" = "#F39B7F",
  "54023_b" = "#4DBBD5", "56784_c" = "#3C5488"
)

color_config <- list(
  cluster_colors = cluster_palette_21,
  condition_colors = condition_colors,
  sample_colors = sample_colors,
  n_clusters = 21
)

# ==============================================================================
# SAMPLE METADATA
# ==============================================================================

samples <- data.table(
  sample_id  = c("18321_a", "24774_d", "54023_b", "56784_c"),
  dir_name   = c("0027420_18321_a", "0027420_24774_d", "0027420_54023_b", "0027420_56784_c"),
  condition  = c("irHepatitis", "irHepatitis", "AIH", "AIH"),
  x_shift    = c(0, 12000, 0, 12000),
  y_shift    = c(0, 0, -8000, -8000)
)

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

read_xenium_h5 <- function(h5_path) {
  barcodes   <- h5read(h5_path, "/matrix/barcodes")
  data       <- h5read(h5_path, "/matrix/data")
  indices    <- h5read(h5_path, "/matrix/indices")
  indptr     <- h5read(h5_path, "/matrix/indptr")
  shape      <- h5read(h5_path, "/matrix/shape")
  gene_names <- h5read(h5_path, "/matrix/features/name")

  sparseMatrix(
    i = indices + 1, p = indptr, x = data,
    dims = shape, dimnames = list(gene_names, barcodes)
  )
}

# ==============================================================================
# PIPELINE START
# ==============================================================================

cat("════════════════════════════════════════════════════════════════════════\n")
cat("   XENIUM IO PANEL — PIPELINE                                            \n")
cat("════════════════════════════════════════════════════════════════════════\n\n")

cat("Parameters:\n")
for (p in names(params)) {
  cat(sprintf("  %s = %s\n", p, params[[p]]))
}
cat("\n")

setwd(root_dir)

# ==============================================================================
# STEP 1: LOAD DATA
# ==============================================================================

cat("═══════════════════════════════════════════════════════════════════════\n")
cat(" STEP 1: Loading Xenium data\n")
cat("═══════════════════════════════════════════════════════════════════════\n")

all_matrices <- list()
all_coords   <- list()
all_meta     <- list()
sample_info  <- list()

for (i in 1:nrow(samples)) {
  sid      <- samples$sample_id[i]
  dir_name <- samples$dir_name[i]
  cond     <- samples$condition[i]
  x_off    <- samples$x_shift[i]
  y_off    <- samples$y_shift[i]

  cat(sprintf("\n[%d/%d] %s (%s)\n", i, nrow(samples), sid, cond))
  sample_path <- file.path(root_dir, dir_name)

  h5_file <- file.path(sample_path, "cell_feature_matrix.h5")
  expr_mat <- read_xenium_h5(h5_file)
  cat(sprintf("      Expression: %d genes × %d cells\n", nrow(expr_mat), ncol(expr_mat)))

  cells_file <- file.path(sample_path, "cells.parquet")
  cells_dt <- as.data.table(read_parquet(cells_file))
  cell_ids <- colnames(expr_mat)
  cells_dt <- cells_dt[match(cell_ids, cell_id)]

  new_ids <- paste0(sid, "_", cell_ids)
  colnames(expr_mat) <- new_ids

  coords_dt <- data.table(
    cell_ID = new_ids,
    sdimx   = cells_dt$x_centroid + x_off,
    sdimy   = -(cells_dt$y_centroid) + y_off
  )

  meta_dt <- data.table(
    cell_ID   = new_ids,
    sample_id = sid,
    condition = cond,
    orig_id   = cell_ids
  )

  all_matrices[[sid]] <- expr_mat
  all_coords[[sid]]   <- coords_dt
  all_meta[[sid]]     <- meta_dt
  sample_info[[sid]]  <- list(path = sample_path, x_shift = x_off, y_shift = y_off)
}

# ==============================================================================
# STEP 2: COMBINE SAMPLES
# ==============================================================================

cat("\n═══════════════════════════════════════════════════════════════════════\n")
cat(" STEP 2: Combining samples\n")
cat("═══════════════════════════════════════════════════════════════════════\n")

combined_coords <- rbindlist(all_coords)
combined_meta   <- rbindlist(all_meta)
all_genes       <- Reduce(union, lapply(all_matrices, rownames))

cat(sprintf("\nTotal cells: %d\n", nrow(combined_coords)))
cat(sprintf("Total genes (union): %d\n", length(all_genes)))

combined_expr <- Matrix(0, nrow = length(all_genes), ncol = nrow(combined_coords),
                        sparse = TRUE, dimnames = list(all_genes, combined_coords$cell_ID))

col_offset <- 0
for (sid in names(all_matrices)) {
  mat <- all_matrices[[sid]]
  n_cells <- ncol(mat)
  row_idx <- match(rownames(mat), all_genes)
  col_idx <- (col_offset + 1):(col_offset + n_cells)
  combined_expr[row_idx, col_idx] <- mat
  col_offset <- col_offset + n_cells
}

rm(all_matrices); gc()

# ==============================================================================
# STEP 3: REMOVE CONTROL PROBES
# ==============================================================================

cat("\n═══════════════════════════════════════════════════════════════════════\n")
cat(" STEP 3: Removing control probes\n")
cat("═══════════════════════════════════════════════════════════════════════\n")

control_patterns <- c("NegControlProbe", "NegControlCodeword", "UnassignedCodeword",
                      "antisense", "BLANK", "Codeword", "Unassigned")
control_genes <- grep(paste(control_patterns, collapse = "|"),
                      rownames(combined_expr), ignore.case = TRUE, value = TRUE)

cat(sprintf("\nControl probes removed: %d\n", length(control_genes)))
combined_expr <- combined_expr[setdiff(rownames(combined_expr), control_genes), ]
cat(sprintf("Genes retained: %d\n", nrow(combined_expr)))

# ==============================================================================
# STEP 4: CREATE GIOTTO OBJECT
# ==============================================================================

cat("\n═══════════════════════════════════════════════════════════════════════\n")
cat(" STEP 4: Creating Giotto object\n")
cat("═══════════════════════════════════════════════════════════════════════\n")

gobj <- createGiottoObject(
  expression   = combined_expr,
  spatial_locs = combined_coords,
  instructions = instrs
)

gobj <- addCellMetadata(gobj, new_metadata = combined_meta,
                        by_column = TRUE, column_cell_ID = "cell_ID")

cat("\n✓ Giotto object created\n")
print(table(pDataDT(gobj)$sample_id))

# ==============================================================================
# STEP 5: ADD DAPI IMAGES
# ==============================================================================

cat("\n═══════════════════════════════════════════════════════════════════════\n")
cat(" STEP 5: Adding DAPI images\n")
cat("═══════════════════════════════════════════════════════════════════════\n")

for (i in 1:nrow(samples)) {

  sid      <- samples$sample_id[i]
  dir_name <- samples$dir_name[i]
  x_off    <- samples$x_shift[i]
  y_off    <- samples$y_shift[i]

  cat(sprintf("\n%s: ", sid))

  tryCatch({
    temp_gobj <- createGiottoXeniumObject(
      xenium_dir = file.path(root_dir, dir_name),
      qv_threshold = params$qv_threshold,
      instructions = instrs
    )

    # DAPI
    dapi_img <- temp_gobj@images[["dapi"]]
    dapi_img@name <- paste0(sid, "_DAPI")
    dapi_img <- spatShift(dapi_img, dx = x_off, dy = y_off)
    gobj <- addGiottoImage(gobj, images = list(dapi_img))

    # Boundary staining
    for (bnd in c("bound1", "bound2", "bound3")) {
      if (bnd %in% names(temp_gobj@images)) {
        bnd_img <- temp_gobj@images[[bnd]]
        bnd_img@name <- paste0(sid, "_", bnd)
        bnd_img <- spatShift(bnd_img, dx = x_off, dy = y_off)
        gobj <- addGiottoImage(gobj, images = list(bnd_img))
      }
    }

    rm(temp_gobj); gc()
    cat("✓\n")

  }, error = function(e) cat(sprintf("✗ %s\n", e$message)))
}

# ==============================================================================
# STEP 6: ADD H&E IMAGES
# ==============================================================================

cat("\n═══════════════════════════════════════════════════════════════════════\n")
cat(" STEP 6: Adding H&E images\n")
cat("═══════════════════════════════════════════════════════════════════════\n")

for (sid in names(sample_info)) {

  info <- sample_info[[sid]]
  sample_path <- info$path
  x_shift <- info$x_shift
  y_shift <- info$y_shift

  he_file <- file.path(sample_path, paste0(sid, "_he_image.ome.tif"))
  if (!file.exists(he_file)) {
    cat(sprintf("\n%s: H&E not found, skipping\n", sid))
    next
  }

  # Pixel size from metadata
  pixel_size <- 0.2125
  exp_file <- file.path(sample_path, "experiment.xenium")
  if (file.exists(exp_file)) {
    tryCatch({
      exp_data <- fromJSON(exp_file)
      if (!is.null(exp_data$pixel_size)) pixel_size <- exp_data$pixel_size
    }, error = function(e) NULL)
  }

  cat(sprintf("\n%s: ", sid))

  tryCatch({
    he_img <- createGiottoLargeImage(
      raster_object = he_file,
      name = paste0(sid, "_HE"),
      negative_y = TRUE
    )

    # Affine transformation
    align_file <- file.path(sample_path, "he_alignment", "matrix.csv")
    if (file.exists(align_file)) {
      affine_raw <- as.matrix(read.csv(align_file, header = FALSE))
      tfm <- matrix(c(1,0,0, 0,1,0, 0,0,1), nrow = 3, byrow = TRUE)
      tfm[1,1] <- affine_raw[1,1] * pixel_size
      tfm[1,2] <- affine_raw[1,2] * pixel_size
      tfm[2,1] <- affine_raw[2,1] * pixel_size
      tfm[2,2] <- affine_raw[2,2] * pixel_size
      tfm[1,3] <- affine_raw[1,3] * pixel_size + x_shift
      tfm[2,3] <- -affine_raw[2,3] * pixel_size + y_shift
      he_img <- affine(he_img, tfm)
    } else {
      he_img <- spatShift(he_img, dx = x_shift, dy = y_shift)
    }

    gobj <- addGiottoImage(gobj, images = list(he_img))
    cat("✓\n")

  }, error = function(e) cat(sprintf("✗ %s\n", e$message)))
}

cat("\nImages loaded:\n")
print(showGiottoImageNames(gobj))

# ==============================================================================
# STEP 7: QUALITY CONTROL AND FILTERING
# ==============================================================================

cat("\n═══════════════════════════════════════════════════════════════════════\n")
cat(" STEP 7: Quality control and filtering\n")
cat("═══════════════════════════════════════════════════════════════════════\n")

gobj <- addStatistics(gobj, expression_values = "raw")
meta <- pDataDT(gobj)

cat(sprintf("\nBefore filtering: %d cells\n", nrow(meta)))
cat(sprintf("  Median transcripts/cell: %.0f\n", median(meta$total_expr)))
cat(sprintf("  Median genes/cell: %.0f\n", median(meta$nr_feats)))

gobj <- filterGiotto(gobj,
                     expression_threshold = 1,
                     feat_det_in_min_cells = params$min_cells_per_gene,
                     min_det_feats_per_cell = params$min_transcripts_per_cell)

gobj <- addStatistics(gobj, expression_values = "raw")
meta <- pDataDT(gobj)

cat(sprintf("\nAfter filtering: %d cells, %d genes\n", nrow(meta), length(featIDs(gobj))))
print(table(meta$sample_id))

saveRDS(gobj, file.path(out_dir, "objects", "gobj_01_filtered.rds"))

# ==============================================================================
# STEP 8: NORMALIZATION
# ==============================================================================

cat("\n═══════════════════════════════════════════════════════════════════════\n")
cat(" STEP 8: Normalization\n")
cat("═══════════════════════════════════════════════════════════════════════\n")

gobj <- normalizeGiotto(gobj, scalefactor = params$scale_factor, verbose = FALSE)
gobj <- addStatistics(gobj, expression_values = "normalized")

cat(sprintf("\n✓ Library size normalization (scale factor = %d)\n", params$scale_factor))

# ==============================================================================
# STEP 9: PCA
# ==============================================================================

cat("\n═══════════════════════════════════════════════════════════════════════\n")
cat(" STEP 9: Principal component analysis\n")
cat("═══════════════════════════════════════════════════════════════════════\n")

gobj <- runPCA(gobj,
               expression_values = "normalized",
               feats_to_use = featIDs(gobj),
               center = TRUE,
               scale_unit = TRUE,
               ncp = params$n_pcs)

cat(sprintf("\n✓ PCA complete (%d components)\n", params$n_pcs))

screePlot(gobj, ncp = params$n_pcs,
          save_plot = TRUE,
          save_param = list(save_name = "01_scree_plot"))

# ==============================================================================
# STEP 10: HARMONY BATCH CORRECTION
# ==============================================================================

cat("\n═══════════════════════════════════════════════════════════════════════\n")
cat(" STEP 10: Harmony batch correction\n")
cat("═══════════════════════════════════════════════════════════════════════\n")

cat(sprintf("  theta = %d (diversity penalty)\n", params$harmony_theta))
cat(sprintf("  sigma = %.2f (soft k-means width)\n", params$harmony_sigma))
cat(sprintf("  lambda = %d (ridge regression penalty)\n", params$harmony_lambda))

gobj <- runGiottoHarmony(gobj,
                         vars_use = "sample_id",
                         do_pca = FALSE,
                         dim_reduction_to_use = "pca",
                         dim_reduction_name = "pca",
                         dimensions_to_use = 1:params$harmony_dims,
                         theta = params$harmony_theta,
                         sigma = params$harmony_sigma,
                         lambda = params$harmony_lambda,
                         name = "harmony")

cat(sprintf("\n✓ Harmony complete (dims 1:%d, batch = sample_id)\n", params$harmony_dims))

# ==============================================================================
# STEP 11: UMAP
# ==============================================================================

cat("\n═══════════════════════════════════════════════════════════════════════\n")
cat(" STEP 11: UMAP embedding\n")
cat("═══════════════════════════════════════════════════════════════════════\n")

gobj <- runUMAP(gobj,
                dim_reduction_to_use = "harmony",
                dim_reduction_name = "harmony",
                dimensions_to_use = 1:params$umap_dims,
                n_neighbors = params$umap_neighbors,
                min_dist = params$umap_min_dist,
                name = "umap")

cat(sprintf("\n✓ UMAP complete (dims 1:%d, n_neighbors=%d, min_dist=%.4f)\n",
            params$umap_dims, params$umap_neighbors, params$umap_min_dist))

# ==============================================================================
# STEP 12: CLUSTERING
# ==============================================================================

cat("\n═══════════════════════════════════════════════════════════════════════\n")
cat(" STEP 12: Leiden clustering\n")
cat("═══════════════════════════════════════════════════════════════════════\n")

gobj <- createNearestNetwork(gobj,
                             dim_reduction_to_use = "harmony",
                             dim_reduction_name = "harmony",
                             dimensions_to_use = 1:params$umap_dims,
                             k = params$snn_k,
                             name = "sNN.harmony")

gobj <- doLeidenCluster(gobj,
                        network_name = "sNN.harmony",
                        resolution = params$leiden_resolution,
                        n_iterations = 200,
                        name = "leiden")

meta <- pDataDT(gobj)
n_clusters <- length(unique(meta$leiden))

cat(sprintf("\n✓ Clustering complete (k=%d, resolution=%.2f)\n",
            params$snn_k, params$leiden_resolution))
cat(sprintf("  Clusters identified: %d\n", n_clusters))

# Use predefined colors
cluster_colors <- cluster_palette_21[as.character(1:min(n_clusters, 21))]
if (n_clusters > 21) {
  extra <- colorRampPalette(c("#FF6B6B", "#4ECDC4", "#45B7D1"))(n_clusters - 21)
  cluster_colors <- c(cluster_colors, setNames(extra, as.character(22:n_clusters)))
}

saveRDS(gobj, file.path(out_dir, "objects", "gobj_02_clustered.rds"))

# ==============================================================================
# STEP 13: VISUALIZATION
# ==============================================================================

cat("\n═══════════════════════════════════════════════════════════════════════\n")
cat(" STEP 13: Generating plots\n")
cat("═══════════════════════════════════════════════════════════════════════\n")

# UMAP - Clusters
plotUMAP(gobj, cell_color = "leiden", point_size = 0.5,
         cell_color_code = cluster_colors,
         title = sprintf("Leiden clusters (n=%d)", n_clusters),
         save_plot = TRUE,
         save_param = list(save_name = "02_UMAP_leiden", base_width = 10))
cat("  ✓ UMAP - clusters\n")

# UMAP - Samples
plotUMAP(gobj, cell_color = "sample_id", point_size = 0.5,
         cell_color_code = sample_colors,
         title = "Sample distribution (Harmony-corrected)",
         save_plot = TRUE,
         save_param = list(save_name = "03_UMAP_samples", base_width = 10))
cat("  ✓ UMAP - samples\n")

# UMAP - Condition
plotUMAP(gobj, cell_color = "condition", point_size = 0.5,
         cell_color_code = condition_colors,
         title = "Condition",
         save_plot = TRUE,
         save_param = list(save_name = "04_UMAP_condition", base_width = 10))
cat("  ✓ UMAP - condition\n")

# Spatial - Clusters
spatPlot2D(gobj, cell_color = "leiden", point_size = 0.3,
           cell_color_code = cluster_colors,
           title = "Spatial distribution - Clusters",
           save_plot = TRUE,
           save_param = list(save_name = "05_spatial_leiden",
                             base_width = 14, base_height = 10))
cat("  ✓ Spatial - clusters\n")

# Spatial - Condition
spatPlot2D(gobj, cell_color = "condition", point_size = 0.3,
           cell_color_code = condition_colors,
           title = "Spatial distribution - Condition",
           save_plot = TRUE,
           save_param = list(save_name = "06_spatial_condition",
                             base_width = 14, base_height = 10))
cat("  ✓ Spatial - condition\n")

# Per-sample H&E overlays
meta <- pDataDT(gobj)
img_names <- names(gobj@images)

for (sid in samples$sample_id) {
  he_name <- paste0(sid, "_HE")
  if (he_name %in% img_names) {
    sample_cells <- meta[sample_id == sid]$cell_ID
    gobj_sub <- subsetGiotto(gobj, cell_ids = sample_cells)

    tryCatch({
      spatPlot2D(gobj_sub,
                 cell_color = "leiden",
                 cell_color_code = cluster_colors,
                 show_image = TRUE,
                 image_name = he_name,
                 point_size = 0.6,
                 point_alpha = 0.8,
                 title = paste0(sid, " - Clusters on H&E"),
                 save_plot = TRUE,
                 save_param = list(save_name = paste0("07_HE_", sid),
                                   base_width = 12, base_height = 10))
      cat(sprintf("  ✓ H&E overlay - %s\n", sid))
    }, error = function(e) cat(sprintf("  ✗ H&E overlay - %s: %s\n", sid, e$message)))
  }
}

# ==============================================================================
# STEP 14: CLUSTER MARKERS
# ==============================================================================

cat("\n═══════════════════════════════════════════════════════════════════════\n")
cat(" STEP 14: Cluster markers (scran method with FDR)\n")
cat("═══════════════════════════════════════════════════════════════════════\n")

markers <- findMarkers_one_vs_all(gobj,
                                  method = "scran",
                                  expression_values = "normalized",
                                  cluster_column = "leiden",
                                  min_feats = 3)

fwrite(markers, file.path(out_dir, "tables", "cluster_markers.csv"))
cat(sprintf("\n✓ %d marker-cluster associations identified\n", nrow(markers)))

# Heatmap
top_genes <- unique(markers[, head(.SD, 5), by = cluster]$feats)
plotMetaDataHeatmap(gobj,
                    expression_values = "normalized",
                    metadata_cols = "leiden",
                    selected_feats = top_genes,
                    save_plot = TRUE,
                    save_param = list(save_name = "08_marker_heatmap",
                                      base_width = 12, base_height = 10))
cat("  ✓ Marker heatmap saved\n")

# ==============================================================================
# STEP 15: SAVE RESULTS
# ==============================================================================

cat("\n═══════════════════════════════════════════════════════════════════════\n")
cat(" STEP 15: Saving results\n")
cat("═══════════════════════════════════════════════════════════════════════\n")

dir.create(file.path(out_dir, "objects"), recursive = TRUE, showWarnings = FALSE)
saveRDS(gobj, file.path(out_dir, "objects", "gobj_final.rds"))
saveRDS(color_config, file.path(out_dir, "objects", "color_config.rds"))

# Save the canonical session snapshot consumed by every downstream script.
# Object name `gobject` matches what 04, 05, 08, 11, 12a etc. expect via
# `load(SESSION_SNAPSHOT)`. Giotto's stochastic functions (runPCA,
# runGiottoHarmony, runUMAP, doLeidenCluster) all default to
# `set_seed = TRUE, seed_number = 1234`, so this snapshot is deterministic
# given identical inputs.
gobject <- gobj
save(gobject, file = SESSION_SNAPSHOT)
cat(sprintf("Wrote %s\n", SESSION_SNAPSHOT))

# Parameters
params_dt <- data.table(parameter = names(params), value = as.character(unlist(params)))
fwrite(params_dt, file.path(out_dir, "tables", "analysis_parameters.csv"))

# Summary
meta <- pDataDT(gobj)
summary_dt <- data.table(
  metric = c("Total cells", "Total genes", "Samples", "Conditions", "Clusters",
             "Median transcripts/cell", "Median genes/cell"),
  value = c(nrow(meta), length(featIDs(gobj)), length(unique(meta$sample_id)),
            length(unique(meta$condition)), n_clusters,
            round(median(meta$total_expr)), round(median(meta$nr_feats)))
)
fwrite(summary_dt, file.path(out_dir, "tables", "summary_statistics.csv"))

# Cluster composition
cluster_comp <- meta[, .N, by = .(leiden, condition, sample_id)]
fwrite(cluster_comp, file.path(out_dir, "tables", "cluster_composition.csv"))

# Color palettes
color_dt <- data.table(cluster = names(cluster_colors), hex = unname(cluster_colors))
fwrite(color_dt, file.path(out_dir, "tables", "cluster_colors.csv"))

# ==============================================================================
# COMPLETE
# ==============================================================================

cat("\n")
cat("════════════════════════════════════════════════════════════════════════\n")
cat("                         PIPELINE COMPLETE                              \n")
cat("════════════════════════════════════════════════════════════════════════\n\n")

print(summary_dt)

cat(sprintf("\nOutput: %s\n", out_dir))
cat("\nKey files:\n")
cat("  objects/gobj_final.rds       - Final Giotto object\n")
cat("  objects/color_config.rds     - Color palettes (use for all downstream)\n")
cat("  tables/cluster_markers.csv   - Markers with FDR\n")
cat("  tables/cluster_colors.csv    - Hex codes for clusters\n")
cat("  plots/07_HE_*.png            - H&E overlay plots\n")

cat("\n════════════════════════════════════════════════════════════════════════\n")
