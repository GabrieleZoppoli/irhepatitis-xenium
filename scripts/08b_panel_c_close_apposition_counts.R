#!/usr/bin/env Rscript
# ============================================================================
# Figure 2 — CD4 close-apposition to hepatocytes (full tissue sections).
#
# For each sample, computes the fraction of CD4 T cells (cluster 7, passing
# Xenium QC) whose centroid lies within 15 µm of any hepatocyte (clusters 1,
# 6, 10, 18) cell-boundary vertex. Aggregates per condition (irHep =
# 18321_a + 24774_d; AIH = 54023_b + 56784_c).
#
# Used in the Figure 2 results paragraph to support the qualitative
# observation that CD4 T cells are more often hepatocyte-apposed in AIH
# than in irHepatitis. Output table written to:
#   data/processed/tables/cd4_close_apposition.csv
# ============================================================================
suppressPackageStartupMessages({
  library(data.table); library(arrow); library(FNN)
})

source(file.path(Sys.getenv("IRHEP_BUNDLE_ROOT", here::here()), "scripts", "config.R"))
source(file.path(SCRIPTS_DIR, "02_cluster_config.R"))

OUT_DIR <- file.path(PROCESSED_DIR, "tables")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
OUT_CSV <- file.path(OUT_DIR, "cd4_close_apposition.csv")

CD4_CLUSTER <- 7L
HEP_CLUSTERS <- c(1L, 6L, 10L, 18L)
APPOSITION_UM <- 15

ok_cell <- function(ca, na, tc) {
  if (is.na(ca) || is.na(na) || na < 8 || is.na(tc) || tc < 10) return(FALSE)
  TRUE
}

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

all_samples <- unique(c(SAMPLES$irHep, SAMPLES$AIH))
per_sample <- list()
for (sid in all_samples) {
  cat(sprintf("Processing %s...\n", sid))
  a <- as.data.table(read_parquet(file.path(sample_raw_dir(sid), "cells.parquet")))
  setnames(a, "cell_id", "cell_id_raw")
  a[, sample_id := sid]
  a[, ok := mapply(ok_cell, cell_area, nucleus_area, total_counts)]
  cs <- merge(coords[sample_id == sid], a[, .(cell_id_raw, sample_id, ok)],
              by = c("cell_id_raw", "sample_id"))
  cs <- cs[ok == TRUE]

  cd4 <- cs[leiden == CD4_CLUSTER]
  hep_ids <- cs[leiden %in% HEP_CLUSTERS]$cell_id_raw

  if (nrow(cd4) == 0 || length(hep_ids) == 0) {
    per_sample[[sid]] <- data.table(sample_id = sid,
                                    n_cd4 = nrow(cd4),
                                    n_close = NA_integer_,
                                    frac_close = NA_real_)
    next
  }

  cb <- as.data.table(read_parquet(file.path(sample_raw_dir(sid), "cell_boundaries.parquet")))
  hep_verts <- cb[cell_id %in% hep_ids, .(vertex_x, vertex_y)]
  if (nrow(hep_verts) == 0) {
    per_sample[[sid]] <- data.table(sample_id = sid,
                                    n_cd4 = nrow(cd4),
                                    n_close = NA_integer_,
                                    frac_close = NA_real_)
    next
  }

  # CD4 centroids are in Giotto spatial coords (sdimx, sdimy); cell-boundary
  # vertices are in raw Xenium coords (vertex_x, vertex_y). Both per-sample
  # systems are co-located up to global shifts that cancel within a sample.
  # For within-sample nearest-neighbour distance, use raw coordinates on
  # both sides — load CD4 centroids from the cells.parquet (x_centroid /
  # y_centroid columns) to match the boundary frame.
  ca <- as.data.table(read_parquet(file.path(sample_raw_dir(sid), "cells.parquet")))
  cd4_xy <- ca[cell_id %in% cd4$cell_id_raw, .(x = x_centroid, y = y_centroid)]
  if (nrow(cd4_xy) == 0) {
    per_sample[[sid]] <- data.table(sample_id = sid,
                                    n_cd4 = nrow(cd4),
                                    n_close = NA_integer_,
                                    frac_close = NA_real_)
    next
  }

  nn <- FNN::get.knnx(as.matrix(hep_verts[, .(vertex_x, vertex_y)]),
                      as.matrix(cd4_xy), k = 1)
  cd4_xy[, min_dist := nn$nn.dist[, 1]]
  n_close <- sum(cd4_xy$min_dist <= APPOSITION_UM)

  per_sample[[sid]] <- data.table(
    sample_id   = sid,
    n_cd4       = nrow(cd4_xy),
    n_close     = n_close,
    frac_close  = n_close / nrow(cd4_xy)
  )
  cat(sprintf("  %s: %d / %d CD4 cells within %g µm of hepatocyte boundary (%.1f%%)\n",
              sid, n_close, nrow(cd4_xy), APPOSITION_UM,
              100 * n_close / nrow(cd4_xy)))
}

per_sample_dt <- rbindlist(per_sample)
per_sample_dt[, condition := fcase(
  sample_id %in% SAMPLES$irHep, "irHep",
  sample_id %in% SAMPLES$AIH,   "AIH"
)]

per_condition_dt <- per_sample_dt[, .(
  n_cd4_total    = sum(n_cd4, na.rm = TRUE),
  n_close_total  = sum(n_close, na.rm = TRUE),
  frac_close     = sum(n_close, na.rm = TRUE) / sum(n_cd4, na.rm = TRUE)
), by = condition]

cat("\n=== Per-sample ===\n")
print(per_sample_dt)
cat("\n=== Per-condition (pooled) ===\n")
print(per_condition_dt)

fwrite(per_sample_dt, OUT_CSV)
cat(sprintf("\nWrote %s\n", OUT_CSV))
cat("\nMain-text figures: irHep = %s; AIH = %s\n" |>
      sprintf(per_condition_dt[condition == "irHep",
                sprintf("%d / %d = %.1f%%", n_close_total, n_cd4_total,
                        100 * frac_close)],
              per_condition_dt[condition == "AIH",
                sprintf("%d / %d = %.1f%%", n_close_total, n_cd4_total,
                        100 * frac_close)]))
