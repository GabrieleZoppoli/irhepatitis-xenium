#!/usr/bin/env Rscript
# ============================================================================
# 01a_add_spatial_networks.R
# Takes the post-clustering session_snapshot.RData and adds two per-sample
# spatial networks (computed within each sample so no cross-sample edges
# leak):
#   - contact_network: Delaunay graph with maximum_distance = 25 µm,
#     minimum_k = 2; used by `04_panel_a_proximity_compute.R` and
#     `05_panel_b_cellchat.R` (CellChat's contact-range neighbourhood).
#   - local_network:   kNN graph (k = 6) with maximum_distance = 40 µm,
#     minimum_k = 0; used by spatial-aware downstream analyses (LIANA+,
#     supplementary morphology).
#
# Output: data/processed/gobj_with_spatial_networks.rds
# Runtime: ~5 minutes on a workstation.
# ============================================================================
suppressPackageStartupMessages({
  library(Giotto); library(GiottoClass); library(data.table)
})

source(file.path(Sys.getenv("IRHEP_BUNDLE_ROOT", here::here()), "scripts", "config.R"))

CONTACT_MAX_UM <- 25
CONTACT_MIN_K  <- 2
LOCAL_K        <- 6
LOCAL_MAX_UM   <- 40
LOCAL_MIN_K    <- 0
OUT_RDS        <- file.path(PROCESSED_DIR, "gobj_with_spatial_networks.rds")

cat("Loading session_snapshot.RData ...\n")
load(SESSION_SNAPSHOT)
stopifnot(exists("gobject"))
cat(sprintf("  Loaded gobject with %d cells, %d features\n",
            ncol(gobject), nrow(gobject)))

meta <- pDataDT(gobject)
all_samples <- sort(unique(meta$sample_id))
cat(sprintf("  Samples: %s\n", paste(all_samples, collapse = ", ")))

# Build per-sample networks and concatenate the edge tables. We rebuild
# the global @spatial_network slot manually rather than letting Giotto
# create networks on the combined object — Delaunay across a multi-sample
# gobject would otherwise produce nonsensical cross-sample edges.

build_sample_network <- function(net_name, method, params, sid) {
  cat(sprintf("  [%s] %s ...\n", sid, net_name))
  cell_ids <- meta[sample_id == sid, cell_ID]
  g_sub <- subsetGiotto(gobject, cell_ids = cell_ids)
  g_sub <- do.call(createSpatialNetwork,
                   c(list(gobject = g_sub, method = method, name = net_name,
                          verbose = FALSE), params))
  obj <- g_sub@spatial_network$cell[[net_name]]
  list(networkDT = as.data.table(obj@networkDT),
       networkDT_before_filter = tryCatch(
         as.data.table(obj@networkDT_before_filter),
         error = function(e) NULL))
}

merge_per_sample_networks <- function(net_name, method, params) {
  cat(sprintf("\n=== %s ===\n", net_name))
  parts <- lapply(all_samples, function(sid)
    build_sample_network(net_name, method, params, sid))
  combined_dt   <- rbindlist(lapply(parts, `[[`, "networkDT"))
  combined_pre  <- tryCatch(
    rbindlist(Filter(Negate(is.null), lapply(parts, `[[`, "networkDT_before_filter"))),
    error = function(e) NULL)
  cat(sprintf("  Total edges: %d (post-filter)\n", nrow(combined_dt)))

  obj <- new("spatialNetworkObj",
             name = net_name, method = method, parameters = params,
             networkDT = combined_dt,
             networkDT_before_filter = combined_pre,
             spat_unit = "cell", provenance = "cell")
  obj
}

contact <- merge_per_sample_networks(
  "contact_network", "Delaunay",
  list(maximum_distance_delaunay = CONTACT_MAX_UM,
       minimum_k = CONTACT_MIN_K))

local <- merge_per_sample_networks(
  "local_network", "kNN",
  list(k = LOCAL_K,
       maximum_distance_knn = LOCAL_MAX_UM,
       minimum_k = LOCAL_MIN_K))

# Attach to the combined gobject under @spatial_network$cell
if (is.null(gobject@spatial_network)) gobject@spatial_network <- list()
if (is.null(gobject@spatial_network$cell)) gobject@spatial_network$cell <- list()
gobject@spatial_network$cell$contact_network <- contact
gobject@spatial_network$cell$local_network   <- local

saveRDS(gobject, OUT_RDS)
cat(sprintf("\nWrote %s\n", OUT_RDS))
