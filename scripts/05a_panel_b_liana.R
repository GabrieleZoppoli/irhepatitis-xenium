#!/usr/bin/env Rscript
# ============================================================================
# LIANA+ multi-method CCC analysis per sample, cross-validating CellChat v2
# findings. Runs 5 methods (NATMI, Connectome, logFC, SingleCellSignalR,
# CellPhoneDB) and aggregates via robust rank aggregation.
# Output: per-sample liana tables + coherent cross-condition summary.
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(liana)
  library(SingleCellExperiment)
  library(future)
})
set.seed(42)

options(future.globals.maxSize = 8 * 1024^3)
# LIANA is not as parallelized as MISTy — running per-sample sequentially is fine,
# inner parallelism is limited. Keep plan minimal.
plan(multisession, workers = 4)

source(file.path(Sys.getenv("IRHEP_BUNDLE_ROOT", here::here()), "scripts", "config.R"))
source(file.path(SCRIPTS_DIR, "02_cluster_config.R"))

out_dir <- file.path(PROCESSED_DIR, "liana", "per_sample")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================================
# 1. Load Giotto object
# ============================================================================
cat("Loading Giotto object...\n")
gobj <- readRDS(file.path(PROCESSED_DIR, "gobj_with_spatial_networks.rds"))
meta <- as.data.table(GiottoClass::pDataDT(gobj))
leiden_col <- grep("leiden", names(meta), value = TRUE)[1]
if (leiden_col != "leiden") meta[, leiden := get(leiden_col)]
meta[, cluster_name := cell_type_labels[as.character(leiden)]]
meta <- meta[!leiden %in% as.integer(excluded_clusters)]
expr_norm <- GiottoClass::getExpression(gobj, values = "normalized", output = "matrix")
expr_norm <- expr_norm[, meta$cell_ID]
expr_raw  <- GiottoClass::getExpression(gobj, values = "raw", output = "matrix")
expr_raw  <- expr_raw[, meta$cell_ID]
cat(sprintf("Loaded raw (%d x %d) + normalized expression matrices\n",
    nrow(expr_raw), ncol(expr_raw)))

samples_info <- data.table(
  sample_id = c("18321_a", "24774_d", "54023_b", "56784_c"),
  condition = c("irHepatitis", "irHepatitis", "AIH", "AIH")
)

# ============================================================================
# 2. Per-sample LIANA run
# ============================================================================
run_liana_one_sample <- function(sid) {
  cat(sprintf("\n========== %s ==========\n", sid))
  cells <- meta[sample_id == sid, cell_ID]
  cat(sprintf("  %d cells\n", length(cells)))

  # Build SingleCellExperiment with BOTH counts and logcounts (LIANA requires both)
  expr_norm_sub <- expr_norm[, cells]
  expr_raw_sub  <- expr_raw[, cells]
  labels_sub <- meta[cell_ID %in% cells]$cluster_name
  sce <- SingleCellExperiment(
    assays = list(counts = expr_raw_sub, logcounts = expr_norm_sub),
    colData = DataFrame(cell_ID = cells, label = labels_sub, row.names = cells)
  )

  cat("  Running liana_wrap (natmi, connectome, logfc, sca, cellphonedb)...\n")
  t0 <- Sys.time()
  lres <- liana_wrap(
    sce,
    method = c("natmi", "connectome", "logfc", "sca", "cellphonedb"),
    resource = c("Consensus"),
    idents_col = "label",
    min_cells = 10,
    verbose = FALSE
  )
  cat(sprintf("  Completed in %.1f min\n", as.numeric(Sys.time() - t0, units = "mins")))

  cat("  Aggregating via robust rank aggregation...\n")
  lagg <- liana_aggregate(lres, verbose = FALSE)

  return(list(per_method = lres, aggregated = lagg))
}

# ============================================================================
# 3. Run for all samples
# ============================================================================
all_results <- list()
for (i in seq_len(nrow(samples_info))) {
  sid  <- samples_info$sample_id[i]
  cond <- samples_info$condition[i]
  res <- run_liana_one_sample(sid)

  saveRDS(res, file.path(out_dir, sprintf("liana_%s.rds", sid)))
  cat(sprintf("  Saved liana_%s.rds\n", sid))

  # Write the aggregated table with sample + condition annotation
  agg <- as.data.table(res$aggregated)
  agg[, sample_id := sid]
  agg[, condition := cond]
  fwrite(agg, file.path(out_dir, sprintf("liana_agg_%s.csv", sid)))
  all_results[[sid]] <- agg
}

# ============================================================================
# 4. Combine across samples and save master table
# ============================================================================
master <- rbindlist(all_results, fill = TRUE)
fwrite(master, file.path(out_dir, "liana_agg_all_samples.csv"))
cat(sprintf("\nWrote master table: %d rows (4 samples)\n", nrow(master)))

plan(sequential)
cat("\nDone.\n")
