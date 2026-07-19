#!/usr/bin/env Rscript
# ============================================================================
# Per-sample CellChat run for Figure 2 Panels B (chord) and C (L–R).
# ============================================================================
#
# PURPOSE:
#   Run CellChat v2 once per sample (4 samples), then aggregate per-condition
#   so we can apply the same per-sample coherence rule we use in Panel A.
#   This gives us:
#     - per-sample × per-LR-pair × per-cell-type-pair communication probabilities
#     - per-condition mean probability + sample-level sig coherence
#     - irHep − AIH differential at the L-R level
#
#   The point of going per-sample is to NOT lose sample-level variability and
#   to mirror Panel A's analytical structure (per-sample → per-condition →
#   coherence rule for highlighting).
#
# UNBIASED:
#   - Full CellChatDB.human (all categories)
#   - All cell-type pairs, no pre-specified targets
#   - The choice of which interactions to highlight in Panel B/C is made
#     downstream from the unbiased output, not pre-specified here.
#
# OUTPUTS (in {PROCESSED_DIR}/cellchat/per_sample/):
#   - cellchat_{sample_id}_ss10.rds          (one per sample, the CellChat objects)
#   - cellchat_per_sample_communications.csv (full table, all 4 samples)
#   - cellchat_per_condition_mean.csv        (per-condition mean probs)
#   - cellchat_differential_unbiased.csv     (irHep − AIH per L-R pair × cell-type pair)
#
# DEPENDENCIES: CellChat (>= 2.2), Giotto (>= 4.2), GiottoClass, data.table
# ============================================================================

suppressPackageStartupMessages({
  library(CellChat)
  library(data.table)
  library(Matrix)
  library(Giotto)
  library(GiottoClass)
  library(future)
})


source(file.path(Sys.getenv("IRHEP_BUNDLE_ROOT", here::here()), "scripts", "config.R"))
set.seed(42)

# Parallel plan for CellChat permutation tests (computeCommunProb uses future).
# 8 workers = ~11% of 72-core machine, leaves headroom for concurrent jobs.
# maxSize bumped to 4 GB because CellChat's exported globals (FunMean, data.use)
# exceed the 500 MB default for the larger AIH samples (52k+ cells).
options(future.globals.maxSize = 4 * 1024^3)
plan(multisession, workers = 8)

# --- Config ---
SPOT_SIZE         <- 10     # tol = 5 µm contact tolerance (CellChat v2 recommended)
INTERACTION_RANGE <- 250    # µm — max diffusion range for secreted signaling
SCALE_DISTANCE    <- 6.5    # CellChat v2 default; controls rate of distance decay
CONTACT_RANGE     <- 25     # µm — biologically honest for Cell-Cell Contact
                             # pathways: ~1–2 cell diameters reflects true adjacency.
MIN_CELLS         <- 10

working_dir <- RAW_DIR
out_dir     <- file.path(PROCESSED_DIR, "cellchat", "per_sample")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

source(file.path(SCRIPTS_DIR, "02_cluster_config.R"))

samples_info <- data.table(
  sample_id = c("18321_a", "24774_d", "54023_b", "56784_c"),
  condition = c("irHepatitis", "irHepatitis", "AIH", "AIH")
)

# ============================================================================
# 1. Load Giotto object + extract data once
# ============================================================================
cat("Loading Giotto object...\n")
gobj <- readRDS(file.path(PROCESSED_DIR, "gobj_with_spatial_networks.rds"))

meta <- as.data.table(GiottoClass::pDataDT(gobj))
leiden_col <- grep("leiden", names(meta), value = TRUE)[1]
if (leiden_col != "leiden") meta[, leiden := get(leiden_col)]
meta[, cluster_name := cell_type_labels[as.character(leiden)]]
meta <- meta[!leiden %in% as.integer(excluded_clusters)]

if (!"condition" %in% names(meta)) {
  meta <- merge(meta, samples_info[, .(sample_id, condition)],
                by = "sample_id", all.x = TRUE)
}

# Spatial coordinates
spat_locs <- as.data.table(GiottoClass::getSpatialLocations(gobj, output = "data.table"))
meta <- merge(meta, spat_locs[, .(cell_ID, sdimx, sdimy)],
              by = "cell_ID", all.x = TRUE)

# Expression matrix (genes × cells, normalized)
cat("Extracting normalized expression matrix...\n")
expr_mat <- GiottoClass::getExpression(gobj, values = "normalized", output = "matrix")
expr_mat <- expr_mat[, meta$cell_ID]

cat(sprintf("  %d cells (after exclusion) × %d genes\n",
            nrow(meta), nrow(expr_mat)))
cat(sprintf("  Cell-type breakdown:\n"))
print(meta[, .N, by = cluster_name][order(-N)])

# ============================================================================
# 2. Per-sample CellChat helper
# ============================================================================
run_cellchat_one_sample <- function(sid, cond) {
  cat(sprintf("\n========== %s (%s) ==========\n", sid, cond))

  cells <- meta[sample_id == sid, cell_ID]
  cat(sprintf("  %d cells\n", length(cells)))

  meta_sub <- meta[cell_ID %in% cells]
  data.input <- expr_mat[, cells]

  cellchat_meta <- data.frame(
    labels  = meta_sub$cluster_name,
    samples = meta_sub$sample_id,
    row.names = meta_sub$cell_ID,
    stringsAsFactors = FALSE
  )

  spatial.locs <- as.data.frame(meta_sub[, .(sdimx, sdimy)])
  rownames(spatial.locs) <- meta_sub$cell_ID

  spatial.factors <- data.frame(ratio = 1, tol = SPOT_SIZE / 2)

  cat("  Creating CellChat object...\n")
  cellchat <- createCellChat(
    object          = data.input,
    meta            = cellchat_meta,
    group.by        = "labels",
    datatype        = "spatial",
    coordinates     = spatial.locs,
    spatial.factors = spatial.factors
  )

  cellchat@DB <- CellChatDB.human

  cat("  Pre-processing (subset, OE genes, OE interactions)...\n")
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat,
                                                variable.both = FALSE)

  cat("  Computing communication probabilities (spatial)...\n")
  cellchat <- computeCommunProb(
    cellchat,
    type              = "truncatedMean",
    trim              = 0.1,
    distance.use      = TRUE,
    interaction.range = INTERACTION_RANGE,
    scale.distance    = SCALE_DISTANCE,
    contact.dependent = TRUE,
    contact.range     = CONTACT_RANGE
  )

  cellchat <- filterCommunication(cellchat, min.cells = MIN_CELLS)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)

  comm <- subsetCommunication(cellchat)
  cat(sprintf("  → %d L-R interactions detected\n", nrow(comm)))

  return(cellchat)
}

# ============================================================================
# 3. Run all 4 samples
# ============================================================================
all_comms <- list()

for (i in seq_len(nrow(samples_info))) {
  sid  <- samples_info$sample_id[i]
  cond <- samples_info$condition[i]

  rds_path <- file.path(out_dir, sprintf("cellchat_%s_ss%d.rds", sid, SPOT_SIZE))
  cellchat_obj <- run_cellchat_one_sample(sid, cond)
  saveRDS(cellchat_obj, rds_path)
  cat(sprintf("  Saved: cellchat_%s_ss%d.rds\n", sid, SPOT_SIZE))

  comm <- as.data.table(subsetCommunication(cellchat_obj))
  comm[, sample_id := sid]
  comm[, condition := cond]
  all_comms[[sid]] <- comm
}

# ============================================================================
# 4. Combine per-sample tables
# ============================================================================
all_dt <- rbindlist(all_comms, fill = TRUE)
fwrite(all_dt, file.path(out_dir, "cellchat_per_sample_communications.csv"))
cat(sprintf("\nWrote per-sample communications: %d rows\n", nrow(all_dt)))

# ============================================================================
# 5. Per-condition aggregation (mean prob across samples)
# ============================================================================
cat("\n=== Per-condition aggregation ===\n")

# Identify the right key columns. CellChat returns: source, target, ligand,
# receptor, interaction_name, interaction_name_2, pathway_name, prob, pval
key_cols <- c("source", "target", "ligand", "receptor",
              "interaction_name", "interaction_name_2", "pathway_name")

# Sanity check
missing_keys <- setdiff(key_cols, names(all_dt))
if (length(missing_keys) > 0) {
  cat("WARNING: missing key columns: ",
      paste(missing_keys, collapse = ", "), "\n")
  cat("Available columns: ", paste(names(all_dt), collapse = ", "), "\n")
}

agg_dt <- all_dt[, .(
  mean_prob   = mean(prob, na.rm = TRUE),
  n_samples   = .N,
  n_sig       = sum(pval < 0.05, na.rm = TRUE),
  any_sig     = any(pval < 0.05, na.rm = TRUE),
  all_sig     = all(pval < 0.05, na.rm = TRUE)
), by = c("condition", key_cols)]

fwrite(agg_dt, file.path(out_dir, "cellchat_per_condition_mean.csv"))
cat(sprintf("Wrote per-condition aggregated: %d rows\n", nrow(agg_dt)))

# ============================================================================
# 6. Differential: irHep − AIH (full outer join)
# ============================================================================
cat("\n=== Computing differential (irHep − AIH) ===\n")

irh <- agg_dt[condition == "irHepatitis"]
aih <- agg_dt[condition == "AIH"]
setnames(irh,
         c("mean_prob", "n_samples", "n_sig", "any_sig", "all_sig"),
         c("irHep_mean_prob", "irHep_n_samples",
           "irHep_n_sig", "irHep_any_sig", "irHep_all_sig"))
setnames(aih,
         c("mean_prob", "n_samples", "n_sig", "any_sig", "all_sig"),
         c("AIH_mean_prob", "AIH_n_samples",
           "AIH_n_sig", "AIH_any_sig", "AIH_all_sig"))
irh[, condition := NULL]
aih[, condition := NULL]

diff_dt <- merge(irh, aih, by = key_cols, all = TRUE)
diff_dt[is.na(irHep_mean_prob), irHep_mean_prob := 0]
diff_dt[is.na(AIH_mean_prob),   AIH_mean_prob   := 0]
for (c in c("irHep_n_sig", "AIH_n_sig", "irHep_n_samples", "AIH_n_samples")) {
  diff_dt[is.na(get(c)), (c) := 0L]
}
for (c in c("irHep_any_sig", "AIH_any_sig", "irHep_all_sig", "AIH_all_sig")) {
  diff_dt[is.na(get(c)), (c) := FALSE]
}

# Linear differential and fold-change with epsilon
EPS <- 1e-6
diff_dt[, delta_prob := irHep_mean_prob - AIH_mean_prob]
diff_dt[, log2_fc    := log2((irHep_mean_prob + EPS) / (AIH_mean_prob + EPS))]
diff_dt[, abs_delta  := abs(delta_prob)]

# Coherent rule analogue (Panel A philosophy):
#   "robust in condition X" = at least 1 sample of X has pval<0.05 (any direction)
#                             AND no sample of X is sig in opposite direction.
# For CellChat there isn't a 'direction' the way there is for proximity
# (each L-R pair only has one direction: ligand → receptor). So 'robust in
# condition' simplifies to: any_sig in that condition (i.e. at least 1 sample
# at p<0.05). The user can adjust this threshold downstream.
diff_dt[, robust_count := as.integer(irHep_any_sig) + as.integer(AIH_any_sig)]

# Sort by |delta_prob| (absolute differential)
setorder(diff_dt, -abs_delta)

fwrite(diff_dt, file.path(out_dir, "cellchat_differential_unbiased.csv"))
cat(sprintf("Wrote differential: %d L-R interactions\n", nrow(diff_dt)))

# ============================================================================
# 7. Console summary
# ============================================================================
cat("\n=== Summary ===\n")
cat(sprintf("  Per-sample CellChat objects: %d (in %s)\n",
            nrow(samples_info), out_dir))
cat(sprintf("  Total per-sample interactions: %d\n", nrow(all_dt)))
cat(sprintf("  Per-condition aggregated rows: %d\n", nrow(agg_dt)))
cat(sprintf("  Differential rows (full outer): %d\n", nrow(diff_dt)))

cat("\n--- Top 20 by |delta_prob| ---\n")
top20 <- head(diff_dt, 20)
for (i in seq_len(nrow(top20))) {
  r <- top20[i]
  cat(sprintf("  %-30s -> %-30s  %-22s  Δ=%+.5f  irHep=%.5f  AIH=%.5f\n",
              substr(r$source, 1, 30),
              substr(r$target, 1, 30),
              substr(r$interaction_name_2, 1, 22),
              r$delta_prob,
              r$irHep_mean_prob,
              r$AIH_mean_prob))
}

cat("\n--- Top 10 unique pathways by total |delta_prob| ---\n")
pw_dt <- diff_dt[, .(total_abs_delta = sum(abs_delta, na.rm = TRUE),
                     n_pairs         = .N),
                  by = pathway_name]
setorder(pw_dt, -total_abs_delta)
print(head(pw_dt, 10))

cat("\n================================================================\n")
cat("CellChat per-sample COMPLETE\n")
cat("================================================================\n")
cat(sprintf("Output dir: %s\n", out_dir))
