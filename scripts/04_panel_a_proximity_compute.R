#!/usr/bin/env Rscript
# ============================================================================
# Figure 2 Panel A — spatial proximity enrichment.
#
# Computes cell-type spatial proximity enrichment per sample via
# Giotto::cellProximityEnrichment (Dries et al., Genome Biology 2021,
# s13059-021-02286-2). The enrichment statistic is
#   log2((observed + 1) / (expected + 1))
# where "expected" comes from Monte Carlo simulations (random cluster-label
# assignment). cellProximityEnrichment operates on one spatial network at a
# time and permutes labels globally within that network; to preserve
# sample-level composition we subset the Giotto object per sample and run
# the function on each subset. Per-sample enrichments are then averaged
# within condition and differenced (irHep − AIH) for the Panel A
# differential heatmap.
#
# Outputs (in data/processed/spatial_graph_analysis_giotto/):
#   - giotto_proximity_per_sample.csv     (raw cellProximityEnrichment output)
#   - giotto_proximity_per_condition.csv  (condition means)
#   - giotto_proximity_differential.csv   (irHep − AIH)
#   - giotto_proximity_summary.txt        (top hits)
#
# Dependencies: Giotto (>= 4.2), GiottoClass, data.table.
# ============================================================================

suppressPackageStartupMessages({
  library(Giotto)
  library(GiottoClass)
  library(data.table)
})

set.seed(42)

# --- Config ---
N_SIMS <- 1000
source(file.path(Sys.getenv("IRHEP_BUNDLE_ROOT", here::here()), "scripts", "config.R"))
source(file.path(SCRIPTS_DIR, "02_cluster_config.R"))

out_dir <- file.path(PROCESSED_DIR, "spatial_graph_analysis_giotto")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

samples_info <- data.table(
  sample_id = c("18321_a", "24774_d", "54023_b", "56784_c"),
  condition = c("irHepatitis", "irHepatitis", "AIH", "AIH")
)

# ============================================================================
# 1. Load Giotto object; add cluster_name column; identify cells per sample
# ============================================================================
cat("Loading Giotto object...\n")
gobj <- readRDS(file.path(PROCESSED_DIR, "gobj_with_spatial_networks.rds"))

meta <- pDataDT(gobj)
leiden_col <- grep("leiden", names(meta), value = TRUE)[1]
if (leiden_col != "leiden") meta[, leiden := get(leiden_col)]
meta[, cluster_name := cell_type_labels[as.character(leiden)]]

excluded_cell_ids <- meta[leiden %in% as.integer(excluded_clusters), cell_ID]
cat(sprintf("  Total cells: %d\n", nrow(meta)))
cat(sprintf("  Excluding clusters %s: %d cells dropped\n",
            paste(excluded_clusters, collapse = ", "), length(excluded_cell_ids)))

# Add cluster_name to Giotto metadata if not already present
if (!"cluster_name" %in% names(pDataDT(gobj))) {
  gobj <- addCellMetadata(
    gobj,
    new_metadata = data.frame(cell_ID = meta$cell_ID,
                              cluster_name = meta$cluster_name),
    by_column = TRUE,
    column_cell_ID = "cell_ID"
  )
  cat("  Added cluster_name column to Giotto metadata\n")
} else {
  cat("  cluster_name already present in Giotto metadata\n")
}

# ============================================================================
# 2. Run cellProximityEnrichment per sample
# ============================================================================
cat("\n=== Running cellProximityEnrichment per sample ===\n")

per_sample_results <- list()

for (i in seq_len(nrow(samples_info))) {
  sid  <- samples_info$sample_id[i]
  cond <- samples_info$condition[i]
  cat(sprintf("\n--- Sample: %s (%s) ---\n", sid, cond))

  # Cells in this sample, after exclusion
  sample_cells <- meta[sample_id == sid &
                       !leiden %in% as.integer(excluded_clusters), cell_ID]
  cat(sprintf("  Cells after exclusion: %d\n", length(sample_cells)))

  # Subset Giotto object to this sample
  gobj_sample <- tryCatch(
    subsetGiotto(gobj, cell_ids = sample_cells),
    error = function(e) {
      cat("  WARNING: subsetGiotto failed: ", conditionMessage(e), "\n")
      cat("  Trying subset() generic...\n")
      tryCatch(
        subset(gobj, cell_ids = sample_cells),
        error = function(e2) {
          stop(sprintf("Could not subset Giotto object for sample %s: %s",
                       sid, conditionMessage(e2)))
        }
      )
    }
  )

  n_subset <- ncol(gobj_sample)
  cat(sprintf("  Subsetted Giotto object: %d cells\n", n_subset))

  # Run cellProximityEnrichment
  cat(sprintf("  Running cellProximityEnrichment (%d simulations)...\n", N_SIMS))
  cPG <- cellProximityEnrichment(
    gobject               = gobj_sample,
    cluster_column        = "cluster_name",
    spatial_network_name  = "contact_network",
    number_of_simulations = N_SIMS,
    adjust_method         = "BH",
    set_seed              = TRUE,
    seed_number           = 42L + i
  )

  # Inspect + extract results
  cat(sprintf("  Output class: %s\n", paste(class(cPG), collapse = ", ")))
  if (is.list(cPG)) {
    cat(sprintf("  List names: %s\n", paste(names(cPG), collapse = ", ")))
  }

  er <- NULL
  if (is.list(cPG) && "enrichm_res" %in% names(cPG)) {
    er <- as.data.table(cPG$enrichm_res)
  } else if (is.data.frame(cPG)) {
    er <- as.data.table(cPG)
  } else if (is.list(cPG) && length(cPG) == 1 &&
             (is.data.frame(cPG[[1]]) || inherits(cPG[[1]], "data.table"))) {
    er <- as.data.table(cPG[[1]])
  } else {
    cat("  !! Unexpected output structure — inspecting further\n")
    print(str(cPG))
  }

  if (!is.null(er)) {
    er[, sample_id := sid]
    er[, condition := cond]
    per_sample_results[[sid]] <- er
    cat(sprintf("  Returned %d rows, columns: %s\n",
                nrow(er), paste(names(er), collapse = ", ")))
  }
}

# ============================================================================
# 3. Combine per-sample results
# ============================================================================
per_sample_dt <- rbindlist(per_sample_results, fill = TRUE)
cat(sprintf("\nCombined per-sample: %d rows\n", nrow(per_sample_dt)))

# Extract from/to labels if unified_int is the only pair column
if ("unified_int" %in% names(per_sample_dt) &&
    !("from_label" %in% names(per_sample_dt))) {
  # Giotto uses "A--B" format with double hyphen
  per_sample_dt[, c("from_label", "to_label") :=
                  tstrsplit(unified_int, "--", fixed = TRUE)]
  cat("  Extracted from_label/to_label from unified_int\n")
}

fwrite(per_sample_dt, file.path(out_dir, "giotto_proximity_per_sample.csv"))
cat(sprintf("Wrote: giotto_proximity_per_sample.csv\n"))

# ============================================================================
# 4. Aggregate per-condition (mean of per-sample enrichments)
# ============================================================================
cat("\n=== Aggregating per-condition ===\n")

# Identify the enrichment column — Giotto usually calls it 'enrichm'
enrichm_col <- NULL
for (candidate in c("enrichm", "enrichment", "log2_enrichm")) {
  if (candidate %in% names(per_sample_dt)) {
    enrichm_col <- candidate
    break
  }
}

if (is.null(enrichm_col)) {
  cat("  !! Could not identify enrichment column. Available columns:\n")
  print(names(per_sample_dt))
  stop("Cannot continue without enrichment column")
}
cat(sprintf("  Using enrichment column: %s\n", enrichm_col))

# Build aggregation formula using .SDcols
agg_cols <- c(enrichm_col)
if ("original" %in% names(per_sample_dt))     agg_cols <- c(agg_cols, "original")
if ("sim_mean" %in% names(per_sample_dt))     agg_cols <- c(agg_cols, "sim_mean")
if ("p_higher_orig" %in% names(per_sample_dt)) agg_cols <- c(agg_cols, "p_higher_orig")
if ("p_lower_orig" %in% names(per_sample_dt))  agg_cols <- c(agg_cols, "p_lower_orig")
if ("PI_value" %in% names(per_sample_dt))     agg_cols <- c(agg_cols, "PI_value")

per_cond <- per_sample_dt[,
  lapply(.SD, mean, na.rm = TRUE),
  by = .(condition, from_label, to_label),
  .SDcols = agg_cols
]
setnames(per_cond, agg_cols, paste0("mean_", agg_cols))

n_per_pair <- per_sample_dt[, .(n_samples = .N),
                            by = .(condition, from_label, to_label)]
per_cond <- merge(per_cond, n_per_pair,
                  by = c("condition", "from_label", "to_label"))

fwrite(per_cond, file.path(out_dir, "giotto_proximity_per_condition.csv"))
cat(sprintf("  Wrote: giotto_proximity_per_condition.csv (%d rows)\n", nrow(per_cond)))

# ============================================================================
# 5. Differential: irHep mean - AIH mean
# ============================================================================
cat("\n=== Computing differential (irHep - AIH) ===\n")

mean_enrichm_col <- paste0("mean_", enrichm_col)

irh <- per_cond[condition == "irHepatitis",
                .(from_label, to_label,
                  irHep_enrichm = get(mean_enrichm_col))]
aih <- per_cond[condition == "AIH",
                .(from_label, to_label,
                  AIH_enrichm = get(mean_enrichm_col))]

diff_dt <- merge(irh, aih, by = c("from_label", "to_label"), all = TRUE)
diff_dt[is.na(irHep_enrichm), irHep_enrichm := 0]
diff_dt[is.na(AIH_enrichm),   AIH_enrichm := 0]
diff_dt[, delta_enrichm := irHep_enrichm - AIH_enrichm]
setorder(diff_dt, -delta_enrichm)

fwrite(diff_dt, file.path(out_dir, "giotto_proximity_differential.csv"))
cat(sprintf("  Wrote: giotto_proximity_differential.csv (%d rows)\n", nrow(diff_dt)))

# ============================================================================
# 6. Cross-check against the supplementary permutation analysis
# ============================================================================
cat("\n=== Cross-check vs supplementary permutation analysis ===\n")

supp_file <- file.path(PROCESSED_DIR, "spatial_graph_analysis",
                      "proximity_permutation_differential.csv")

if (file.exists(supp_file)) {
  supp <- fread(supp_file)
  supp[, pair_lo := pmin(cluster_from_name, cluster_to_name)]
  supp[, pair_hi := pmax(cluster_from_name, cluster_to_name)]

  diff_dt[, pair_lo := pmin(from_label, to_label)]
  diff_dt[, pair_hi := pmax(from_label, to_label)]

  # Deduplicate diff_dt to unique symmetric pairs (take mean)
  diff_sym <- diff_dt[, .(
    delta_enrichm = mean(delta_enrichm),
    irHep_enrichm = mean(irHep_enrichm),
    AIH_enrichm   = mean(AIH_enrichm)
  ), by = .(pair_lo, pair_hi)]

  merged <- merge(
    diff_sym[, .(pair_lo, pair_hi,
                 giotto_delta = delta_enrichm,
                 giotto_irHep = irHep_enrichm,
                 giotto_AIH   = AIH_enrichm)],
    supp[, .(pair_lo, pair_hi,
            custom_diff = diff,
            custom_irHep = irHep_enrich,
            custom_AIH   = AIH_enrich)],
    by = c("pair_lo", "pair_hi")
  )

  sp_diff  <- cor(merged$giotto_delta, merged$custom_diff,
                  method = "spearman", use = "complete.obs")
  pe_diff  <- cor(merged$giotto_delta, merged$custom_diff,
                  method = "pearson",  use = "complete.obs")
  sp_irh   <- cor(merged$giotto_irHep, merged$custom_irHep,
                  method = "spearman", use = "complete.obs")
  sp_aih   <- cor(merged$giotto_AIH,   merged$custom_AIH,
                  method = "spearman", use = "complete.obs")

  cat(sprintf("  Pairs in common: %d\n", nrow(merged)))
  cat(sprintf("  Spearman (differential): %.3f\n", sp_diff))
  cat(sprintf("  Pearson  (differential): %.3f\n", sp_diff))
  cat(sprintf("  Spearman (irHep enrich): %.3f\n", sp_irh))
  cat(sprintf("  Spearman (AIH enrich):   %.3f\n", sp_aih))

  # Top-10 overlap
  g_top <- head(merged[order(-giotto_delta)], 10)
  c_top <- head(merged[order(-custom_diff)], 10)
  g_top_k <- paste(g_top$pair_lo, g_top$pair_hi, sep = "||")
  c_top_k <- paste(c_top$pair_lo, c_top$pair_hi, sep = "||")
  cat(sprintf("  Top 10 overlap (irHep-enriched): %d / 10\n",
              length(intersect(g_top_k, c_top_k))))

  g_bot <- head(merged[order(giotto_delta)], 10)
  c_bot <- head(merged[order(custom_diff)], 10)
  g_bot_k <- paste(g_bot$pair_lo, g_bot$pair_hi, sep = "||")
  c_bot_k <- paste(c_bot$pair_lo, c_bot$pair_hi, sep = "||")
  cat(sprintf("  Top 10 overlap (AIH-enriched):   %d / 10\n",
              length(intersect(g_bot_k, c_bot_k))))

  sink(file.path(out_dir, "giotto_proximity_validation.txt"))
  cat("================================================================\n")
  cat("Cross-check: Giotto::cellProximityEnrichment vs supplementary permutation analysis\n")
  cat(sprintf("Date: %s\n", Sys.Date()))
  cat("================================================================\n\n")
  cat(sprintf("Pairs in common: %d\n", nrow(merged)))
  cat(sprintf("Spearman (differential): %.4f\n", sp_diff))
  cat(sprintf("Pearson  (differential): %.4f\n", pe_diff))
  cat(sprintf("Spearman (irHep enrich): %.4f\n", sp_irh))
  cat(sprintf("Spearman (AIH enrich):   %.4f\n", sp_aih))
  cat("\n--- Top 15 Giotto irHep-enriched side-by-side ---\n")
  show <- head(merged[order(-giotto_delta)], 15)
  cat(sprintf("  %-33s <-> %-33s  giotto  custom\n", "from", "to"))
  for (i in seq_len(nrow(show))) {
    r <- show[i]
    cat(sprintf("  %-33s <-> %-33s  %+6.3f  %+6.3f\n",
                substr(r$pair_lo, 1, 33), substr(r$pair_hi, 1, 33),
                r$giotto_delta, r$custom_diff))
  }
  sink()
  cat("  Wrote: giotto_proximity_validation.txt\n")
} else {
  cat(sprintf("  Custom pipeline output not found at %s\n", supp_file))
}

# ============================================================================
# 7. Summary
# ============================================================================
sink(file.path(out_dir, "giotto_proximity_summary.txt"))
cat("================================================================\n")
cat("SPATIAL PROXIMITY via Giotto::cellProximityEnrichment\n")
cat(sprintf("Date: %s\n", Sys.Date()))
cat(sprintf("Tool: Giotto v%s\n", as.character(packageVersion("Giotto"))))
cat("Reference: Dries et al., Genome Biology 2021 (s13059-021-02286-2)\n")
cat(sprintf("Simulations: %d per sample\n", N_SIMS))
cat("Strategy: subset Giotto per sample -> run cellProximityEnrichment -> aggregate\n")
cat("Spatial network: contact_network\n")
cat(sprintf("Excluded clusters: %s\n", paste(excluded_clusters, collapse = ", ")))
cat("================================================================\n\n")

cat("--- Top 20 irHep-enriched pairs (delta_enrichm) ---\n")
for (i in 1:min(20, nrow(diff_dt))) {
  r <- diff_dt[i]
  cat(sprintf("  %-33s <-> %-33s  delta=%+6.3f  irHep=%+6.3f  AIH=%+6.3f\n",
              substr(r$from_label, 1, 33), substr(r$to_label, 1, 33),
              r$delta_enrichm, r$irHep_enrichm, r$AIH_enrichm))
}
cat("\n--- Top 20 AIH-enriched pairs (delta_enrichm) ---\n")
aih_top <- head(diff_dt[order(delta_enrichm)], 20)
for (i in seq_len(nrow(aih_top))) {
  r <- aih_top[i]
  cat(sprintf("  %-33s <-> %-33s  delta=%+6.3f  irHep=%+6.3f  AIH=%+6.3f\n",
              substr(r$from_label, 1, 33), substr(r$to_label, 1, 33),
              r$delta_enrichm, r$irHep_enrichm, r$AIH_enrichm))
}
sink()
cat("  Wrote: giotto_proximity_summary.txt\n")

cat("\n================================================================\n")
cat("Giotto proximity analysis COMPLETE\n")
cat("================================================================\n")
cat(sprintf("Output dir: %s\n", out_dir))
