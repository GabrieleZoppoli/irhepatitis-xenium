#!/usr/bin/env Rscript
# ============================================================================
# Aggregate per-sample COMMOT output and compute condition-level differential
# flow for each (L-R, source_ct, target_ct) triple.
#
#   1. Reads per-sample commot_results.csv (one per sample).
#   2. Separates:
#        (a) individual L-R pairs (e.g. "CCL5-CCR1")
#        (b) pathway-level aggregates (e.g. "CCL", "MHC-I") — reported separately
#        (c) "total-total" summary row — DROPPED (not meaningful for differential).
#   3. Tracks THREE layers of missingness:
#        - "pair_tested_in_sample": the L-R pair passed COMMOT's min_cell_pct=0.05
#          filter in that sample (presence in the per-sample CSV anywhere).
#        - "cluster_present_in_sample": source or target cluster had >= 10 cells
#          (COMMOT's hard-coded skip). Below threshold => triple silently dropped
#          in the per-sample CSV; must be NA, not zero.
#          Example: "Ig+ hepatocyte" has only 2 cells in 24774_d (<10), so all
#          triples involving it are absent in that sample's CSV.
#        - "triple_absent_given_tested_and_present": if pair was tested AND both
#          endpoints have >=10 cells, then missing triple = true zero flow.
#   4. Computes coherence under two regimes:
#        - STRICT (primary, publishable): L-R pair must be tested in ALL 4
#          samples; then apply min(enriched) > max(other) rule on mean_flow.
#        - PRESENCE (exploratory): L-R pair tested in both samples of one
#          condition but in 0/2 of the other. These are reported with a
#          "presence-driven" flag.
#
# R-only; no Python dependencies post-COMMOT. Conda env: xenium_env.
# Rationale for R: user preference + consistent with downstream plotting and
# cross-reference with CellChat outputs (cellchat_coherent_interactions.csv).
#
# Output (under data/processed/commot/):
#   - commot_all_samples.csv              — concatenation, with pair/pathway flag
#   - commot_per_condition_lr.csv         — wide per-triple, L-R pairs only
#   - commot_per_condition_pathway.csv    — wide per-triple, pathway aggregates
#   - commot_differential_strict.csv      — coherent L-R pairs tested in all 4
#   - commot_differential_presence.csv    — presence/absence L-R pairs
#   - commot_differential_pathway.csv     — coherent pathway-level aggregates
#   - commot_aggregate_summary.txt        — counts + top findings
#
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

source(file.path(Sys.getenv("IRHEP_BUNDLE_ROOT", here::here()), "scripts", "config.R"))
COMMOT_DIR  <- file.path(PROCESSED_DIR, "commot")
PER_SAMPLE  <- file.path(COMMOT_DIR, "per_sample")
CELLS_INPUT <- file.path(COMMOT_DIR, "cells_input")
OUT_DIR     <- COMMOT_DIR
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

SAMPLES <- data.table(
  sample_id = c("18321_a", "24774_d", "54023_b", "56784_c"),
  condition = c("irHepatitis", "irHepatitis", "AIH", "AIH")
)

ir_samples  <- SAMPLES[condition == "irHepatitis", sample_id]
aih_samples <- SAMPLES[condition == "AIH",         sample_id]

MIN_CELL_PER_CLUSTER <- 10  # Must match COMMOT per-sample script (src_idx/tgt_idx skip)

cluster_counts <- rbindlist(lapply(SAMPLES$sample_id, function(sid) {
  cfile <- file.path(CELLS_INPUT, paste0(sid, ".csv"))
  if (!file.exists(cfile)) stop(sprintf("Missing cells.csv for %s", sid))
  cells <- fread(cfile)
  cc <- cells[, .(n_cells = .N), by = cluster_name]
  cc[, sample_id := sid]
  cc
}))
# Wide form: cluster × sample
cluster_wide <- dcast(cluster_counts, cluster_name ~ sample_id,
                      value.var = "n_cells", fill = 0)
for (sid in SAMPLES$sample_id) {
  if (!sid %in% names(cluster_wide)) cluster_wide[, (sid) := 0]
}
setnames(cluster_wide, old = SAMPLES$sample_id,
         new = paste0("ncell_", SAMPLES$sample_id))

# Boolean: cluster_present_in_sample (>= MIN_CELL_PER_CLUSTER)
for (sid in SAMPLES$sample_id) {
  ncell_col   <- paste0("ncell_",   sid)
  present_col <- paste0("present_", sid)
  cluster_wide[, (present_col) := get(ncell_col) >= MIN_CELL_PER_CLUSTER]
}

cat("\nClusters with <", MIN_CELL_PER_CLUSTER, "cells in any sample (silently dropped):\n", sep = "")
flagged <- cluster_wide[
  apply(cluster_wide[, paste0("present_", SAMPLES$sample_id), with = FALSE], 1, function(x) !all(x))
]
print(flagged)

per_sample <- rbindlist(lapply(SAMPLES$sample_id, function(sid) {
  f <- file.path(PER_SAMPLE, sid, "commot_results.csv")
  if (!file.exists(f)) stop(sprintf("Missing: %s", f))
  dt <- fread(f)
  if (!"sample_id" %in% names(dt)) dt[, sample_id := sid]
  dt
}), fill = TRUE)

# Drop total-total summary rows — they are not meaningful for differential analysis.
per_sample <- per_sample[lr_id != "total-total"]

# Classify: L-R pair (ligand-receptor notation contains "-") vs pathway (no "-").
# COMMOT naming convention: "CCL5-CCR1" = pair; "CCL" = pathway aggregate.
# Heteromeric receptors use "_" (e.g. "IL15-IL15RA_IL2RB_IL2RG") — still a pair.
per_sample[, kind := fifelse(grepl("-", lr_id), "pair", "pathway")]

fwrite(per_sample, file.path(OUT_DIR, "commot_all_samples.csv"))
cat(sprintf("Loaded %d rows from %d samples (after dropping total-total)\n",
            nrow(per_sample), uniqueN(per_sample$sample_id)))
cat(sprintf("  L-R pair rows:    %d\n", sum(per_sample$kind == "pair")))
cat(sprintf("  Pathway rows:     %d\n", sum(per_sample$kind == "pathway")))

tested <- unique(per_sample[, .(lr_id, sample_id, kind)])
tested[, pair_tested := TRUE]

tested_wide <- dcast(tested, lr_id + kind ~ sample_id,
                     value.var = "pair_tested", fill = FALSE)
# Ensure all 4 sample columns exist
for (sid in SAMPLES$sample_id) {
  if (!sid %in% names(tested_wide)) tested_wide[, (sid) := FALSE]
}
setnames(tested_wide,
         old = SAMPLES$sample_id,
         new = paste0("tested_", SAMPLES$sample_id))

tested_wide[, n_ir_tested  := rowSums(.SD),
            .SDcols = paste0("tested_", ir_samples)]
tested_wide[, n_aih_tested := rowSums(.SD),
            .SDcols = paste0("tested_", aih_samples)]
tested_wide[, comparability := fcase(
  n_ir_tested == 2 & n_aih_tested == 2, "all4",
  n_ir_tested == 2 & n_aih_tested == 0, "ir_only",
  n_ir_tested == 0 & n_aih_tested == 2, "aih_only",
  default = "partial"
)]

cat(sprintf("\nComparability of L-R pairs across samples:\n"))
print(tested_wide[, .N, by = .(kind, comparability)])

pivot_wide <- function(dt_sub) {
  wide <- dcast(dt_sub,
                lr_id + source + target ~ sample_id,
                value.var = "mean_flow", fill = 0)
  for (sid in SAMPLES$sample_id) {
    if (!sid %in% names(wide)) wide[, (sid) := 0]
  }
  # Attach pair-level comparability info
  wide <- merge(wide, tested_wide[, .(lr_id, comparability,
                                      n_ir_tested, n_aih_tested)],
                by = "lr_id", all.x = TRUE)
  # For sample columns where pair was NOT tested, rewrite 0 → NA.
  # (When pair_tested_in_sample is FALSE, structural missing — not real zero.)
  for (sid in SAMPLES$sample_id) {
    col_name    <- sid
    tested_col  <- paste0("tested_", sid)
    if (tested_col %in% names(tested_wide)) {
      not_tested <- wide$lr_id %in% tested_wide[get(tested_col) == FALSE, lr_id]
      wide[not_tested, (col_name) := NA_real_]
    }
  }
  # Cluster-level missingness: if source OR target cluster has <MIN_CELL_PER_CLUSTER
  # in that sample, COMMOT silently dropped it — mask as NA.
  # Build lookup: (cluster, sample) → present?
  for (sid in SAMPLES$sample_id) {
    present_col <- paste0("present_", sid)
    if (!present_col %in% names(cluster_wide)) next
    absent_clusters <- cluster_wide[get(present_col) == FALSE, cluster_name]
    if (length(absent_clusters) > 0) {
      mask <- wide$source %in% absent_clusters | wide$target %in% absent_clusters
      wide[mask, (sid) := NA_real_]
    }
  }
  wide
}

wide_pair    <- pivot_wide(per_sample[kind == "pair"])
wide_pathway <- pivot_wide(per_sample[kind == "pathway"])

compute_coherence <- function(wide) {
  # Count non-NA samples per triple per side.
  wide[, n_ir_notNA := rowSums(!is.na(.SD)),
       .SDcols = ir_samples]
  wide[, n_aih_notNA := rowSums(!is.na(.SD)),
       .SDcols = aih_samples]
  # mean_flow per condition (ignores NA)
  wide[, irHep_mean := rowMeans(.SD, na.rm = TRUE),
       .SDcols = ir_samples]
  wide[, AIH_mean   := rowMeans(.SD, na.rm = TRUE),
       .SDcols = aih_samples]
  # min/max per condition, na.rm=FALSE: any NA propagates, so we only trust
  # these values when n_*_notNA == 2.
  wide[, irHep_min := apply(.SD, 1, function(x) min(x)),  .SDcols = ir_samples]
  wide[, irHep_max := apply(.SD, 1, function(x) max(x)),  .SDcols = ir_samples]
  wide[, AIH_min   := apply(.SD, 1, function(x) min(x)),  .SDcols = aih_samples]
  wide[, AIH_max   := apply(.SD, 1, function(x) max(x)),  .SDcols = aih_samples]
  wide
}

wide_pair    <- compute_coherence(wide_pair)
wide_pathway <- compute_coherence(wide_pathway)

# Pseudocount: 1% of global median non-zero mean_flow, across pairs only
nz_pair <- per_sample[kind == "pair" & mean_flow > 0, mean_flow]
PSEUDO_PAIR <- if (length(nz_pair) > 0) 0.01 * median(nz_pair) else 1e-9
nz_pw <- per_sample[kind == "pathway" & mean_flow > 0, mean_flow]
PSEUDO_PW <- if (length(nz_pw) > 0) 0.01 * median(nz_pw) else 1e-9
cat(sprintf("\nPseudocount (L-R pairs):   %.3e\n", PSEUDO_PAIR))
cat(sprintf("Pseudocount (pathway agg): %.3e\n", PSEUDO_PW))

wide_pair[, log2FC := log2((irHep_mean + PSEUDO_PAIR) /
                           (AIH_mean   + PSEUDO_PAIR))]
wide_pathway[, log2FC := log2((irHep_mean + PSEUDO_PW) /
                              (AIH_mean   + PSEUDO_PW))]

# STRICT coherence: only for pairs tested in all 4 samples.
# min_enriched > max_other, both sides must have all 4 NA-free values.
strict_coherent <- function(wide) {
  # Strict eligibility requires:
  #   (a) L-R pair tested in all 4 samples (comparability == "all4")
  #   (b) Both endpoints have n_cells >= MIN_CELL_PER_CLUSTER in all 4 samples
  #       → n_ir_notNA == 2 AND n_aih_notNA == 2 (no NA masks from step 3)
  wide[, is_strict_eligible := comparability == "all4" &
                               n_ir_notNA  == 2 &
                               n_aih_notNA == 2]
  wide[, coherent_strict_irHep := is_strict_eligible & irHep_min > AIH_max]
  wide[, coherent_strict_AIH   := is_strict_eligible & AIH_min   > irHep_max]
  wide[, coherent_strict := coherent_strict_irHep | coherent_strict_AIH]
  wide[, enriched_strict := fifelse(coherent_strict_irHep, "irHepatitis",
                              fifelse(coherent_strict_AIH, "AIH",
                                      NA_character_))]
  wide
}

wide_pair    <- strict_coherent(wide_pair)
wide_pathway <- strict_coherent(wide_pathway)

# PRESENCE coherence: L-R pair tested in both samples of one condition,
# but in 0/2 of the other. Only applicable to pair-level (not pathway).
# ir_only → presence-irHep; aih_only → presence-AIH.
wide_pair[, enriched_presence := fcase(
  comparability == "ir_only"  & irHep_min > 0, "irHepatitis",
  comparability == "aih_only" & AIH_min   > 0, "AIH",
  default = NA_character_
)]

# Save wide tables (per-triple)
fwrite(wide_pair,    file.path(OUT_DIR, "commot_per_condition_lr.csv"))
fwrite(wide_pathway, file.path(OUT_DIR, "commot_per_condition_pathway.csv"))

diff_strict  <- wide_pair[coherent_strict == TRUE][order(-abs(log2FC))]
diff_presc   <- wide_pair[!is.na(enriched_presence)][order(-abs(log2FC))]
diff_pw      <- wide_pathway[coherent_strict == TRUE][order(-abs(log2FC))]

fwrite(diff_strict, file.path(OUT_DIR, "commot_differential_strict.csv"))
fwrite(diff_presc,  file.path(OUT_DIR, "commot_differential_presence.csv"))
fwrite(diff_pw,     file.path(OUT_DIR, "commot_differential_pathway.csv"))

sink(file.path(OUT_DIR, "commot_aggregate_summary.txt"))

cat("================================================================\n")
cat("  COMMOT aggregate summary\n")
cat("  Generated:", format(Sys.time()), "\n")
cat("================================================================\n\n")

cat("Samples:\n")
print(SAMPLES)
cat("\n")

cat("L-R pair tested counts (from comparability table):\n")
print(tested_wide[kind == "pair",
                  .N,
                  keyby = .(n_ir_tested, n_aih_tested, comparability)])
cat("\n")

cat(sprintf("Pairs tested in all 4 samples (strict-eligible): %d\n",
            uniqueN(tested_wide[kind == "pair" & comparability == "all4", lr_id])))
cat(sprintf("Pairs tested in both AIH but 0 irHep:          %d\n",
            uniqueN(tested_wide[kind == "pair" & comparability == "aih_only", lr_id])))
cat(sprintf("Pairs tested in both irHep but 0 AIH:          %d\n",
            uniqueN(tested_wide[kind == "pair" & comparability == "ir_only", lr_id])))
cat("\n")

cat(sprintf("STRICT differential (L-R pair × triple, all 4 tested):\n"))
cat(sprintf("  Coherent total:       %d\n", nrow(diff_strict)))
cat(sprintf("  irHepatitis-enriched: %d\n", sum(diff_strict$enriched_strict == "irHepatitis")))
cat(sprintf("  AIH-enriched:         %d\n", sum(diff_strict$enriched_strict == "AIH")))
cat("\n")
cat(sprintf("PRESENCE differential (L-R pair only tested one side):\n"))
cat(sprintf("  Coherent total:       %d\n", nrow(diff_presc)))
cat(sprintf("  irHepatitis-present:  %d\n", sum(diff_presc$enriched_presence == "irHepatitis")))
cat(sprintf("  AIH-present:          %d\n", sum(diff_presc$enriched_presence == "AIH")))
cat("\n")
cat(sprintf("PATHWAY differential (pathway aggregates, all 4 tested):\n"))
cat(sprintf("  Coherent total:       %d\n", nrow(diff_pw)))
cat(sprintf("  irHepatitis-enriched: %d\n", sum(diff_pw$enriched_strict == "irHepatitis")))
cat(sprintf("  AIH-enriched:         %d\n", sum(diff_pw$enriched_strict == "AIH")))

print_top <- function(dt, side, side_col, n = 20) {
  cat(sprintf("\nTop %d %s-enriched (by |log2FC|):\n", n, side))
  show <- head(dt[get(side_col) == side,
                  .(lr_id, source, target, log2FC, irHep_mean, AIH_mean,
                    comparability)], n)
  print(show)
}

print_top(diff_strict, "irHepatitis", "enriched_strict", 20)
print_top(diff_strict, "AIH",          "enriched_strict", 20)
print_top(diff_pw,     "irHepatitis", "enriched_strict", 10)
print_top(diff_pw,     "AIH",          "enriched_strict", 10)

cat(sprintf("\nPresence-driven irHepatitis pairs (tested only irHep side, top 15):\n"))
print(head(diff_presc[enriched_presence == "irHepatitis",
                      .(lr_id, source, target, irHep_mean, AIH_mean)], 15))
cat(sprintf("\nPresence-driven AIH pairs (tested only AIH side, top 15):\n"))
print(head(diff_presc[enriched_presence == "AIH",
                      .(lr_id, source, target, irHep_mean, AIH_mean)], 15))

sink()

cat("\n=== Output files ===\n")
cat(sprintf("  %s\n", file.path(OUT_DIR, "commot_all_samples.csv")))
cat(sprintf("  %s\n", file.path(OUT_DIR, "commot_per_condition_lr.csv")))
cat(sprintf("  %s\n", file.path(OUT_DIR, "commot_per_condition_pathway.csv")))
cat(sprintf("  %s\n", file.path(OUT_DIR, "commot_differential_strict.csv")))
cat(sprintf("  %s\n", file.path(OUT_DIR, "commot_differential_presence.csv")))
cat(sprintf("  %s\n", file.path(OUT_DIR, "commot_differential_pathway.csv")))
cat(sprintf("  %s\n", file.path(OUT_DIR, "commot_aggregate_summary.txt")))

cat("\nDone.\n")
