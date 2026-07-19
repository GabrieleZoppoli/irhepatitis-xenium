#!/usr/bin/env Rscript
# ============================================================================
# Supplementary permutation analysis for spatial proximity enrichment.
#
# Computes empirical p-values for cell–cell proximity enrichment by permuting
# cluster labels while keeping the spatial network fixed. Tests both overall
# (all samples) and per-condition enrichment.
#
# Method:
#   1. Extract spatial network edges from Giotto object.
#   2. For N_PERM iterations, shuffle cluster labels (within each sample to
#      preserve sample-level composition), recount edges per pair.
#   3. Enrichment = log2((observed + 1) / (expected + 1)).
#   4. Empirical p-value = fraction of permutations with
#      |enrichment| >= observed.
#   5. FDR correction (Benjamini–Hochberg) across all pairs.
#
# Outputs (in data/processed/spatial_graph_analysis/):
#   - proximity_permutation_results.csv         (overall)
#   - proximity_permutation_per_condition.csv   (irHep, AIH separately)
#   - proximity_permutation_differential.csv
#   - proximity_permutation_summary.txt         (human-readable summary)
# ============================================================================

suppressPackageStartupMessages({
  library(Giotto); library(GiottoClass); library(data.table)
})


source(file.path(Sys.getenv("IRHEP_BUNDLE_ROOT", here::here()), "scripts", "config.R"))
set.seed(42)

# --- Config ---
N_PERM <- 1000
out_dir <- file.path(PROCESSED_DIR, "spatial_graph_analysis")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# --- Load cluster config ---
source(file.path(SCRIPTS_DIR, "02_cluster_config.R"))

samples_info <- data.table(
  sample_id = c("18321_a", "24774_d", "54023_b", "56784_c"),
  condition = c("irHepatitis", "irHepatitis", "AIH", "AIH")
)

# --- Load Giotto object + spatial network ---
cat("Loading Giotto object...\n")
gobj <- readRDS(file.path(PROCESSED_DIR, "gobj_with_spatial_networks.rds"))

meta <- pDataDT(gobj)
leiden_col <- grep("leiden", names(meta), value = TRUE)[1]
if (leiden_col != "leiden") meta[, leiden := get(leiden_col)]
meta[, cluster_name := cell_type_labels[as.character(leiden)]]
meta <- meta[!leiden %in% as.integer(excluded_clusters)]
if (!"condition" %in% names(meta)) {
  meta <- merge(meta, samples_info[, .(sample_id, condition)],
                by = "sample_id", all.x = TRUE)
}
cat(sprintf("  Conditions: %s\n", paste(unique(meta$condition), collapse = ", ")))

cat("Extracting spatial network...\n")
spat_net <- getSpatialNetwork(gobj, name = "contact_network", output = "networkDT")
valid_cells <- meta$cell_ID
spat_net <- spat_net[from %in% valid_cells & to %in% valid_cells]

cat(sprintf("  %s cells, %s edges\n", format(nrow(meta), big.mark = ","),
            format(nrow(spat_net), big.mark = ",")))

cell_cluster <- meta[, .(cell_ID, leiden, sample_id, condition)]
setkey(cell_cluster, cell_ID)

# --- Helper: count edges per cluster pair from a label vector ---
count_edges <- function(net, labels) {
  net_cl <- data.table(
    cl_from = labels[net$from],
    cl_to   = labels[net$to]
  )
  net_cl[, `:=`(
    c1 = pmin(cl_from, cl_to),
    c2 = pmax(cl_from, cl_to)
  )]
  counts <- net_cl[, .N, by = .(c1, c2)]
  counts
}

# --- Helper: compute enrichment matrix from edge counts ---
compute_enrichment <- function(edge_counts, cluster_sizes, total_edges) {
  total_cells <- sum(cluster_sizes)
  clusters <- sort(as.integer(names(cluster_sizes)))
  n_cl <- length(clusters)

  pairs <- CJ(c1 = clusters, c2 = clusters)
  pairs <- pairs[c1 <= c2]
  pairs[, key_str := paste(c1, c2, sep = "_")]

  edge_counts[, key_str := paste(c1, c2, sep = "_")]
  pairs <- merge(pairs, edge_counts[, .(key_str, N)],
                 by = "key_str", all.x = TRUE)
  pairs[is.na(N), N := 0]

  pairs[, `:=`(
    n1 = cluster_sizes[as.character(c1)],
    n2 = cluster_sizes[as.character(c2)]
  )]
  pairs[, expected := ifelse(c1 == c2,
    (n1 / total_cells) * (n2 / total_cells) * total_edges,
    (n1 / total_cells) * (n2 / total_cells) * total_edges * 2
  )]
  pairs[, enrichment := log2((N + 1) / (expected + 1))]
  pairs
}

# ===========================================================================
# 1. OVERALL PROXIMITY (all samples combined)
# ===========================================================================
cat("\n=== Overall proximity permutation test ===\n")

labels_obs <- setNames(cell_cluster$leiden, cell_cluster$cell_ID)
cluster_sizes <- table(labels_obs)
total_edges <- nrow(spat_net)

obs_counts <- count_edges(spat_net, labels_obs)
obs_enrich <- compute_enrichment(obs_counts, cluster_sizes, total_edges)

cat(sprintf("Running %d permutations...\n", N_PERM))
n_pairs <- nrow(obs_enrich)
perm_extreme <- rep(0L, n_pairs)

cells_by_sample <- split(cell_cluster$cell_ID, cell_cluster$sample_id)
labels_by_sample <- split(cell_cluster$leiden, cell_cluster$sample_id)

pb_interval <- max(1, N_PERM %/% 20)

for (p in seq_len(N_PERM)) {
  if (p %% pb_interval == 0) cat(sprintf("  permutation %d/%d\n", p, N_PERM))

  perm_labels <- labels_obs
  for (sid in names(cells_by_sample)) {
    cids <- cells_by_sample[[sid]]
    perm_labels[cids] <- sample(labels_by_sample[[sid]])
  }

  perm_counts <- count_edges(spat_net, perm_labels)
  perm_enrich <- compute_enrichment(perm_counts, cluster_sizes, total_edges)

  perm_vals <- perm_enrich$enrichment[match(obs_enrich$key_str,
                                             perm_enrich$key_str)]
  perm_extreme <- perm_extreme + (abs(perm_vals) >= abs(obs_enrich$enrichment))
}

obs_enrich[, p_value := (perm_extreme + 1) / (N_PERM + 1)]
obs_enrich[, fdr := p.adjust(p_value, method = "BH")]
obs_enrich[, cluster_from_name := cell_type_labels[as.character(c1)]]
obs_enrich[, cluster_to_name := cell_type_labels[as.character(c2)]]

fwrite(obs_enrich[, .(c1, c2, cluster_from_name, cluster_to_name,
                       observed = N, expected = round(expected, 1),
                       enrichment = round(enrichment, 4),
                       p_value = round(p_value, 5), fdr = round(fdr, 5))],
       file.path(out_dir, "proximity_permutation_results.csv"))

n_sig <- sum(obs_enrich$fdr < 0.05)
n_sig_pos <- sum(obs_enrich$fdr < 0.05 & obs_enrich$enrichment > 0)
n_sig_neg <- sum(obs_enrich$fdr < 0.05 & obs_enrich$enrichment < 0)
cat(sprintf("\nOverall results: %d pairs tested, %d significant (FDR < 0.05)\n",
            n_pairs, n_sig))
cat(sprintf("  %d co-localized (positive), %d avoiding (negative)\n",
            n_sig_pos, n_sig_neg))


# ===========================================================================
# 2. PER-CONDITION PROXIMITY
# ===========================================================================
cat("\n=== Per-condition proximity permutation test ===\n")

results_cond <- list()

for (cond in c("irHepatitis", "AIH")) {
  cat(sprintf("\n--- %s ---\n", cond))

  cond_cells <- cell_cluster[condition == cond]
  cond_cell_ids <- cond_cells$cell_ID

  net_cond <- spat_net[from %in% cond_cell_ids & to %in% cond_cell_ids]

  labels_cond <- setNames(cond_cells$leiden, cond_cells$cell_ID)
  sizes_cond <- table(labels_cond)
  total_edges_cond <- nrow(net_cond)

  cat(sprintf("  %s cells, %s edges\n",
              format(length(cond_cell_ids), big.mark = ","),
              format(total_edges_cond, big.mark = ",")))

  obs_counts_c <- count_edges(net_cond, labels_cond)
  obs_enrich_c <- compute_enrichment(obs_counts_c, sizes_cond, total_edges_cond)

  cat(sprintf("  Running %d permutations...\n", N_PERM))
  perm_extreme_c <- rep(0L, nrow(obs_enrich_c))

  cond_by_sample <- split(cond_cells$cell_ID, cond_cells$sample_id)
  labels_cond_by_sample <- split(cond_cells$leiden, cond_cells$sample_id)

  for (p in seq_len(N_PERM)) {
    if (p %% pb_interval == 0) cat(sprintf("    permutation %d/%d\n", p, N_PERM))

    perm_lab <- labels_cond
    for (sid in names(cond_by_sample)) {
      cids <- cond_by_sample[[sid]]
      perm_lab[cids] <- sample(labels_cond_by_sample[[sid]])
    }

    pc <- count_edges(net_cond, perm_lab)
    pe <- compute_enrichment(pc, sizes_cond, total_edges_cond)

    perm_vals_c <- pe$enrichment[match(obs_enrich_c$key_str, pe$key_str)]
    perm_extreme_c <- perm_extreme_c +
      (abs(perm_vals_c) >= abs(obs_enrich_c$enrichment))
  }

  obs_enrich_c[, p_value := (perm_extreme_c + 1) / (N_PERM + 1)]
  obs_enrich_c[, fdr := p.adjust(p_value, method = "BH")]
  obs_enrich_c[, condition := cond]
  obs_enrich_c[, cluster_from_name := cell_type_labels[as.character(c1)]]
  obs_enrich_c[, cluster_to_name := cell_type_labels[as.character(c2)]]

  n_sig_c <- sum(obs_enrich_c$fdr < 0.05)
  n_pos_c <- sum(obs_enrich_c$fdr < 0.05 & obs_enrich_c$enrichment > 0)
  n_neg_c <- sum(obs_enrich_c$fdr < 0.05 & obs_enrich_c$enrichment < 0)
  cat(sprintf("  %s: %d pairs, %d significant (FDR<0.05): %d co-loc, %d avoid\n",
              cond, nrow(obs_enrich_c), n_sig_c, n_pos_c, n_neg_c))

  results_cond[[cond]] <- obs_enrich_c
}

cond_combined <- rbindlist(results_cond)
fwrite(cond_combined[, .(condition, c1, c2, cluster_from_name, cluster_to_name,
                          observed = N, expected = round(expected, 1),
                          enrichment = round(enrichment, 4),
                          p_value = round(p_value, 5), fdr = round(fdr, 5))],
       file.path(out_dir, "proximity_permutation_per_condition.csv"))

# --- Differential enrichment with significance ---
irh <- results_cond[["irHepatitis"]][, .(key_str, c1, c2,
  irHep_enrich = enrichment, irHep_fdr = fdr)]
aih <- results_cond[["AIH"]][, .(key_str,
  AIH_enrich = enrichment, AIH_fdr = fdr)]
diff_dt <- merge(irh, aih, by = "key_str", all = TRUE)
diff_dt[is.na(irHep_enrich), irHep_enrich := 0]
diff_dt[is.na(AIH_enrich), AIH_enrich := 0]
diff_dt[, diff := irHep_enrich - AIH_enrich]
diff_dt[, cluster_from_name := cell_type_labels[as.character(c1)]]
diff_dt[, cluster_to_name := cell_type_labels[as.character(c2)]]
diff_dt[, sig_either := (irHep_fdr < 0.05) | (AIH_fdr < 0.05)]
diff_dt[, sig_both := (irHep_fdr < 0.05) & (AIH_fdr < 0.05)]

fwrite(diff_dt[, .(c1, c2, cluster_from_name, cluster_to_name,
                    irHep_enrich = round(irHep_enrich, 4),
                    AIH_enrich = round(AIH_enrich, 4),
                    diff = round(diff, 4),
                    irHep_fdr = round(irHep_fdr, 5),
                    AIH_fdr = round(AIH_fdr, 5),
                    sig_either, sig_both)],
       file.path(out_dir, "proximity_permutation_differential.csv"))


# ===========================================================================
# 3. SUMMARY
# ===========================================================================

sig_coloc <- obs_enrich[fdr < 0.05 & enrichment > 0][order(-enrichment)]
sig_avoid <- obs_enrich[fdr < 0.05 & enrichment < 0][order(enrichment)]

irh_sig <- results_cond[["irHepatitis"]][fdr < 0.05 & enrichment > 0][
  order(-enrichment)]
aih_sig <- results_cond[["AIH"]][fdr < 0.05 & enrichment > 0][
  order(-enrichment)]

irh_specific <- diff_dt[irHep_fdr < 0.05 & irHep_enrich > 0 &
                         (AIH_fdr >= 0.05 | AIH_enrich <= 0)][order(-diff)]
aih_specific <- diff_dt[AIH_fdr < 0.05 & AIH_enrich > 0 &
                         (irHep_fdr >= 0.05 | irHep_enrich <= 0)][order(diff)]

sink(file.path(out_dir, "proximity_permutation_summary.txt"))
cat("================================================================\n")
cat("SPATIAL PROXIMITY PERMUTATION TEST RESULTS\n")
cat(sprintf("Date: %s\n", Sys.Date()))
cat(sprintf("Permutations: %d (within-sample shuffling)\n", N_PERM))
cat(sprintf("Cells: %s | Edges: %s\n",
            format(nrow(meta), big.mark = ","),
            format(nrow(spat_net), big.mark = ",")))
cat("Enrichment = log2((observed+1)/(expected+1))\n")
cat("P-values: two-sided empirical | FDR: Benjamini-Hochberg\n")
cat("================================================================\n\n")

cat("--- OVERALL (all samples combined) ---\n")
cat(sprintf("Pairs tested: %d\n", n_pairs))
cat(sprintf("Significant (FDR<0.05): %d (%d co-localized, %d avoiding)\n\n",
            n_sig, n_sig_pos, n_sig_neg))

cat("Top 20 significant CO-LOCALIZATIONS (overall):\n")
cat(sprintf("  %-35s  enrich    FDR\n", "Pair"))
cat(paste(rep("-", 60), collapse = ""), "\n")
for (i in 1:min(20, nrow(sig_coloc))) {
  r <- sig_coloc[i]
  pair_str <- paste0(r$cluster_from_name, " \u2014 ", r$cluster_to_name)
  cat(sprintf("  %-35s  %+.3f  %.4f\n", pair_str, r$enrichment, r$fdr))
}

cat("\n\n--- PER-CONDITION ---\n")
for (cond in c("irHepatitis", "AIH")) {
  rc <- results_cond[[cond]]
  cat(sprintf("\n%s: %d cells\n", cond,
              nrow(cell_cluster[condition == cond])))
  cat(sprintf("  Significant (FDR<0.05): %d (%d co-loc, %d avoid)\n",
              sum(rc$fdr < 0.05),
              sum(rc$fdr < 0.05 & rc$enrichment > 0),
              sum(rc$fdr < 0.05 & rc$enrichment < 0)))
}

cat("\n\nirHepatitis-SPECIFIC co-localizations (sig in irHep, not in AIH):\n")
cat(sprintf("  %-35s  irHep    AIH     diff\n", "Pair"))
cat(paste(rep("-", 70), collapse = ""), "\n")
for (i in 1:min(20, nrow(irh_specific))) {
  r <- irh_specific[i]
  pair_str <- paste0(r$cluster_from_name, " \u2014 ", r$cluster_to_name)
  cat(sprintf("  %-35s  %+.3f  %+.3f  %+.3f\n",
              pair_str, r$irHep_enrich, r$AIH_enrich, r$diff))
}

cat("\n\nAIH-SPECIFIC co-localizations (sig in AIH, not in irHep):\n")
cat(sprintf("  %-35s  irHep    AIH     diff\n", "Pair"))
cat(paste(rep("-", 70), collapse = ""), "\n")
for (i in 1:min(20, nrow(aih_specific))) {
  r <- aih_specific[i]
  pair_str <- paste0(r$cluster_from_name, " \u2014 ", r$cluster_to_name)
  cat(sprintf("  %-35s  %+.3f  %+.3f  %+.3f\n",
              pair_str, r$irHep_enrich, r$AIH_enrich, r$diff))
}

sink()

cat("\n================================================================\n")
cat("PERMUTATION TEST COMPLETE\n")
cat("================================================================\n")
cat(sprintf("Results: %s\n", out_dir))
cat(sprintf("  Overall: proximity_permutation_results.csv\n"))
cat(sprintf("  Per-condition: proximity_permutation_per_condition.csv\n"))
cat(sprintf("  Differential: proximity_permutation_differential.csv\n"))
cat(sprintf("  Summary: proximity_permutation_summary.txt\n"))
