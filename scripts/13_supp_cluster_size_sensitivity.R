#!/usr/bin/env Rscript
# ============================================================================
# 13_supp_cluster_size_sensitivity.R  (Supplementary Table S9)
# Cluster-size sensitivity check for the Figure 2B CCC analysis.
# Takes the top-10 coherent flows per
# direction (by |log2FC|) from the canonical Figure 2B file, and report for
# each whether both endpoint clusters retain at least 30 and at least 50
# cells in every sample. This documents the robustness of the manuscript's
# top-ranked claims to a more conservative cluster-size threshold than the
# n >= 10 filter used in the published figure.
# Inputs:
#   data/processed/tables/cluster_composition.csv (per-cluster per-sample N)
#   data/processed/clustering/cluster_characterization.csv (Leiden -> name)
#   data/processed/cellchat/cellchat_coherent_interactions.csv (canonical)
# Output:
#   tables/supplementary/Table_S9_cluster_size_sensitivity.csv
# ============================================================================
suppressPackageStartupMessages({
  library(data.table)
})

source(file.path(Sys.getenv("IRHEP_BUNDLE_ROOT", here::here()),
                 "scripts", "config.R"))

cluster_comp <- fread(file.path(PROCESSED_DIR, "tables",
                                "cluster_composition.csv"))
cluster_char <- fread(file.path(PROCESSED_DIR, "clustering",
                                "cluster_characterization.csv"))
canonical <- fread(CELLCHAT_COHERENT)

# Per-cluster minimum per-sample cell count, by cluster_name (matches the
# names CellChat uses in canonical$source / $target).
id2name <- setNames(cluster_char$cluster_name, as.character(cluster_char$leiden))
excluded <- c("12", "15", "21")
cluster_comp[, leiden_chr := as.character(leiden)]
cluster_comp <- cluster_comp[!(leiden_chr %in% excluded)]
cluster_comp[, cluster_name := id2name[leiden_chr]]
min_per_cluster <- cluster_comp[, .(min_n_per_sample = min(N)),
                                 by = cluster_name]

# Annotate canonical flows with the worst (smallest) per-sample cell count
# of source and target clusters.
src_min <- min_per_cluster[, .(cluster_name, src_min = min_n_per_sample)]
tgt_min <- min_per_cluster[, .(cluster_name, tgt_min = min_n_per_sample)]
canonical <- merge(canonical, src_min, by.x = "source",
                   by.y = "cluster_name", all.x = TRUE)
canonical <- merge(canonical, tgt_min, by.x = "target",
                   by.y = "cluster_name", all.x = TRUE)
canonical[, triple_min_per_sample := pmin(src_min, tgt_min, na.rm = TRUE)]
canonical[, preserved_at_n30 := triple_min_per_sample >= 30]
canonical[, preserved_at_n50 := triple_min_per_sample >= 50]

# Top-10 per direction by |log2FC| (the manuscript's headline flows).
top10 <- rbindlist(lapply(c("irHep", "AIH"), function(d) {
  sub <- canonical[direction == d]
  sub <- sub[order(-abs(log2FC))]
  sub <- head(sub, 10)
  sub[, rank := seq_len(.N)]
  sub
}))

out <- top10[, .(direction, rank, source, target, pathway_label, log2FC,
                  triple_min_per_sample,
                  preserved_at_n30, preserved_at_n50)]

out_dir <- file.path(BUNDLE_ROOT, "tables", "supplementary")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
fwrite(out, file.path(out_dir, "Table_S9_cluster_size_sensitivity.csv"))
cat(sprintf("Wrote Table_S8 (%d rows).\n", nrow(out)))

summary_dt <- out[, .(n_preserved_at_n30 = sum(preserved_at_n30),
                      n_preserved_at_n50 = sum(preserved_at_n50)),
                  by = direction]
print(summary_dt)
