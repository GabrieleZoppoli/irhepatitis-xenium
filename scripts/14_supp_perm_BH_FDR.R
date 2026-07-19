#!/usr/bin/env Rscript
# ============================================================================
# 14_supp_perm_BH_FDR.R  (Supplementary Table S10)
# Multi-test correction across CellChat-coherent triples.
# The canonical coherent rule is a between-sample agreement
# filter, not a formal hypothesis test. We supply a per-triple
# Benjamini-Hochberg-adjusted p-value derived by combining the per-sample
# CellChat permutation-based p-values across the two samples of each
# condition via Fisher's method, then BH-adjusting across the full set of
# coherent triples.
# Inputs:
#   data/processed/cellchat/per_sample/
#       cellchat_per_sample_communications.csv (per-sample prob and pval)
#       cellchat_differential_unbiased.csv (coherent set source)
# Output:
#   tables/supplementary/Table_S10_perm_BH_FDR.csv
# ============================================================================
suppressPackageStartupMessages({
  library(data.table)
})

source(file.path(Sys.getenv("IRHEP_BUNDLE_ROOT", here::here()),
                 "scripts", "config.R"))

CCDIR      <- file.path(PROCESSED_DIR, "cellchat", "per_sample")
per_sample <- fread(file.path(CCDIR, "cellchat_per_sample_communications.csv"))
unbiased   <- fread(file.path(CCDIR, "cellchat_differential_unbiased.csv"))

multi_receptor_families <- c("MHC-I", "MHC-II", "CXCL", "CCL")
make_pathway_label <- function(pathway_name, receptor) {
  ifelse(pathway_name %in% multi_receptor_families,
         paste0(pathway_name, "_", receptor),
         pathway_name)
}
unbiased[, pathway_label := make_pathway_label(pathway_name, receptor)]
unbiased[, direction := ifelse(log2_fc > 0, "irHep", "AIH")]

# Canonical coherent rule (matches Methods text and the 47-row file used
# by Figure 2B): in the enriched condition, at least one sample sig; in
# the opposite condition, none sig.
unbiased[, coherent := (
  (direction == "irHep" & irHep_n_sig >= 1 & AIH_n_sig == 0) |
  (direction == "AIH"   & AIH_n_sig   >= 1 & irHep_n_sig == 0)
)]
coh <- unbiased[coherent == TRUE]
coh[, abs_lfc := abs(log2_fc)]
coh <- coh[, .SD[which.max(abs_lfc)],
           by = .(source, target, pathway_label)]

# Per-condition p-value via Fisher's combined chi-square (k = 2 samples).
# Replace per-sample pval = 0 with 1/(B+1) where B = 100 is the CellChat
# default permutation budget (avoids -Inf log on combination).
B <- 100
pseudo <- 1 / (B + 1)
per_sample[, pval_floor := pmax(pval, pseudo)]

# Aggregate per (source, target, ligand, receptor, condition) to
# Fisher-combined p across two samples.
combine_fisher <- function(p) {
  k <- length(p)
  if (k == 0) return(NA_real_)
  X <- -2 * sum(log(p))
  pchisq(X, df = 2 * k, lower.tail = FALSE)
}
per_cond <- per_sample[, .(
  comb_p = combine_fisher(pval_floor),
  n_samples = .N
), by = .(source, target, ligand, receptor, pathway_name, condition)]

per_cond_wide <- dcast(per_cond,
                       source + target + ligand + receptor + pathway_name ~ condition,
                       value.var = "comb_p")
setnames(per_cond_wide, c("irHepatitis", "AIH"),
         c("comb_p_irHep", "comb_p_AIH"), skip_absent = TRUE)

# Merge into coherent set; the "enriched" condition's combined p drives
# the BH adjustment.
joined <- merge(coh, per_cond_wide,
                by = c("source", "target", "ligand", "receptor", "pathway_name"),
                all.x = TRUE)
joined[, enriched_p := ifelse(direction == "irHep",
                              comb_p_irHep, comb_p_AIH)]
joined[, BH_FDR := p.adjust(enriched_p, method = "BH")]

# Schema for Table S9.
out <- joined[, .(
  source, target, pathway_label, direction,
  log2FC = log2_fc,
  p_perm = enriched_p,
  BH_FDR
)]
out <- out[order(direction, -abs(log2FC))]

out_dir <- file.path(BUNDLE_ROOT, "tables", "supplementary")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
fwrite(out, file.path(out_dir, "Table_S10_perm_BH_FDR.csv"))
cat(sprintf("Wrote Table_S9 (%d rows).\n", nrow(out)))

cat(sprintf("Top-10 SR-irHep flows passing BH FDR < 0.10: %d\n",
            sum(out[direction == "irHep"][order(-abs(log2FC))][1:10]$BH_FDR < 0.10,
                na.rm = TRUE)))
cat(sprintf("Top-10 AIH flows passing BH FDR < 0.10: %d\n",
            sum(out[direction == "AIH"][order(-abs(log2FC))][1:10]$BH_FDR < 0.10,
                na.rm = TRUE)))
