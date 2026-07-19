#!/usr/bin/env Rscript
# ============================================================================
# Cross-reference COMMOT strict coherent findings with CellChat Panel B.
#
# Inputs:
#   - figure_2_panels/cellchat_coherent_interactions.csv   (CellChat Panel B)
#   - data/processed/commot/commot_differential_strict.csv     (COMMOT L-R)
#   - data/processed/commot/commot_differential_pathway.csv    (COMMOT pw)
#
# Matching strategy:
#   CellChat uses a sub-pathway label (pathway × receptor), e.g. "MHC-I_CD8A".
#   COMMOT outputs (a) individual L-R pair (e.g. "HLA-B-CD8A") and (b) pathway
#   aggregate (e.g. "MHC-I"). We build two lookup tables:
#     (i)  pathway_label -> set of COMMOT L-R pair ids
#     (ii) pathway_label -> COMMOT pathway aggregate id
#   A CellChat (source, target, pathway_label, direction) is "confirmed" if any
#   matching COMMOT (lr_id, source, target) triple has the same direction.
#
# Output:
#   - data/processed/commot/commot_cellchat_agreement.csv
#   - data/processed/commot/commot_cellchat_agreement.txt (summary)
#
# R-only; uses data.table. Conda env: xenium_env.
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

source(file.path(Sys.getenv("IRHEP_BUNDLE_ROOT", here::here()), "scripts", "config.R"))
COMMOT_DIR <- file.path(PROCESSED_DIR, "commot")
CC_FILE  <- if (file.exists(file.path(PROCESSED_DIR, "cellchat", "cellchat_coherent_interactions.csv"))) {
  file.path(PROCESSED_DIR, "cellchat", "cellchat_coherent_interactions.csv")
} else {
  file.path(PROCESSED_DIR, "figure_2_panels", "cellchat_coherent_interactions.csv")
}
COM_LR   <- file.path(COMMOT_DIR, "commot_differential_strict.csv")
COM_PW   <- file.path(COMMOT_DIR, "commot_differential_pathway.csv")
OUT_DIR  <- COMMOT_DIR

cc  <- fread(CC_FILE)
com_lr <- fread(COM_LR)
com_pw <- fread(COM_PW)

cat(sprintf("CellChat coherent interactions:  %d\n", nrow(cc)))
cat(sprintf("COMMOT strict coherent triples:   %d\n", nrow(com_lr)))
cat(sprintf("COMMOT pathway coherent triples:  %d\n", nrow(com_pw)))

com_pair_ids <- unique(com_lr$lr_id)
com_pw_ids   <- unique(com_pw$lr_id)

cat(sprintf("\nCOMMOT unique L-R pairs in strict set: %s\n",
            paste(com_pair_ids, collapse = ", ")))
cat(sprintf("COMMOT unique pathway aggregates in strict set: %s\n",
            paste(com_pw_ids, collapse = ", ")))

# Manual mapping based on CellChatDB conventions.
# CellChat sub-pathway label "CCL_CCR1" means: pathway=CCL, restricted to CCR1 receptors.
# COMMOT L-R pair "CCL5-CCR1" is one specific match; the pathway-level "CCL" is less specific.
map_cellchat_to_commot <- function(pathway_label) {
  # Split sub-pathway: "CCL_CCR1" -> pathway="CCL", receptor_keyword="CCR1"
  parts <- strsplit(pathway_label, "_", fixed = TRUE)[[1]]
  pathway <- parts[1]
  rec_key <- if (length(parts) > 1) parts[2] else NA

  # A) match L-R pairs where ligand matches pathway name as prefix and receptor matches rec_key
  #    Handles "CCL_CCR1" → "CCL.*-CCR1", "CXCL_CXCR3" → "CXCL.*-CXCR3"
  #    For "MHC-I_CD8A" the pathway has a hyphen; treat as ligand-family "HLA" + receptor "CD8A"
  #    (special case: CellChat MHC-I_CD8A maps to COMMOT HLA-B-CD8A)
  regex_lr <- if (pathway == "MHC-I") {
                paste0("^HLA.*-", rec_key, "$")
              } else if (pathway %in% c("CCL", "CXCL")) {
                paste0("^", pathway, "[0-9]+[A-Z]*-", rec_key, "$")
              } else if (pathway == "CSF") {
                "^CSF1-CSF1R$"
              } else if (pathway == "VEGF") {
                "^VEGFA-FLT1$"
              } else if (pathway == "FN1") {
                "^FN1-"
              } else if (pathway == "FASLG") {
                "^FASLG-FAS$"
              } else if (pathway == "ICOS") {
                "^ICOSLG-"
              } else if (pathway == "CD40") {
                "^CD40LG?-"
              } else if (pathway == "CD80") {
                "^CD80-"
              } else if (pathway == "CD86") {
                "^CD86-"
              } else if (pathway == "CEACAM") {
                "^CEACAM"
              } else {
                paste0("^", pathway)
              }
  matching_pairs <- grep(regex_lr, com_pair_ids, value = TRUE)
  # B) pathway-level aggregate: exact match if pathway is in com_pw_ids
  pw_match <- if (pathway %in% com_pw_ids) pathway else NA_character_
  list(pairs = matching_pairs, pathway = pw_match)
}

cc[, `:=`(commot_pairs_match = NA_character_,
          commot_pathway_match = NA_character_)]
for (i in seq_len(nrow(cc))) {
  m <- map_cellchat_to_commot(cc$pathway_label[i])
  cc[i, commot_pairs_match   := paste(m$pairs, collapse = "|")]
  cc[i, commot_pathway_match := m$pathway]
}

com_dir <- function(x) fifelse(x == "irHepatitis", "irHep", "AIH")
com_lr[, dir_short  := com_dir(enriched_strict)]
com_pw[, dir_short  := com_dir(enriched_strict)]

cc[, commot_pair_confirmed    := FALSE]
cc[, commot_pathway_confirmed := FALSE]
cc[, commot_pair_matches_any  := FALSE]
cc[, commot_pair_coherent_opposite := FALSE]

for (i in seq_len(nrow(cc))) {
  src <- cc$source[i];  tgt <- cc$target[i];  dir <- cc$direction[i]
  pairs_str <- cc$commot_pairs_match[i]
  pathway   <- cc$commot_pathway_match[i]
  if (!is.na(pairs_str) && nzchar(pairs_str)) {
    pairs <- strsplit(pairs_str, "|", fixed = TRUE)[[1]]
    # Any matching triple (same direction)
    hit_same <- com_lr[lr_id %in% pairs & source == src & target == tgt &
                       dir_short == dir]
    hit_opp  <- com_lr[lr_id %in% pairs & source == src & target == tgt &
                       dir_short != dir]
    hit_any  <- com_lr[lr_id %in% pairs & source == src & target == tgt]
    cc[i, commot_pair_confirmed       := nrow(hit_same) > 0]
    cc[i, commot_pair_coherent_opposite := nrow(hit_opp) > 0]
    cc[i, commot_pair_matches_any     := nrow(hit_any) > 0]
  }
  if (!is.na(pathway) && nzchar(pathway)) {
    hit_pw <- com_pw[lr_id == pathway & source == src & target == tgt &
                     dir_short == dir]
    cc[i, commot_pathway_confirmed := nrow(hit_pw) > 0]
  }
}

cc[, commot_any_confirmed := commot_pair_confirmed | commot_pathway_confirmed]

fwrite(cc, file.path(OUT_DIR, "commot_cellchat_agreement.csv"))

sink(file.path(OUT_DIR, "commot_cellchat_agreement.txt"))
cat("================================================================\n")
cat("  COMMOT ↔ CellChat cross-reference (Panel B)\n")
cat("  Generated:", format(Sys.time()), "\n")
cat("================================================================\n\n")

cat(sprintf("CellChat coherent interactions:        %d\n", nrow(cc)))
cat(sprintf("  (%d AIH-enriched, %d irHep-enriched)\n",
            sum(cc$direction == "AIH"), sum(cc$direction == "irHep")))
cat(sprintf("COMMOT strict coherent L-R triples:    %d\n", nrow(com_lr)))
cat(sprintf("COMMOT pathway coherent triples:       %d\n", nrow(com_pw)))
cat("\n")

cat("Per-direction agreement (CellChat → confirmed by matching COMMOT triple):\n")
print(cc[, .(
  N                         = .N,
  any_confirmed             = sum(commot_any_confirmed),
  pair_confirmed            = sum(commot_pair_confirmed),
  pathway_confirmed         = sum(commot_pathway_confirmed),
  pair_coherent_opposite    = sum(commot_pair_coherent_opposite),
  pair_matches_any_direction= sum(commot_pair_matches_any)
), by = direction])

cat("\n")
cat("Per-pathway agreement (CellChat → COMMOT any-match):\n")
print(cc[, .(
  N = .N,
  confirmed = sum(commot_any_confirmed)
), keyby = pathway_label])

cat("\nCellChat coherent interactions NOT confirmed by COMMOT:\n")
print(cc[commot_any_confirmed == FALSE,
         .(source, target, pathway_label, direction, log2FC,
           commot_pairs_match, commot_pathway_match)])

cat("\nCellChat coherent interactions CONFIRMED by COMMOT:\n")
print(cc[commot_any_confirmed == TRUE,
         .(source, target, pathway_label, direction, log2FC)])

sink()

cat(sprintf("\nWrote: %s\n", file.path(OUT_DIR, "commot_cellchat_agreement.csv")))
cat(sprintf("Wrote: %s\n", file.path(OUT_DIR, "commot_cellchat_agreement.txt")))
cat("\nDone.\n")
