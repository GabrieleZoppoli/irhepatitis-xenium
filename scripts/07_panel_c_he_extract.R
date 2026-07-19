#!/usr/bin/env Rscript
# ============================================================================
# Writes roi_config.json for the 4 Panel C supplementary ROIs and calls
# the shared Python H&E/DAPI extractor (07a_panel_c_he_extract.py).
# ============================================================================
suppressPackageStartupMessages({ library(jsonlite) })

source(file.path(Sys.getenv("IRHEP_BUNDLE_ROOT", here::here()), "scripts", "config.R"))
BASE <- PROCESSED_DIR
WD <- RAW_DIR
OUT  <- file.path(BASE, "figure_2_panels", "supp_panelC_crops")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# 4 ROIs -- all 250x250 um for denser tissue context.
# AIH hep ROI 104 updated to (17500, -8900) -- selected by scan_aih_hep_roi.R
# as the 250 um window in 56784_c that maximizes hepatocyte-CD4 direct contact
# while keeping cholangiocytes near zero (so it is a CLEAN hepatocyte region):
#   n_hep=184, n_chol=1, n_CD4=20, n_CD8=91,
#   n_hep_with_CD4_contact=23 (up from 18 in the prior ROI),
#   n_hep_with_CD8_contact=89.
# The prior ROI had 47 CD4 in window but only 18 touched hepatocytes; the new
# ROI has fewer CD4 but more of them are hepatocyte-engaged, which is the
# biological point of this panel.
#
# AIH bile-duct ROI 103 updated 2026-04-27 from 56784_c (13575, -9905) to
# 54023_b (4350, -9850). Reason: prior ROI scored only 6 on the Chol+CD8T2/T3
# hotspot metric (sqrt(n_chol * n_cd8t23) approximated by chol_with_cd23 *
# cd23_with_chol mass), while sample 54023_b has the top-scoring AIH bin at
# (4350, -9850), score 49 (n_chol=220, n_cd8t23=20, n_cd8t1=116). Pairing the
# top irHep hotspot (score 64, sample 18321_a) with a low-density AIH bin
# manufactures visual contrast that does not reflect Panel A/B per-sample
# means. Switching ROI 103 to the symmetric top-scoring AIH bin enforces
# the same selection rule across both conditions. See scan_aih_chol_cd8t23.csv
# and scan_irhep_chol_cd8t23.csv for the full ranked tables.
rois <- list(
  # ROW 1 -- Bile-duct axis (Chol + CD8T2/T3). Top-1 Chol+CD8T2/T3 hotspot per
  # condition (sqrt(n_Chol * n_CD8T2/T3) score). Per-condition rule applied
  # symmetrically. See scan_{irhep,aih}_chol_cd8t23.csv.
  list(cluster_id = 101, sample_id = "18321_a",
       roi_x = 1473, roi_y = -2249, roi_size = 250),     # irHep bile duct
  list(cluster_id = 103, sample_id = "54023_b",
       roi_x = 4350, roi_y = -9850, roi_size = 250),     # AIH bile duct
  # ROW 2 -- Hepatocyte axis (Hep + CD4T). Selected by the INTERSPERSION
  # metric: fraction of CD4 cells located within 25 um of a hepatocyte
  # (i.e. CD4 cells AT THE PARENCHYMAL INTERFACE, not in a separate aggregate).
  # The AIH narrative requires CD4 cells INTERSPERSED between and within
  # hepatocytes -- a single tight lymphoid follicle adjacent to but apart
  # from Hep would not display the AIH-CD4 spatial signature. Selection rule:
  #   - irHep: window with similar Hep+CD4 cell counts to the AIH window,
  #     but very LOW cd4_at_interface_frac (CD4 dispersed away from Hep)
  #   - AIH:   window with similar CD4 count, but very HIGH
  #     cd4_at_interface_frac (CD4 truly interspersed in parenchyma)
  # See scan_{irhep,aih}_hep_cd4_intersperse.csv.
  #   irHep ROI 102: 18321_a (848, -2924), Hep=118, CD4=50,
  #     CD4-at-Hep-interface = 3/50 (6.0%), Hep-with-CD4 = 4/118 (3.4%)
  #     -- Hep parenchyma with CD4 in a separate clump, almost no CD4 touching Hep.
  #   AIH ROI 104:   56784_c (16100, -8480), Hep=137, CD4=51,
  #     CD4-at-Hep-interface = 35/51 (68.6%), Hep-with-CD4 = 51/137 (37.2%)
  #     -- Hep parenchyma with CD4 cells interspersed throughout, virtually
  #     every CD4 cell touching at least one Hep.
  # Matched CD4 abundance (50 vs 51) + 11x interspersion differential.
  list(cluster_id = 102, sample_id = "18321_a",
       roi_x = 848, roi_y = -2924, roi_size = 250),     # irHep hepatocyte
  list(cluster_id = 104, sample_id = "56784_c",
       roi_x = 16400, roi_y = -8630, roi_size = 250)    # AIH hepatocyte
)

config_json <- list(
  working_dir = WD,
  output_dir  = OUT,
  pixel_size  = 0.2125,
  rois        = rois
)
json_path <- file.path(OUT, "roi_config.json")
write_json(config_json, json_path, auto_unbox = TRUE, pretty = TRUE)
cat(sprintf("Wrote %s\n", json_path))

# Run the shared Python extractor (uses cluster_id numbers in output filenames)
py_script <- file.path(SCRIPTS_DIR, "07a_panel_c_he_extract.py")
cat(sprintf("Running: python3 %s %s\n", py_script, json_path))
status <- system(sprintf("python3 '%s' '%s'", py_script, json_path))
cat(sprintf("Python exit: %d\n", status))

# List files produced
cat("\nFiles produced:\n")
for (f in list.files(OUT, full.names = TRUE)) cat("  ", f, "\n")
