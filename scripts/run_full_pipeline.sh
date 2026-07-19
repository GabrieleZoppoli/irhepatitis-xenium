#!/usr/bin/env bash
# ============================================================================
# run_full_pipeline.sh — sequential end-to-end re-run of the irHepatitis
# Xenium publication pipeline. Multi-day runtime (dominated by COMMOT,
# ~19 h on the largest sample). Usage:
#   IRHEP_BUNDLE_ROOT=/path/to/bundle bash run_full_pipeline.sh
# Logs are written to data/processed/logs/.
# ============================================================================
set -euo pipefail

if [ -z "${IRHEP_BUNDLE_ROOT:-}" ]; then
  IRHEP_BUNDLE_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
  export IRHEP_BUNDLE_ROOT
fi
echo "Bundle root: $IRHEP_BUNDLE_ROOT"

LOGDIR="$IRHEP_BUNDLE_ROOT/data/processed/logs"
mkdir -p "$LOGDIR"

run_R () {
  local script="$1"; local label="$2"
  echo "[$(date +'%F %T')] -> $label  ($script)"
  Rscript "$IRHEP_BUNDLE_ROOT/scripts/$script" \
    > "$LOGDIR/${label}.log" 2>&1
  echo "[$(date +'%F %T')] <- $label  done"
}

run_py () {
  local script="$1"; local label="$2"
  echo "[$(date +'%F %T')] -> $label  ($script)"
  python3 "$IRHEP_BUNDLE_ROOT/scripts/$script" \
    > "$LOGDIR/${label}.log" 2>&1
  echo "[$(date +'%F %T')] <- $label  done"
}

# Stage 1 — Giotto loading, clustering, spatial networks (1–2 h)
# 02_cluster_config.R is sourced by downstream scripts; no need to run separately.
run_R 01_giotto_normalize_cluster.R         stage_01_giotto
run_R 01a_add_spatial_networks.R            stage_01a_spatial_networks

# Stage 2 — Figure 1 (minutes)
run_R 03_figure_1_main.R                    stage_03_figure_1

# Stage 3 — Panel A: proximity (~30 min)
run_R 04_panel_a_proximity_compute.R        stage_04_panel_a_compute
run_R 04a_panel_a_proximity_heatmap.R       stage_04a_panel_a_heatmap
run_R 04b_panel_a_panel_only.R              stage_04b_panel_a_panel
run_R 11_supp_proximity_permutation.R       stage_11_supp_prox_permut

# Stage 4 — Panel B: CellChat (~40 min)
run_R 05_panel_b_cellchat.R                 stage_05_cellchat

# Stage 5 — Panel B: LIANA cross-validation (~5 min)
run_R 05a_panel_b_liana.R                   stage_05a_liana

# Stage 6 — Panel B: COMMOT (~19 h on 54023_b alone) — Python
run_py 05b_panel_b_commot.py                stage_05b_commot

# Stage 7 — COMMOT aggregator + cross-reference (minutes)
run_R 05c_panel_b_commot_aggregate.R        stage_05c_commot_aggregate
run_R 05d_panel_b_commot_crossref.R         stage_05d_commot_crossref

# Stage 8 — Panel B Sankey (<1 min)
run_R 06_panel_b_sankey.R                   stage_06_sankey

# Stage 9 — Panel C: H&E extract + ROI build + score distribution (minutes)
run_R 07_panel_c_he_extract.R               stage_07_he_extract
run_R 08_panel_c_rois.R                     stage_08_panel_c
run_R 08a_panel_c_score_distribution.R      stage_08a_score_dist
run_R 08b_panel_c_close_apposition_counts.R stage_08b_close_apposition

# Stage 10 — Figure 2 assembly (<1 min)
run_R 09_figure_2_assembled.R               stage_09_assembled

# Stage 11 — Supplementary lens ROIs (Figure 1 supplement)
run_R 10_supp_he_lens_rois.R                stage_10_supp_lens

# Stage 12 — Additional supplementary tables and figures
#   13:  Table S9   — cluster-size sensitivity at n>=30/n>=50
#   14:  Table S10  — Fisher-combined permutation p-values + BH FDR
#   15:  Table S11 (exhaustion/TRM markers), Tables S12/S12b/S12c
#        (clustering stability ARI + silhouette + correspondence),
#        Figure S3 (clustering stability)
#   15a: Table S11b + Figure S4 — multi-marker TRM signature for CD8T1
#   15b: Tables S11c/S11c2/S11c3 + Figure S5 — cholangiocyte IFN-response
#        composite (STAT1/IRF1/IDO1/ISG15/MX1) vs CXCL9/10/11 co-expression
run_R 13_supp_cluster_size_sensitivity.R    stage_13_supp_cluster_size
run_R 14_supp_perm_BH_FDR.R                 stage_14_supp_perm_BH_FDR
run_R 15_supp_giotto_analyses.R             stage_15_supp_giotto_analyses
run_R 15a_supp_cd8t1_trm_signature.R        stage_15a_supp_trm_signature
run_R 15b_supp_ifn_response_score.R         stage_15b_supp_ifn_response

echo "[$(date +'%F %T')] FULL PIPELINE COMPLETE"
echo "Logs: $LOGDIR"
echo "Figures: $IRHEP_BUNDLE_ROOT/figures/"
