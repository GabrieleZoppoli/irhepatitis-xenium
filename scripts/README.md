# Scripts — pipeline overview

All scripts source `config.R` (R) or `config.py` (Python) for paths, sample IDs, and per-sample shifts. Override the bundle root by setting `IRHEP_BUNDLE_ROOT` before launching:

```bash
export IRHEP_BUNDLE_ROOT=/path/to/irhepatitis-xenium
```

## Pipeline order

| # | Script | Stage | Inputs | Outputs |
|---|---|---|---|---|
| 01 | `01_giotto_normalize_cluster.R` | QC-filter, log-normalisation, PCA + Harmony batch correction, UMAP and SNN-Leiden clustering across all 4 Xenium-Analyzer per-sample outputs (Xenium-segmented cells + decoded transcripts) | `data/raw/0027420_*` | `data/processed/session_snapshot.RData` |
| 02 | `02_cluster_config.R` | Cluster naming, palette, lineage groupings (sourced by all downstream scripts) | `session_snapshot.RData` | `cluster_config.rds` |
| 03 | `03_figure_1_main.R` | Figure 1 build (cellular landscape) | `session_snapshot.RData`, `cluster_config.rds` | `figures/main/Figure_1.pdf` |
| 04 | `04_panel_a_proximity_compute.R` | Per-sample `cellProximityEnrichment` (Giotto, 1000 permutations); writes `giotto_proximity_per_sample.csv` | `gobj_with_spatial_networks.rds` | proximity CSV (input to 04a/04b) |
| 04a | `04a_panel_a_proximity_heatmap.R` | Coherent-rule heatmap + lollipop variants (strict / relaxed); top-10 differential pairs | `giotto_proximity_per_sample.csv` | exploration heatmaps + differential table |
| 04b | `04b_panel_a_panel_only.R` | Stand-alone Panel A render for the composite | `giotto_proximity_per_sample.csv` | `fig2_a_heatmap_only_panel.{pdf,png}` + legend block |
| 05 | `05_panel_b_cellchat.R` | CellChat v2 per-sample (contact-range 25 µm), ligand-split MHC-I/CXCL/CCL families, coherent rule | `session_snapshot.RData` | per-sample CellChat objects, coherent flows CSV |
| 05a | `05a_panel_b_liana.R` | LIANA+ 5-method consensus per-sample (NATMI, Connectome, logFC, SCA, CellPhoneDB) | `session_snapshot.RData` | LIANA coherent CSV |
| 05b | `05b_panel_b_commot.py` | COMMOT spatial optimal-transport CCC per sample (Python) — long-running (≈19 h on largest sample) | per-sample expression + spatial coords | per-sample COMMOT CSV |
| 05c | `05c_panel_b_commot_aggregate.R` | Aggregate per-sample COMMOT into condition-level differentials with explicit 3-layer NA handling (pair-not-tested / cluster-too-small / genuine-zero) | COMMOT per-sample CSVs | `commot_differential_*.csv` |
| 05d | `05d_panel_b_commot_crossref.R` | CellChat ↔ COMMOT cross-reference; agreement table | CellChat coherent + COMMOT differentials | `commot_cellchat_agreement.csv` |
| 06 | `06_panel_b_sankey.R` | Render Panel B Sankey (top-10 coherent flows per direction by |log₂FC|) | `cellchat_coherent_interactions.csv` | Panel B PDF/PNG |
| 07 | `07_panel_c_he_extract.R` | Wrapper that writes ROI config JSON and calls Python H&E crop pipeline | ROI coords (in script) + raw Xenium H&E pyramids | per-ROI H&E PNGs |
| 07a | `07a_panel_c_he_extract.py` | Affine-warp Xenium-aligned H&E into Giotto coordinate space; outputs 600×600 px crops at 0.85 µm/px | H&E pyramids + `he_alignment/matrix.csv` | per-ROI H&E PNGs |
| 08 | `08_panel_c_rois.R` | Panel C build: 4 representative ROIs (2 anatomy rows × 2 condition columns) with row-specific selective highlighting | `session_snapshot.RData` + cluster boundaries + H&E crops | Panel C PDF/PNG |
| 08a | `08a_panel_c_score_distribution.R` | ROI-selection diagnostic: distribution of all sliding-bin Chol+CD8T2/T3 scores per sample with the chosen ROIs marked | `scan_*_chol_cd8t23.csv` | `figure_S_panelC_score_distribution.pdf` |
| 09 | `09_figure_2_assembled.R` | Compose Figure 2 from Panel A, B, C panels | per-panel PNGs | `figures/main/Figure_2.pdf` |
| 10 | `10_supp_he_lens_rois.R` | Supplementary 600 µm "lens" ROI panels for Figure 1 (full cluster palette) | `session_snapshot.RData` + H&E crops | Supp lens ROI panels |
| 10a | `10a_supp_he_lens_extract.py` | H&E crop extractor for Fig 1 lens panels | H&E pyramids | per-sample H&E crops |
| 11 | `11_supp_proximity_permutation.R` | Per-sample proximity raw enrichment + permutation diagnostics (supplement to Panel A) | `session_snapshot.RData` | per-sample proximity tables |
| 13 | `13_supp_cluster_size_sensitivity.R` | Cluster-size sensitivity of the Figure 2B CCC flows at n≥30 / n≥50 (Table S9) | `cellchat_coherent_interactions.csv` | `Table_S9_*.csv` |
| 14 | `14_supp_perm_BH_FDR.R` | Fisher-combined per-sample CellChat permutation p-values, BH-adjusted (Table S10) | per-sample CellChat p-values | `Table_S10_*.csv` |
| 15 | `15_supp_giotto_analyses.R` | Exhaustion/TRM marker expression and clustering-stability analyses (Tables S11, S12/S12b/S12c; Fig S3) | `session_snapshot.RData` | supplementary tables + Fig S3 |
| 15a | `15a_supp_cd8t1_trm_signature.R` | Multi-marker TRM signature across CD8 sub-clusters (Table S11b, Fig S4) | `session_snapshot.RData` | `Table_S11b_*.csv`, Fig S4 |
| 15b | `15b_supp_ifn_response_score.R` | Cholangiocyte IFN-response composite score (Tables S11c/S11c2/S11c3, Fig S5) | `session_snapshot.RData` | `Table_S11c*_*.csv`, Fig S5 |

## Conventions

- **Path resolution.** `BUNDLE_ROOT` is detected via `IRHEP_BUNDLE_ROOT` env var → `here::here()` → fallback. Override at runtime if needed.
- **Sample IDs.** `18321_a` and `24774_d` are irHepatitis; `54023_b` and `56784_c` are AIH. Per-sample x/y shifts for composite-tissue plots are in `config.R::SAMPLE_SHIFTS`.
- **Cluster names.** Single source of truth in `02_cluster_config.R`. Excludes cluster 21 (n=4) and clusters 12, 15 (segmentation fragments).
- **CCC method choice.** CellChat v2 with `contact.range = 25 µm` is primary (Panel B). LIANA+ and COMMOT are independent cross-validation pipelines whose outputs feed Supplementary Tables S1–S2; CellChat-LIANA-COMMOT triple-confirmed flows are reported in the manuscript text.
- **Random seeds.** Set in `01_giotto_normalize_cluster.R` (Leiden seed) and in proximity permutations; otherwise R/Python defaults.

## Reproducibility caveats

- COMMOT runtime scales with sample size; sample 54023_b (52 k cells, 33 ligand–receptor pairs) takes ≈19 h on a 32-core 376 GB workstation.
- The `IRHEP_BUNDLE_ROOT` env var must be set before any script runs; otherwise `here::here()` will resolve to the user's current directory and paths will be wrong.
- The H&E crop pipeline (`07a` / `10a`) requires a Python environment with `tifffile`, `scipy.ndimage`, `PIL`. The conda spec is at `env/commot_env.yml` (also satisfies COMMOT requirements).
