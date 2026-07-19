# irhepatitis-xenium

Analysis code and processed data for *A bile-duct CD8/CXCR3 axis defines steroid-resistant checkpoint-inhibitor hepatitis*, a 10× Xenium spatial-transcriptomics study comparing steroid-resistant immune-related hepatitis to treatment-naive autoimmune hepatitis.

## Repository layout

```
irhepatitis-xenium/
├── manuscript/                     # Brief Communication + supplementary
├── methods/                        # methods text
├── figures/
│   ├── main/                       # Figure 1, Figure 2 (PDF + PNG)
│   └── supplementary/              # Supp Fig S1–S5
├── tables/supplementary/           # Tables S1–S12
├── scripts/
│   ├── config.R / config.py        # central paths, sample IDs
│   ├── env/{xenium_env,commot_env}.yml
│   ├── 01_ … 15_ …                 # numbered pipeline scripts
│   ├── format_docx.py              # docx formatter
│   └── run_full_pipeline.sh        # sequential launcher
└── data/
    ├── data_architecture.md
    ├── raw/                        # symlinks; raw data being deposited at GEO (ongoing)
    └── processed/                  # tracked CSVs / panel sources; session snapshot at Zenodo
```

## Quick start — rebuild the published figures (~5 minutes)

Requires the bundle plus `session_snapshot.RData` from the Zenodo deposit (see *Data availability* below).

```bash
conda activate xenium_env       # R 4.5.2 + Giotto, CellChat, LIANA+, ggalluvial, magick, FNN
export IRHEP_BUNDLE_ROOT=$(pwd)
cd scripts
Rscript 04b_panel_a_panel_only.R          # Figure 2 Panel A
Rscript 06_panel_b_sankey.R               # Figure 2 Panel B
Rscript 08_panel_c_rois.R                 # Figure 2 Panel C
Rscript 09_figure_2_assembled.R           # Figure 2 composite
Rscript 03_figure_1_main.R                # Figure 1
```

## Full pipeline reproduction (~24–36 h)

End-to-end re-run from raw Xenium output. Dominant cost: COMMOT spatial optimal transport (~19 h on the largest sample).

```bash
bash scripts/run_full_pipeline.sh
```

| Stage | Script(s) | Approx. runtime |
|---|---|---|
| 1. Giotto normalisation, clustering, spatial networks | `01_giotto_normalize_cluster.R`, `01a_add_spatial_networks.R`, `02_cluster_config.R` | 1–2 h |
| 2. Figure 1 build + H&E overviews | `03_figure_1_main.R`, `03b_sample_he_overviews.py` | minutes |
| 3. Panel A proximity (1,000 perms) | `04_panel_a_proximity_compute.R`, `04a_panel_a_proximity_heatmap.R`, `04b_panel_a_panel_only.R` | 30 min |
| 4. Panel B CellChat per-sample | `05_panel_b_cellchat.R` | 40 min |
| 5. Panel B LIANA+ per-sample | `05a_panel_b_liana.R` | 5 min |
| 6. Panel B COMMOT per-sample | `05b_panel_b_commot.py` | **~19 h on 54023_b** |
| 7. COMMOT aggregator + cross-reference | `05c_panel_b_commot_aggregate.R`, `05d_panel_b_commot_crossref.R` | 5 min |
| 8. Panel B Sankey | `06_panel_b_sankey.R` | <1 min |
| 9. Panel C H&E extract | `07_panel_c_he_extract.R`, `07a_panel_c_he_extract.py` | 5 min |
| 10. Panel C ROI build + auxiliaries | `08_panel_c_rois.R`, `08a_panel_c_score_distribution.R`, `08b_panel_c_close_apposition_counts.R` | <1 min |
| 11. Figure 2 assembly | `09_figure_2_assembled.R` | <1 min |
| 12. Supplementary analyses | `10_*`, `11_*`, `13_*`–`15_*` (H&E lens ROIs, proximity permutation, cluster-size sensitivity, BH FDR, Giotto/TRM/IFN extras) | ~1 h |

## Software

Two conda environments under `scripts/env/`:
- `xenium_env.yml` — R 4.5.2 + Giotto 4.2.2 + GiottoClass + CellChat 2.x + LIANA+ + ggalluvial + ggplot2 + cowplot + magick + FNN + arrow + data.table.
- `commot_env.yml` — Python 3.10+ + commot 0.0.3 + anndata + tifffile + scipy.ndimage + PIL + scanpy.

## Data availability

| Tier | Content | Location | Status |
|---|---|---|---|
| Raw | Xenium output per sample (TIFF pyramids, transcripts, cell boundaries, gene panel) — ~13 GB/sample × 4 samples | GEO accession | deposition ongoing |
| Processed (large) | Giotto session snapshot (~580 MB), per-sample CellChat / LIANA+ / COMMOT outputs | Zenodo DOI | deposition ongoing |
| Processed (small) | Aggregated CSV tables, H&E overview PNGs, supplementary figure source PDFs | this repo (`data/processed/`, `figures/`) | live |

`.gitignore` excludes `data/raw/`, `data/processed/session_snapshot.RData`, and other large `.rds`/`.RData` binaries.

## License

- Code: MIT (`LICENSE`)
- Manuscript text and figures: CC-BY-4.0 (`LICENSE-CONTENT`)
