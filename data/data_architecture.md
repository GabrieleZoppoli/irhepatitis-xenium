# Data architecture

## Raw data (`data/raw/`)

Symlinked from the working directory `data/10x data/first run - IO panel/data/catalyst_release_CAT097_GZ_UG_Sep27/`. Total raw size ~63 GB (~13 GB per sample, dominated by H&E and DAPI TIFF pyramids). Symlinks preserve the layout; for archival the four sample directories are being deposited at the Gene Expression Omnibus (GEO), with the accession cited in the manuscript once assigned.

Four samples, two per condition:

| Sample | Condition | Cells (n) |
|---|---|---|
| `0027420_18321_a` | irHepatitis (checkpoint-inhibitor induced) | ~20,600 |
| `0027420_24774_d` | irHepatitis | ~18,600 |
| `0027420_54023_b` | autoimmune hepatitis (AIH) | ~52,300 |
| `0027420_56784_c` | AIH | ~54,500 |

### Per-sample directory contents (10x Xenium standard output)

```
0027420_<sample>/
├── 0027420_<sample>.ome.tif         # DAPI pyramid (cell-segmentation reference)
├── <sample>_he_image.ome.tif        # post-Xenium H&E pyramid
├── <sample>_he_alignment_files.zip  # vendor H&E↔Xenium alignment artefacts
├── he_alignment/
│   ├── matrix.csv                   # 3×3 affine, level-0 H&E pixels → Xenium pixels
│   └── keypoints.csv
├── morphology.ome.tif               # cellular morphology imaging
├── morphology_focus/                # focus-stack metadata
├── analysis/                        # 10x default cluster output
├── analysis.zarr.zip
├── analysis_summary.html            # 10x QC summary
├── aux_outputs/
├── cell_boundaries.parquet          # cell polygon vertices (used in Panel C)
├── cell_boundaries.csv.gz
├── nucleus_boundaries.parquet
├── nucleus_boundaries.csv.gz
├── cell_feature_matrix.h5           # gene × cell count matrix
├── cell_feature_matrix.zarr.zip
├── cell_feature_matrix/             # MEX-format export
├── cells.parquet                    # per-cell metadata (centroid, area, segmentation_method, etc.)
├── cells.csv.gz
├── cells.zarr.zip
├── transcripts.parquet              # per-transcript record (x, y, z, feature_name, cell_id, qv)
├── transcripts.zarr.zip
├── gene_panel.json                  # 380-gene IO panel definition
├── experiment.xenium                # vendor experiment manifest
├── metrics_summary.csv              # vendor QC metrics
└── giotto_results/                  # local reproducibility checkpoints
```

### 10x Xenium IO Panel

380-gene curated immuno-oncology panel; full gene list in each sample's `gene_panel.json`. Approximately 1,000 gene-panel control probes are included as negative-control codewords for QC.

## Processed data (`data/processed/`)

Outputs of analysis pipeline that downstream scripts depend on.

| File | Source script | Purpose |
|---|---|---|
| `session_snapshot.RData` | `scripts/02_giotto_processing_clustering.R` | Full Giotto object: log-norm, HVG, PCA, UMAP, Leiden clusters (k=20). Loaded by every Fig 2 script. **629 MB.** |
| `cluster_config.rds` | `scripts/03_cluster_config.R` | Stable cluster identity + colour palette mapping (cluster id → name). |
| `clustering/cluster_markers.csv` | Giotto `findMarkers_one_vs_all` | Per-cluster top differentially expressed genes (logFC, FDR). |
| `clustering/cluster_characterization.csv` | `02_giotto_processing_clustering.R` | Per-cluster cell counts per sample × condition. |
| `cellchat/fig2_b_coherent_interactions_v2.csv` | `04_panel_b_cellchat_per_sample.R` | 47 coherent ligand-receptor flows (CellChat v2, 25 µm contact range, ligand-split families). |
| `cellchat/fig2_b_ss10_vs_v2_comparison.csv` | sensitivity analysis | CellChat 100 µm vs 25 µm contact-range comparison. |
| `liana/fig2_b_liana_coherent.csv` | `04b_panel_b_liana.R` | LIANA+ 5-method consensus coherent flows (cross-method validation). |
| `commot/commot_*.csv` | `04c_panel_b_commot_per_sample.py` + aggregator | COMMOT (spatial-aware optimal-transport CCC) per-sample and differential outputs; cross-reference with CellChat. |
| `scan_tables/scan_*.csv` | `06_panel_c_roi_scans.R` | Sliding-window ranked tables used to select Panel C ROIs (anatomy × T-cell density per sample). |

## Software environments

- **R**: `xenium_env` conda env with R 4.5.2, Giotto 4.2.2, GiottoClass, CellChat 2.x, LIANA+, ggalluvial, ggplot2, cowplot, magick, FNN, arrow, data.table.
- **Python**: `commot_env` conda env (spec at `scripts/env/commot_env.yml`) for COMMOT and image extraction (tifffile, scipy.ndimage, PIL, commot).

## Reproduction order

1. Raw Xenium output (`data/raw/`) → `02_giotto_processing_clustering.R` → `session_snapshot.RData`
2. `session_snapshot.RData` → Panel A (proximity), Panel B (CellChat / LIANA / COMMOT), Panel C (ROIs)
3. Final figures assembled by `09_figure_2_assembled.R`

See `scripts/README.md` for the full script-by-script pipeline.
