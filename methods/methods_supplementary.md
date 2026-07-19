# Supplementary methods

## Cohort framing

The cohort comprises four FFPE liver biopsies: two from patients with steroid-resistant ICI hepatitis (SR-irHep; sample identifiers 18321_a, 24774_d) following PD-1 blockade and second-line corticosteroid failure, and two from patients with treatment-naïve autoimmune hepatitis (AIH) type I (sample identifiers 54023_b, 56784_c). Throughout the manuscript SR-irHep is the asymmetric subject and AIH is the autoimmune comparator. Pre-ICI baseline liver biopsies were not available, a structural limitation of irAE studies of this nature.

## Cross-method validation of cell–cell communication

CellChat v2 inferences (Figure 2B) were cross-validated by two independent computational approaches operating on the same per-sample data.

**LIANA+ consensus.** Per-sample expression-only communication scoring was performed with LIANA+ (v1.x) using five orthogonal scoring methods — NATMI, Connectome, logFC, SingleCellSignalR (SCA) and CellPhoneDB — against the LIANA Consensus database with Robust Rank Aggregation. LIANA+ is not spatially aware and tests whether each CellChat-coherent flow is recoverable from gene-expression statistics alone. Of the 13 testable CellChat-coherent flows (ligand and receptor both with non-zero expression in the relevant clusters), 8 were confirmed by LIANA+ (Robust Rank Aggregation rank in the top decile of the corresponding cell-pair distribution).

**COMMOT spatial optimal transport.** Per-sample spatial-aware CCC was computed with COMMOT (v0.0.3) using the CellChatDB ligand–receptor database, species = "human", cell-cell distance threshold `dis_thr` = 250 µm and heteromeric receptor handling = TRUE. COMMOT was run on each sample individually (the largest sample, 54023_b, with 52,266 cells × 33 testable L–R pairs, required ≈ 19 h on a 32-core / 376 GB workstation). Per-sample COMMOT scores were aggregated to per-condition differentials with explicit three-layer handling of structural missingness:
1. L–R pair not testable in a sample because the per-cell expression filter (`min_cell_pct` = 0.05) failed in that sample → NA, not zero.
2. Source or target cluster present at fewer than 10 cells in a sample → NA.
3. Pair tested, both endpoints present, OT flow returned zero → 0.

Eleven L–R pairs were strictly comparable across all four samples; twelve additional pairs were AIH-only-testable and reported in Supplementary Table S5 with the presence-based comparison flagged.

**Cross-method agreement.** Of 17 L–R × direction triples that were testable by both CellChat and COMMOT, 9 agreed in direction; the 8 disagreements clustered on FN1 and CSF1 axes, where the longer COMMOT distance scale (250 µm diffusion) captures matrix-deposition and cytokine-diffusion patterns that lie outside the 25 µm CellChat contact range. The triple-method-confirmed core (CellChat ∩ LIANA+ ∩ COMMOT, same direction) comprised 2 of 47 CellChat-coherent triples: CCL3/4/5 → CCR1 from VascEC and SinEC to CD8T2 in SR-irHep, and CCL3/4/5 → CCR1 from SinEC to IgM⁺ B cell in AIH. Supplementary Fig. S1 visualises COMMOT spatial signalling vectors over tissue for the headline axes.

## Per-sample quality control

Per-sample QC metrics — total cells after QC, median nuclear and cell area (µm²), and nucleated fraction — are reported in Supplementary Table S2; the filter parameters used at clustering ingestion (≥ 4 transcripts per cell; gene retained if detected in ≥ 3 cells; QV ≥ 20) are listed alongside. Per-cluster per-sample cell counts (Supplementary Table S4) confirm that all 18 active clusters contain ≥ 28 cells per sample, with the single exception of Ig⁺ hepatocytes (cluster 20: 12, 2, 243 and 478 cells in samples 18321_a, 24774_d, 54023_b and 56784_c respectively); this cluster was excluded from the headline Figure 2B display by the minimum-cluster filter and is interpreted separately in Supplementary Note 1.

## Region-of-interest selection

The four ROIs displayed in Figure 2C were selected by pathology review of the four tissue sections as representative of each condition × anatomy axis. Bile-duct ROIs (Figure 2C, top row) were centred on a duct surrounded by lymphocytic infiltrate in SR-irHep (sample 18321_a) and a duct in AIH (sample 54023_b). Parenchymal ROIs (Figure 2C, bottom row) were chosen as representative hepatocyte-dense fields in the same samples (a SR-irHep field with CD4⁺ cells largely outside the hepatocyte plates and an AIH field with CD4⁺ cells in close apposition to hepatocytes). Panel C is illustrative of the cohort-level findings at single-cell tissue scale and is not used as inferential evidence on its own.

## Statistical considerations at n = 2 per condition

The cohort size is the structural limit of this study. The coherent rule used for both proximity (Figure 2A) and CCC (Figure 2B) is a relaxed alternative to requiring both small-n samples to clear an FDR threshold independently, while still preventing within-condition direction conflicts from being counted as hits. Conclusions in the main text are therefore framed as spatially descriptive of the observed cohort and warrant validation in an independent series before generalisation. Figure 2C is illustrative of the cohort-level findings at single-cell tissue scale and is not used as inferential evidence on its own. Beyond the coherent rule, a label-permutation analysis on the Delaunay proximity step (1,000 cluster-label permutations across the four samples; Supplementary Table S7) confirms that the leading differential pairs do not arise by chance.

## Reproducibility

A complete, sequential end-to-end pipeline launcher is provided at `scripts/run_full_pipeline.sh` of the analysis repository. The pipeline is dominated by COMMOT runtime; an end-to-end re-run from raw Xenium output requires approximately one to two days of compute on a modern workstation. All scripts source a centralised path-configuration file (`config.R` for R, `config.py` for Python) that resolves data and output paths from the `IRHEP_BUNDLE_ROOT` environment variable, allowing the bundle to be relocated without script edits. Random seeds are fixed for the Leiden clustering (`01_giotto_normalize_cluster.R`), the proximity-permutation step (`04_panel_a_proximity_compute.R`) and ROI arrow jitter (`08_panel_c_rois.R`); all other steps are deterministic.

## Choice of normalisation, integration and clustering parameters

The choices encoded in `01_giotto_normalize_cluster.R` were arrived at after exploration along several axes; the values shipped are the ones whose downstream cluster stability and biological plausibility were highest, but the published bundle contains only the canonical pipeline:

- **Framework.** Earlier iterations were prototyped in Seurat (versions 4 and 5). The pipeline was migrated to Giotto for native handling of Xenium spatial metadata and cell boundaries; downstream analyses (CellChat with `contact.range`, Giotto proximity enrichment) require Giotto's cell-spatial representation.
- **Integration / batch correction.** Several Giotto-supported strategies were compared (no integration; CCA-style anchor methods; Harmony at multiple parameter regimes). Harmony was retained at `theta = 2`, `sigma = 0.1`, `lambda = 1`; lower `theta` left visible per-sample sub-clusters of common lineages, higher `theta` over-mixed biologically distinct hepatocyte subtypes.
- **Filtering thresholds.** Cell-level QC retained cells with at least 4 transcripts and genes detected in at least 3 cells; `qv ≥ 20` is the Xenium recommended decoding-quality threshold. The chosen `min_transcripts_per_cell = 4` sits just above the dominant low-end mode of the per-cell transcript histogram pooled across the four samples; tighter cutoffs discarded a substantial fraction of biologically informative cholangiocytes and CD4⁺ T cells, while looser cutoffs admitted segmentation-fragment noise.
- **Dim-reduction and UMAP.** 50 PCs were retained on the basis of the elbow of the per-component variance plot; downstream UMAP used the first 15 PCs (`umap_dims = 15`, `umap_neighbors = 15`, `umap_min_dist = 0.0005`). Larger `min_dist` collapsed the proliferating CD8 subcluster into the cytotoxic CD8 subcluster on the UMAP, even though they were distinct under Leiden.
- **Leiden resolution.** Resolutions 0.3, 0.5, 0.7 and 1.0 were inspected. 0.5 produced 21 primary clusters whose marker profiles best matched canonical liver and immune lineages without over-splitting hepatocytes or merging the three CD8 subclusters; the formal stability analysis comparing 0.5 against 0.3, 0.7 and 1.0 is reported in the Supplementary Methods §1.7 of the supplementary materials.

The trial-and-error scripts that produced these comparisons are retained in a separate analysis archive and are not included in the publication bundle to keep the published pipeline concise; their findings inform the parameter choices documented above.

## Cluster naming and lineage assignment

Canonical v3 cluster names, marker genes and lineage assignments are defined in `scripts/02_cluster_config.R` and reproduced in Supplementary Table S3. Cluster names follow the convention: lineage prefix + ordinal subcluster (e.g. "Hepatocyte 1–4", "CD8+ T cell 1–3", "Macrophage 1–2"); cluster colours are lineage-coherent (hepatocytes = warm reds/oranges; T cells = blues; B/plasma = purples/magentas; macrophages = greens; cholangiocytes = gold; sinusoidal and vascular endothelial = teal/brown; stellate = steel blue). The cluster characterisation rationale (canonical markers per cluster vs. published single-cell liver atlases) is documented in `methods/cell_type_classification.md`.

## Supplementary note 1: Ig⁺ hepatocyte cluster

Cluster 20 carries hepatocyte (ARG1, SDC1, A2M, FN1) and plasma-cell (IGHG1, IGHGP, IGKC, JCHAIN) markers and is highly imbalanced across samples (12, 2, 243 and 478 cells in samples 18321_a, 24774_d, 54023_b and 56784_c respectively). The cluster is interpretable as either (i) a hepatocyte–plasmablast doublet artefact in close spatial proximity, or (ii) a *bona fide* hepatocyte subset with low-level immunoglobulin co-expression in a chronic-inflammation context. The cohort size and per-sample imbalance preclude resolution between these possibilities; the cluster is therefore excluded from the headline Figure 2B display by the minimum-cluster filter (n ≥ 10 cells per sample) and is reported separately in Supplementary Table S4 with the caveat made explicit. The phenotype is referenced briefly in the main-text Discussion only as a follow-up direction for a deeper-transcriptome panel or scRNA-seq of dissociated cells.
