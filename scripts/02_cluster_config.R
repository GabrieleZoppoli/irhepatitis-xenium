#!/usr/bin/env Rscript
# ============================================================================
# 02_cluster_config.R — Central cluster naming, colour palette, and lineage
# groupings. Sourced by every downstream analysis and figure script. Provides
# `cluster_config` (named list), `cell_type_labels`, `cluster_colors`,
# `lineage_groups`, `cluster_order`, `condition_colors`, `sample_colors`.
# Excludes cluster 21 (n=4 cells) and clusters 12/15 (segmentation fragments)
# from all downstream analyses.
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

source(file.path(Sys.getenv("IRHEP_BUNDLE_ROOT", here::here()),
                 "scripts", "config.R"))

table_dir <- file.path(PROCESSED_DIR, "clustering")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)


# ============================================================================
# 1. CLUSTER NAME MAPPING (cluster number -> display name)
# ============================================================================
# Cluster 21 excluded (4 cells, all irHep — removed from analysis)
# Clusters 12, 15 excluded based on morphometric assessment:
#   - Cluster 12: segmentation fragments (median area 20 um2, sole marker IGKC)
#   - Cluster 15: hepatocytes in poorly-preserved tissue (markers ARG1/SDC1/MARCO)

cell_type_labels <- c(
  "1"  = "Hepatocyte 1",
  "2"  = "CD8+ T cell 1",
  "3"  = "Macrophage 1",
  "4"  = "Stellate cell",
  "5"  = "Macrophage 2",
  "6"  = "Hepatocyte 2",
  "7"  = "CD4+ T cell",
  "8"  = "Cholangiocyte",
  "9"  = "IgG+ plasmablast / plasma cell",
  "10" = "Hepatocyte 3",
  "11" = "Sinusoidal endothelial cell",
  "12" = "Uncharacterized 1",
  "13" = "IgG+ B cell",
  "14" = "Vascular endothelial cell",
  "15" = "Uncharacterized 2",
  "16" = "IgM+ B cell",
  "17" = "CD8+ T cell 2",
  "18" = "Hepatocyte 4",
  "19" = "CD8+ T cell 3",
  "20" = "Ig+ hepatocyte"
)

# Clusters to exclude
excluded_clusters <- c("21", "12", "15")

# Sequential display numbers for publication (original Leiden → 1-18)
# Internal data keeps original Leiden numbers; display_number is for figures/text only
active_cluster_nums <- sort(as.integer(setdiff(names(cell_type_labels), excluded_clusters)))
display_number <- setNames(as.character(seq_along(active_cluster_nums)),
                           as.character(active_cluster_nums))
# Reverse lookup: display number → original Leiden number
leiden_from_display <- setNames(names(display_number), display_number)


# ============================================================================
# 2. LINEAGE GROUPINGS
# ============================================================================

lineage_groups <- list(
  Hepatocyte           = c("1", "6", "10", "18", "20"),
  `T cell`             = c("2", "17", "19", "7"),
  `B / Plasma cell`    = c("9", "13", "16"),
  Myeloid              = c("3", "5"),
  `Stromal/Vascular`   = c("4", "8", "11", "14")
)

# Lineage label per cluster (for easy lookup)
cluster_lineage <- setNames(
  rep(names(lineage_groups), lengths(lineage_groups)),
  unlist(lineage_groups)
)


# ============================================================================
# 3. DISPLAY ORDER (for figure axes/legends)
# ============================================================================
# Grouped by lineage: Hepatocytes -> T cells -> B cells -> Macrophages ->
# Stromal/Vascular

cluster_order <- unlist(lineage_groups, use.names = FALSE)
cluster_name_order <- cell_type_labels[cluster_order]


# ============================================================================
# 4. COLOR PALETTE
# ============================================================================
# Lineage-coherent palette with clear inter-lineage distinction.
# Hepatocytes = warm reds/oranges, T cells = blues, B cells = purples,
# Macrophages = greens, Stromal/Vascular = mixed (cholangiocytes = gold),
# Uncharacterized = grey.

cluster_colors_by_number <- c(
  # Hepatocytes — warm reds/oranges
  "1"  = "#E64B35",   # Hepatocyte 1 (brick red)
  "6"  = "#FF8C00",   # Hepatocyte 2 (dark orange)
  "10" = "#D4A373",   # Hepatocyte 3 (tan)
  "18" = "#BC4749",   # Hepatocyte 4 (crimson)
  # T cells — blues
  "2"  = "#4DBBD5",   # CD8+ T cell 1 (sky blue)
  "17" = "#1976D2",   # CD8+ T cell 2 (medium blue)
  "19" = "#0D47A1",   # CD8+ T cell 3 (deep blue)
  "7"  = "#7BC8A4",   # CD4+ T cell (mint — distinct from CD8s)
  # B / Plasma cells — purples/magentas
  "9"  = "#9C27B0",   # IgG+ plasmablast / plasma cell (purple)
  "13" = "#E040FB",   # IgG+ B cell (fuchsia)
  "20" = "#E8967A",   # Ig+ hepatocyte (salmon — hepatocyte family)
  "16" = "#5C6BC0",   # IgM+ B cell (indigo)
  # Macrophages — greens
  "3"  = "#00A087",   # Macrophage 1 (teal-green)
  "5"  = "#2E7D32",   # Macrophage 2 (forest green)
  # Stromal / Vascular — distinct per type
  "4"  = "#3C5488",   # Stellate cell (steel blue)
  "8"  = "#FFB300",   # Cholangiocyte (amber/gold — distinct from hepatocytes)
  "11" = "#26A69A",   # Sinusoidal endothelial cell (teal)
  "14" = "#8D6E63",   # Vascular endothelial cell (brown)
  # Uncharacterized — greys
  "12" = "#D0D0D0",   # Uncharacterized 1 (EXCLUDED — segmentation artifacts)
  "15" = "#D0D0D0"    # Uncharacterized 2 (EXCLUDED — poor tissue quality)
)

# Named by cell type label (for ggplot2 scale_*_manual)
cluster_colors <- setNames(
  cluster_colors_by_number[names(cell_type_labels)],
  cell_type_labels
)

# Condition and sample colors
condition_colors <- c("irHepatitis" = "#E64B35", "AIH" = "#4DBBD5")
sample_colors <- c(
  "18321_a" = "#E64B35", "24774_d" = "#F39B7F",
  "54023_b" = "#4DBBD5", "56784_c" = "#3C5488"
)


# ============================================================================
# 5. TOP MARKERS PER CLUSTER (for characterization table)
# ============================================================================

top_markers_list <- list(
  "1"  = "ARG1, CCL16, SDC1, FN1, VEGFA",
  "2"  = "KLRK1, CD8A, CD2, IL2RB, GZMA",
  "3"  = "CD163, CSF1R, MPEG1, FCGR3A, C1QB",
  "4"  = "DCN, IGFBP7, CXCL12, LUM, PDGFRA",
  "5"  = "S100A9, ITGAX, MPEG1, FCGR2A, CD163",
  "6"  = "ARG1, SDC1, C1S, APOE, FN1",
  "7"  = "IL7R, CD2, DGKA, CD3E, LTB",
  "8"  = "CXCL6, EPCAM, SOX9, CXCL1, VEGFA",
  "9"  = "IGHGP, IGHG2, IGHG3, IGHG4, JCHAIN",
  "10" = "ARG1, SDC1, APOE, C1S, CCND1",
  "11" = "FCGR2B, FLT1, CCL14, IL1R1, NOTCH1",
  "12" = "IGKC",
  "13" = "IGHGP, IGHG2, IGHG1, IGKC, IGHG3",
  "14" = "SPARCL1, IGFBP7, FLT1, PLVAP, RGS5",
  "15" = "ARG1, SDC1, C1S, MARCO, APOE",
  "16" = "MS4A1, IGHM, BANK1, TNFRSF13C, FCMR",
  "17" = "CD8A, KLRK1, CD2, CXCL10, CD3E",
  "18" = "MKI67, FN1, CDK1, SDC1, CCL16",
  "19" = "MKI67, CDK1, STMN1, CD8A, TUBA1B",
  "20" = "IGHG1, IGHGP, IGHG2, IGHG3, IGKC"
)


# ============================================================================
# 6. BUILD CLUSTER CHARACTERIZATION TABLE
# ============================================================================
cat("=== Building cluster characterization table ===\n")

comp_file <- file.path(PROCESSED_DIR, "tables", "cluster_composition.csv")
if (file.exists(comp_file)) {
  comp <- fread(file = comp_file)
  comp <- comp[!leiden %in% as.integer(excluded_clusters)]

  # Per-cluster totals
  cluster_totals <- comp[, .(n_cells = sum(N)), by = leiden]

  # Per-cluster per-condition
  cluster_cond <- comp[, .(n = sum(N)), by = .(leiden, condition)]
  cluster_cond_wide <- dcast(cluster_cond, leiden ~ condition, value.var = "n", fill = 0)

  # Ensure both condition columns exist
  if (!"irHepatitis" %in% names(cluster_cond_wide))
    cluster_cond_wide[, irHepatitis := 0L]
  if (!"AIH" %in% names(cluster_cond_wide))
    cluster_cond_wide[, AIH := 0L]

  char_dt <- merge(cluster_totals, cluster_cond_wide, by = "leiden")

  # Add names, lineage, markers
  char_dt[, `:=`(
    cluster_name = cell_type_labels[as.character(leiden)],
    lineage      = cluster_lineage[as.character(leiden)],
    top_markers  = sapply(as.character(leiden), function(x) top_markers_list[[x]])
  )]

  # Proportions
  total_irHep <- char_dt[, sum(irHepatitis)]
  total_AIH   <- char_dt[, sum(AIH)]
  char_dt[, `:=`(
    irHep_pct = round(irHepatitis / total_irHep * 100, 1),
    AIH_pct   = round(AIH / total_AIH * 100, 1)
  )]

  # Reorder columns
  setcolorder(char_dt, c("leiden", "cluster_name", "lineage", "n_cells",
                          "irHepatitis", "AIH", "irHep_pct", "AIH_pct",
                          "top_markers"))

  # Sort by display order
  char_dt[, display_order := match(as.character(leiden), cluster_order)]
  setorder(char_dt, display_order)
  char_dt[, display_order := NULL]

  # Save
  fwrite(char_dt, file.path(table_dir, "cluster_characterization.csv"))
  cat(sprintf("  Saved cluster_characterization.csv: %d clusters, %s total cells\n",
              nrow(char_dt), format(char_dt[, sum(n_cells)], big.mark = ",")))
  cat(sprintf("  irHep: %s cells | AIH: %s cells\n",
              format(total_irHep, big.mark = ","),
              format(total_AIH, big.mark = ",")))
} else {
  cat("  WARNING: cluster_composition.csv not found, skipping characterization table\n")
}


# ============================================================================
# 7. SAVE CONFIG RDS
# ============================================================================

cluster_config <- list(
  cell_type_labels     = cell_type_labels,
  excluded_clusters    = excluded_clusters,
  cluster_colors       = cluster_colors,
  cluster_colors_by_number = cluster_colors_by_number,
  cluster_order        = cluster_order,
  cluster_name_order   = cluster_name_order,
  lineage_groups       = lineage_groups,
  cluster_lineage      = cluster_lineage,
  condition_colors     = condition_colors,
  sample_colors        = sample_colors,
  top_markers          = top_markers_list,
  n_clusters           = length(cell_type_labels)
)

saveRDS(cluster_config, file.path(PROCESSED_DIR, "cluster_config.rds"))
cat(sprintf("  Saved cluster_config.rds (%d clusters)\n", length(cell_type_labels)))

cat("\n02_cluster_config.R complete.\n")
