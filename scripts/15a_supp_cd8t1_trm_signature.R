#!/usr/bin/env Rscript
# ============================================================================
# Multi-marker TRM signature for CD8T1 vs other CD8 sub-clusters.
#
# Computes per-cluster expression and a composite TRM score across the five
# canonical tissue-resident-memory markers available on the Xenium 380-gene
# Human Immuno-Oncology Panel: ITGAE (CD103), ITGA1 (CD49a), ZNF683 (HOBIT),
# CXCR6 and RUNX3. CD69, the most commonly cited TRM marker, is not on the
# panel.
#
# Output:
#   tables/supplementary/Table_S11b_cd8t1_trm_signature.csv
#   figures/supplementary/Figure_S4_cd8t1_trm_signature.pdf
# ============================================================================
set.seed(1234)

suppressPackageStartupMessages({
  library(Giotto)
  library(GiottoClass)
  library(data.table)
  library(Matrix)
  library(ggplot2)
  library(patchwork)
})

source(file.path(Sys.getenv("IRHEP_BUNDLE_ROOT", here::here()),
                 "scripts", "config.R"))

OUT_TABLES <- file.path(BUNDLE_ROOT, "tables", "supplementary")
OUT_FIGS   <- file.path(BUNDLE_ROOT, "figures", "supplementary")
dir.create(OUT_TABLES, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_FIGS,   recursive = TRUE, showWarnings = FALSE)

cat("Loading Giotto object ...\n")
gobj <- readRDS(file.path(PROCESSED_DIR, "objects", "gobj_final.rds"))
meta <- pDataDT(gobj)
expr <- getExpression(gobj, values = "normalized", output = "exprObj")
expr_mat <- expr[]

cluster_char <- fread(file.path(PROCESSED_DIR, "clustering",
                                "cluster_characterization.csv"))
id2name <- setNames(cluster_char$cluster_name, as.character(cluster_char$leiden))
excluded <- c("12", "15", "21")
meta[, leiden_chr := as.character(leiden)]
meta[, cluster_name := id2name[leiden_chr]]
active_meta <- meta[!(leiden_chr %in% excluded)]
clusters_active <- sort(unique(active_meta$cluster_name))

trm_markers <- c("ITGAE", "ITGA1", "ZNF683", "CXCR6", "RUNX3")
present <- intersect(trm_markers, rownames(expr_mat))
absent  <- setdiff(trm_markers, rownames(expr_mat))
cat(sprintf("TRM markers used: %s\n", paste(present, collapse = ", ")))
if (length(absent) > 0) {
  cat(sprintf("TRM markers absent from panel: %s\n",
              paste(absent, collapse = ", ")))
}

# Per-cell composite TRM score: mean of log-normalised expression across
# the available TRM markers.
cell_ids <- active_meta$cell_ID
m_sub <- expr_mat[present, cell_ids, drop = FALSE]
trm_score <- colMeans(as.matrix(m_sub))
active_meta[, trm_score := trm_score]

# Per-cluster TRM score summary
cluster_summary <- active_meta[, .(
  n_cells        = .N,
  mean_trm_score = mean(trm_score),
  median_trm_score = median(trm_score),
  pct_high       = 100 * mean(trm_score > 0.5)
), by = cluster_name][order(-mean_trm_score)]

# Per-cluster per-marker breakdown
build_marker_row <- function(marker, this_cluster) {
  cells <- active_meta[cluster_name == this_cluster, cell_ID]
  vals  <- expr_mat[marker, cells]
  data.table(
    cluster   = this_cluster,
    marker    = marker,
    mean_logexpr  = mean(vals),
    pct_cells_pos = 100 * mean(vals > 0)
  )
}
per_marker <- rbindlist(lapply(present, function(m) {
  rbindlist(lapply(clusters_active, function(c) build_marker_row(m, c)))
}))

# Combine into a tidy table for output
out_table <- merge(per_marker, cluster_summary[, .(cluster_name, n_cells,
                                                    mean_trm_score,
                                                    pct_high)],
                   by.x = "cluster", by.y = "cluster_name")
out_table[, marker_class := "TRM"]
setcolorder(out_table, c("cluster", "n_cells", "marker", "marker_class",
                          "mean_logexpr", "pct_cells_pos",
                          "mean_trm_score", "pct_high"))
setorder(out_table, -mean_trm_score, cluster, marker)

fwrite(out_table, file.path(OUT_TABLES,
                            "Table_S11b_cd8t1_trm_signature.csv"))
cat(sprintf("  Wrote Table_S11b_cd8t1_trm_signature.csv (%d rows).\n",
            nrow(out_table)))

# Print quick summary
cat("\nPer-cluster mean TRM composite score (top 12):\n")
print(cluster_summary[1:12])

# CD8 cluster comparison
cat("\nCD8 sub-cluster comparison:\n")
print(cluster_summary[grepl("CD8", cluster_name)])

# ============================================================================
# Figure: per-cluster TRM score barplot + marker-wise heatmap
# ============================================================================
p_bar <- ggplot(cluster_summary,
                aes(x = reorder(cluster_name, mean_trm_score),
                    y = mean_trm_score)) +
  geom_col(aes(fill = grepl("CD8", cluster_name)), width = 0.7) +
  scale_fill_manual(values = c(`TRUE` = "#3B7DBA", `FALSE` = "grey70"),
                    guide = "none") +
  coord_flip() +
  labs(x = NULL, y = "Mean composite TRM score",
       title = "TRM composite score (ITGAE + ITGA1 + ZNF683 + CXCR6 + RUNX3)") +
  theme_minimal(base_size = 9) +
  theme(plot.title = element_text(size = 9, face = "bold"))

hm_data <- copy(per_marker)
hm_data[, cluster := factor(cluster, levels = cluster_summary$cluster_name)]
hm_data[, marker := factor(marker, levels = present)]

p_hm <- ggplot(hm_data, aes(x = marker, y = cluster, fill = mean_logexpr)) +
  geom_tile(colour = "grey92", linewidth = 0.2) +
  scale_fill_gradient(low = "white", high = "#B2182B",
                      name = "mean\nlog-expr") +
  labs(x = NULL, y = NULL,
       title = "Per-marker mean log-expression") +
  theme_minimal(base_size = 9) +
  theme(plot.title = element_text(size = 9, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1))

fig <- p_bar | p_hm

ggsave(file.path(OUT_FIGS, "Figure_S4_cd8t1_trm_signature.pdf"),
       fig, width = 9, height = 5, device = cairo_pdf, bg = "white")
ggsave(file.path(OUT_FIGS, "Figure_S4_cd8t1_trm_signature.png"),
       fig, width = 9, height = 5, dpi = 320, bg = "white")
cat(sprintf("  Wrote Figure_S4_cd8t1_trm_signature.{pdf,png}\n"))
