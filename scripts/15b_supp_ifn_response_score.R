#!/usr/bin/env Rscript
# ============================================================================
# Cholangiocyte IFN-response composite score (Supplementary).
#
# Computes a per-cell IFN-response composite score across the five canonical
# IFN-stimulated genes available on the 380-gene Xenium Immuno-Oncology
# Panel (STAT1, IRF1, IDO1, ISG15, MX1; IFNG itself is NOT on the panel),
# then characterises this score in cholangiocytes specifically:
#   (i)   per-cluster mean IFN-response score (all active clusters);
#   (ii)  cholangiocyte IFN-response distribution by condition
#         (SR-irHep vs AIH);
#   (iii) per-cholangiocyte-cell correlation of IFN-response score with the
#         composite CXCL9/10/11 score, the chemokines whose CXCR3-receptor
#         carriers are the apposed CD8 cells of Figure 2A.
#
# Rationale: addresses the "putatively IFN-gamma-driven" hedge attached to
# the CXCR3 chemokine axis in the main text (see Tokunaga et al., Cancer
# Treat Rev 2018). The IFNG ligand itself is not captured on the panel, so
# the source-identification question cannot be answered here; what can be
# answered is whether the IFN-response transcriptional programme is
# elevated in SR-irHep cholangiocytes and co-localises at the single-cell
# level with the CXCL9/10/11 chemokines.
#
# Output:
#   tables/supplementary/Table_S11c_ifn_response_score.csv
#   figures/supplementary/Figure_S5_ifn_response_score.pdf
#   figures/supplementary/Figure_S5_ifn_response_score.png
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

# Attach condition (irHep / AIH) from sample_id if missing, then
# normalise label variants ("irHepatitis", "ICI") to "irHep" for consistency
if (!"condition" %in% names(meta)) {
  samples_info <- data.table(
    sample_id = c("18321_a", "24774_d", "54023_b", "56784_c"),
    condition = c("irHep",   "irHep",   "AIH",     "AIH")
  )
  meta <- merge(meta, samples_info, by = "sample_id", all.x = TRUE)
}
meta[condition %in% c("irHepatitis", "ICI"), condition := "irHep"]

active_meta <- meta[!(leiden_chr %in% excluded)]
clusters_active <- sort(unique(active_meta$cluster_name))

ifn_genes <- c("STAT1", "IRF1", "IDO1", "ISG15", "MX1")
cxcl_genes <- c("CXCL9", "CXCL10", "CXCL11")

stopifnot(all(ifn_genes  %in% rownames(expr_mat)))
stopifnot(all(cxcl_genes %in% rownames(expr_mat)))
cat(sprintf("IFN-response score genes: %s\n", paste(ifn_genes,  collapse = ", ")))
cat(sprintf("CXCL9-11 chemokine genes: %s\n", paste(cxcl_genes, collapse = ", ")))

# ----------------------------------------------------------------------------
# Per-cell composite scores: mean of log-normalised expression
# ----------------------------------------------------------------------------
cell_ids <- active_meta$cell_ID
ifn_score  <- colMeans(as.matrix(expr_mat[ifn_genes,  cell_ids, drop = FALSE]))
cxcl_score <- colMeans(as.matrix(expr_mat[cxcl_genes, cell_ids, drop = FALSE]))
active_meta[, ifn_score  := ifn_score]
active_meta[, cxcl_score := cxcl_score]

# ----------------------------------------------------------------------------
# Per-cluster mean IFN-response score (all active clusters)
# ----------------------------------------------------------------------------
cluster_summary <- active_meta[, .(
  n_cells          = .N,
  mean_ifn_score   = mean(ifn_score),
  median_ifn_score = median(ifn_score),
  pct_high         = 100 * mean(ifn_score > 0.5)
), by = cluster_name][order(-mean_ifn_score)]

# Per-cluster per-marker breakdown
build_marker_row <- function(marker, this_cluster) {
  cells <- active_meta[cluster_name == this_cluster, cell_ID]
  vals  <- expr_mat[marker, cells]
  data.table(
    cluster       = this_cluster,
    marker        = marker,
    marker_class  = "IFN-response",
    mean_logexpr  = mean(vals),
    pct_cells_pos = 100 * mean(vals > 0)
  )
}
per_marker <- rbindlist(lapply(ifn_genes, function(m) {
  rbindlist(lapply(clusters_active, function(c) build_marker_row(m, c)))
}))

# ----------------------------------------------------------------------------
# Cholangiocyte-focused analyses
# ----------------------------------------------------------------------------
chol_cluster_name <- clusters_active[grepl("Chol|Cholangiocyte", clusters_active,
                                            ignore.case = TRUE)]
cat(sprintf("Cholangiocyte cluster identified as: %s\n",
            paste(chol_cluster_name, collapse = ", ")))
stopifnot(length(chol_cluster_name) >= 1)

chol_meta <- active_meta[cluster_name %in% chol_cluster_name]
cat(sprintf("Cholangiocyte cells: %d total (SR-irHep %d, AIH %d)\n",
            nrow(chol_meta),
            nrow(chol_meta[condition == "irHep"]),
            nrow(chol_meta[condition == "AIH"])))

# Per-condition IFN-response score summary
chol_per_cond <- chol_meta[, .(
  n_cells          = .N,
  mean_ifn_score   = mean(ifn_score),
  median_ifn_score = median(ifn_score),
  mean_cxcl_score  = mean(cxcl_score),
  pct_ifn_high     = 100 * mean(ifn_score  > 0.5),
  pct_cxcl_high    = 100 * mean(cxcl_score > 0.5)
), by = condition]

# Per-cell IFN-response vs CXCL9/10/11 Spearman correlation in cholangiocytes
chol_corr <- chol_meta[, .(
  n_cells        = .N,
  spearman_rho   = cor(ifn_score, cxcl_score, method = "spearman"),
  pearson_r      = cor(ifn_score, cxcl_score, method = "pearson")
), by = condition]
cat("\nCholangiocyte IFN-response vs CXCL9-11 correlation (per condition):\n")
print(chol_corr)

# Wilcoxon test: cholangiocyte IFN-response between conditions
wt <- wilcox.test(
  ifn_score ~ condition,
  data = chol_meta,
  exact = FALSE
)
cat(sprintf("\nCholangiocyte IFN-response score, SR-irHep vs AIH: Wilcoxon p = %.3e\n",
            wt$p.value))

# Per-condition cholangiocyte summary table (printed for log + appended to CSV)
cat("\nCholangiocyte per-condition summary (IFN-response + CXCL9/10/11):\n")
print(chol_per_cond)

# Per-sample cholangiocyte breakdown (coherent-rule check)
chol_per_sample <- chol_meta[, .(
  n_cells          = .N,
  mean_ifn_score   = mean(ifn_score),
  mean_cxcl_score  = mean(cxcl_score),
  spearman_rho     = cor(ifn_score, cxcl_score, method = "spearman")
), by = .(sample_id, condition)][order(condition, sample_id)]
cat("\nCholangiocyte per-sample breakdown (coherent-rule check):\n")
print(chol_per_sample)

# ----------------------------------------------------------------------------
# Assemble output table
# ----------------------------------------------------------------------------
out_table <- merge(per_marker,
                   cluster_summary[, .(cluster_name, n_cells, mean_ifn_score,
                                       pct_high)],
                   by.x = "cluster", by.y = "cluster_name")
setcolorder(out_table, c("cluster", "n_cells", "marker", "marker_class",
                          "mean_logexpr", "pct_cells_pos",
                          "mean_ifn_score", "pct_high"))
setorder(out_table, -mean_ifn_score, cluster, marker)
fwrite(out_table, file.path(OUT_TABLES,
                            "Table_S11c_ifn_response_score.csv"))
cat(sprintf("\n  Wrote Table_S11c_ifn_response_score.csv (%d rows).\n",
            nrow(out_table)))

# Cholangiocyte-specific per-condition + per-sample summary (companion CSV)
fwrite(chol_per_cond, file.path(OUT_TABLES,
                                "Table_S11c2_ifn_response_cholangiocyte_per_condition.csv"))
fwrite(chol_per_sample, file.path(OUT_TABLES,
                                "Table_S11c3_ifn_response_cholangiocyte_per_sample.csv"))
cat("  Wrote Table_S11c2 (per-condition) + Table_S11c3 (per-sample).\n")

# Per-cluster top-12 console summary
cat("\nPer-cluster mean IFN-response composite score (top 12):\n")
print(cluster_summary[1:12])

# ----------------------------------------------------------------------------
# Figure S5: three-panel layout
#   (A) per-cluster IFN-response barplot (cholangiocyte highlighted)
#   (B) cholangiocyte IFN-response distribution by condition (violin)
#   (C) cholangiocyte per-cell scatter IFN-response vs CXCL9-11, per condition
# ----------------------------------------------------------------------------
p_bar <- ggplot(cluster_summary,
                aes(x = reorder(cluster_name, mean_ifn_score),
                    y = mean_ifn_score)) +
  geom_col(aes(fill = cluster_name %in% chol_cluster_name), width = 0.7) +
  scale_fill_manual(values = c(`TRUE` = "#D4A017", `FALSE` = "grey70"),
                    guide = "none") +
  coord_flip() +
  labs(x = NULL, y = "Mean IFN-response composite score",
       title = "(A) Per-cluster IFN-response\n(STAT1 + IRF1 + IDO1 + ISG15 + MX1)") +
  theme_minimal(base_size = 9) +
  theme(plot.title = element_text(size = 9, face = "bold"))

chol_meta[, condition := factor(condition, levels = c("AIH", "irHep"),
                                 labels = c("AIH", "SR-irHep"))]
p_violin <- ggplot(chol_meta, aes(x = condition, y = ifn_score,
                                   fill = condition)) +
  geom_violin(width = 0.7, alpha = 0.85, colour = "grey20") +
  geom_boxplot(width = 0.18, outlier.size = 0.2,
               fill = "white", alpha = 0.75) +
  scale_fill_manual(values = c("AIH" = "#3B7DBA", "SR-irHep" = "#B2182B"),
                    guide = "none") +
  labs(x = NULL, y = "IFN-response score",
       title = sprintf("(B) Cholangiocyte\nIFN-response by condition\n(Wilcoxon p = %.2e)",
                       wt$p.value)) +
  theme_minimal(base_size = 9) +
  theme(plot.title = element_text(size = 9, face = "bold"))

p_scatter <- ggplot(chol_meta,
                    aes(x = ifn_score, y = cxcl_score, colour = condition)) +
  geom_point(alpha = 0.4, size = 0.4) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.5) +
  scale_colour_manual(values = c("AIH" = "#3B7DBA", "SR-irHep" = "#B2182B")) +
  labs(x = "IFN-response score",
       y = "CXCL9/10/11 score",
       colour = NULL,
       title = "(C) Cholangiocyte per-cell co-expression\nIFN-response vs CXCL9/10/11") +
  theme_minimal(base_size = 9) +
  theme(plot.title = element_text(size = 9, face = "bold"),
        legend.position = "top")

fig <- p_bar | p_violin | p_scatter
ggsave(file.path(OUT_FIGS, "Figure_S5_ifn_response_score.pdf"),
       fig, width = 13, height = 5.2, device = cairo_pdf, bg = "white")
ggsave(file.path(OUT_FIGS, "Figure_S5_ifn_response_score.png"),
       fig, width = 13, height = 5.2, dpi = 320, bg = "white")
cat("  Wrote Figure_S5_ifn_response_score.{pdf,png}\n")
