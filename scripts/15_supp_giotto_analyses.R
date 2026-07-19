#!/usr/bin/env Rscript
# ============================================================================
# 15_supp_giotto_analyses.R
# Two Giotto-object-dependent supplementary analyses:
#
#   Table S10  :  Per-cluster expression of canonical CD8 exhaustion and
#                 tissue-resident-memory markers (those present on the
#                 Xenium 380-gene Immuno-Oncology Panel)
#   Fig S3     :  Clustering stability across Leiden resolutions
#                 (ARI vs chosen + per-cluster silhouette + res-0.3->0.5 splits)
#
# CD69 was audited as absent from the panel and is excluded from Table S10.
# FOXP3/Treg spatial analysis is not included here and is not reported.
# ============================================================================
set.seed(1234)  # belt-and-braces; Giotto stochastic fns also default to 1234

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
t0 <- Sys.time()
gobj <- readRDS(file.path(PROCESSED_DIR, "objects", "gobj_final.rds"))
cat(sprintf("Loaded in %.1f s\n", as.numeric(difftime(Sys.time(), t0, units = "s"))))

meta <- pDataDT(gobj)
expr <- getExpression(gobj, values = "normalized", output = "exprObj")
expr_mat <- expr[]  # dgCMatrix: genes x cells

cluster_char <- fread(file.path(PROCESSED_DIR, "clustering",
                                "cluster_characterization.csv"))
id2name <- setNames(cluster_char$cluster_name, as.character(cluster_char$leiden))
excluded <- c("12", "15", "21")

meta[, leiden_chr := as.character(leiden)]
meta[, cluster_name := id2name[leiden_chr]]
active_meta <- meta[!(leiden_chr %in% excluded)]
clusters_active <- sort(unique(active_meta$cluster_name))

# ============================================================================
# TABLE S10 -- per-cluster exhaustion / TRM marker expression
# ============================================================================
cat("\n[1/2] Table S10 -- per-cluster exhaustion / TRM marker expression\n")

exhaust_markers <- c("PDCD1", "LAG3", "TIGIT", "HAVCR2", "TOX", "EOMES",
                     "ENTPD1", "PDCD1LG2")
trm_markers <- c("ITGAE", "CXCR6", "RUNX3")  # CD69 absent from panel
all_markers <- c(exhaust_markers, trm_markers)

markers_present <- intersect(all_markers, rownames(expr_mat))
markers_absent  <- setdiff(all_markers, rownames(expr_mat))
if (length(markers_absent) > 0) {
  cat(sprintf("  Markers absent from panel: %s\n",
              paste(markers_absent, collapse = ", ")))
}

# Per-cluster mean log-expression and pct cells positive.
# IMPORTANT: avoid name collision between the iteration variable and
# the data.table column also called `cluster_name`.
build_marker_row <- function(marker, this_cluster) {
  cells <- active_meta[cluster_name == this_cluster, cell_ID]
  vals  <- expr_mat[marker, cells]
  data.table(
    cluster   = this_cluster,
    marker    = marker,
    mean_logexpr_per_cluster = mean(vals),
    pct_cells_pos = 100 * mean(vals > 0)
  )
}

table_s10 <- rbindlist(lapply(markers_present, function(m) {
  rbindlist(lapply(clusters_active, function(c) build_marker_row(m, c)))
}))
table_s10[, marker_class := ifelse(marker %in% exhaust_markers,
                                    "exhaustion", "TRM")]
setcolorder(table_s10, c("cluster", "marker", "marker_class",
                          "mean_logexpr_per_cluster", "pct_cells_pos"))

fwrite(table_s10, file.path(OUT_TABLES, "Table_S11_exhaustion_trm_markers.csv"))
cat(sprintf("  Wrote Table_S10 (%d rows).\n", nrow(table_s10)))

# Short readable summary: top-3 clusters by mean expression per marker.
top3_summary <- table_s10[
  order(-mean_logexpr_per_cluster),
  .SD[1:3], by = marker
]
cat("\n  Top-3 expressing clusters per marker:\n")
print(top3_summary)

# ============================================================================
# FIG S3 -- Clustering stability across Leiden resolutions
# ============================================================================
cat("\n[2/2] Fig S3 -- Clustering stability across Leiden resolutions\n")

if (!"sNN.harmony" %in% list_nearest_networks(gobj)$name) {
  gobj <- createNearestNetwork(gobj,
                               dim_reduction_to_use = "harmony",
                               dim_reduction_name = "harmony",
                               dimensions_to_use = 1:15,
                               k = 15, name = "sNN.harmony")
}

resolutions <- c(0.3, 0.5, 0.7, 1.0)
for (r in resolutions) {
  col_name <- sprintf("leiden_r%s", gsub("\\.", "", as.character(r)))
  if (col_name %in% names(pDataDT(gobj))) next
  gobj <- doLeidenCluster(gobj,
                          network_name = "sNN.harmony",
                          resolution = r,
                          n_iterations = 200,
                          name = col_name)
}
m <- pDataDT(gobj)
keep_cols <- c("cell_ID", "leiden",
               sprintf("leiden_r%s", gsub("\\.", "",
                                          as.character(resolutions))))
m <- m[, keep_cols, with = FALSE]

ari <- function(a, b) {
  tab <- table(a, b)
  n  <- sum(tab)
  ai <- rowSums(tab); bj <- colSums(tab)
  sum_choose <- function(x) sum(choose(x, 2))
  idx <- sum_choose(as.vector(tab))
  expected <- sum_choose(ai) * sum_choose(bj) / choose(n, 2)
  max_idx  <- 0.5 * (sum_choose(ai) + sum_choose(bj))
  if (max_idx == expected) return(0)
  (idx - expected) / (max_idx - expected)
}

ari_dt <- data.table(
  resolution = c(0.3, 0.5, 0.7, 1.0),
  ari_vs_chosen = c(
    ari(m$leiden_r03, m$leiden),
    1.0,
    ari(m$leiden_r07, m$leiden),
    ari(m$leiden_r1,  m$leiden)
  ),
  n_clusters = c(
    length(unique(m$leiden_r03)),
    length(unique(m$leiden)),
    length(unique(m$leiden_r07)),
    length(unique(m$leiden_r1))
  )
)
cat("ARI vs chosen (resolution 0.5):\n")
print(ari_dt)
fwrite(ari_dt, file.path(OUT_TABLES, "Table_S12_clustering_stability_ARI.csv"))

pA <- ggplot(ari_dt, aes(x = factor(resolution), y = ari_vs_chosen,
                          fill = factor(resolution))) +
  geom_col(width = 0.6, alpha = 0.7) +
  geom_text(aes(label = sprintf("%.2f\n(k=%d)", ari_vs_chosen, n_clusters)),
            vjust = -0.3, size = 2.6) +
  scale_y_continuous(limits = c(0, 1.05),
                     expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Leiden resolution", y = "ARI vs chosen (0.5)",
       title = "A. Adjusted Rand index vs chosen resolution") +
  theme_minimal(base_size = 9) +
  theme(legend.position = "none")

# Per-cell silhouette over the Harmony embedding at resolution 0.5.
suppressPackageStartupMessages(library(cluster))
embed <- getDimReduction(gobj, reduction_method = "harmony",
                         name = "harmony", output = "dimObj")
embed_mat <- embed[]
if (nrow(embed_mat) > 30000) {
  set.seed(1234)
  s_idx <- sample(seq_len(nrow(embed_mat)), 30000)
} else {
  s_idx <- seq_len(nrow(embed_mat))
}
embed_sub <- embed_mat[s_idx, 1:15]
lab_sub <- as.integer(m$leiden[s_idx])
sil <- silhouette(lab_sub, dist(embed_sub))
sil_dt <- data.table(cluster = lab_sub, silhouette = sil[, 3])
med_sil <- sil_dt[, .(median_silhouette = median(silhouette),
                      n_cells_sampled = .N), by = cluster][order(cluster)]
fwrite(med_sil,
       file.path(OUT_TABLES, "Table_S12b_per_cluster_silhouette.csv"))

pB <- ggplot(med_sil,
             aes(x = factor(cluster), y = median_silhouette,
                 fill = median_silhouette > 0)) +
  geom_col(width = 0.6, alpha = 0.7) +
  geom_hline(yintercept = 0, colour = "grey50", linewidth = 0.4) +
  scale_fill_manual(values = c("FALSE" = "#E76F51", "TRUE" = "#2A9D8F"),
                    guide = "none") +
  labs(x = "Leiden cluster (res = 0.5)",
       y = "Median per-cell silhouette",
       title = "B. Median silhouette per cluster (n=30,000 sampled)") +
  theme_minimal(base_size = 9)

cross_dt <- m[, .N, by = .(leiden_r03, leiden, leiden_r07)]
fwrite(cross_dt, file.path(OUT_TABLES,
                           "Table_S12c_cluster_correspondence.csv"))
pC_dat <- cross_dt[, .(n_pairs_03_to_05 = uniqueN(leiden)),
                    by = leiden_r03][order(leiden_r03)]
pC <- ggplot(pC_dat,
             aes(x = factor(leiden_r03), y = n_pairs_03_to_05)) +
  geom_col(width = 0.6, alpha = 0.7, fill = "#3C5488") +
  labs(x = "Leiden cluster at resolution 0.3",
       y = "Number of res-0.5 sub-clusters this maps to",
       title = "C. Resolution 0.3 -> 0.5 splits") +
  theme_minimal(base_size = 9)

fig_s3 <- (pA | pB) / pC +
  plot_layout(heights = c(1.4, 1)) +
  plot_annotation(title = "Supplementary Figure S3. Clustering stability across resolutions")

ggsave(file.path(OUT_FIGS, "Figure_S3_clustering_stability.pdf"),
       fig_s3, width = 10, height = 7)
ggsave(file.path(OUT_FIGS, "Figure_S3_clustering_stability.png"),
       fig_s3, width = 10, height = 7, dpi = 200)
cat("  Wrote Figure_S3_clustering_stability.{pdf,png}\n")

cat("\nDone.\n")
