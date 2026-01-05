rm(list = ls())
.rs.restartR(clean = T)
.libPaths()
library(reticulate)
reticulate::use_condaenv(condaenv = "ting-crc-to-liver-mets", required = T)

library(Matrix)
library(BPCells)
library(anndataR)
library(nebula)
library(tinyplot)
library(ggplot2)
library(dplyr)

# dir.create(path = "11_plots_outs")
outs <- "11_plots_outs"

## ITER1 HVGs

nfit1 <- qs2::qs_read("05_hvg_selection_iter1_outs/hvg_nfit.qs2")
str(nfit1)

rhdf5::h5ls(file = "crc-to-liver-mets.h5ad")
hvgs <- rhdf5::h5read(file = "crc-to-liver-mets.h5ad", name = "var") |> 
  bind_cols()

ggplot() + 
  geom_point(data = hvgs, mapping = aes(x = log_cpc, y = overdispersion, color = highly_variable), shape = 16) + 
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black")) + 
  ggrepel::geom_text_repel(data = hvgs[hvgs$overdispersion > 25,], 
                           mapping = aes(x = log_cpc, y = overdispersion, label = `_index`)) +
  theme_bw()
# ggsave(filename = file.path(outs, "overdispersions_iter1.pdf"), height = 6, width = 6.5)

fitted <- loess(hvgs$log_cpc_se ~ hvgs$log_cpc, span = 0.3)
ggplot() + 
  geom_point(data = hvgs, mapping = aes(x = log_cpc, y = log_cpc_se), shape = 16) +
  geom_line(mapping = aes(x = fitted$x[order(fitted$x)], y = fitted$fitted[order(fitted$x)]), color = "red", linewidth = 2, linetype = "dashed") +
  theme_bw()
# ggsave(filename = file.path(outs, "mean-se_trend_iter1.pdf"), height = 6, width = 6.5)

## ITER1 scVI
scvi <- reticulate::import("scvi")
ad <- reticulate::import("anndata")
adata <- ad$read_h5ad("crc-to-liver-mets.h5ad")
model1 <- scvi$model$SCVI$load(adata = adata, dir_path = "06_scVI_iter1_outs/scVI_model1")
train_test_results = model1$history["elbo_train"]
train_test_results["elbo_validation"] = model1$history["elbo_validation"]
train_test_results <- dplyr::bind_cols(train_test_results)

# pdf(file = file.path(outs, "iter1_elbow.pdf"), width = 6, height = 6)
plot(x = rownames(train_test_results), y = train_test_results$elbo_train, col = "dodgerblue", xlab = "epoch", ylab = "ELBO statistic")
points(x = rownames(train_test_results), y = train_test_results$elbo_validation, col = "orange")
# dev.off()

## ITER1 clusters
UMAP <- adata$obsm["X_umap"]
CLUSTERS <- adata$obs$initial_cluster
# pdf(file = file.path(outs, "iter1_clusters.pdf"), width = 8, height = 6)
plot_embedding(CLUSTERS, UMAP, rasterize = T)
# dev.off()

## ITER1 outliers
OUTLIERS <- adata$obs$counts_outlier
OUTLIERS_SCUTTLE <- adata$obs$counts_outlier_scuttle
# pdf(file = file.path(outs, "iter1_outliers_original.pdf"), width = 7, height = 6)
plot_embedding(OUTLIERS, UMAP, rasterize = T, labels_discrete = F)
# dev.off()
# pdf(file = file.path(outs, "iter1_outliers_scuttle.pdf"), width = 7, height = 6)
plot_embedding(OUTLIERS_SCUTTLE, UMAP, rasterize = T, labels_discrete = F)
# dev.off()

# pdf(file = file.path(outs, "iter1_outliers_original_xy.pdf"), width = 8, height = 6)
tinyplot::plt(y_centroid ~ x_centroid | counts_outlier, 
              data = adata$obs |> dplyr::filter(core_global == "mets2_B1"), 
              pal = "Classic Tableau",
              pch = 16, cex = 0.25, legend = legend(pt.cex = 2), asp = 1)
# dev.off()
# pdf(file = file.path(outs, "iter1_outliers_scuttle_xy.pdf"), width = 8, height = 6)
tinyplot::plt(y_centroid ~ x_centroid | counts_outlier_scuttle, 
              data = adata$obs |> dplyr::filter(core_global == "mets2_B1"), 
              pal = "Classic Tableau",
              pch = 16, cex = 0.25, legend = legend(pt.cex = 2), asp = 1)
# dev.off()

## ITER2 clusters
adata <- ad$read_h5ad("crc-to-liver-mets_filtered.h5ad")
UMAP <- adata$obsm["X_umap"]
CLUSTERS <- adata$obs$final_cluster
# pdf(file = file.path(outs, "iter2_clusters.pdf"), width = 8, height = 6)
plot_embedding(CLUSTERS, UMAP, rasterize = T)
# dev.off()

## ITER2 markers
norm <- Matrix::t(adata$layers["lognorm"]) |> as("CsparseMatrix")
rownames(norm) <- (adata$var |> rownames())
# pdf(file = file.path(outs, "iter2_markers_p1.pdf"), width = 10, height = 5.5)
plot_embedding(
  source = norm,
  embedding = UMAP,
  features = c("EPCAM", "EEF1G",
               "FN1",
               "CD68", 
               "PECAM1", 
               "SERPINA1"),
  rasterize = T, 
  colors_continuous = viridis::viridis(n = 71)
)
# dev.off()
# pdf(file = file.path(outs, "iter2_markers_p2.pdf"), width = 10, height = 5.5)
plot_embedding(
  source = norm,
  embedding = UMAP,
  features = c("CD4", "CD8A", "TCF7", "CD3E"), 
  rasterize = T,
  colors_continuous = viridis::viridis(n = 71)
)
# dev.off()

## ITER2 clusters in XY
# pdf(file = file.path(outs, "iter2_clusters_xy_mets2_B1.pdf"), width = 8, height = 6)
tinyplot::plt(y_centroid ~ x_centroid | final_cluster, 
              data = adata$obs |> dplyr::filter(core_global == "mets2_B1"), 
              pal = "Polychrome 36",
              pch = 16, cex = 0.25, legend = legend(pt.cex = 2), asp = 1)
# dev.off()
# pdf(file = file.path(outs, "iter2_clusters_xy_prim2_C1.pdf"), width = 8, height = 6)
tinyplot::plt(y_centroid ~ x_centroid | final_cluster, 
              data = adata$obs |> dplyr::filter(core_global == "prim2_C1"), 
              pal = "Polychrome 36",
              pch = 16, cex = 0.25, legend = legend(pt.cex = 2), asp = 1)
# dev.off()

## ITER2 segmentation example
fp <- "../../ting-crc-to-liver-mets-data/prim2/stow/cell_boundaries.parquet"
seg <- arrow::read_parquet(file = fp)
tmp <- inner_join(x = adata$obs |> dplyr::filter(core_global == "prim2_C1"), y = seg, by = "cell_id")
ggplot() + 
  geom_polygon(data = tmp, mapping = aes(x = vertex_x, y = vertex_y, group = cell_id, fill = final_cluster), color = "white", linewidth = 0.05) + 
  coord_fixed() + 
  #ggthemes::scale_fill_tableau("Classic 20") +
  ggthemes::scale_fill_stata() +
  scale_x_continuous(limits = c(1500, 2750), expand = c(0, 0)) +
  scale_y_continuous(limits = c(7500, 9000), expand = c(0, 0)) +
  theme_void() + 
  Seurat::NoLegend()
# ggsave(filename = file.path(outs, "iter2_clusters_xy_prim2_C1_zoomed_seg.svg"), height = 8, width = 8)

