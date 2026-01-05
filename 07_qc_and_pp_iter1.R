rm(list = ls())
.rs.restartR(clean = T)
.libPaths()

library(BPCells)
library(Matrix)

## QC & PP ITER 1 --------------------------------------------------------------
# We want to identify the noise clusters and filter them out.
# Next, we will re-run the pipeline without these cells.

# Reading the data
adata = anndataR::read_h5ad(path = "crc-to-liver-mets.h5ad", mode = "r+")

# UMAP using the scVI dimensions
UM <- uwot::umap(X = adata$obsm$X_scVI, n_neighbors = 30, min_dist = 0.1, metric = "cosine", nn_method = "annoy", spread = 1, 
                 fast_sgd = T, n_sgd_threads = 40, verbose = T)

# Quick viz
umemb <- UM |> as.data.frame()
colnames(umemb) <- c("umap_1", "umap_2")
umemb$core_global <- adata$obs$core_global
umemb$sample_type <- adata$obs$sample_type

tinyplot::plt(umap_2 ~ umap_1 | core_global, data = umemb, pal = "Polychrome 36", pch = ".", legend = legend(pt.cex = 9), asp = 1)
tinyplot::plt(umap_2 ~ umap_1 | sample_type, data = umemb, pch = ".", legend = legend(pt.cex = 9), asp = 1)

# Adding to the anndata and saving
adata$obsm$X_umap <- UM
anndataR::write_h5ad(object = adata, path = "crc-to-liver-mets.h5ad", mode = "w")

# Normalization (out-of-memory)
cts <- adata$layers$counts |>  as("CsparseMatrix")
libsizes <- Matrix::rowSums(cts)

# This all happens outside of the R session
norm <- write_matrix_dir(mat = cts, dir = "norm_tmp", overwrite = T)
norm <- convert_matrix_type(matrix = norm, type = "float")
norm <- BPCells::multiply_rows(mat = norm, vec = 1/libsizes)
norm <- log1p(norm*1000)/log(2) # Change of base rule!
colnames(norm) <- adata$var_names
rownames(norm) <- adata$obs_names
norm <- write_matrix_dir(mat = norm, dir = "norm_tmp", overwrite = T)

# Then I bring the data back into R and convert it to the expected format
norm <- BPCells::open_matrix_dir(dir = "norm_tmp")
norm <- BPCells::write_matrix_memory(mat = norm, compress = F)
norm <- as(norm, "dgCMatrix")

# NOTE: We can do this very efficiently in-memory as follows: 
# scaling_factor <- 1000
# norm_factors <- Matrix::Diagonal(x = scaling_factor/adata$obs$transcript_counts, names=rownames(adata$layers$counts))
# norm <- ((norm_factors %*% adata$layers$counts) |> log1p())/log(2)

# I stole this trick from NanoString and added the change of base part. Knowing
# how to use BP cells is good for other things though, such as z-score scaling.

# Adding to the anndata and saving
adata$layers$lognorm <- norm
anndataR::write_h5ad(object = adata, path = "crc-to-liver-mets.h5ad", mode = "w")

# Saving space
remove(norm)
remove(cts)
unlink("norm_tmp", recursive = T)

# Checking out the results with the normalized data
plot_embedding(
  source = Matrix::t(adata$layers$lognorm) |> as("CsparseMatrix") |> magrittr::set_rownames(value = adata$var_names),
  embedding = adata$obsm$X_umap,
  features = c("EPCAM", "EEF1G",
               "FN1",
               "CD68", 
               "PECAM1", 
               "SERPINA1"),
  rasterize = T, 
  colors_continuous = viridis::viridis(n = 71)
)
plot_embedding(
  source = Matrix::t(adata$layers$lognorm) |> as("CsparseMatrix") |> magrittr::set_rownames(value = adata$var_names),
  embedding = adata$obsm$X_umap,
  features = c("CD4", "CD8A", "TCF7", "CD3E"), 
  rasterize = T,
  colors_continuous = viridis::viridis(n = 71)
)

# Louvain clustering
clusts <- knn_hnsw(adata$obsm$X_scVI, k = 30, metric = "cosine", ef = 300) |> # Find approximate nearest neighbors
  knn_to_snn_graph() |> # Convert to a SNN graph
  cluster_graph_louvain(resolution = 0.5) # Perform graph-based clustering

# dir.create("07_qc_and_pp_iter1_outs")
# pdf(file = "07_qc_and_pp_iter1_outs/umap_by_clusters.pdf", width = 6, height = 6)
plot_embedding(clusts, adata$obsm$X_umap, rasterize = T)
# dev.off()

adata$obs$initial_cluster <- clusts

# Finding outliers
(adata$obs$transcript_counts < 15) |> sum()
(adata$obs$transcript_counts > 2000) |> sum()
(adata$obs$nucleus_count > 2) |> sum()
(adata$obs$cell_area < 10) |> sum()
(adata$obs$cell_area > 400) |> sum()

adata$obs$transcript_counts |> hist(breaks = 200)
abline(v = c(15, 2000), col = "red")

adata$obs$nucleus_count |> hist(breaks = 20)
abline(v = 2, col = "red")

adata$obs$cell_area |> hist(breaks = 200)
abline(v = c(10, 400), col = "red")

adata$obs$counts_outlier <- ifelse((adata$obs$transcript_counts < 15) | (adata$obs$transcript_counts > 2000), yes = T, no = F)
adata$obs$nuclei_outlier <- ifelse(adata$obs$nucleus_count > 2, yes = T, no = F)
adata$obs$area_outlier <- ifelse((adata$obs$cell_area < 10) | (adata$obs$cell_area > 400), yes = T, no = F)

# Saving the anndata object again
anndataR::write_h5ad(object = adata, path = "crc-to-liver-mets.h5ad", mode = "w")

# Investigating for QC
plot_embedding(adata$obs$counts_outlier, adata$obsm$X_umap, rasterize = T, labels_discrete = F) # Seems too strict!
plot_embedding(adata$obs$nuclei_outlier, adata$obsm$X_umap, rasterize = T, labels_discrete = F) # Seems okay
plot_embedding(adata$obs$area_outlier, adata$obsm$X_umap, rasterize = T, labels_discrete = F) # May be too strict...

(table(clusts) < 100) |> which()

umemb <- adata$obsm$X_umap |> as.data.frame()
colnames(umemb) <- c("umap_1", "umap_2")
umemb$transcript_counts <- adata$obs$transcript_counts

tinyplot::plt(y_centroid ~ x_centroid | initial_cluster, 
              data = adata$obs |> dplyr::filter(slide == "prim2"), 
              pal = "Polychrome 36",
              pch = ".", legend = legend(pt.cex = 9), asp = 1)
tinyplot::plt(y_centroid ~ x_centroid | initial_cluster, 
              data = adata$obs |> dplyr::filter(core_global == "prim2_C1"), 
              pal = "Polychrome 36",
              pch = 16, cex = 0.25, legend = legend(pt.cex = 2), asp = 1)
tinyplot::plt(y_centroid ~ x_centroid | initial_cluster, 
              data = adata$obs |> dplyr::filter(slide == "mets2"), 
              pal = "Polychrome 36",
              pch = ".", legend = legend(pt.cex = 9), asp = 1)
tinyplot::plt(y_centroid ~ x_centroid | initial_cluster, 
              data = adata$obs |> dplyr::filter(core_global == "mets2_F2"), 
              pal = "Polychrome 36",
              pch = 16, cex = 0.25, legend = legend(pt.cex = 2), asp = 1)
tinyplot::plt(y_centroid ~ x_centroid | initial_cluster, 
              data = adata$obs |> dplyr::filter(core_global == "mets2_C2"), 
              pal = "Polychrome 36",
              pch = 16, cex = 0.25, legend = legend(pt.cex = 2), asp = 1)

# Should we be thinking about the thresholds this way?
par(mar = c(8, 5, 4, 1) + 0.1, las = 2)
boxplot(adata$obs$transcript_counts ~ adata$obs$core_global, 
        pch = 16, cex = 0.5, xlab = NA, ylab = "RNA counts",
        col = palette.colors(palette = "Polychrome 36", n = dplyr::n_distinct(adata$obs$core_global)))
outliers <- table(adata$obs$core_global, adata$obs$counts_outlier) 
sweep(x = outliers, STATS = rowSums(outliers), MARGIN = 1, FUN = "/") |> plot(main = "counts outliers")

# B1 from mets2 has very low counts. Many of the cells would get filtered out... 
# However, many of these cells appear to be real cells that were clustered 
# properly!
tinyplot::plt(y_centroid ~ x_centroid | initial_cluster, 
              data = adata$obs |> dplyr::filter(core_global == "mets2_B1"), 
              pal = "Polychrome 36",
              pch = 16, cex = 0.25, legend = legend(pt.cex = 2), asp = 1)
tinyplot::plt(y_centroid ~ x_centroid | counts_outlier, 
              data = adata$obs |> dplyr::filter(core_global == "mets2_B1"), 
              pal = "Classic Tableau",
              pch = 16, cex = 0.25, legend = legend(pt.cex = 2), asp = 1)

# In most scRNA-seq pipelines nowadays, QC thresholds are set for each sample
# before they are integrated together. Let's apply this principle here for each
# core.

# Using scuttle to get the metrics
library(scuttle)
cts <- adata$layers$counts |> Matrix::t()
dimnames(cts) <- list(adata$var_names, adata$obs_names)
qc_metrics <- scuttle::perCellQCMetrics(x = cts)
summary(qc_metrics$sum)
summary(qc_metrics$detected)

# Finding core-wise (batch-wise) outliers for counts, which would be affected by
# the batch.
all_outliers <- isOutlier(qc_metrics$sum, type = "both", log = T, batch = adata$obs$core_global, 2.5)
summary(all_outliers)
attributes(all_outliers)$thresholds

# Adding to the AnnData
adata$obs$counts_outlier_scuttle <- as.vector(all_outliers)
tinyplot::plt(y_centroid ~ x_centroid | counts_outlier_scuttle, 
              data = adata$obs |> dplyr::filter(core_global == "mets2_B1"), 
              pal = "Classic Tableau",
              pch = 16, cex = 0.25, legend = legend(pt.cex = 2), asp = 1)

# Looks a bit better
outliers_scuttle <- table(adata$obs$core_global, adata$obs$counts_outlier_scuttle) 
sweep(x = outliers_scuttle, STATS = rowSums(outliers_scuttle), MARGIN = 1, FUN = "/") |> plot(main = "counts outliers")

# Revising our area filtering
areas <- adata$obs$cell_area
par(mar=c(4, 4, 4, 4))
hist(x = log(areas), breaks = 50)
ol <- isOutlier(nc, log = F, type = "both", nmads = 3)
(th <- attr(ol, "threshold")[1:2])
abline(v = th, col="blue")

# Adding to the AnnData and saving again
adata$obs$area_outlier_scuttle <- as.vector(ol)
anndataR::write_h5ad(object = adata, path = "crc-to-liver-mets.h5ad", mode = "w")

# Examining the embedding again: 
plot_embedding(adata$obs$counts_outlier_scuttle, adata$obsm$X_umap, rasterize = T, labels_discrete = F)
plot_embedding(adata$obs$area_outlier_scuttle, adata$obsm$X_umap, rasterize = T, labels_discrete = F)

# Checking cluster-wise
props <- adata$obs |> dplyr::group_by(initial_cluster) |> 
  dplyr::summarise(prop_ol = mean(counts_outlier_scuttle)) |>
  dplyr::arrange(desc(prop_ol)) |>  data.frame()
Ns <- adata$obs |> dplyr::group_by(initial_cluster) |> 
  dplyr::tally() |>
  dplyr::arrange(desc(n)) |>  data.frame()

dplyr::inner_join(x = props, y = Ns, by = "initial_cluster")
#    initial_cluster     prop_ol      n
# 1               12 0.978309648   2674
# 2               15 0.677966102     59
# 3               19 0.584000000    250
# 4               14 0.541666667     96
# 5               21 0.515306122    196
# 6               22 0.500000000      2
# 7               20 0.384615385     26
# 8               13 0.351648352    182
# 9                9 0.346948141   2179
# 10              17 0.266666667     15
# 11              18 0.226415094     53
# 12               1 0.070967886  89618
# 13               3 0.042908845 103615
# 14              16 0.031670282   2305
# 15              11 0.028414243  16541
# 16               5 0.025536393  85525
# 17               2 0.025183600  69172
# 18               6 0.016064420 296556
# 19               7 0.015815602 195187
# 20               8 0.014278260  42162
# 21              10 0.008847949  16501
# 22               4 0.007359273  12773

# We need to check: 
# -- 12 --> not real (cells with 0 counts)
# -- 15 --> not real (only marker is LINE1orf2)
# -- 19 --> not real (only marker is HSATII and there are very few)
# -- 14 --> not real (same as above)
# -- 21 --> not real (no real markers)
# -- 22 --> not real (2 cells)
# -- 20 --> not real (26 cells)
# -- 13 --> possibly real TRM CD8 T cells... but the markers do not really agree 
# -- 9 --> not real (decent number of cells, but no true markers and looks spatially random)
# -- 17 --> not real (15 cells)
# -- 18 --> not real (53 cells)

library(presto)
markers <- presto::wilcoxauc(X = Matrix::t(adata$layers$lognorm) |> as("CsparseMatrix") |> magrittr::set_rownames(value = adata$var_names), 
                             y = clusts)

# 12
markers |> dplyr::filter(group == "12", logFC > 0)
plot_embedding((adata$obs$transcript_counts == 0), adata$obsm$X_umap, rasterize = T, labels_discrete = F)

# 21
markers |> dplyr::filter(group == "21", logFC > 0)

# 15
markers |> dplyr::filter(group == "15", logFC > 0)

# 19
markers |> dplyr::filter(group == "19", logFC > 0)
table(adata$obs$core_global, adata$obs$initial_cluster == "19")
tmp <- adata$obs |> dplyr::filter(core_global == "mets2_B1")
tmp$c19 <- ifelse(tmp$initial_cluster == "19", yes = "c19", no = "!c19")
tinyplot::plt(y_centroid ~ x_centroid | c19, 
              data = tmp,
              pch = 16, cex = 0.25, legend = legend(pt.cex = 2), asp = 1)

# 14
markers |> dplyr::filter(group == "14", logFC > 0) |> dplyr::arrange(desc(logFC))
table(adata$obs$core_global, adata$obs$initial_cluster == "14")
tmp <- adata$obs |> dplyr::filter(core_global == "mets2_B1")
tmp$c14 <- ifelse(tmp$initial_cluster == "14", yes = "c14", no = "!c14")
tinyplot::plt(y_centroid ~ x_centroid | c14, 
              data = tmp,
              pch = 16, cex = 0.25, legend = legend(pt.cex = 2), asp = 1)

# 13
markers |> dplyr::filter(group == "13", logFC > 0) |> dplyr::arrange(desc(logFC)) 
table(adata$obs$core_global, adata$obs$initial_cluster == "13")
tmp <- adata$obs |> dplyr::filter(core_global == "mets2_B3")
tmp$c13 <- ifelse(tmp$initial_cluster == "13", yes = "c13", no = "!c13")
tinyplot::plt(y_centroid ~ x_centroid | c13, 
              data = tmp,
              pch = 16, cex = 0.25, legend = legend(pt.cex = 2), asp = 1)
plot_embedding(ifelse(adata$obs$initial_cluster == "13", yes = "c13", no = "!c13"), embedding = adata$obsm$X_umap, rasterize = T)

# 9
markers |> dplyr::filter(group == "9", logFC > 0) |> dplyr::arrange(desc(logFC)) 
table(adata$obs$core_global, adata$obs$initial_cluster == "9")
tmp <- adata$obs |> dplyr::filter(core_global == "prim2_F3")
tmp$c9 <- ifelse(tmp$initial_cluster == "9", yes = "c9", no = "!c9")
tinyplot::plt(y_centroid ~ x_centroid | c9, 
              data = tmp,
              pch = 16, cex = 0.25, legend = legend(pt.cex = 2), asp = 1)
plot_embedding(ifelse(adata$obs$initial_cluster == "9", yes = "c9", no = "!c9"), embedding = adata$obsm$X_umap, rasterize = T)

# I do not see any systematic trends on the area outlier axis
props <- adata$obs |> dplyr::group_by(initial_cluster) |> 
  dplyr::summarise(prop_ol = mean(area_outlier_scuttle)) |>
  dplyr::arrange(desc(prop_ol)) |>  data.frame()
dplyr::inner_join(x = props, y = Ns, by = "initial_cluster") 

# Same for nuclei
props <- adata$obs |> dplyr::group_by(initial_cluster) |> 
  dplyr::summarise(prop_ol = mean(nuclei_outlier)) |>
  dplyr::arrange(desc(prop_ol)) |>  data.frame()
dplyr::inner_join(x = props, y = Ns, by = "initial_cluster") 

# Final filtering
cellstokeep <- !(adata$obs$counts_outlier_scuttle | 
                   adata$obs$area_outlier_scuttle | 
                   adata$obs$nuclei_outlier | 
                   (adata$obs$initial_cluster %in% c(12, 22))) # Obvious noise
mean(cellstokeep) # 0.9631768

reticulate::py_require("anndata")
ad <- reticulate::import("anndata")
adata <- ad$read_h5ad("crc-to-liver-mets.h5ad")
adata_filtered <- adata[cellstokeep, ]
adata_filtered$write_h5ad("crc-to-liver-mets_filtered.h5ad")


# Note: some other cool things we could do related to data storage: 
# BPCells::write_matrix_hdf5(mat = norm, path = "test.hdf5", group = "matrix")
# BPCells::write_matrix_anndata_hdf5(mat = norm, "crc-to-liver-mets.h5ad", group = "/layers/lognorm")
# rhdf5::h5ls(file = "test.hdf5")
# rhdf5::h5ls("crc-to-liver-mets.h5ad")
# test <- rhdf5::h5read(file = "test.hdf5", name = "matrix")

