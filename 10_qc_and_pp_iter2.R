rm(list = ls())
.rs.restartR(clean = T)
.libPaths()

library(BPCells)
library(Matrix)

## QC & PP ITER 2 --------------------------------------------------------------
# We are now re-running the pipeline without the filtered cells.

# Reading the data
adata = anndataR::read_h5ad(path = "crc-to-liver-mets_filtered.h5ad", mode = "r+")

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
anndataR::write_h5ad(object = adata, path = "crc-to-liver-mets_filtered.h5ad", mode = "w")

# Normalization
# I stole this trick from NanoString and added the change of base part.
cts <- adata$layers$counts |>  as("CsparseMatrix")
scaling_factor <- 1000
norm_factors <- Matrix::Diagonal(x = scaling_factor/adata$obs$transcript_counts, names=rownames(adata$layers$counts))
norm <- ((norm_factors %*% adata$layers$counts) |> log1p())/log(2)

# Adding to the anndata and saving
adata$layers$lognorm <- norm
anndataR::write_h5ad(object = adata, path = "crc-to-liver-mets_filtered.h5ad", mode = "w")

# Saving space
remove(norm)
remove(cts)

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
  features = c("CD4", "CD8A", "CD3E", "CD19"), 
  rasterize = T,
  colors_continuous = viridis::viridis(n = 71)
)

# Louvain clustering
clusts <- knn_hnsw(adata$obsm$X_scVI, k = 30, metric = "cosine", ef = 1000) |> # Find approximate nearest neighbors
  knn_to_snn_graph() |> # Convert to a SNN graph
  cluster_graph_louvain(resolution = 0.5) # Perform graph-based clustering

dir.create("10_qc_and_pp_iter2_outs")
pdf(file = "10_qc_and_pp_iter2_outs/umap_by_clusters.pdf", width = 6, height = 6)
plot_embedding(clusts, adata$obsm$X_umap, rasterize = T)
dev.off()

# Adding to the AnnData and re-saving
adata$obs$final_cluster <- clusts
anndataR::write_h5ad(object = adata, path = "crc-to-liver-mets_filtered.h5ad", mode = "w")

