rm(list = ls())
.rs.restartR(clean = T)
.libPaths()

## counts matrix ##

# Can use h5 then construct the matrix
rhdf5::h5ls("../../ting-crc-to-liver-mets-data/prim1/stow/cell_feature_matrix.h5")
cts <- rhdf5::h5read(file = "../../ting-crc-to-liver-mets-data/prim1/stow/cell_feature_matrix.h5", name = "matrix")
str(cts)

library(Matrix)
m <- Matrix::sparseMatrix(x = cts$data, 
                          dims = cts$shape, 
                          i = (cts$indices+1), 
                          p = cts$indptr, 
                          repr = "C", 
                          dimnames = list(cts$features$name, cts$barcodes)
                          )
rna <- m[which((cts$features$feature_type == "Gene Expression")),]
colSums(rna)[1:25]

# Can also read the MatrixMarket format
testMM <- Seurat::Read10X(data.dir = "../../ting-crc-to-liver-mets-data/prim1/stow/cell_feature_matrix")

# These functions work for this task too
# Matrix::readMM()
# DropletUtils::read10xCounts()


## metadata ##

# csv is easy, but parquet works too.
meta <- data.table::fread("../../ting-crc-to-liver-mets-data/prim1/stow/cells.csv.gz")
meta$transcript_counts[1:25]

all(meta$transcript_counts[1:25] == colSums(rna)[1:25]) # TRUE

mindray_discrete <- c("#10BBB9", "#233987", "#030306", "#C7360B", "#F6BF15")
mindray_palette <- colorRampPalette(mindray_discrete)(256)
image(matrix(1:256, ncol=1), col = mindray_palette, axes = FALSE, main = "Mindray Doppler Palette")

par(mar = c(0, 0, 0, 0))
plot(x = meta$x_centroid, y = -meta$y_centroid, pch = ".", asp = 1, 
     col = mindray_palette[
       pmin(256, 1 + round(255 * meta$transcript_counts / quantile(meta$transcript_counts, 0.90)))
     ])


## segmentation ##

# same as with the metadata
seg <- arrow::read_parquet("../../ting-crc-to-liver-mets-data/prim1/stow/cell_boundaries.parquet")

library(ggplot2)
seg[1:68000,] |> 
  ggplot() + 
  geom_polygon(mapping = aes(x = vertex_x, y = vertex_y, group = cell_id))


## transcripts ##

# Many tools are used for reading large datasets like this. I used to use apache
# arrow combined with dplyr's collect() function.

# DuckDB is more robust and reliable. Plus, it is just a good tool to become 
# familiar with.

# AI helped with this: 
library(DBI)
library(duckdb)

con <- dbConnect(duckdb(), dbdir = ":memory:")
on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)

# Works with a single file path or a folder/glob of files
path <- "../../ting-crc-to-liver-mets-data/prim1/stow/transcripts.parquet"

tx <- dbGetQuery(con, sprintf("
  SELECT
    x_location, y_location, qv, is_gene,
    CAST(cell_id AS VARCHAR)           AS cell_id,
    CAST(feature_name AS VARCHAR)      AS feature_name,
    CAST(fov_name AS VARCHAR)          AS fov_name,
    CAST(codeword_category AS VARCHAR) AS codeword_category
  FROM read_parquet('%s')
  WHERE codeword_category = 'custom_gene'
    AND feature_name IN ('HSATII', 'LINE1orf1', 'HERVK')
", path))

dbDisconnect(con, shutdown = T) # Playing it safe

head(tx) # Looks good

# Same plotting tricks:
# png(filename = "tx_test_plot.png", width = 10, height = 12, units = "in", res = 96)
plot(x = tx$x_location, y = -tx$y_location, pch = ".", asp = 1,
     col = dplyr::case_when(tx$feature_name == "HSATII" ~ "red", 
                            tx$feature_name == "LINE1orf1" ~ "darkblue", 
                            tx$feature_name == "HERVK" ~ "green3")
     )
legend(x = "topright", pch = 16, 
       col = c("red", "darkblue", "green3"), 
       legend = c("HSATII", "LINE1orf1", "HERVK"), pt.cex = 2)
# dev.off()

# png(filename = "tx_test_plot2.png", width = 10, height = 12, units = "in", res = 150)
scattermore::scattermoreplot(x = tx$x_location, y = -tx$y_location, asp = 1, pch = ".",
                             col = dplyr::case_when(tx$feature_name == "HSATII" ~ "red", 
                                                    tx$feature_name == "LINE1orf1" ~ "darkblue", 
                                                    tx$feature_name == "HERVK" ~ "green3")
                             )
legend(x = "topright", pch = 16, 
       col = c("red", "darkblue", "green3"), 
       legend = c("HSATII", "LINE1orf1", "HERVK"), pt.cex = 2)
# dev.off()

# This function helps too: 
aspectratio <- (range(tx$y_location)[2] - range(tx$y_location)[1]) / (range(tx$x_location)[2] - range(tx$x_location)[1])
scattermore::scattermore(xy = as.matrix(tx[,1:2]), size = c(1028, 1028)) |> plot(asp = aspectratio) 

par(mar = c(5, 4, 4, 2) + 0.1)


## loom ##

# After much tinkering, it seems that loomR is no longer maintained, so it 
# is no longer really compatible with scverse libraries in Python.
# It can still be useful for storing the data and query certain parts of it: 

# Creating a loom file
gene_mask <- cts$features$feature_type == "Gene Expression"
genemeta <- lapply(cts$features, function(x) x[gene_mask])[2:5]
loomR::create(filename = "prim1.loom", 
              data = rna, gene.attrs = genemeta, cell.attrs = meta,
              do.transpose = T, overwrite = T)

# Connecting
lfile <- loomR::connect(filename = "prim1.loom", mode = "r+")

# Practice slicing
lfile[["row_attrs/id"]][1:5]
lfile[["col_attrs/cell_names"]][1:5]
lfile[["matrix"]][1:5, 1:5]
lfile$row.attrs$gene_names[1:5]
which(lfile$row.attrs$gene_names[] == "LINE1orf1")

# Some lazy computation
library(loomR)
nUMI_map <- lfile$map(FUN = rowSums, MARGIN = 2, chunk.size = 500, dataset.use = "matrix", 
                      display.progress = F)

all(nUMI_map == meta$transcript_counts)  # TRUE

# Getting the sparse matrix again: 
mat <- lfile[["matrix"]][,] |> as("CsparseMatrix") # Not ideal memory usage
remove(mat)

# Closing
lfile$close_all()

# The up-to-date loom is really only loompy. Here we implement this (AI helped):
library(reticulate)
reticulate::py_install("loompy")
loompy <- import("loompy")
np <- import("numpy")

# Convert R metadata lists to numpy byte arrays (fixed-length ASCII)
meta$segmentation_method <- gsub(pattern = "µ", replacement = "u", x = meta$segmentation_method)
to_fixed_ascii <- function(x) np$array(enc2utf8(as.character(x)), dtype="S")

row_attrs <- dict()
for (nm in names(genemeta)) {
  row_attrs[[nm]] <- to_fixed_ascii(genemeta[[nm]])
}

col_attrs <- dict()
for (nm in colnames(meta)) {
  if (is.character(meta[[nm]])) {
    col_attrs[[nm]] <- to_fixed_ascii(meta[[nm]])
  }
  else {
    col_attrs[[nm]] <- meta[[nm]]
  }
}

rna_t <- as(rna, "TsparseMatrix")
scipy_sparse <- import("scipy.sparse")

coo <- scipy_sparse$coo_matrix(
  reticulate::tuple(
    np$array(rna_t@x),
    reticulate::tuple(np$array(rna_t@i), np$array(rna_t@j))
  ),
  shape = reticulate::tuple(nrow(rna_t), ncol(rna_t))
)

loompy$create(
  "prim1_fixed_sparse.loom",
  coo,
  row_attrs = row_attrs,
  col_attrs = col_attrs
)

# Reading as sparse
loomfile <- anndata::read_loom(filename = "prim1_fixed_sparse.loom", sparse = T)

# Or, this works too
ad <- import("anndata")
loomfile <- ad$io$read_loom("prim1_fixed_sparse.loom", sparse = T)
loomfile$obs

# I checked and I can open this in Python.


## AnnData ##

# Can just write the loom file:
anndata::write_h5ad(anndata = loomfile, filename = "prim1.h5ad")

# Can build the anndata in a couple ways without loom: 
metalist <- as.list(meta)
prim1ad <- anndata::AnnData(X = Matrix::t(rna), obs = metalist, var = genemeta)
prim1ad
# anndata::write_h5ad()

genemetadf <- do.call(what = cbind, args = genemeta) |> as.data.frame()
prim1ad <- anndataR::AnnData(X = Matrix::t(rna), obs = meta, var = genemetadf)
prim1ad
# anndataR::write_h5ad()

# Could also just read the data into Python and make an AnnData object there. 

## qs2 ##

# Good for saving seurat objects and R data files... it's super fast.
rownames(meta) <- meta$cell_id
prim1list <- list("cts" = rna, "meta" = meta, "genes" = genemetadf)

qs2::qs_save(prim1list, file = "prim1.qs2")


## DELETING THESE SAVED OBJECTS ##
# We do not need them anymore.
file.remove(c("prim1.qs2", "prim1.h5ad", "prim1_fixed_sparse.loom", "prim1.loom"))

