rm(list = ls())
.rs.restartR(clean = T)
.libPaths()

library(anndataR)
library(nebula)

adata <- read_h5ad(path = "crc-to-liver-mets_filtered.h5ad", mode = "r+")
adata$obs |> dplyr::glimpse() # Good
adata$var |> dplyr::glimpse() # We will overwrite this below
adata$X |> abind::acorn() # Good

adata$obs$counts_outlier_scuttle |> table() # Just making sure we actually filtered

## HVG selection Round 2 -------------------------------------------------------
# -- At this point, we have done QC via iteration 1.
# -- Again, we will select 2500 genes out of the 5000.
# -- Again, we will use NEBULA to identify the genes.
# -- Afterwards, we will move on to finishing our processing.

# Setting up the data 
cts <- adata$X
dimnames(cts) <- list(adata$obs_names, adata$var_names)
cts <- Matrix::t(cts) |> as("CsparseMatrix")

meta <- adata$obs
meta$core <- as.character(meta$core) |> as.factor()
idx <- meta |> dplyr::arrange(core) |> rownames()
cts <- cts[,idx]
meta <- meta[idx,]
mm <- model.matrix(~1, data = meta)
eff <- meta$transcript_counts

# Fitting the nebula model
# nfit <- nebula::nebula(count = cts, id = meta$core, pred = mm, offset = eff,
#                        covariance = F, cpc = 0.001, ncore = 6)
# dir.create("08_hvg_selection_iter2_outs")
# qs2::qs_save(object = nfit, file = "08_hvg_selection_iter2_outs/hvg_nfit.qs2")
nfit <- qs2::qs_read(file = "08_hvg_selection_iter2_outs/hvg_nfit.qs2")

# Plotting the results
tokeep <- nfit$convergence >= -10
hvg_data <- data.frame(mean_expr = nfit$summary$`logFC_(Intercept)`,
                       se_expr = nfit$summary$`se_(Intercept)`,
                       overdispersion = nfit$overdispersion$Cell, 
                       gene = nfit$summary$gene)
hvg_data <- hvg_data[tokeep,]
hvgs <- dplyr::arrange(hvg_data, desc(overdispersion))[1:2500, "gene"]

# pdf("08_hvg_selection_iter2_outs/2500_HVGs.pdf", width = 6, height = 6)
plot(x = hvg_data$mean_expr, hvg_data$overdispersion, pch = 16, cex = 0.5,
     col = ifelse(hvg_data$gene %in% hvgs, yes = "red", no = "black"))
# dev.off()

# pdf("08_hvg_selection_iter2_outs/mean-variance_trend.pdf", width = 6, height = 6)
plot(x = hvg_data$mean_expr, y = sqrt(hvg_data$se_expr), pch = 16, cex = 0.5)
fitted <- loess(sqrt(hvg_data$se_expr) ~ hvg_data$mean_expr, span = 0.3)
lines(x = fitted$x[order(fitted$x)], y = fitted$fitted[order(fitted$x)], col = "red", lwd = 4)
# dev.off()

# Adding these results to the AnnData
adata$var[["overdispersion"]] <- NA
adata$var[hvg_data$gene, "overdispersion"] <- hvg_data$overdispersion
adata$var[["log_cpc"]] <- NA
adata$var[hvg_data$gene, "log_cpc"] <- hvg_data$mean_expr
adata$var[["log_cpc_se"]] <- NA
adata$var[hvg_data$gene, "log_cpc_se"] <- hvg_data$se_expr

# Adding ensembl ids
ens <- rhdf5::h5read(file = "prim2.h5ad", name = "/var") |> dplyr::bind_rows()
adata$var[["ensembl_id"]] <- NA
adata$var[ens$name, "ensembl_id"] <- ens$id

# Re-saving
anndataR::write_h5ad(object = adata, path = "crc-to-liver-mets_filtered.h5ad", mode = "w")

