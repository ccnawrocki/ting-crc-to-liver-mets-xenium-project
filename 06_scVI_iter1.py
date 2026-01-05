##  Creating a conda environment -----------------------------------------------
# I did the following in the terminal:

# $ conda create -n ting-crc-to-liver-mets python=3.11
# $ conda activate ting-crc-to-liver-mets
# $ pip install -U scvi-tools
# $ pip install rpy2
# $ conda env export > ting-crc-to-liver-mets.yml

# As long as we install libraries in the conda environment, we can import them
# within the project, after using the reticulate call below.

## Using this environment in our project ---------------------------------------
# I did the following in the console, using the R kernel: 

# rm(list = ls())
# .rs.restartR(clean = T)
# reticulate::conda_list() # Shows all the environments and their paths
# reticulate::use_condaenv("ting-crc-to-liver-mets", conda = "/opt/homebrew/Caskroom/miniforge/base/bin/conda")

## Setup -----------------------------------------------------------------------
# Actual Python code begins here:

import anndata as ad
import scvi # scVI for deep learning methods
import scanpy as sc # scanpy for processing functions
import rpy2.robjects as ro # rpy2 lets us make R function calls
import pandas as pd
import torch
import matplotlib.pyplot as plt

adata = sc.read_h5ad("crc-to-liver-mets.h5ad")


##### ITERATION 1 --------------------------------------------------------------
# -- We will keep all samples together, rather than processing primaries and 
# mets separately. 
# -- We will integrate the data, using core as the batch variable.
# -- We will use scVI & scanpy for processing steps.
# -- Next, we will identify noise clusters and filter them out.

## scVI modeling
adata.layers["counts"] = adata.X.copy().astype(int)
adata.layers["counts"]

hvgs = adata.var.sort_values("overdispersion", ascending = False).iloc[0:2500].index.tolist()
adata.var["highly_variable"] = adata.var_names.isin(hvgs)
adata.var["highly_variable"].sum()

adata.obs.loc[:,"core_global"] = adata.obs.slide.str.cat(adata.obs.core, sep="_").astype("category")

# scvi.model.SCVI.setup_anndata(
#   adata = adata,
#   layer = "counts",
#   batch_key = "core_global"
# )
# 
# model1 = scvi.model.SCVI(
#     adata,
#     gene_likelihood="nb" 
# )
# 
# # Claude helped me with the nitty-gritty for speeding up the training: 
# model1.train(
#     check_val_every_n_epoch = 5,  # Validate less frequently (was 1)
#     plan_kwargs={
#         "lr": 1e-3,  # Slightly higher learning rate can speed convergence
#     },
#     batch_size=2048,
#     max_epochs = 250,
#     early_stopping = True,
#     early_stopping_patience = 20,
#     early_stopping_monitor="elbo_validation", 
#     accelerator = "mps" # This is the apple silicon GPU. Use "gpu" for CUDA.
# ) 
# # RUNTIME: 1:42:36

# !mkdir 06_scVI_iter1_outs
# model1.save("06_scVI_iter1_outs/scVI_model1", overwrite=True)
model1 = scvi.model.SCVI.load(adata = adata, dir_path = "06_scVI_iter1_outs/scVI_model1")

train_test_results = model1.history["elbo_train"]
train_test_results["elbo_validation"] = model1.history["elbo_validation"]
train_test_results.plot(logy=True)
plt.show()

SCVI_LATENT_KEY = "X_scVI"
latent = model1.get_latent_representation()
adata.obsm[SCVI_LATENT_KEY] = latent

adata.write_h5ad("crc-to-liver-mets.h5ad")
globals().clear()

