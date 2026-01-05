# In the console, using the R kernel: 
# reticulate::py_require("anndata")

import anndata as ad
import glob

!ls *.h5ad
adata_paths = glob.glob("*.h5ad")

adata = {}
for p in adata_paths:
  nm = p.replace(".h5ad", "")
  adata[nm] = ad.read_h5ad(p)

adata = ad.concat(adata, index_unique = "_")
adata.shape
# (1106043, 5027)

# I manually screened cores that were not suitable for this study. 
# I remove them here: 
tokeep = ~(((adata.obs["slide"] == "mets2") & (adata.obs["core"].isin(["D3", "C1"]))) | ((adata.obs["slide"] == "prim2") & (adata.obs["core"].isin(["F2", "B1", "E1", "B2", "C2", "E3"]))))
adata = adata[tokeep,]
adata.shape
# (935687, 5027)

adata.write_h5ad("crc-to-liver-mets.h5ad")

