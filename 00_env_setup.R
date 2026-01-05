renv::init(bare = T)
.libPaths()

Sys.which("clang") # "/usr/bin/clang" 

renv::install("bioc::rhdf5")
renv::install("Matrix")

renv::install(packages = list("bioc::SingleCellExperiment", "Seurat"))
renv::install("bnprks/BPCells/r")

renv::install("pak")
renv::install("remotes")
remotes::install_github("scverse/anndataR", ref="084a187") # This is an older version that works with R 4.4.2

renv::install("korsunskylab/spatula")
renv::install("immunogenomics/presto@glmm")
renv::install("immunogenomics/singlecellmethods")

renv::install(packages = list("Bioc::limma", "Bioc::edgeR", "Bioc::DESeq2", "Bioc::glmGamPoi"))

renv::install("bioc::sparseMatrixStats")
install.packages("irlba", type = "source")

renv::install(packages = list("lme4", "lmerTest"))

install.packages("https://cran.r-project.org/bin/macosx/big-sur-arm64/contrib/4.5/Rfast_2.1.5.2.tgz", repos = NULL, type = "binary")

renv::snapshot()
.rs.restartR()

