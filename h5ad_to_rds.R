#!/xxx/bin/Rscript

## written by H.J.Kim
## Jan 2022
## file format convertion from h5ad to seurat
## ".X" of AnnData of h5ad is converted to Seurat_obj[["RNA"]]@counts
## NOTE : cells about 1,000,000 cannot be converted to Seurat object by this way 

setwd("./")

library(Seurat)
library(Scater)
library(dplyr)
library(tidyr)
library(zellkonverter)

h5ad_path <- "/xxx/"
h5ad_file <- c("data.h5ad")

S_name <- paste0(h5ad_path, h5ad_file)
S <- readH5AD(S_name)
S <- as.Seurat(S, counts = "X", data = "X")
S <- NormalizeData(S, normalization.method = "LogNormalize", scale.factor = 10000)
S <- FindVariableFeatures(S, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(S)
S <- ScaleData(S, features = all.genes)
saveRDS(S, paste0(h5ad_file, ".to.rds"))

