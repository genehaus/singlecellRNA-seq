
## written by Hyojin Kim 
## March 2022

## R version 4.1.0 (2021-05-18)
## Seurat version 4.3.0
## SingleCellExperiment version 1.16.0
## zellkonverter version 1.4.0
## ggplot2 version 3.4.2

library(Seurat)
library(dplyr)
library(ggplot2)
library(SingleCellExperiment)
library(tibble)
library(zellkonverter)
library(future)
set.seed(42)
args <- commandArgs(trailingOnly=TRUE)
plan("multicore", workers = 10)


## set up
## ---------------------- ## 
indir = args[1] 
file = args[2] 

sc = readRDS(file = paste0(indir, file))
sc.sce <- as.SingleCellExperiment(sc)
writeH5AD(sc.sce, X_name = "counts", file = paste0(indir, file, ".to.h5ad"))
