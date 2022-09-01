### 02 - Preprocessing single-cell RNA-seq
## Author: Diego Mañanes
## Date: 22/09/01

## loading packages
library("Seurat")
library("dplyr")
library("Matrix")
library("ggplot2")
library("scater")

## paths
projectPath <- here::here()
dataPath <- here::here("data")
metadataPath <- here::here("metadata")
sourceCode <- here::here("src")

## loading data
samplesMetadata <- readRDS(
  file = file.path(metadataPath, "samples_metadata_def.rds")
)
cellsMetadata <- readRDS(
  file = file.path(metadataPath, "cells_metadata_def.rds")
)
rawCounts <- readRDS(
  file = file.path(dataPath, "raw_counts_sparse.rds")
)
genesMetadata <- readRDS(
  file = file.path(metadataPath, "geneAnnotations.rds")
)
rawCounts <- rawCounts[genesMetadata$external_gene_name,]

## filtering low quality cells
sce <- SingleCellExperiment(
  assays = list(counts = rawCounts[genesMetadata$external_gene_name, ]),
  colData = cellsMetadata,
  rowData = genesMetadata
)
zeroExpr <- sum(rowSums(counts(sce)) == 0)
sce <- sce[rowSums(counts(sce) > 0) > 0, ]

MTGeneNames <- genesMetadata$external_gene_name[
  which(genesMetadata$chromosome_name == "MT")]
MTGeneNames <- intersect(MTGeneNames, rownames(sce))
sce <- addPerCellQC(
  x = sce, subsets = list(MT = MTGeneNames)
)
## min counts
minCounts <- 400
maxCounts <- 3000
filter_by_total_counts <- sce$sum > minCounts & sce$sum < maxCounts
## min detected genes
minGenes <- 200
maxGenes <- 1250
filter_by_total_genes <- sce$detected > minGenes & sce$detected < maxGenes
## min perc mitochondrial genes
maxMTperc <- 25
filter_by_mt_perc <- sce$subsets_MT_percent < maxMTperc
## manual filtering
sce$manual_filter <- (
  filter_by_total_counts & filter_by_total_genes & filter_by_mt_perc
)
validCells <- colnames(sce[, sce$manual_filter])
sum(sce$manual_filter)
filterCounts <- rawCounts[, colnames(rawCounts) %in% validCells]
seuratSC.Raw <- CreateSeuratObject(
  counts = filterCounts, 
  project = "HCL_Macrophages", 
  names.field = 1, 
  min.cells = 50,
  min.features = 200,
  names.delim = "-",
  meta.data = as.data.frame(colData(sce)[validCells, ])
)

## log normalize and save
seuratSC.Norm <- seuratSC.Raw %>% NormalizeData(
  normalization.method = "LogNormalize",
  vscale.factor = 10000,
  verbose = TRUE
)
saveRDS(seuratSC.Norm, file.path(dataPath, "seuratSC.Norm.rds"))
