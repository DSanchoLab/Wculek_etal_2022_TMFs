### 01 - Prepare data HCL
## Author: Diego Mañanes
## Date: 22/09/01

## loading packages
library("dplyr")
library("Matrix")
library("data.table")

## paths
projectPath <- here::here()
dataPath <- here::here("data")
metadataPath <- here::here("metadata")
sourceCode <- here::here("src")
dataFigsharePath <- here::here(dataPath, "dge_raw_data") 
dataFigsharePathR <- here::here(dataPath, "dge")
metadataFigsharePath <- here::here(metadataPath, "annotation_rmbatch_data_revised417") 

## loading metadata
sampleMetadata <- readRDS(
  file = file.path(metadataPath, "samples_metadata_def.rds")
)

cellsMetadata <- readRDS(
  file = file.path(metadataPath, "cells_metadata_def.rds")
)

## downloading data from figshare
# options(timeout = max(1000, getOption("timeout")))
urlMetadata <- "https://figshare.com/ndownloader/files/22447898"
destFileMetadata <- file.path(metadataFigsharePath)
download.file(urlMetadata, paste0(destFileMetadata, ".zip"))
unzip(
  paste0(destFileMetadata, ".zip"), exdir = metadataPath
)
unlink(paste0(destFileMetadata, ".zip"))
urlData <- "https://figshare.com/ndownloader/files/23062979"
destFileData <- file.path(dataFigsharePath)
download.file(urlData, paste0(destFileData, ".tar.gz"))
untar(paste0(destFileData, ".tar.gz"), exdir = dataPath)
unlink(paste0(destFileData, ".tar.gz"))

## take matrices
listMatrices <- list()
## just in case
# if (!file.exists(file.path(dataPath, "macrophages_samples"))) 
#   dir.create(file.path(dataPath, "macrophages_samples"))

for (i in seq(nrow(sampleMetadata)))  {
  rowSel <- sampleMetadata[i, ]
  # read dataset
  dataset <- as.matrix(
    read.delim(
      file = gzfile(file.path(dataFigsharePathR, rowSel$file_name_figshare)), 
      sep = " "
    )
  )
  # read metadata
  metadata <- fread(
    file = file.path(metadataFigsharePath, rowSel$file_name_metadata), 
    sep = ",", data.table = FALSE
  )
  ## filtering
  filtMetadata <- metadata %>% filter(
    grepl(pattern = "macro|kuppfer", x = Celltype, ignore.case = TRUE)
  )
  filtDataset <- dataset[, filtMetadata$Cell_id]
  fileNew <- paste0(
    gsub(".txt.gz", "", rowSel$file_name_figshare), "_macrophages", ".tsv"
  )
  ## just in case
  # fwrite(
  #   data.frame(filtDataset),
  #   file = file.path(dataPath, "macrophages_samples", fileNew),
  #   sep = "\t", row.names = TRUE, col.names = TRUE
  # )
  listMatrices[[i]] <- filtDataset
} 
## generate raw count matrix
rawCountsMacro <- Reduce(
  function(a, b) {
    newM <- merge(a, b, by = "row.names", all = TRUE)
    rn <- newM[, 1] 
    newM <- newM[,-1]
    row.names(newM) <- rn 
    return(newM) 
  }, 
  listMatrices[-1], 
  init = listMatrices[[1]]
)
rawCountsMacro[is.na(rawCountsMacro)] <- 0

## plain text
if (!file.exists(file.path(dataPath, "macrophages_samples", 
                           "all_raw_counts_macrophages.tsv"))) {
  fwrite(
    rawCountsMacro,
    file = file.path(dataPath, "macrophages_samples", "all_raw_counts_macrophages.tsv"),
    sep = "\t", row.names = TRUE, col.names = TRUE
  )  
}
## rds object
if (!file.exists(file.path(objectsPath, "raw_counts_sparse.rds"))) {
  rawCountsMacroS <- Matrix(as.matrix(rawCountsMacro), sparse = TRUE)
  saveRDS(
    object = rawCountsMacroS, 
    file = file.path(dataPath, "raw_counts_sparse.rds")
  )  
}
