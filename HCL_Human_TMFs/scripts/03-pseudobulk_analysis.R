### 03 - Pseudo-bulk analysis
## Author: Diego Mañanes
## Date: 22/09/01

library("Seurat")
library("dplyr")
library("ggplot2")
library("fgsea")
library("org.Hs.eg.db")
library("GO.db")
library("msigdbr")
library("Matrix.utils")
library("edgeR")
library("ggpubr")
library("gage")
library("gageData")
library("ComplexHeatmap")

## paths
projectPath <- here::here()
dataPath <- here::here("data")
metadataPath <- here::here("metadata")
sourceCode <- here::here("src")

## helper functions
source(file.path(sourceCode, "helperFunctions.R"))

## loading data
samplesMetadata <- readRDS(
  file = file.path(metadataPath, "samples_metadata_def.rds")
)
cellsMetadata <- readRDS(
  file = file.path(metadataPath, "cells_metadata_def.rds")
)
genesMetadata <- readRDS(
  file = file.path(metadataPath, "geneAnnotations.rds")
)
seuratSC.Sc.Norm <- readRDS(
  file = file.path(dataPath, "seuratSC.Norm.rds")
)

############# taking genesets #############
# GO - Biological Process
listBPGO <- getGOgenes(
  OrgDb = org.Hs.eg.db, selgo = "All", keytype = "SYMBOL", ont = "BP"
)
# KEGG pathways
data(kegg.gs)
listKeggF <- lapply(
  X = seq_along(kegg.gs), FUN = function(x) {
    x1 <- unique(genesMetadata$external_gene_name[genesMetadata$entrezgene %in% 
                                                    kegg.gs[[x]]])
    return(unique(x1[x1 %in% rownames(seuratSC.Sc.Norm)]))
  }
) %>% setNames(names(kegg.gs))

############# diagnostic plots #############
## formatting names 
seuratSC.Sc.Norm@meta.data <- seuratSC.Sc.Norm@meta.data %>% mutate(
  Sample = gsub(pattern = "Adult", replacement = "", x = Sample),
  Sample = case_when(
    Sample == "AscendingColon" ~ "Ascending colon",
    Sample == "SigmoidColon" ~ "Sigmoid colon",
    Sample == "SpleenParenchyma" ~ "Spleen",
    Sample == "TransverseColon" ~ "Transverse colon",
    TRUE ~ Sample
  )
)  %>% filter(Sample != "Cerebellum")
seuratSC.Sc.Norm <- subset(seuratSC.Sc.Norm, subset = Sample != "Cerebellum")

dfNumCells <- seuratSC.Sc.Norm@meta.data %>% 
  group_by(CellType = Sample) %>% 
  summarize(Number = n()) %>% arrange(desc(Number)) %>% 
  mutate(CellType = factor(CellType, levels = CellType))

## number of cells per tissue (raw, without cutoff)
ggplot(data = dfNumCells, mapping = aes(x = CellType, y = Number, fill = CellType)) + 
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = color.list()) +
  geom_text(aes(label = Number), vjust = -0.3, color = "black", size = 2.5) +
  theme_minimal() + 
  ggtitle("Number of cells per tissue") + 
  geom_hline(yintercept = 400, color = "red", alpha = 0.5) + ggtheme_custom + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

aggUMIsExpressTissue <- Matrix.utils::aggregate.Matrix(
  t(seuratSC.Sc.Norm@assays$RNA@counts),
  groupings = factor(seuratSC.Sc.Norm@meta.data$Sample),
  fun = "sum"
)
dfPlot <- data.frame(
  Tissue = rownames(aggUMIsExpressTissue),
  TotalGenes = apply(aggUMIsExpressTissue, 1, function(x) sum(x != 0)),
  TotalCounts = rowSums(aggUMIsExpressTissue),
  NumCells = table(seuratSC.Sc.Norm@meta.data$Sample),
  Sparsity = apply(aggUMIsExpressTissue, 1, function(x) sum(x == 0) / length(x))
)
ggplot(dfPlot, aes(x = TotalGenes, y = NumCells.Freq, color = Tissue)) +
  scale_color_manual(values = color.list()) +
  geom_point(size = 4) + ggtitle("Number of detected genes vs number of cells") +
  geom_hline(yintercept = 400, color = "red") +
  ylab("Number of cells") + xlab("Number of detected genes") +
  theme_minimal() + ggtheme_custom

maxMT <- 15
filter_by_MT_features <- seuratSC.Sc.Norm$subsets_MT_percent <= maxMT
seuratSC.Sc.Norm$filter_by_MT_features <- filter_by_MT_features
rawCountsFilt <- seuratSC.Sc.Norm@assays$RNA@counts[ , seuratSC.Sc.Norm$subsets_MT_percent <= maxMT]
seuratSC.Sc.Norm <- seuratSC.Sc.Norm[, seuratSC.Sc.Norm$subsets_MT_percent <= maxMT]

## setting colors
colorsTissues <-  c(
  "Adipose" = "#1f78b4", 
  "Esophagus" = "#fdbf6f", 
  "Gallbladder" = "#a6cee3", 
  "Heart" = "#b3d88b",
  "Kidney" = "#e3de5b",
  "Liver" = "#f58020", 
  "Lung" = "#e31a1c", 
  "Omentum" = "#cab3d6", 
  "Pleura" = "#33a248", 
  "Spleen" = "#6b3f98"
)

totalCountsPerTissue <- aggregate(
  x = seuratSC.Sc.Norm@meta.data$nCount_RNA, 
  by = list(seuratSC.Sc.Norm@meta.data$Sample), 
  FUN = sum
)
cpmsCounts <- edgeR::cpm(as.matrix(seuratSC.Sc.Norm@assays$RNA@counts))
totalCPMsPerTissue <- aggregate(
  x = colSums(seuratSC.Sc.Norm@assays$RNA@data), 
  by = list(seuratSC.Sc.Norm@meta.data$Sample), 
  FUN = mean
)
rownames(totalCPMsPerTissue) <- totalCPMsPerTissue$Group.1
rownames(totalCountsPerTissue) <- totalCountsPerTissue$Group.1

dfNumCells <- data.frame(
  Number = as.vector(
    sort(table(seuratSC.Sc.Norm@meta.data$Sample), decreasing = TRUE)
  ),
  CellType = factor(
    names(
      sort(table(seuratSC.Sc.Norm@meta.data$Sample), decreasing = TRUE)
    ),
    levels = names(
      sort(table(seuratSC.Sc.Norm@meta.data$Sample), decreasing = TRUE)
    )
  )
) %>% filter(Number >= 400, CellType != "Cerebellum") %>% mutate(
  TotalCounts = totalCountsPerTissue[as.vector(CellType), 2],
  AvgCountsPerCell = round(TotalCounts / Number, 2)
)
cellMatrixFilf <- seuratSC.Sc.Norm@assays$RNA@counts[
  , seuratSC.Sc.Norm@meta.data$Sample %in% dfNumCells$CellType
]
cellsMetadataFilf <- seuratSC.Sc.Norm@meta.data[
  seuratSC.Sc.Norm@meta.data$Sample %in% dfNumCells$CellType,
]

## number of TMFs per tissue
dfNumCells$CellType <- factor(
  dfNumCells$CellType,
  levels = with(dfNumCells, CellType[order(dfNumCells$Number, decreasing = T)])
)

ggplot(dfNumCells, aes(x = CellType, y = Number, fill = CellType)) + 
  geom_bar(stat = "identity", color = "black") + 
  scale_fill_manual(values = colorsTissues) + 
  geom_text(aes(label = Number), vjust = -0.3, color = "black", size = 2.5) +
  theme_minimal() + 
  ggtitle("Total number of cells") + ggtheme_custom + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1), 
    legend.title = element_text(face = "bold", size = 10)
  ) 

############# generation of pseudo-bulk samples #############
samples1 <- dfNumCells %>% filter(Number < 600)
samples2 <- dfNumCells %>% filter(Number > 600 & Number < 1000)
samples3 <- dfNumCells %>% filter(Number > 1000)

set.seed(123)
## setting number of cells per pseudo-bulk sample
listNumCells1 <- splitNumCells(samples1, numRep = 3)
names(listNumCells1) <- samples1$CellType
listNumCells2 <- splitNumCells(samples2, numRep = 5)
names(listNumCells2) <- samples2$CellType
listNumCells3 <- splitNumCells(samples3, numRep = 10)
names(listNumCells3) <- samples3$CellType

## selecting randomly cells for each replicate
cellsGroup1 <- splitCells(cellsMetadataFilf, listNumCells1)
cellsGroup2 <- splitCells(cellsMetadataFilf, listNumCells2)
cellsGroup3 <- splitCells(cellsMetadataFilf, listNumCells3)

## generating final expression matrix
cpmsCounts <- edgeR::cpm(seuratSC.Sc.Norm@assays$RNA@counts)
groupsList <- list(cellsGroup1, cellsGroup2, cellsGroup3)
listTotalDef <- lapply(
  X = groupsList, FUN = function(group) {
    generateFinalMatrixCPMs(
      groupCells = group, totalCounts = cpmsCounts
    )
  }
)
cpmCountsGroupsF <- do.call(cbind, args = listTotalDef)

## generating samples metadata
pseudoSamplesMetadata <- data.frame(
  Sample = colnames(cpmCountsGroupsF),
  Tissue = gsub(
    pattern = "_\\d*", replacement = "", x = colnames(cpmCountsGroupsF)
  )
) %>% mutate(
  Group = case_when(
    Tissue %in% samples1$CellType ~ "Group1",
    Tissue %in% samples2$CellType ~ "Group2",
    Tissue %in% samples3$CellType ~ "Group3",
    TRUE ~ "NA"
  )
)

dfMetrics <- data.frame(
  TotalGenes = apply(cpmCountsGroupsF, 2, function(x) sum(x != 0)),
  TotalCounts = colSums(cpmCountsGroupsF),
  Sparsity = round(apply(cpmCountsGroupsF, 2, function(x) sum(x == 0) / length(x)), 2),
  Sample = pseudoSamplesMetadata$Sample,
  Tissue = pseudoSamplesMetadata$Tissue,
  Group = pseudoSamplesMetadata$Group
)

## number of detected genes per pseudo-bulk sample
ggplot(dfMetrics, aes(x = Tissue, y = TotalGenes, fill = Tissue)) +
  geom_boxplot() +
  geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 0.5) +
  scale_fill_manual(values = colorsTissues) +
  theme_minimal() + ggtitle("Total genes detected") + 
  ggtheme_custom + theme(axis.text.x = element_text(angle = 90))

## remove genes without expression in any sample and set cutoffs to remove low 
# expressed genes
cpmCountsGroupsFOri <- cpmCountsGroupsF
cpmCountsGroupsF <- cpmCountsGroupsFOri[rowSums(cpmCountsGroupsFOri >= 10) >= 5, ]

## log2(CPM + 1)
logcpmTMMGroupsf <- log2(cpmCountsGroupsF + 1)

## similarity between pseudo-bulk samples
sampleDists <- dist(t(logcpmTMMGroupsf), method = "euclidean") 
col_fun <- circlize::colorRamp2(
  c(0, mean(sampleDists), max(sampleDists)), c("red", "white", "blue")
)
ha <- HeatmapAnnotation(
  Tissue = pseudoSamplesMetadata[["Tissue"]],
  col = list(Tissue = colorsTissues)
)

heatmapDis <- ComplexHeatmap::Heatmap(
  1 - as.matrix(sampleDists)/max(as.matrix(sampleDists)), 
  row_names_gp = gpar(fontsize = 7.5),
  column_names_gp = gpar(fontsize = 7.5), 
  name = "Similarity",
  column_title_gp = gpar(fontface = "bold"), 
  column_title = "Euclidean similarity between pseudo-bulk samples",
  column_names_rot = 45,
  show_row_dend = FALSE,
  row_title_side = "left",
  row_names_side = "left", 
  top_annotation = ha,
  show_column_names = FALSE,
  use_raster = TRUE, 
  raster_quality = 5
)
heatmapDis

############# analysis based on principal component analysis #############
PCAall <- prcomp(
  t(logcpmTMMGroupsf)[, which(apply(t(logcpmTMMGroupsf), 2, var) != 0)], 
  scale. = TRUE
)

variance <- round(factoextra::get_eigenvalue(PCAall)[c(1, 2), 2], 1)
p <- ggplot(data.frame(PCAall[["x"]]), aes(x = PC1, y = PC2)) +
  geom_point(
    aes(fill = as.factor(pseudoSamplesMetadata$Tissue)), 
    pch = 21, size = 3*1.5, color = "black"
  ) +
  scale_fill_manual(name = "Tissue", values = colorsTissues) + 
  xlab(paste0("PC1 (", variance[1], "%)")) + 
  ylab(paste0("PC2 (", variance[2], "%)")) + 
  geom_vline(xintercept = 0, size = 0.4, linetype = "dashed") + 
  geom_hline(yintercept = 0, size = 0.4, linetype = "dashed")  + 
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  theme_minimal() + 
  theme(
    plot.title = element_text(face = "bold", hjust = 0, size = 12),
    legend.title = element_text(face = "bold", hjust = 0.5, size = 11),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    text = element_text(family = "Helvetica"), aspect.ratio = 1
  ) + ggtitle("PCA human macrophages")
p


loadingVectors <- PCAall$rotation[, 1:3]
rankPC1 <- sort(loadingVectors[, 1], decreasing = TRUE)
rankPC2 <- sort(loadingVectors[, 2], decreasing = TRUE)
rankPC3 <- sort(loadingVectors[, 3], decreasing = TRUE)

## FGSEA on loading vectors with GO terms
fgseaPCAGOBP <- lapply(
  X = colnames(loadingVectors),
  FUN = function(x) {
    genes <- sort(loadingVectors[, x], decreasing = TRUE)
    fgseaMultilevel(
      listBPGO, stats = genes, minSize = 15, maxSize = 500
    ) 
  }
)
# PC1
topPathwaysUpGO1 <- fgseaPCAGOBP[[1]][NES > 0][head(order(padj), n = 10), pathway]
topPathwaysDownGO1 <- fgseaPCAGOBP[[1]][NES < 0][head(order(padj), n = 10), pathway]
topPathwaysGO1 <- c(topPathwaysUpGO1, rev(topPathwaysDownGO1))
grid::grid.newpage()
grDevices::dev.interactive()
plotGseaTable(
  listBPGO[topPathwaysGO1], rankPC1, fgseaPCAGOBP[[1]], 
  gseaParam = 0.5, render = TRUE
)
# PC2
topPathwaysUpGO2 <- fgseaPCAGOBP[[2]][NES > 0][head(order(padj), n = 10), pathway]
topPathwaysDownGO2 <- fgseaPCAGOBP[[2]][NES < 0][head(order(padj), n = 10), pathway]
topPathwaysGO2 <- c(topPathwaysUpGO2, rev(topPathwaysDownGO2))
grid::grid.newpage()
grDevices::dev.interactive()
plotGseaTable(
  listBPGO[topPathwaysGO2], rankPC2, fgseaPCAGOBP[[2]], 
  gseaParam = 0.5
)
# PC3
topPathwaysUpGO3 <- fgseaPCAGOBP[[3]][NES > 0][head(order(padj), n = 10), pathway]
topPathwaysDownGO3 <- fgseaPCAGOBP[[3]][NES < 0][head(order(padj), n = 10), pathway]
topPathwaysGO3 <- c(topPathwaysUpGO3, rev(topPathwaysDownGO3))
grid::grid.newpage()
grDevices::dev.interactive()
plotGseaTable(
  listBPGO[topPathwaysGO3], rankPC3, fgseaPCAGOBP[[3]], 
  gseaParam = 0.5
)

## FGSEA on loading vectors with KEGG ddbb
fgseaPCAKEGG <- lapply(
  X = colnames(loadingVectors),
  FUN = function(x) {
    genes <- sort(loadingVectors[, x], decreasing = TRUE)
    fgseaMultilevel(
      listKeggF, stats = genes, minSize = 15, maxSize = 500
    ) 
  }
)
# PC1
topPathwaysUpKEGG1 <- fgseaPCAKEGG[[1]][NES > 0][head(order(padj), n = 10), pathway]
topPathwaysDownKEGG1 <- fgseaPCAKEGG[[1]][NES < 0][head(order(padj), n = 10), pathway]
topPathwaysKEGG1 <- c(topPathwaysUpKEGG1, rev(topPathwaysDownKEGG1))
grid::grid.newpage()
grDevices::dev.interactive()
plotGseaTable(
  listKeggF[topPathwaysKEGG1], rankPC1, fgseaPCAKEGG[[1]], 
  gseaParam = 0.5, render = TRUE
)
# PC2
topPathwaysUpKEGG2 <- fgseaPCAKEGG[[2]][NES > 0][head(order(padj), n = 10), pathway]
topPathwaysDownKEGG2 <- fgseaPCAKEGG[[2]][NES < 0][head(order(padj), n = 10), pathway]
topPathwaysKEGG2 <- c(topPathwaysUpKEGG2, rev(topPathwaysDownKEGG2))
grid::grid.newpage()
grDevices::dev.interactive()
plotGseaTable(
  listKeggF[topPathwaysKEGG2], rankPC2, fgseaPCAKEGG[[2]], 
  gseaParam = 0.5, render = TRUE
)
# PC3
topPathwaysUpKEGG3 <- fgseaPCAKEGG[[3]][NES > 0][head(order(padj), n = 10), pathway]
topPathwaysDownKEGG3 <- fgseaPCAKEGG[[3]][NES < 0][head(order(padj), n = 10), pathway]
topPathwaysKEGG3 <- c(topPathwaysUpKEGG3, rev(topPathwaysDownKEGG3))
grid::grid.newpage()
grDevices::dev.interactive()
plotGseaTable(
  listKeggF[topPathwaysKEGG3], rankPC3, fgseaPCAKEGG[[3]], 
  gseaParam = 0.5, render = TRUE
)

## distribution of OXPHOS genes on PC1's loading vector 
plotEnrichment(
  listKeggF[["hsa00190 Oxidative phosphorylation"]], rankPC1) + 
  labs(title = "Oxidative phosphorylation (PC1)")+ 
  theme(plot.title = element_text(face = "bold")) + 
  ylab("Enrichment score") + 
  xlab("Rank based on PC1 scores")

oxphosFilt <- with(
  listKeggF, 
  `hsa00190 Oxidative phosphorylation`[`hsa00190 Oxidative phosphorylation` %in% 
                                       rownames(logcpmTMMGroupsf)]
)

`selGenes <- PCAall$rotation[, 1]
dfBarplotPCA <- data.frame(
  Genes = factor(names(selGenes), levels = names(selGenes)),
  LoadingVectors = selGenes,
  Pathway = names(selGenes) %in% oxphosFilt
)
dfRug <- data.frame(
  Rug = selGenes[names(selGenes) %in% oxphosFilt],
  Pathway = TRUE
)
dfRugB <- data.frame(
  RugB = selGenes[!names(selGenes) %in% oxphosFilt],
  Pathway = TRUE
)
pp <- ggplot(data = dfBarplotPCA, aes(x = LoadingVectors, fill = Pathway)) + 
  geom_density(alpha = 0.6, position = 'identity', color = "black") + 
  geom_rug(data = dfRug, mapping = aes(x = Rug), color = "#ed1c24", alpha = 0.8) + 
  theme_classic() + 
  ggtitle(paste0("Distribution of PC1 scores")) + 
  guides(fill = guide_legend(title = "Genes")) + 
  theme(
    legend.title = element_text(face = "bold", hjust = 0.5), 
    plot.title = element_text(face = "bold")
  ) + 
  scale_fill_manual(
    labels = c(
      paste0("Rest of genes (", sum(!dfBarplotPCA$Pathway), " genes)"), 
      paste0("OXPHOS genes (", sum(dfBarplotPCA$Pathway), " genes)")
    ), values = c("#818285", "#ed1c24")
  ) + xlab("PCA scores") + ylab("Density") 

ksTest <- ks.test(
  x = selGenes[names(selGenes) %in% oxphosFilt],
  y = selGenes[!names(selGenes) %in% oxphosFilt]
)
textKS <- paste(
  ksTest$method, 
  paste(names(ksTest$statistic), "statistic =", round(ksTest$statistic, 2)), 
  paste("p-value <", ifelse(ksTest$p.value == 0, "2.2e-16", ksTest$p.value)), 
  sep = "\n"
)

pp + annotate("text", x = 0.015, y = 45, label = textKS)

## enrichment analysis: corr. between pathway activities
selGOTerms <- read.csv(file.path(metadataPath, "selected_go_terms.csv"), header = FALSE)
listIntPathways <- append(
  getGOgenes(org.Hs.eg.db, selGOTerms[[2]], keytype = "SYMBOL"),
  listKeggF
)
cn <- colnames(logcpmTMMGroupsf)
gageListNoRef <- list()
reference <- seq(ncol(logcpmTMMGroupsf))
for (i in unique(pseudoSamplesMetadata$Tissue)) {
  sample <- grep(i, pseudoSamplesMetadata$Tissue, fixed = TRUE)
  gageListNoRef[[i]] <- gage(
    exprs = logcpmTMMGroupsf,
    use.fold = TRUE,
    gsets = listIntPathways,
    ref = reference,
    samp = sample, 
    rank.test = TRUE,
    compare = 'as.group'
  )
}
orderPathways <- rownames(gageListNoRef$Kidney$stats)
matrixNoRefGAGE <- as.data.frame(
  t(
    do.call(
      what = cbind, 
      args = lapply(
        X = seq_along(gageListNoRef), 
        FUN = function(x) {
          res <- gageListNoRef[[x]][["stats"]][orderPathways, "stat.mean", drop = FALSE]
          colnames(res) <- names(gageListNoRef)[x]
          return(res)
        }
      )
    )
  )
)
matrixNoRefGAGE <- matrixNoRefGAGE[, sapply(
  matrixNoRefGAGE, function(x) if(any(is.na(x))) return(FALSE) else TRUE
)] %>% mutate(Tissue = rownames(matrixNoRefGAGE))

## corr. plots
color.line <- rgb(43, 43, 43, max = 255, alpha = 60, names = "blue50")
for (i in c(
  "GO:0032367: intracellular cholesterol transport", 
  "GO:0009062: fatty acid catabolic process",
  "GO:0044242: cellular lipid catabolic process"
)) {
  print(
    ggscatter(
      matrixNoRefGAGE, x = "hsa00190 Oxidative phosphorylation", 
      y = i, 
      fill = "Tissue", size = 3*1.5, shape = 21, palette = colorsTissues,
      add = "none", conf.int = FALSE, 
      cor.coef = TRUE, cor.method = "pearson", 
      title = "T-score using log-FC (pseudo-bulk)"
    ) + geom_smooth(method = "lm", colour = color.line, se = FALSE) + 
      theme_classic() + 
      theme(
        plot.title = element_text(face = "bold", hjust = 0, size = 12),
        legend.title = element_text(face = "bold", size = 11),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        text = element_text(family = "Helvetica"),
        aspect.ratio=1
      ) + coord_fixed() + 
      scale_x_continuous(limits = symmetric_limits) + 
      scale_y_continuous(limits = symmetric_limits)
  )
}
`