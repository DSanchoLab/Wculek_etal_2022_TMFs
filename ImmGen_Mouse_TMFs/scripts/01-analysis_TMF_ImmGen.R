### 01 - Analysis of TMF bulk RNA-seq data from ImmGen project
## Author: Diego Mañanes
## Date: 22/09/01

## loading packages
library("dplyr")
library("data.table")
library("edgeR")
library("org.Mm.eg.db")
library("ggplot2")
library("fgsea")
library("GO.db")
library("gage")
library("gageData")
library("ggpubr")
library("ComplexHeatmap")

## paths
projectPath <- here::here()
dataPath <- here::here("data")
metadataPath <- here::here("metadata")
sourceCode <- here::here("src")

## helper functions
source(file.path(sourceCode, "helperFunctions.R"))

############# loading data #############
## GSE109125
rawCountsGSE109125 <- data.table::fread(
  file.path(dataPath, "GSE109125_Gene_count_table_GENCODE_vM25.csv.gz")
) %>% as.data.frame()
geneSymbol <- rawCountsGSE109125[[1]]
index <- lapply( 
  X = list("MF.RP", "MF.PC", "MF.Alv.Lu", "MF.microglia.CNS", "MF.226.II.480lo.PC"), 
  FUN = function (x, y) grep(pattern = x, x = colnames(y), ignore.case = T),
  y = rawCountsGSE109125
) %>% unlist()
rawCountsGSE109125 <- rawCountsGSE109125[, index]
rownames(rawCountsGSE109125) <- toupper(geneSymbol)
samplesMetadataGSE109125 <- data.table::fread(
  file = file.path(metadataPath, "GSE109125_SraRunTable.txt")
)
colnames(samplesMetadataGSE109125) <- tolower(
  gsub(pattern = " ", replacement = "_", x = colnames(samplesMetadataGSE109125))
)
## GSE122108
rawCountsGSE122108 <- data.table::fread(
  file.path(dataPath, "GSE122108_Gene_count_table_GENCODE_vM25.csv.gz")
)%>% as.data.frame()
geneSymbol <- rawCountsGSE122108[[1]]
index <- lapply(
  X = list("MF.KC.Clec4FpTim4p64p.Lv", "MF.6Gn480hi.Kd", "MF.480p64pMerTKp.WAT",
           "MF.480p64pMerTKp.MAT", "MF.480p64pLyve1p.MS", 
           "MF.64p.Th"), 
  FUN = function (x, y) grep(pattern = x, x = colnames(y), ignore.case = T),
  y = rawCountsGSE122108
) %>% unlist()
rawCountsGSE122108 <- rawCountsGSE122108[, index]
rownames(rawCountsGSE122108) <- toupper(geneSymbol)
samplesMetadataGSE122108 <- data.table::fread(
  file = file.path(metadataPath, "GSE122108_SraRunTable.txt")
)
colnames(samplesMetadataGSE122108) <- tolower(
  gsub(pattern = " ", replacement = "_", x = colnames(samplesMetadataGSE122108))
)
## build samples metadata
vec <- readLines(file.path(metadataPath, "series_samples_mod.txt"))
dfIDs <- data.frame(
  ID = strsplit(x = vec[1], split = ",")[[1]],
  Name = strsplit(x = vec[2], split = ", ")[[1]] %>% 
    gsub(pattern = "_RNA-seq ", replacement = "", x = .) %>% 
    gsub(pattern = " ", replacement = "", x = .)
)
dfIDs <- dfIDs[-1, ]
namesGSE109125 <- dfIDs[dfIDs$Name %in% colnames(rawCountsGSE109125), ]
namesGSE122108 <- dfIDs[dfIDs$Name %in% colnames(rawCountsGSE122108), ]
samplesMetadata <- rbind(
  samplesMetadataGSE122108[samplesMetadataGSE122108$`geo_accession_(exp)` %in% 
                             namesGSE122108$ID, ],
  samplesMetadataGSE109125[samplesMetadataGSE109125$`geo_accession_(exp)` %in% 
                             namesGSE109125$ID, ], 
  fill = TRUE
)
samplesMetadata <- merge(
  samplesMetadata, dfIDs, by.x = "geo_accession_(exp)", by.y = "ID"
)
samplesMetadata$Groups <- gsub(
  pattern = "#\\d|\\.\\d$", replacement = "", x = samplesMetadata$Name, perl = TRUE
)
samplesMetadata$TissueDef <- ifelse(
  test = is.na(samplesMetadata$tissue), 
  yes = samplesMetadata$source_name, 
  no = samplesMetadata$tissue
)
## build final count matrix
rawCounts <- cbind(
  rawCountsGSE109125, rawCountsGSE122108[rownames(rawCountsGSE109125),]
)
rawCounts <- rawCounts[, samplesMetadata$Name]
## build genes metadata
gtfMouse <- rtracklayer::readGFF(
  file.path(metadataPath, "gencode.vM25.annotation.gtf.gz")
)
gtfMouse$gene_name <- toupper(gtfMouse$gene_name)
entrez <- read.table(file.path(metadataPath, "gencode.vM25.metadata.EntrezGene"))

colnames(entrez) <- c("transcript_id", "entrez_id")
gtfMouse <- merge(
  x = gtfMouse, y = entrez, by.x = "transcript_id", by.y = "transcript_id"
)

gtfMouse <- gtfMouse[gtfMouse$gene_name %in% rownames(rawCounts), ]
gtfMouse <- gtfMouse[!duplicated(gtfMouse$gene_name), ]
rownames(gtfMouse) <- gtfMouse$gene_name
## remove unmatched genes
rawCounts <- rawCounts[gtfMouse$gene_name, ]
genesMetadata <- gtfMouse[rownames(rawCounts), ]

## Formatting names and setting colors
samplesMetadata$TissueDef[samplesMetadata$TissueDef == "Inguinal fat"] <- "Inguinal Fat"

colorsTissues <- c(
  "Brain" = "#b3d88b", 
  "Kidney" = "#e3de5b", 
  "Peritoneal Cavity" = "#fdbf6f", 
  "Inguinal Fat" = "#1b79b6", 
  "Mesenteric Fat" = "#a6cfe5", 
  "Liver" = "#f58020", 
  "Lung" = "#e41e26", 
  "Mesenteric Sheet" = "#cab3d6", 
  "Spleen" = "#6b3f98",
  "Thymus" = "#33a248" 
)

## normalization 
rawCountsF <- rawCounts[rowSums(rawCounts) > 1, ]
keep.exprs <- filterByExpr(rawCountsF, group = samplesMetadata$Groups)
rawCountsF2 <- rawCountsF[keep.exprs, ]
DGEDataF <- DGEList(counts = rawCountsF2, group = samplesMetadata$Groups)
DGEDataF <- calcNormFactors(DGEDataF, method = "TMM")
lcpmTMM <- cpm(DGEDataF, log = TRUE)
genesMetadataF <- genesMetadata[rownames(lcpmTMM), ]

############# taking genesets #############
## GO - Biological Process
listBPGO <- getGOgenes(
  OrgDb = org.Mm.eg.db, selgo = "All", keytype = "SYMBOL", ont = "BP"
)
## KEGG pathways
data(kegg.sets.mm)
listKegg <- lapply(
  X = seq_along(kegg.sets.mm), FUN = function(x) {
    genesMetadataF$gene_name[genesMetadataF$entrez_id %in% kegg.sets.mm[[x]]]
  }
) %>% setNames(names(kegg.sets.mm))

############# analysis based on PCA #############
PCAall <- prcomp(t(lcpmTMM), center = TRUE, scale. = TRUE)

variance <- round(factoextra::get_eigenvalue(PCAall)[c(1, 2), 2], 1)
p <- ggplot(data.frame(PCAall[["x"]]), aes(x = PC1, y = PC2)) +
  geom_point(
    aes(fill = as.factor(samplesMetadata$TissueDef)), 
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
    text = element_text(family = "Helvetica"), 
    aspect.ratio = 1
  ) + scale_x_continuous(limits = symmetric_limits) + 
  scale_y_continuous(limits = symmetric_limits) + 
  ggtitle("PCA human macrophages")
p

## FGSEA on loading vectors with GO terms
loadingVectors <- PCAall$rotation[, 1:3]
rankPC1 <- sort(loadingVectors[, 1], decreasing = TRUE)
rankPC2 <- sort(loadingVectors[, 2], decreasing = TRUE)
rankPC3 <- sort(loadingVectors[, 3], decreasing = TRUE)
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

## FGSEA on loading vectors with KEGG pathways
fgseaPCAKEGG <- lapply(
  X = colnames(loadingVectors),
  FUN = function(x) {
    genes <- sort(loadingVectors[, x], decreasing = TRUE)
    fgseaMultilevel(
      listKegg, stats = genes, minSize = 15, maxSize = 500
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
  listKegg[topPathwaysKEGG1], rankPC1, fgseaPCAKEGG[[1]], 
  gseaParam = 0.5, render = TRUE
)

# PC2
topPathwaysUpKEGG2 <- fgseaPCAKEGG[[2]][NES > 0][head(order(padj), n = 10), pathway]
topPathwaysDownKEGG2 <- fgseaPCAKEGG[[2]][NES < 0][head(order(padj), n = 10), pathway]
topPathwaysKEGG2 <- c(topPathwaysUpKEGG2, rev(topPathwaysDownKEGG2))

grid::grid.newpage()
grDevices::dev.interactive()
plotGseaTable(
  listKegg[topPathwaysKEGG2], rankPC2, fgseaPCAKEGG[[2]], 
  gseaParam = 0.5, render = TRUE
)

# PC3
topPathwaysUpKEGG3 <- fgseaPCAKEGG[[3]][NES > 0][head(order(padj), n = 10), pathway]
topPathwaysDownKEGG3 <- fgseaPCAKEGG[[3]][NES < 0][head(order(padj), n = 10), pathway]
topPathwaysKEGG3 <- c(topPathwaysUpKEGG3, rev(topPathwaysDownKEGG3))

grid::grid.newpage()
grDevices::dev.interactive()
plotGseaTable(
  listKegg[topPathwaysKEGG3], rankPC3, fgseaPCAKEGG[[3]], 
  gseaParam = 0.5, render = TRUE
)

## distribution of OXPHOS genes on PC1's loading vector 
plotEnrichment(
  listKegg$`mmu00190 Oxidative phosphorylation`, rankPC1) + 
  labs(title = "Oxidative phosphorylation (PC1)")+ 
  theme(plot.title = element_text(face = "bold")) + 
  ylab("Enrichment score") + 
  xlab("Rank based on PC1 scores")
plotEnrichment(
  listKegg$`mmu00190 Oxidative phosphorylation`, rankPC2) + 
  labs(title = "Oxidative phosphorylation (PC2)")+ 
  theme(plot.title = element_text(face = "bold")) + 
  ylab("Enrichment score") + 
  xlab("Rank based on PC2 scores")


oxphosFilt <- with(
  listKegg, 
  `mmu00190 Oxidative phosphorylation`[`mmu00190 Oxidative phosphorylation` %in% 
                                         rownames(lcpmTMM)]
)
selGenes <- PCAall$rotation[, 1]
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

pp + annotate("text", x = -0.011, y = 42, label = textKS)

## enrichment analysis: corr. between pathway activities
selGOTerms <- read.csv(file.path(metadataPath, "selected_go_terms.csv"), header = FALSE)
listIntPathways2 <- append(
  getGOgenes(org.Mm.eg.db, selGOTerms[[2]], keytype = "SYMBOL"),
  listKegg[c("mmu00190 Oxidative phosphorylation", 
              "mmu00010 Glycolysis / Gluconeogenesis", 
              "mmu00071 Fatty acid metabolism")]
)
cn <- colnames(lcpmTMM)
gageListNoRef <- list()
reference <- seq(ncol(lcpmTMM))
for (i in unique(samplesMetadata$TissueDef)) {
  sample <- grep(i, samplesMetadata$TissueDef, fixed = TRUE)
  gageListNoRef[[i]] <- gage(
    exprs = lcpmTMM,
    use.fold = TRUE,
    gsets = listIntPathways2,
    ref = reference,
    samp = sample,
    compare = 'as.group'
  )
}
orderPathways <- rownames(gageListNoRef$Lung$stats)
orderPathways <- orderPathways[!is.na(orderPathways)]
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
)]
matrixNoRefGAGE$Tissue <- rownames(matrixNoRefGAGE)
color.line <- rgb(43, 43, 43, max = 255, alpha = 60, names = "blue50")
for (i in c(
  "GO:0009062: fatty acid catabolic process", 
  "GO:0032376: positive regulation of cholesterol transport"
)) {
  p1 <- ggscatter(
    matrixNoRefGAGE, x = "mmu00190 Oxidative phosphorylation", 
    y = i, 
    fill = "Tissue", size = 3*1.5, shape = 21, palette = colorsTissues,
    add = "none", conf.int = TRUE, 
    cor.coef = TRUE, cor.method = "pearson", 
    title = "T-score using log-FC (pseudo-bulk)"
  ) + geom_smooth(method = "lm", colour = color.line, se = FALSE) + 
    # xlab("Oxidative phosphorylation") +
    # ylab(sub(pattern = "GO\\:\\d{7} ", replacement = "", x = i)) +
    theme_classic2() + 
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
  print(p1)
}

############# differential expression analysis with limma #############

mDesign <- model.matrix(~ 0 + samplesMetadata$TissueDef)
colnames(mDesign) <- gsub(
  pattern = "samplesMetadata\\$TissueDef", 
  replacement = "", 
  x = colnames(mDesign)
) %>% gsub(pattern = " ", replacement = "_", x = .)
rownames(mDesign) <- samplesMetadata$Name
DGEfiltDataVoom <- voom(DGEDataF, mDesign, plot = FALSE)
fit <- lmFit(DGEfiltDataVoom, mDesign)
resListNoInt <- list()
tableTotal <- list()
## diff expression analysis
for (group in colnames(mDesign)) {
  contrast <- paste0(
    group, " - ((", paste(colnames(mDesign), collapse = " + "), ") / 10)"
  )
  comp <- makeContrasts(contrast, levels = mDesign)
  resListNoInt[[group]] <- contrasts.fit(fit, contrast = comp)
  resListNoInt[[group]] <- eBayes(resListNoInt[[group]])
  tableTotal[[group]] <- topTable(
    resListNoInt[[group]], sort.by = "P", n = Inf
  )
}
resMTgenesLimma <- lapply(
  X = tableTotal, 
  FUN = function(x) x[grep(pattern = "MT-", x = rownames(x)), ]
)
mtDif <- lapply(
  X = resMTgenesLimma, 
  FUN = function(x) {
    dexp <- rownames(x %>% filter(adj.P.Val <= 0.05))
    if (identical(dexp, character(0))) {
      return(NULL)
    } else {
      dexp
    }
  }
)

############# heatmaps #############

## TFAM
intGenesAgg <- aggregate(
  x = t(lcpmTMM["TFAM", , drop = FALSE]), 
  by = list(samplesMetadata$TissueDef), FUN = mean
)
tfamExprs <- lcpmTMM["TFAM", , drop = FALSE]
tfamRefMean <- rowMeans(tfamExprs)
tfamFC <- sapply(
  X = unique(samplesMetadata$TissueDef), FUN = function(tissue) {
    samp.index <- which(samplesMetadata$TissueDef == tissue)
    rowMeans(tfamExprs[, samp.index, drop = FALSE]) - tfamRefMean  
  }
)
tfamFCmatrix <- t(matrix(tfamFC, dimnames = list(names(tfamFC), "TFAM")))
heatmapFinal <- ComplexHeatmap::Heatmap(
  tfamFCmatrix, 
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 8), name = "LogFC",
  column_title_gp = gpar(fontface = "bold"), 
  column_title = "Mouse TMFs",
  column_names_rot = 90,
  cluster_columns = TRUE, 
  show_row_dend = FALSE,
  row_names_side = "left",
  border_gp = gpar(col = "black"),
  height = nrow(tfamFCmatrix)*unit(8, "mm"),
  heatmap_legend_param = list(legend_height = unit(3, "cm"))
)
heatmapFinal

## MT genes
mtGenes <- grep(pattern = "MT-", x = rownames(lcpmTMM), value = TRUE)
mtExprs <- lcpmTMM[mtGenes, ]
mtRefMean <- rowMeans(mtExprs)
mtgenesFC <- sapply(
  X = unique(samplesMetadata$TissueDef), FUN = function(tissue) {
    samp.index <- which(samplesMetadata$TissueDef == tissue)
    rowMeans(mtExprs[, samp.index]) - mtRefMean  
  }
)
orderSorted <- sort(colMeans(mtgenesFC), decreasing = FALSE)
mtgenesFC <- mtgenesFC[, names(orderSorted)]

col_fun <- circlize::colorRamp2(c(-4, 0, 4), c("blue", "white", "red"))

columnAvg2 <- HeatmapAnnotation(
  Average = anno_simple(
    colMeans(mtgenesFC), col = col_fun, border = TRUE
  ),
  Average2 = anno_text(
    as.character(round(colMeans(mtgenesFC), 3)), 
    gp = gpar(fontsize = 8, col = c(rep("white", 3), rep("black", 7))), 
    rot = 0, just = 0.5,
    offset = 2.4
  ), 
  ColNames = anno_text(
    colnames(mtgenesFC), 
    rot = 45, 
    offset = unit(1.1, "npc"), 
    just = "right",
    gp = gpar(fontsize = 8)
  ),
  col = list(Average = col_fun),
  annotation_name_gp = gpar(fontsize = 9, fontface = "bold")
)
ComplexHeatmap::Heatmap(
  mtgenesFC, 
  show_column_names = FALSE,
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 8, just = "centre"), 
  name = "LogFC",
  column_title_gp = gpar(fontface = "bold"), 
  column_title = "Detected mitochondrial genes across tissues",
  column_names_rot = 0,
  border = TRUE,
  cluster_columns = FALSE, 
  show_row_dend = FALSE,
  heatmap_width = unit(150, "mm"),
  bottom_annotation = columnAvg2,
  cell_fun = function(j, i, x, y, width, height, fill) {
    if (
      rownames(mtgenesFC)[i] %in% mtDif[[gsub(
        pattern = " ", replacement = "_", x = colnames(mtgenesFC)[j]
      )]]
    ) {
      diff <- "*"
      if(mtgenesFC[i, j] < -0.5) {
        grid.text(
          paste0(sprintf("%.2f", mtgenesFC[i, j]), diff), x, y, 
          gp = gpar(fontsize = 8, col = "white", fontface = "bold")
        )  
      } else {
        grid.text(
          paste0(sprintf("%.2f", mtgenesFC[i, j]), diff), x, y, 
          gp = gpar(fontsize = 8, col = "black", fontface = "bold")
        )
      }
    } else {
      diff <- ""
      if(mtgenesFC[i, j] < -0.5) {
        grid.text(
          paste0(sprintf("%.2f", mtgenesFC[i, j]), diff), x, y, 
          gp = gpar(fontsize = 8, col = "white")
        )  
      } else {
        grid.text(
          paste0(sprintf("%.2f", mtgenesFC[i, j]), diff), x, y, 
          gp = gpar(fontsize = 8, col = "black")
        )
      }
    }
  }
)
