################################################################################
############################### Helper functions ###############################
################################################################################
# Author: Diego Mañanes
# Date: 22/09/22
# Description: helper functions for transcriptomics analysis
###############################################################################

## list of custom colors (using RColorBrewer)
color.list <- function() {
  color.list.2 <- c(
    RColorBrewer::brewer.pal(12, "Paired"), "#d45b91", "#374738",
    RColorBrewer::brewer.pal(12, "Set3"),
    RColorBrewer::brewer.pal(8, "Pastel2"),
    "#333333", "#5D5D5D",
    "#888888", "#B3B3B3"
  )
  color.list.2[11] <- "#e3dc5b"
  color.list.2[15] <- "#60c4b4"
  return(color.list.2)
}


## get genes of GO terms using AnnotationDbi R package
getGOgenes <- function(
  OrgDb, 
  selgo, 
  keytype,
  ont = "All"
) {
  kt <- keytypes(OrgDb)
  if (!keytype %in% kt) stop("keytype is not supported...")
  goterms <- Ontology(GO.db::GOTERM)
  gotermsNames <- Term(GO.db::GOTERM)
  if (selgo != "All") {
    goterms <- goterms[names(goterms) %in% selgo]
  }
  if (ont != "All") {
    goterms <- goterms[goterms %in% ont]
  }
  go2gene <- suppressMessages(
    mapIds(
      OrgDb, keys = names(goterms), column = keytype,
      keytype = "GOALL", multiVals = 'list'
    )
  )
  goAnno <- stack(go2gene)
  colnames(goAnno) <- c(keytype, "GOALL")
  goAnno <- unique(goAnno[!is.na(goAnno[,1]), ])
  goAnno$ONTOLOGYALL <- goterms[goAnno$GOALL]
  bpGOF <- suppressWarnings(
    lapply(go2gene, function(x) if (!is.na(x)) return(unique(toupper(x))))
  )
  bpGOF[sapply(bpGOF, is.null)] <- NULL
  names(bpGOF) <- paste(names(bpGOF), gotermsNames[names(bpGOF)], sep = ": ")
  return(bpGOF)
}

## ggplot theme for plots
ggtheme_custom <- theme(
  plot.title = element_text(face = "bold", hjust = 0),
  legend.title = element_text(face = "bold", hjust = 0.5)
) 

## set of pseudo-bulk functions (DISCLAIMER: not intended to be general)
# limit of cells
setLimit <- function(vec, limit) {
  subs <- sum(vec) - limit
  if (subs > 0) {
    for (i in seq(abs(subs))) {
      loc <- sample(x = seq_along(vec), size = 1)
      vec[loc] <- vec[loc] - 1
    }
  } else if (subs < 0) {
    for (i in seq(abs(subs))) {
      loc <- sample(x = seq_along(vec), size = 1)
      vec[loc] <- vec[loc] + 1
    }
  }
  return(vec)
}
# split a group of cells into a desired number 
splitNumCells <- function(totalNum, numRep, col = "Number") {
  lapply(
    seq_along(ceiling(totalNum[, col] / numRep)), 
    FUN = function(x, num.rep) {
      return(
        setLimit(
          vec = rep(ceiling(totalNum[, col] / 3)[x], num.rep), 
          limit = totalNum[x, col]
        )
      )
    }, num.rep = numRep
  )
}
# selection of the number of cells provided randomly
splitCells <- function(totalCells, listNum) {
  listSel <- lapply(
    X = seq_along(listNum), FUN = function(x, totalCells) {
      cells <- totalCells %>% filter(Sample %in% names(listNum[x]))
      listSel <- list()
      n <- 1
      for (i in listNum[[x]]) {
        sel <- sample(x = rownames(cells), size = i)
        listSel[[n]] <- sel
        cells <- cells[-which(rownames(cells) %in% sel), ]
        n <- n + 1
      }
      return(listSel)
    }, totalCells = totalCells
  )
  names(listSel) <- names(listNum)
  return(listSel)
}
# generate the final expression matrix (raw counts)
generateFinalMatrix <- function(groupCells, totalCounts) {
  finalList <- lapply(
    X = seq_along(groupCells), 
    FUN = function(x) {
      finalMatrix <- matrix(
        NA, nrow = nrow(totalCounts), ncol = length(groupCells[[x]])
      )
      for (i in seq_along(groupCells[[x]])) {
        finalMatrix[, i] <- rowSums(totalCounts[, groupCells[[x]][[i]]])
      }
      colnames(finalMatrix) <- paste(
        names(groupCells)[x], seq(length(groupCells[[x]])), sep = "_"
      )
      rownames(finalMatrix) <- rownames(totalCounts)
      return(finalMatrix)
    }
  )
  return(do.call(cbind, args = finalList))
}
# generate the final expression matrix (avg logCPMs)
generateFinalMatrixCPMs <- function(groupCells, totalCounts) {
  finalList <- lapply(
    X = seq_along(groupCells), 
    FUN = function(x) {
      finalMatrix <- matrix(
        NA, nrow = nrow(totalCounts), ncol = length(groupCells[[x]])
      )
      for (i in seq_along(groupCells[[x]])) {
        finalMatrix[, i] <- rowMeans(totalCounts[, groupCells[[x]][[i]]])
      }
      colnames(finalMatrix) <- paste(
        names(groupCells)[x], seq(length(groupCells[[x]])), sep = "_"
      )
      rownames(finalMatrix) <- rownames(totalCounts)
      return(finalMatrix)
    }
  )
  return(do.call(cbind, args = finalList))
}

## symmetric axes for ggplot
symmetric_limits <- function(x) {
  max <- max(abs(x))
  c(-max, max)
}

## density PCA plots
densityPCAscores <- function(pcaObj, PC, pathway) {
  selGenes <- pcaObj$rotation[, PC]
  dfBarplotPCA <- data.frame(
    Genes = factor(names(selGenes), levels = names(selGenes)),
    LoadingVectors = selGenes,
    Pathway = names(selGenes) %in% pathway
  )
  p <- ggplot(data = dfBarplotPCA, aes(x = LoadingVectors, fill = Pathway)) + 
    geom_density(alpha = 0.6, position = 'identity', color = "white") + 
    theme_classic() + 
    ggtitle(paste0("Loading scores (", PC, ")")) + 
    theme(plot.title = element_text(face = "bold"))
  return(p)
}



