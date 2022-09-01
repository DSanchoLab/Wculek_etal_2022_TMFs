################################################################################
############################### Helper functions ###############################
################################################################################
# Author: Diego Mañanes
# Date: 22/09/22
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

symmetric_limits <- function(x) {
  max <- max(abs(x))
  c(-max, max)
}
