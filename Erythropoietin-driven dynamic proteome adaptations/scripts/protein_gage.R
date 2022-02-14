
library(gage)
library(GO.db)
library(org.Mm.eg.db)
library(openxlsx)
library(pheatmap)
library(viridis)
library(msigdbr)

# BUILD DB
m_df <- msigdbr(species = "Mus musculus")
m_df$gs_subcat <- gsub("\\:", "_", m_df$gs_subcat)
m_df$gs_subcat <- paste0(m_df$gs_cat, "_", m_df$gs_subcat)
m_df$gs_subcat <- gsub("_$", "", m_df$gs_subcat)

subcat.unique <- sort(unique(m_df$gs_subcat))

dbList <- lapply(subcat.unique, function(i){
  m_df.sub <- m_df[m_df$gs_subcat == i, ]
  terms <- unique(m_df.sub$gs_name)
  db <- mclapply(terms, function(j){
    return(m_df.sub$entrez_gene[m_df.sub$gs_name == j])
  }, mc.cores = 8)
  names(db) <- terms
  return(db)
})
names(dbList) <- subcat.unique


############################################
###                                      ###
###               FUNCTION               ###
###                                      ###
############################################

entrez2symbol <- function(entrez)
{
  symbol <- mget(as.character(entrez), org.Mm.egSYMBOL, ifnotfound=NA)
  symbol <- unlist(lapply(symbol, function(i) return(i[1])))
  return(symbol)
}


symbol2entrez <- function(symbol)
{
  entrez <- mget(as.character(symbol), org.Mm.egSYMBOL2EG, ifnotfound=NA)
  entrez <- unlist(lapply(entrez, function(i) return(i[1])))
  entrez <- unique(entrez[!is.na(entrez)])
  return(entrez)
}

symbol2entrezV2 <- function(symbol)
{
  entrez <- mget(as.character(symbol), org.Mm.egSYMBOL2EG, ifnotfound=NA)
  entrez <- unlist(lapply(entrez, function(i) return(i[1])))
  #entrez <- unique(entrez[!is.na(entrez)])
  return(entrez)
}


gagePipeline <- function(test, ref, db, up = NULL, down = NULL)
{
  gageMat <- cbind(ref, test)
  gos <- gage(gageMat, gsets = db, ref = seq(1:ncol(ref)),
              same.dir = T, compare = "unpaired", rank = T,
              saaTest = gs.tTest, set.size = c(5, 500))
  significant.groups <- sigGeneSet(gos, cutoff = 1, qpval = c("q.val"))
  
  # Greater
  greater <- significant.groups$greater
  greater.entrez <- lapply(rownames(greater), function(i) return(intersect(db[[i]], up)))
  greater.symbol <- lapply(greater.entrez, entrez2symbol)
  greater.nb <- lapply(greater.entrez, length)
  
  greater.entrez <- unlist(lapply(greater.entrez, toString))
  greater.symbol <- unlist(lapply(greater.symbol, toString))
  greater.nb <- unlist(greater.nb)
  greater <- cbind(greater, nb = greater.nb, entrez = greater.entrez, symbol = greater.symbol)
  
  # Less
  less <- significant.groups$less
  less.entrez <- lapply(rownames(less), function(i) return(intersect(db[[i]], down)))
  less.symbol <- lapply(less.entrez, entrez2symbol)
  less.nb <- lapply(less.entrez, length)
  
  less.entrez <- unlist(lapply(less.entrez, toString))
  less.symbol <- unlist(lapply(less.symbol, toString))
  less.nb <- unlist(less.nb)
  less <- cbind(less, nb = less.nb, entrez = less.entrez, symbol = less.symbol)
  
  return(list(greater = greater, less = less))
}



my.read_xlsx <- function(inFile)
{
  mysheets <- getSheetNames(inFile)
  mList <- lapply(mysheets, read.xlsx, xlsxFile = inFile)
  names(mList) <- mysheets
  return(mList)
}

getHeatmapMatrix <- function(gsea_list, adjusted = FALSE)
{
  trms <- sort(unique(unlist(lapply(gsea_list, function(i) i[,1]))))
  m <- matrix(NA, nrow = length(trms), ncol = length(gsea_list))
  
  for(i in 1:length(gsea_list))
  {
    gs <- gsea_list[[i]]
    idx <- match(gs[,1], trms)
    if(adjusted)m[idx, i] <- as.numeric(gs$"q.val")		
    else(m[idx, i] <- as.numeric(gs$"p.val"))
  }
  m[is.na(m)] <- 1	
  rownames(m) <- trms
  colnames(m) <- names(gsea_list)
  return(m)
}

gageHeatmap <- function(fh_list, outFile, nb, adjusted = FALSE)
{
  gsMat <- getHeatmapMatrix(fh_list, adjusted)	
  colnames(gsMat) <- names(fh_list)	
  
  # re-order gsMat columns
  idx.up <- grep("UP$", colnames(gsMat))
  idx.down <- grep("DOWN$", colnames(gsMat))
  gsMat <- gsMat[, c(idx.up, idx.down)]
  
  # add annotation
  ann.col <- data.frame(Sign = c(rep("UP", length(idx.up)), rep("DOWN", length(idx.down))))
  rownames(ann.col) <- colnames(gsMat)
  myColor.ann <- list(Sign = c(UP = "red", DOWN = "blue"))
  
  # select top X gene-sets per column
  idxMat <- apply(gsMat, 2, order)
  idxMat <- idxMat[1:nb, ]
  gsMat <- gsMat[unique(as.numeric(idxMat)), ]	
  
  gsMat <- -log10(gsMat)	
  gsMat[gsMat > 5] <- 5 # set the limit to 5
  
  paletteLength <- 10
  myColor <- magma(paletteLength)
  #myColor <- viridis(paletteLength)
  #myColor <- colorRampPalette(c("white", "red"))(paletteLength)
  myBreaks <- c(seq(0, max(gsMat), length.out=ceiling(paletteLength)))
  
  doClust <- ifelse(nrow(gsMat) > 1, TRUE, FALSE)
  pheatmap(gsMat, color = myColor, breaks = myBreaks, filename = outFile,
           cluster_cols = FALSE, cluster_rows = doClust,
           annotation_col = ann.col,
           annotation_colors = myColor.ann,
           cellwidth = 10, cellheight = 10, fontsize_row = 8, fontsize_col = 8,
           width = 10, height = 10)
}



############################################
###                                      ###
###                 MAIN                 ###
###                                      ###
############################################


dataDir <- file.path("../data/")
limmaDir <- file.path("../protein_limma/")
gageDir <- file.path("../protein_gage/")
dir.create(gageDir)

##########################################################################################
# CFU

# EXPRESSION
setwd(dataDri)
expMat <- read.xlsx("EPO_proteins_log2.xlsx", sheet = "CFU",rowNames = TRUE)
expMat.entrez <- symbol2entrezV2(rownames(expMat))

toKeep <- !is.na(expMat.entrez) & !(duplicated(expMat.entrez))
expMat <- expMat[toKeep,]
rownames(expMat) <- expMat.entrez[toKeep]

cond <- c("WT", "WT", "WT", "KO", "KO", "KO")



lapply(1:length(dbList), function(i){
  
  mydb <- dbList[[i]]
  mydb.name <- names(dbList)[i]
  
  #####################################
  # WT vs. KO
  
  setwd(limmaDir)
  limmaMat <- read.xlsx("CFU_WT-KO_limma.xlsx", sheet=1)
  limmaFC <- limmaMat$logFC
  limmaUP <- limmaMat$entrez[limmaFC > 0]
  limmaDOWN <- limmaMat$entrez[limmaFC < 0]
  
  testMat <- expMat[,cond=="WT"]
  refMat <- expMat[,cond=="KO", drop = FALSE]
  gage.res <- gagePipeline(testMat, refMat, mydb, limmaUP, limmaDOWN)
  
  setwd(file.path(gageDir))
  write.xlsx(list(UP = gage.res$greater, DOWN = gage.res$less),
             paste("CFU_WT-KO_", mydb.name, "_gage.xlsx", sep = ""), row.names = TRUE)
  
})

##########################################################################################
# placenta

# EXPRESSION
setwd(dataDir)
expMat <- read.xlsx("EPO_proteins_log2.xlsx", sheet = "placenta", rowNames = TRUE)
expMat.entrez <- symbol2entrezV2(rownames(expMat))

toKeep <- !is.na(expMat.entrez) & !(duplicated(expMat.entrez))
expMat <- expMat[toKeep,]
rownames(expMat) <- expMat.entrez[toKeep]

cond <- c("WT", "WT", "WT", "KO", "KO")

lapply(1:length(dbList), function(i){
  
  mydb <- dbList[[i]]
  mydb.name <- names(dbList)[i]
  
  #####################################
  # WT vs. KO
  
  setwd(limmaDir)
  limmaMat <- read.xlsx("placenta_WT-KO_limma.xlsx", sheet=1)
  limmaFC <- limmaMat$logFC
  limmaUP <- limmaMat$entrez[limmaFC > 0]
  limmaDOWN <- limmaMat$entrez[limmaFC < 0]
  
  testMat <- expMat[,cond=="WT"]
  refMat <- expMat[,cond=="KO", drop = FALSE]
  gage.res <- gagePipeline(testMat, refMat, mydb, limmaUP, limmaDOWN)
  
  setwd(file.path(gageDir))
  write.xlsx(list(UP = gage.res$greater, DOWN = gage.res$less),
             paste("placenta_WT-KO_", mydb.name, "_gage.xlsx", sep = ""), row.names = TRUE)
  
})

##########################################################################################
# liver

# EXPRESSION
setwd(dataDir)
expMat <- read.xlsx("EPO_proteins_log2.xlsx", sheet = "liver", rowNames = TRUE)
expMat.entrez <- symbol2entrezV2(rownames(expMat))

toKeep <- !is.na(expMat.entrez) & !(duplicated(expMat.entrez))
expMat <- expMat[toKeep,]
rownames(expMat) <- expMat.entrez[toKeep]

cond <- c("WT", "WT", "WT", "WT", "KO", "KO", "KO")

lapply(1:length(dbList), function(i){
  
  mydb <- dbList[[i]]
  mydb.name <- names(dbList)[i]
  
  #####################################
  # WT vs. KO
  
  setwd(limmaDir)
  limmaMat <- read.xlsx("liver_WT-KO_limma.xlsx", sheet=1)
  limmaFC <- limmaMat$logFC
  limmaUP <- limmaMat$entrez[limmaFC > 0]
  limmaDOWN <- limmaMat$entrez[limmaFC < 0]
  
  testMat <- expMat[,cond=="WT"]
  refMat <- expMat[,cond=="KO", drop = FALSE]
  gage.res <- gagePipeline(testMat, refMat, mydb, limmaUP, limmaDOWN)
  
  setwd(file.path(gageDir))
  write.xlsx(list(UP = gage.res$greater, DOWN = gage.res$less),
             paste("liver_WT-KO_", mydb.name, "_gage.xlsx", sep = ""), row.names = TRUE)
  
})


##########################################################################################
# HEATMAPS

setwd(gageDir)
lapply(1:length(dbList), function(i){
  mydb.name <- names(dbList)[i]
  
  # load fisher results
  fhFiles <- c(paste0("CFU_WT-KO_", mydb.name, "_gage.xlsx"),
               paste0("liver_WT-KO_", mydb.name, "_gage.xlsx"),
               paste0("placenta_WT-KO_", mydb.name, "_gage.xlsx")
  )
  fhList <- lapply(fhFiles, my.read_xlsx)
  names(fhList) <- gsub(paste0("_", mydb.name, "_gage.xlsx"), "", fhFiles)
  fhList <- do.call(c, fhList)
  
  # plot heatmap
  gageHeatmap(fhList, paste0("WT-KO_", mydb.name, "_gage_heatmap.pdf"), nb = 10, adjusted = FALSE)
  
})



