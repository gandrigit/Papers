
library(org.Hs.eg.db)
library(openxlsx)
library(pheatmap)
library(GO.db)
library(fgsea)
library(GeneAnswers)
library(viridis)
library(gridExtra)
library(msigdbr)
library(ggplot2)
library(ungeviz)
library(limma)
library(tidyr)


source("tools.r")

# BUILD DB
m_df <- msigdbr(species = "Homo sapiens")
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

my.read_xlsx <- function(inFile)
{
  mysheets <- getSheetNames(inFile)
  mList <- lapply(mysheets, read.xlsx, xlsxFile = inFile)
  names(mList) <- mysheets
  return(mList)
}


fgsea_pipeline <- function(rank, pathways, minSize=15, maxSize=500, nperm=10000){
  fg <- fgsea(rank, pathways = pathways, minSize=minSize, maxSize=maxSize, nperm=nperm)
  fg <- as.data.frame(fg)
  fg.entrez <- fg$leadingEdge
  fg.symbol <- lapply(fg.entrez, entrez2symbol)
  fg.symbol <- unlist(lapply(fg.symbol, function(i) paste(i, collapse = ",")))
  avg.delta <- unlist(lapply(fg.entrez, function(i) mean(rank[i])))
  fg <- data.frame(fg, leadingEdge.symbol = fg.symbol, avgRank = avg.delta)
  fg <- fg[order(fg$pval), ]
  return(fg)
}


symbol2entrez <- function(symbol)
{
  entrez <- mget(as.character(symbol), org.Hs.egSYMBOL2EG, ifnotfound=NA)
  entrez <- unique(unlist(lapply(entrez, function(i) return(i[1]))))
  #entrez <- entrez[!is.na(entrez)]
  return(entrez)
}

entrez2symbol <- function(entrez)
{
  symbol <- mget(as.character(entrez), org.Hs.egSYMBOL, ifnotfound=NA)
  symbol <- unique(unlist(lapply(symbol, function(i) return(i[1]))))
  return(symbol)
}

ensembl2entrez <- function(ensembl)
{
  entrez <- mget(as.character(ensembl), org.Hs.egENSEMBL2EG, ifnotfound=NA)
  entrez <- lapply(entrez, function(i) return(i[1]))
  return(unlist(entrez))
}

refseq2entrez <- function(refseq){
  entrez <- mget(as.character(refseq), org.Hs.egREFSEQ2EG, ifnotfound=NA)
  entrez <- unique(unlist(lapply(entrez, function(i) return(i[1]))))
  return(entrez)
}






getHeatmapMatrix <- function(gsea_list, scoreName)
{
  trms <- sort(unique(unlist(lapply(gsea_list, function(i) i[,1]))))
  m <- matrix(NA, nrow = length(trms), ncol = length(gsea_list))
  
  for(i in 1:length(gsea_list))
  {
    gs <- gsea_list[[i]]
    idx <- match(unlist(gs[,1]), trms)
    m[idx, i] <- as.numeric(gs[, scoreName])
  }
  m[is.na(m)] <- 0	
  rownames(m) <- trms
  colnames(m) <- names(gsea_list)
  return(m)
}

getHeatmapMatrix.pv <- function(gsea_list, adjusted = FALSE)
{
  trms <- sort(unique(unlist(lapply(gsea_list, function(i) i[,1]))))
  m <- matrix(NA, nrow = length(trms), ncol = length(gsea_list))
  
  for(i in 1:length(gsea_list))
  {
    gs <- gsea_list[[i]]
    idx <- match(unlist(gs[,1]), trms)
    if(adjusted)m[idx, i] <- as.numeric(gs$"pval")		
    else(m[idx, i] <- as.numeric(gs$"padj"))
  }
  m[is.na(m)] <- 1	
  rownames(m) <- trms
  colnames(m) <- names(gsea_list)
  return(m)
}



	

############################################
###                                      ###
###                 MAIN                 ###
###                                      ###
############################################


# PARAMETERS
fgseaDir <- file.path("../TCGA_fgsea/")
dir.create(fgseaDir, recursive = TRUE, showWarnings = FALSE)

##################################
# Rank genes based on limma output

####################################
# LOAD TCGA COAD DATA
setwd(file.path("../data"))
rna <- read.xlsx("rna_COAD_allGenes.xlsx", sheet = 1, rowNames = TRUE)
prot <- read.xlsx("prot_COAD_allGenes.xlsx", sheet = 1, rowNames = TRUE)
prot[prot == 0] <- NA

colnames(prot) <- substring(colnames(prot), 1, 16)

idx <- sapply(colnames(prot), grep, x = colnames(rna))
idx <- sapply(idx, function(i){
  if(length(i) == 0) return(NA)
  return(i[1])
})

prot <- prot[, !is.na(idx)]
idx <- idx[!is.na(idx)]
rna <- rna[, idx]
colnames(rna) <- colnames(prot)

rna.tcga <- rna
prot.tcga <- prot

colnames(rna.tcga) <- substring(colnames(rna.tcga), 1, 12)
colnames(prot.tcga) <- substring(colnames(prot.tcga), 1, 12)


####################################
# LOAD HEALTHY COLON DATA
setwd(file.path("~/Research/Sajib/P2_half_lives/data"))
rna <- read.xlsx("rna_colon_allGenes.xlsx", sheet = 1, rowNames = TRUE)
prot <- read.xlsx("prot_colon_allGenes.xlsx", sheet = 1, rowNames = TRUE)
prot[prot == 0] <- NA

rna.ctl <- rna
prot.ctl <- prot

colnames(rna.ctl) <- substring(gsub("TCGA-COAD.", "", colnames(rna.ctl)), 1, 12)


####################################
# INTERSECT TCGA AND CTL GENES

commonGenes.rna <- intersect(rownames(rna.tcga), rownames(rna.ctl))
rna.tcga <- rna.tcga[commonGenes.rna, ]
rna.ctl <- rna.ctl[commonGenes.rna, ]

commonGenes.prot <- intersect(rownames(prot.tcga), rownames(prot.ctl))
prot.tcga <- prot.tcga[commonGenes.prot, ]
prot.ctl <- prot.ctl[commonGenes.prot, ]


# SUBSTRACT AVG CTL INTENSITY
expr <- rna.tcga - rowMeans(rna.ctl, na.rm = TRUE)

# ADD GENE META DATA
symbol <- rownames(expr)
entrezID <- symbol2entrez(symbol)

expr <- expr[!is.na(entrezID), ]
symbol <- symbol[!is.na(entrezID)]
entrezID <- entrezID[!is.na(entrezID)]

# ENSEMBL TO ENTREZ (IQR)
expr.IQR <- apply(expr, 1, IQR, na.rm = TRUE)
expr <- lapply(unique(entrezID), function(i){
  symbol.current <- symbol[entrezID == i]
  if(length(symbol.current) == 1) expr[symbol.current,] else{
    expr[names(which.max(expr.IQR[symbol.current])), ]
  }
})
expr <- do.call(rbind, expr)
rownames(expr) <- unique(entrezID)

# CPM MATRIX TO RANKED GENE LIST
rankedGenesList <- lapply(1:ncol(expr), function(i){
  r <- as.numeric(expr[, i])
  names(r) <- rownames(expr)
  return(sort(r))
})
names(rankedGenesList) <- colnames(expr)


#########################
# Perform fgsea analysis

dbList <- dbList[c("H")]

# FGSEA
setwd(fgseaDir)
lapply(1:length(dbList), function(i){
  mydb <- dbList[[i]]
  mydb.name <- names(dbList)[i]
  
  plotDir <- fgseaDir   
  
  dir.create(file.path(plotDir, mydb.name), showWarnings = FALSE)
  setwd(file.path(plotDir, mydb.name))   
  
  #########################
  # Perform fgsea analysis
  
  fgseaResList <- lapply(rankedGenesList,
                         fgsea, pathways = mydb, minSize=15, maxSize=1000, nperm=10000
  )
  fgseaResList <- lapply(fgseaResList, function(i) i[order(i$pval),])	
  
  # Add symbols
  fgseaResList <- lapply(fgseaResList, as.data.frame)
  
  fgseaResList <- lapply(fgseaResList, function(k){
    k.symbol <- lapply(k$leadingEdge, entrez2symbol)
    k.symbol <- lapply(k.symbol, function(j) unique(j[!is.na(j)]))
    k.symbol <- lapply(k.symbol, toString)
    k$leadingEdge.symbol <- unlist(k.symbol)
    k
  })
  
  fgseaResList <- lapply(fgseaResList, function(j){
    j$NES[is.na(j$NES)] <- 0
    return(j)
  })
  
  # Divide UP and DOWN
  fgseaResList.UP <- lapply(fgseaResList, function(j) j[j$NES > 0, ])
  names(fgseaResList.UP) <- paste0(names(fgseaResList), ".UP")
  fgseaResList.DOWN <- lapply(fgseaResList, function(j) j[j$NES < 0, ])
  names(fgseaResList.DOWN) <- paste0(names(fgseaResList), ".DOWN")
  fgseaResList <- c(fgseaResList.UP, fgseaResList.DOWN)
  
  # Save fgsea output   
  lapply(1:length(rankedGenesList), function(j)
    write.xlsx(list(UP =  fgseaResList.UP[[j]], DOWN = fgseaResList.DOWN[[j]]),
               paste(names(rankedGenesList)[j], "_", mydb.name, "_fgsea.xlsx", sep = ""),
               row.names = FALSE, firstRow = T, headerStyle = createStyle(textDecoration = 'bold'))
  )	
  
})




#######################################
# SCORE MATRICES

setwd(fgseaDir)
lapply(1:length(dbList), function(i){
  mydb.name <- names(dbList)[i]
  
  # load fisher results
  fhFiles <- file.path(mydb.name, paste0(comp, "_", mydb.name, "_fgsea.xlsx"))
  names(fhFiles) <- comp
  fhList <- lapply(fhFiles, my.read_xlsx)
  fhList <- lapply(fhList, function(j) do.call(rbind, j))
  
  # BUILD SCORE MATRIX
  scMat <- getHeatmapMatrix(fhList, "NES")
  
  write.xlsx(scMat, file.path(mydb.name, paste0(mydb.name, "_NES_matrix.xlsx")), row.names = TRUE)
})

