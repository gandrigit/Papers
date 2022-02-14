library(org.Mm.eg.db)
library("annotate")
library("genefilter")
library(sva)
library(labdsv)
library(gage)
library(limma)
library(colorRamps)
library(GeneAnswers)
library(Cairo)

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



################
# function
getGeneSymFromDbTerm <- function(dbTerm, ids)
{
  geneMat <- c()
  for (term in dbTerm){
    geneID <- mydb[[term]]
    # Map Limma results on GeneID
    geneID <- intersect(ids, geneID)
    geneSym <- as.character(mget(geneID,org.Mm.egSYMBOL,ifnotfound=NA))
    geneSym <- unique(geneSym)
    nb <- length(geneSym)
    geneSym <- toString(geneSym)
    geneMat <- rbind(geneMat, c(nb,geneSym))
  }
  geneMat <- as.matrix(geneMat)
  rownames(geneMat) <- dbTerm
  colnames(geneMat) <- c("nb","GeneList")
  return(geneMat)
}

gageAnalysis <- function(co, fileN)
{
  res <- topTable(fit2, coef=co, adjust="BH", number=nrow(mymat), sort.by="logFC")
  if (is.null(res$ID)) res_names <- rownames(res)
  else res_names <- res$ID
  limmaUp <- c()
  limmaDown <- c()
  for (i in 1:nrow(res)) {
    if (res$logFC[i]>0) limmaUp <- c(limmaUp, res_names[i])
    else limmaDown <- c(limmaDown, res_names[i])
  }
  
  for (i in 1:length(cutOffValues)) { 
    significant.groups=sigGeneSet(gos,cutoff=cutOffValues[i], qpval = c("q.val"))
    #path <- paste(gageD, cfuD, cutOffD[i], fileN, sep = "")
    path <- fileN
    
    # greater
    greater <- significant.groups$greater
    if (!nrow(greater)==0)
    {
      # get gene list
      gMat <- getGeneSymFromDbTerm(rownames(greater), limmaUp)
      finalMat <- cbind(greater, gMat)
      write.table(finalMat, paste(path, "_greater.txt", sep = ""), sep = "\t", quote = F)
      
      # export to xls
      finalMat <- as.data.frame(finalMat)
      finalMat <- cbind(rownames(finalMat),finalMat)
      colnames(finalMat)[1] <- db
      try(WriteXLS("finalMat", paste(path, "_greater.xls", sep = "")))
    }
    
    # less
    less <- significant.groups$less
    if (!nrow(less)==0)
    {
      # get gene list
      gMat <- getGeneSymFromDbTerm(rownames(less), limmaDown)
      finalMat <- cbind(less, gMat)
      write.table(finalMat, paste(path, "_less.txt", sep = ""), sep = "\t", quote = F)
      
      # export to xls
      finalMat <- as.data.frame(finalMat)
      finalMat <- cbind(rownames(finalMat),finalMat)
      colnames(finalMat)[1] <- db
      try(WriteXLS("finalMat", paste(path, "_less.xls", sep = "")))
    }
    
    # PLot the results as a network
    try(plotGAGE(greater, paste(path,"_greater.graphml", sep = "")))
    try(plotGAGE(less, paste(path,"_less.graphml", sep = "")))
  }
}

plotGAGE <- function(gene_sets, graph_name)
{
  group.plot(gene_sets, thresh=0, largest=F, groups=mydb, score.cutoff=1,outputFile=graph_name)
}

get_limma_up <- function(limma_res, cutoff)
{
  limmaUp <- c()
  for (i in 1:length(limma_res[,1]))
  {
    if (limma_res[i,1]>cutoff) limmaUp <- c(limmaUp, rownames(limma_res)[i])
  }	
  return(limmaUp)
}

get_limma_down <- function(limma_res, cutoff)
{
  limmaDown <- c()
  for (i in 1:length(limma_res[,1]))
  {
    if ((limma_res[i,1]<0) && (abs(limma_res[i,1]) > cutoff)) limmaDown <- c(limmaDown, rownames(limma_res)[i])
  }	
  return(limmaDown)
}


# Prepare directory
dataD <- file.path("../prediction")
gageD <- file.path("../prediction_gage")
dir.create(gageD)
cutOffValues <- c(0.05)

######################
# PRED PROTEOME

# REMARK:
# first use protein_prediction.R to generate predicted_protein_intensity.txt

# Read proteome data
setwd(dataD)
raw_prot_data <- read.delim("predicted_protein_intensity.txt", sep="\t", header=FALSE)
genes <- as.character(raw_prot_data[c(2:nrow(raw_prot_data)),1])

# Limma proteomics KO Vs WT
mymat <- as.matrix(raw_prot_data[c(2:nrow(raw_prot_data)), c(2:ncol(raw_prot_data))])
mymat <-apply(mymat,2,as.numeric)
rownames(mymat) <- genes
rownames(mymat) <- as.character(mget(rownames(mymat),org.Mm.egSYMBOL2EG,ifnotfound=NA))


# AGAINST 0to1h TIME POINTS
contrast.matrix <- makeContrasts(e2-e1, e3-e1, e4-e1, e5-e1, e6-e1,levels=design)
# fit and get up and down genes
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)


############      
fileName <- "gage_pred_prot_Epo_2to3_Vs_0to1"
gos <- gage(mymat[,c(1:4)],gsets =mydb,ref=c(1:2),same.dir=T,compare = "unpaired",rank=T, saaTest=gs.KSTest)
gageAnalysis(1, fileName)

############      
fileName <- "gage_pred_prot_Epo_4to5_Vs_0to1"
gos <- gage(mymat[,c(1:2, 5:6)],gsets =mydb,ref=c(1:2),same.dir=T,compare = "unpaired",rank=T, saaTest=gs.KSTest)
gageAnalysis(2, fileName)

############      
fileName <- "gage_pred_prot_Epo_6to7_Vs_0to1"
gos <- gage(mymat[,c(1:2, 7:8)],gsets =mydb,ref=c(1:2),same.dir=T,compare = "unpaired",rank=T, saaTest=gs.KSTest)
gageAnalysis(3, fileName)

############      
fileName <- "gage_pred_prot_Epo_8to14_Vs_0to1"
gos <- gage(mymat[,c(1:2, 9:10)],gsets =mydb,ref=c(1:2),same.dir=T,compare = "unpaired",rank=T, saaTest=gs.KSTest)
gageAnalysis(4, fileName)

############      
fileName <- "gage_pred_prot_Epo_19to24_Vs_0to1"
gos <- gage(mymat[,c(1:2, 11:12)],gsets =mydb,ref=c(1:2),same.dir=T,compare = "unpaired",rank=T, saaTest=gs.KSTest)
gageAnalysis(5, fileName)



# AGAINST PREVIOUS TIME POINTS

contrast.matrix <- makeContrasts(e2-e1, e3-e2, e4-e3, e5-e4, e6-e5,levels=design)
# fit and get up and down genes
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)


############      
fileName <- "gage_pred_prot_Epo_2to3_Vs_0to1_v2"
gos <- gage(mymat[,c(1:4)],gsets =mydb,ref=c(1:2),same.dir=T,compare = "unpaired",rank=T, saaTest=gs.KSTest)
gageAnalysis(1, fileName)

############      
fileName <- "gage_pred_prot_Epo_4to5_Vs_2to3"
gos <- gage(mymat[,c(3:6)],gsets =mydb,ref=c(1:2),same.dir=T,compare = "unpaired",rank=T, saaTest=gs.KSTest)
gageAnalysis(2, fileName)

############      
fileName <- "gage_pred_prot_Epo_6to7_Vs_4to5"
gos <- gage(mymat[,c(5:8)],gsets =mydb,ref=c(1:2),same.dir=T,compare = "unpaired",rank=T, saaTest=gs.KSTest)
gageAnalysis(3, fileName)

############      
fileName <- "gage_pred_prot_Epo_8to14_Vs_6to7"
gos <- gage(mymat[,c(7:10)],gsets =mydb,ref=c(1:2),same.dir=T,compare = "unpaired",rank=T, saaTest=gs.KSTest)
gageAnalysis(4, fileName)

############      
fileName <- "gage_pred_prot_Epo_19to24_Vs_8to14"
gos <- gage(mymat[,c(9:12)],gsets =mydb,ref=c(1:2),same.dir=T,compare = "unpaired",rank=T, saaTest=gs.KSTest)
gageAnalysis(5, fileName)

