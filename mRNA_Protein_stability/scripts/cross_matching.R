#options(java.parameters = "-Xmx8g")

library(plyr)
library(ggplot2)
library(openxlsx)
library(GeneAnswers)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(pheatmap)
library(reshape2)
library(limma)
library(openxlsx)
library(plyr)
library(SummarizedExperiment)
library(limma)
library(org.Hs.eg.db)
library(edgeR)
library(pheatmap)
library(Rtsne)
library(survival) #use survfit and Surv function
library(survminer)
library(parallel)
library(viridis)
library(ggpubr)
library(ggplot2)
library(DESeq2)

source("tools.r")

############################################
###                                      ###
###               FUNCTION               ###
###                                      ###
############################################



entrez2symbol <- function(entrez)
{
  symbol <- mget(as.character(entrez), org.Hs.egSYMBOL, ifnotfound=NA)
  symbol <- unlist(lapply(symbol, function(i) return(i[1])))
  return(symbol)
}

ensembl2entrez <- function(ensembl)
{
  entrez <- mget(as.character(ensembl), org.Hs.egENSEMBL2EG, ifnotfound=NA)
  entrez <- lapply(entrez, function(i) return(i[1]))
  return(unlist(entrez))
}

symbol2entrez <- function(symbol)
{
  entrez <- mget(as.character(symbol), org.Hs.egSYMBOL2EG, ifnotfound=NA)
  entrez <- unlist(lapply(entrez, function(i) return(i[1])))
  entrez <- unique(entrez[!is.na(entrez)])
  return(entrez)
}

entrez2name <- function(entrez)
{
  gn <- mget(as.character(entrez), org.Hs.egGENENAME, ifnotfound=NA)
  gn <- unlist(lapply(gn, function(i) return(i[1])))
  return(gn)
}

toNum <- function(x){
  if(is.numeric(x)) return(x)
  return(as.numeric(levels(x))[x])
}


mouseSymbol2HumanSymbol <- function(symbol.mouse)
{
  entrez.mouse <- mget(as.character(symbol.mouse), org.Mm.egSYMBOL2EG, ifnotfound=NA)
  entrez.mouse <- unlist(lapply(entrez.mouse, function(i) return(i[1])))
  
  homo <- getHomoGeneIDs(entrez.mouse, species="mouse", speciesL="human",mappingMethod='direct')	
  
  entrez.human <- rep(NA, length(entrez.mouse))
  entrez.human[match(names(homo), entrez.mouse)] <- homo
  
  symbol.human <- mget(homo, org.Hs.egSYMBOL, ifnotfound=NA)
  symbol.human <- unlist(lapply(symbol.human, function(i) return(i[1])))
  
  symbol.human.complete <- rep(NA, length(entrez.mouse))
  symbol.human.complete[match(names(symbol.human), entrez.human)] <- symbol.human
  
  return(symbol.human.complete)
}



getDensity <- function(x, y)
{
  dc <- densCols(x, y, colramp=colorRampPalette(c("black", "white")))
  dens <- col2rgb(dc)[1,] + 1L
  dens <- dens/max(dens)
  return(dens)
}


densityScatterplot <- function(mat, outFile)
{
  myCol <-  colorRampPalette(c("#000099", "#00FEFF", "#45FE4F", "#FCFF00", "#FF9400", "#FF3100"))(256)
  
  p <- ggplot(mat)
  p <- p + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                              panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
  p <- p + theme(axis.line.x = element_line(color="black", size = 0.5),
                 axis.line.y = element_line(color="black", size = 0.5))	
  p <- p + geom_point(aes(X, Y, col = Density), size = 2, alpha = 0.5) + scale_colour_gradientn(limits=c(0, 1), colours = myCol)
  p <- p + facet_wrap(~ Group, ncol=5)
  
  #p <- p + scale_x_continuous(limits = c(-1,13), breaks=seq(0, 12, 4))
  #p <- p + scale_y_continuous(limits = c(-1,13), breaks=seq(0, 12, 4))
  p <- p + theme(axis.text.x = element_text(size=20),axis.text.y = element_text(size=20))
  #p <- p + xlab("Dataset 1") + ylab("Dataset 2")
  
  pdf(outFile, width = 16, height = 4)
  plot(p)
  dev.off()
}


correlationPipeline <- function(x, y, genes, groups, name)
{
  # init
  corre <- c()
  pv.greater <- c()
  pv.less <- c()
  
  # loop groups
  for(g in groups)
  {
    idx <- which(genes %in% g)
    cor.true <- cor(toNum(x[idx]), toNum(y[idx]),
                    method = "spearman", use = "complete.obs")
    cor.random <- getRandomCor(length(idx), 10000, toNum(x), toNum(y))	
    
    corre <- c(corre, cor.true)
    pv.greater <- c(pv.greater, getEmpiricalPV(cor.random, cor.true, "greater"))
    pv.less <- c(pv.less, getEmpiricalPV(cor.random, cor.true, "less"))	
  }
  qv.greater <- p.adjust(pv.greater, method = "BH")
  qv.less <- p.adjust(pv.less, method = "BH")
  return(data.frame(Sample = rep(name, length(groups)),
                    Group = names(groups),
                    Correlation = corre,
                    PV.greater = pv.greater,
                    FDR.greater = qv.greater,
                    PV.less = pv.less,
                    FDR.less = qv.less))
}

getRandomCor <- function(size, nb, x, y)
{
  randomIdx <- lapply(1:nb, function(i) return(sample(1:length(x), size)))
  corList <- lapply(randomIdx, function(i) return(cor(x[i], y[i], method = "spearman", use = "complete.obs")))
  return(unlist(corList))
}

getEmpiricalPV <- function(x, y, type)
{
  x.ecdf <- stats::ecdf(x)
  if(type == "greater") return(1-x.ecdf(y))
  else if(type == "less") return(x.ecdf(y))
  return(NA)    
}


getDensity <- function(x, y)
{
  dc <- densCols(x, y, colramp=colorRampPalette(c("black", "white")))
  dens <- col2rgb(dc)[1,] + 1L
  dens <- dens/max(dens)
  return(dens)
}


densityScatterplot <- function(mat, outFile)
{
  myCol <-  colorRampPalette(c("#000099", "#00FEFF", "#45FE4F", "#FCFF00", "#FF9400", "#FF3100"))(256)
  
  p <- ggplot(mat)
  p <- p + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                              panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
  p <- p + theme(axis.line.x = element_line(color="black", size = 0.5),
                 axis.line.y = element_line(color="black", size = 0.5))	
  p <- p + geom_point(aes(value, valueY, col = Density), size = 2, alpha = 0.5) + scale_colour_gradientn(limits=c(0, 1), colours = myCol)
  p <- p + facet_wrap(~ key, ncol=7)
  
  p <- p + scale_x_continuous(breaks=seq(-4, 4, 2))
  p <- p + scale_y_continuous(breaks=seq(-4, 4, 2))
  p <- p + theme(axis.text.x = element_text(size=20),axis.text.y = element_text(size=20))
  #p <- p + xlab("Dataset 1") + ylab("Dataset 2")
  
  pdf(outFile, width = 16, height = 4)
  plot(p)
  dev.off()
}

############################################
###                                      ###
###                 MAIN                 ###
###                                      ###
############################################


# LOAD STABILITY GROUPS
setwd(file.path("../data/"))
ann <- read.xlsx("stability_group_human.xlsx", sheet = 1)

####################################
# LOAD HUMAN TISSUE DATA (DS1)
rna1 <- read.xlsx("rna_DS1.xlsx", sheet = 1, rowNames = TRUE)
prot1 <- read.xlsx("protein_DS1.xlsx", sheet = 1, rowNames = TRUE)

colnames(rna1) <- toupper(colnames(rna1))
colnames(prot1) <- toupper(colnames(prot1))

#ann.sub <- ann[match(rownames(rna), ann$Symbol), ]


####################################
# LOAD HUMAN TISSUE DATA (DS2)
rna2 <- read.xlsx("rna_DS2.xlsx", sheet = 1, rowNames = TRUE)
prot2 <- read.xlsx("protein_DS2.xlsx", sheet = 1, rowNames = TRUE)

colnames(rna2) <- toupper(colnames(rna2))
colnames(prot2) <- toupper(colnames(prot2))

####################################
# LOAD HUMAN TISSUE DATA (DS3)
rna3 <- read.xlsx("rna_DS3.xlsx", sheet = 1, rowNames = TRUE)
prot3 <- read.xlsx("protein_DS3.xlsx", sheet = 1, rowNames = TRUE)

colnames(rna3) <- toupper(colnames(rna3))
colnames(prot3) <- toupper(colnames(prot3))



# INTERSECT THE 3 DATASETS
commonGenes <- Intersect(list(rownames(rna1), rownames(rna2), rownames(rna3)))
commonSamples <- Intersect(list(colnames(rna1), colnames(rna2), colnames(rna3)))


rna1 <- rna1[commonGenes, commonSamples]
prot1 <- prot1[commonGenes, commonSamples]

rna2 <- rna2[commonGenes, commonSamples]
prot2 <- prot2[commonGenes, commonSamples]

rna3 <- rna3[commonGenes, commonSamples]
prot3 <- prot3[commonGenes, commonSamples]


ann.sub <- ann[match(commonGenes, ann$Symbol), ]


# CROSS MATCHING ANALYSIS
outDir <- file.path("../output")
dir.create(outDir)
setwd(outDir)

tissueNames <- colnames(rna1)

for(i in tissueNames){
  rna1.current <- as.numeric(rna1[,i])
  prot1.current <- as.numeric(prot1[,i])
  
  rna2.current <- as.numeric(rna2[,i])
  prot2.current <- as.numeric(prot2[,i])
  
  rna3.current <- as.numeric(rna3[,i])
  prot3.current <- as.numeric(prot3[,i])
  
  genes.current <- ann.sub$Symbol
  
  # remove NA
  rna1.NA <- is.na(rna1.current)
  prot1.NA <- is.na(prot1.current)
  rna2.NA <- is.na(rna2.current)
  prot2.NA <- is.na(prot2.current)
  rna3.NA <- is.na(rna3.current)
  prot3.NA <- is.na(prot3.current)
  
  tokeep <- !(rna1.NA | prot1.NA | rna2.NA | prot2.NA | rna3.NA | prot3.NA)
  
  rna1.current <- rna1.current[tokeep]
  prot1.current <- prot1.current[tokeep]
  rna2.current <- rna2.current[tokeep]
  prot2.current <- prot2.current[tokeep]
  rna3.current <- rna3.current[tokeep]
  prot3.current <- prot3.current[tokeep]
  
  genes.current <- genes.current[tokeep]
  
  ann.current <- ann[match(genes.current, ann$Symbol), ]
  
  groups <- unique(ann.current$Group)
  
  # Bootstrap PV
  groupList <- lapply(groups, function(g) ann.current$Symbol[ann.current$Group == g])
  names(groupList) <- groups
  
  corPV.rna1.prot1 <- correlationPipeline(x = rna1.current, y = prot1.current,
                                          genes = genes.current,
                                          groups = groupList, name = i)
  corPV.rna1.prot1$Cross.matching <- "rna1.prot1"
  
  
  corPV.rna1.prot2 <- correlationPipeline(x = rna1.current, y = prot2.current,
                                          genes = genes.current,
                                          groups = groupList, name = i)
  corPV.rna1.prot2$Cross.matching <- "rna1.prot2"
  
  
  corPV.rna1.prot3 <- correlationPipeline(x = rna1.current, y = prot3.current,
                                          genes = genes.current,
                                          groups = groupList, name = i)
  corPV.rna1.prot3$Cross.matching <- "rna1.prot3"
  
  
  corPV.rna2.prot1 <- correlationPipeline(x = rna2.current, y = prot1.current,
                                          genes = genes.current,
                                          groups = groupList, name = i)
  corPV.rna2.prot1$Cross.matching <- "rna2.prot1"
  
  
  corPV.rna2.prot2 <- correlationPipeline(x = rna2.current, y = prot2.current,
                                          genes = genes.current,
                                          groups = groupList, name = i)
  corPV.rna2.prot2$Cross.matching <- "rna2.prot2"
  
  corPV.rna2.prot3 <- correlationPipeline(x = rna2.current, y = prot3.current,
                                          genes = genes.current,
                                          groups = groupList, name = i)
  corPV.rna2.prot3$Cross.matching <- "rna2.prot3"
  
  corPV.rna3.prot1 <- correlationPipeline(x = rna3.current, y = prot1.current,
                                          genes = genes.current,
                                          groups = groupList, name = i)
  corPV.rna3.prot1$Cross.matching <- "rna3.prot1"
  
  
  corPV.rna3.prot2 <- correlationPipeline(x = rna3.current, y = prot2.current,
                                          genes = genes.current,
                                          groups = groupList, name = i)
  corPV.rna3.prot2$Cross.matching <- "rna3.prot2"
  
  corPV.rna3.prot3 <- correlationPipeline(x = rna3.current, y = prot3.current,
                                          genes = genes.current,
                                          groups = groupList, name = i)
  corPV.rna3.prot3$Cross.matching <- "rna3.prot3"
  
  
  corPV <- do.call(rbind, list(corPV.rna1.prot1, corPV.rna1.prot2, corPV.rna1.prot3,
                               corPV.rna2.prot1, corPV.rna2.prot2, corPV.rna2.prot3,
                               corPV.rna3.prot1, corPV.rna3.prot2, corPV.rna3.prot3))
  
  
  write.xlsx(corPV, paste0(i, "_crossmatching_correlation_bootstrap_PV.xlsx"), row.names = FALSE)  
  
}

##################################################
# PLOT
setwd(outDir)
corFiles <- list.files(pattern = "_crossmatching_correlation_bootstrap_PV.xlsx")
corList <- lapply(corFiles, read.xlsx, sheet = 1)
corMat <- do.call(rbind, corList)

corMat$Sample <- factor(corMat$Sample, levels = rev(sort(unique(corMat$Sample))))

# SET LIMIT TO 10^-4
corMat$PV.greater[corMat $PV.greater < 0.0001] <- 0.0001
corMat$PV.less[corMat$PV.less < 0.0001] <- 0.0001

score <- c()
for(i in 1:nrow(corMat))
{
  if(corMat$PV.greater[i] < corMat$PV.less[i]) score <- c(score, -log10(corMat$PV.greater[i]))
  else(score <- c(score, log10(corMat$PV.less[i])))
}

corMat$Score <- score	

# without grey zone
mycol <- c("rS.pS" = rgb(237, 189, 79, maxColorValue = 255),
           "rU.pS" = rgb(230, 121, 123, maxColorValue = 255),
           "rS.pU" = rgb(105, 153, 196, maxColorValue = 255),
           "rU.pU" = rgb(105, 48, 48, maxColorValue = 255))

p <- ggplot(data=corMat[corMat$Group != "grey.zone", ], aes(x=Sample, y=Score))
p <- p + coord_flip()
p <- p + theme_bw(base_size = 16) + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                          panel.grid.minor = element_blank())
p <- p + theme(axis.line.x = element_line(color="black", size = 0.5),
               axis.line.y = element_line(color="black", size = 0.5))	
p <- p + geom_hline(yintercept = 0, colour="black", linetype = "longdash", alpha = 0.75)
#p <- p + geom_boxplot()
#p <- p + geom_point(aes(y=Score, group=Group, colour = Crossmatch, shape = Group),
#	position = position_dodge(width=0.75), size = 5, alpha = 0.5)
p <- p + geom_point(aes(y=Score, group=Cross.matching, colour = Group, shape = Cross.matching),
                    position = position_dodge(width=0.75), size = 5, alpha = 0.5)	
p <- p + scale_shape_manual(values = c(15, 3, 16, 4, 17, 5, 18, 6, 7))
p <- p + xlab("")
p <- p + scale_color_manual(values = mycol)

#pdf("crossmatching_woGZ.pdf", width = 8, height = 7)
pdf("figure2D.pdf", width = 8, height = 7)
plot(p)
dev.off()


toxlsx <- corMat[corMat$Group != "grey.zone", ]
write.xlsx(toxlsx, "figure2D_data.xlsx")



#################################################
# Correlate RNA1 vs. RNA2 vs. RNA3

rna3[rna3 == 0] <- NA

# RNA1 vs. RNA2
ggmat.d1 <- gather(data.frame(scale(rna1)))
ggmat.d2 <- gather(data.frame(scale(rna2)))
ggmat.d3 <- gather(data.frame(scale(rna3)))

ggmat <- cbind(ggmat.d1, valueY = ggmat.d2$value)

densList <- lapply(unique(ggmat$key), function(i)
  getDensity(ggmat$value[ggmat$key == i], ggmat$valueY[ggmat$key == i]))
ggmat$Density <- unlist(densList)	

densityScatterplot(ggmat, "rna1_vs_rna2_scatterplot.pdf")
write.xlsx(ggmat, "SupFigure4A_data.xlsx")

corMat <- lapply(unique(ggmat$key), function(i){
  ct <- cor.test(ggmat$value[ggmat$key == i], ggmat$valueY[ggmat$key == i])
  return(c(ct$estimate["cor"], ct$p.value))
})
corMat <- do.call(rbind, corMat)
rownames(corMat) <- unique(ggmat$key)
colnames(corMat) <- c("R", "PV")
write.xlsx(corMat, "rna1_vs_rna2_cor_pv.xlsx", row.names = TRUE)


# RNA1 vs. RNA3
ggmat <- cbind(ggmat.d1, valueY = ggmat.d3$value)

densList <- lapply(unique(ggmat$key), function(i)
  getDensity(ggmat$value[ggmat$key == i], ggmat$valueY[ggmat$key == i]))
ggmat$Density <- unlist(densList)	

densityScatterplot(ggmat, "rna1_vs_rna3_scatterplot.pdf")
write.xlsx(ggmat, "SupFigure4B_data.xlsx")

corMat <- lapply(unique(ggmat$key), function(i){
  ct <- cor.test(ggmat$value[ggmat$key == i], ggmat$valueY[ggmat$key == i])
  return(c(ct$estimate["cor"], ct$p.value))
})
corMat <- do.call(rbind, corMat)
rownames(corMat) <- unique(ggmat$key)
colnames(corMat) <- c("R", "PV")
write.xlsx(corMat, "rna1_vs_rna3_cor_pv.xlsx", row.names = TRUE)


# RNA2 vs. RNA3
ggmat <- cbind(ggmat.d2, valueY = ggmat.d3$value)

densList <- lapply(unique(ggmat$key), function(i)
  getDensity(ggmat$value[ggmat$key == i], ggmat$valueY[ggmat$key == i]))
ggmat$Density <- unlist(densList)	

densityScatterplot(ggmat, "rna2_vs_rna3_scatterplot.pdf")
write.xlsx(ggmat, "SupFigure4C_data.xlsx")


corMat <- lapply(unique(ggmat$key), function(i){
  ct <- cor.test(ggmat$value[ggmat$key == i], ggmat$valueY[ggmat$key == i])
  return(c(ct$estimate["cor"], ct$p.value))
})
corMat <- do.call(rbind, corMat)
rownames(corMat) <- unique(ggmat$key)
colnames(corMat) <- c("R", "PV")
write.xlsx(corMat, "rna2_vs_rna3_cor_pv.xlsx", row.names = TRUE)


# PROT1 vs. PROT2
prot3[prot3 == 0] <- NA

ggmat.d1 <- gather(data.frame(scale(prot1)))
ggmat.d2 <- gather(data.frame(scale(prot2)))
ggmat.d3 <- gather(data.frame(scale(prot3)))

ggmat <- cbind(ggmat.d1, valueY = ggmat.d2$value)

densList <- lapply(unique(ggmat$key), function(i)
  getDensity(ggmat$value[ggmat$key == i], ggmat$valueY[ggmat$key == i]))
ggmat$Density <- unlist(densList)	

densityScatterplot(ggmat, "prot1_vs_prot2_scatterplot.pdf")
write.xlsx(ggmat, "SupFigure4D_data.xlsx")


corMat <- lapply(unique(ggmat$key), function(i){
  ct <- cor.test(ggmat$value[ggmat$key == i], ggmat$valueY[ggmat$key == i])
  return(c(ct$estimate["cor"], ct$p.value))
})
corMat <- do.call(rbind, corMat)
rownames(corMat) <- unique(ggmat$key)
colnames(corMat) <- c("R", "PV")
write.xlsx(corMat, "prot1_vs_prot2_cor_pv.xlsx", row.names = TRUE)


# PROT1 vs. PROT3
ggmat <- cbind(ggmat.d1, valueY = ggmat.d3$value)

densList <- lapply(unique(ggmat$key), function(i)
  getDensity(ggmat$value[ggmat$key == i], ggmat$valueY[ggmat$key == i]))
ggmat$Density <- unlist(densList)	

densityScatterplot(ggmat, "prot1_vs_prot3_scatterplot.pdf")
write.xlsx(ggmat, "SupFigure4E_data.xlsx")


corMat <- lapply(unique(ggmat$key), function(i){
  ct <- cor.test(ggmat$value[ggmat$key == i], ggmat$valueY[ggmat$key == i])
  return(c(ct$estimate["cor"], ct$p.value))
})
corMat <- do.call(rbind, corMat)
rownames(corMat) <- unique(ggmat$key)
colnames(corMat) <- c("R", "PV")
write.xlsx(corMat, "prot1_vs_prot3_cor_pv.xlsx", row.names = TRUE)


# PROT2 vs. PROT3
ggmat <- cbind(ggmat.d2, valueY = ggmat.d3$value)

densList <- lapply(unique(ggmat$key), function(i)
  getDensity(ggmat$value[ggmat$key == i], ggmat$valueY[ggmat$key == i]))
ggmat$Density <- unlist(densList)	

densityScatterplot(ggmat, "prot2_vs_prot3_scatterplot.pdf")
write.xlsx(ggmat, "SupFigure4F_data.xlsx")


corMat <- lapply(unique(ggmat$key), function(i){
  ct <- cor.test(ggmat$value[ggmat$key == i], ggmat$valueY[ggmat$key == i])
  return(c(ct$estimate["cor"], ct$p.value))
})
corMat <- do.call(rbind, corMat)
rownames(corMat) <- unique(ggmat$key)
colnames(corMat) <- c("R", "PV")
write.xlsx(corMat, "prot2_vs_prot3_cor_pv.xlsx", row.names = TRUE)



