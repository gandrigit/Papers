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
library(tidyr)
library(ggrepel)

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


deltaPipeline <- function(x1, y1, x2, y2, genes, groups, name)
{
  # init
  cor1 <- c()
  cor2 <- c()
  delta <- c()
  delta.rd <- c()
  pv.greater <- c()
  pv.less <- c()
  x1.mean <- c()
  y1.mean <- c()
  x2.mean <- c()
  y2.mean <- c()
  
  # loop groups
  for(g in groups)
  {
    idx <- which(genes %in% g)
    cor1.true <- cor(toNum(x1[idx]), toNum(y1[idx]),
                     method = "spearman", use = "complete.obs")
    cor2.true <- cor(toNum(x2[idx]), toNum(y2[idx]),
                     method = "spearman", use = "complete.obs")
    delta.true <- cor1.true - cor2.true
    
    delta.random <- getRandomDelta(length(idx), 1000, toNum(x1), toNum(y1), toNum(x2), toNum(y2))	
    
    cor1 <- c(cor1, cor1.true)
    cor2 <- c(cor2, cor2.true)
    delta <- c(delta, delta.true)
    delta.rd <- c(delta.rd, median(delta.random))
    
    x1.mean <- c(x1.mean, mean(toNum(x1[idx])))
    y1.mean <- c(y1.mean, mean(toNum(y1[idx])))
    x2.mean <- c(x2.mean, mean(toNum(x2[idx])))
    y2.mean <- c(y2.mean, mean(toNum(y2[idx])))
    
    pv.greater <- c(pv.greater, getEmpiricalPV(delta.random, delta.true, "greater"))
    pv.less <- c(pv.less, getEmpiricalPV(delta.random, delta.true, "less"))	
  }
  return(data.frame(Sample = rep(name, length(groups)),
                    Group = names(groups),
                    X1.mean = x1.mean,
                    Y1.mean = y1.mean,
                    X2.mean = x2.mean,
                    Y2.mean = y2.mean,
                    Correlation1 = cor1,
                    Correlation2 = cor2,
                    Delta = delta, 
                    Delta.random = delta.rd, 
                    PV.greater = pv.greater,
                    PV.less = pv.less))
}


getRandomCor <- function(size, nb, x, y)
{
  randomIdx <- lapply(1:nb, function(i) return(sample(1:length(x), size)))
  corList <- lapply(randomIdx, function(i) return(cor(x[i], y[i], method = "spearman", use = "complete.obs")))
  return(unlist(corList))
}

getRandomDelta <- function(size, nb, x1, y1, x2, y2)
{
  randomIdx <- lapply(1:nb, function(i) return(sample(1:length(x1), size)))
  corList1 <- sapply(randomIdx, function(i) return(cor(x1[i], y1[i], method = "spearman", use = "complete.obs")))
  corList2 <- sapply(randomIdx, function(i) return(cor(x2[i], y2[i], method = "spearman", use = "complete.obs")))
  return(corList1 - corList2)
}

getEmpiricalPV <- function(x, y, type)
{
  x.ecdf <- stats::ecdf(x)
  if(type == "greater") return(1-x.ecdf(y))
  else if(type == "less") return(x.ecdf(y))
  return(NA)    
}

getNormPV <- function(x, v, lower = TRUE)
{
  return(pnorm(x, mean(v), sd(v), lower.tail = lower))
}



############################################
###                                      ###
###                 MAIN                 ###
###                                      ###
############################################


# LOAD STABILITY GROUPS
setwd(file.path("../data"))
ann <- read.xlsx("stability_group_human.xlsx", sheet = 1)


####################################
# LOAD TCGA COAD DATA
setwd(file.path("../data"))
rna <- read.xlsx("rna_COAD.xlsx", sheet = 1, rowNames = TRUE)
prot <- read.xlsx("prot_COAD.xlsx", sheet = 1, rowNames = TRUE)
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

# common genes
commonGenes <- intersect(rownames(rna), rownames(prot))
rna.tcga <- rna[commonGenes, ]
prot.tcga <- prot[commonGenes, ]

colnames(rna.tcga) <- substring(colnames(rna.tcga), 1, 12)
colnames(prot.tcga) <- substring(colnames(prot.tcga), 1, 12)

ann.sub <- ann[match(commonGenes, ann$Symbol), ]

####################################
# LOAD TCGA SAMPLE ANNOTATION

clinMat <- read.xlsx("TCGA_COAD_clinical.xlsx", sheet = 1)

# Intersect RNA samples
clinMat <- clinMat[match(colnames(rna.tcga), clinMat$bcr_patient_barcode), ]

# DIVIDE EARLY / LATE STAGES
stage <- rep(NA, ncol(rna.tcga))
stage[clinMat$tumor_stage == "stage i" |
        clinMat$tumor_stage == "stage ii" |
        clinMat$tumor_stage == "stage iia" |
        clinMat$tumor_stage == "stage iib"] <- "early"
stage[clinMat$tumor_stage == "stage iii" |
        clinMat$tumor_stage == "stage iiia" |
        clinMat$tumor_stage == "stage iiib" |
        clinMat$tumor_stage == "stage iiic" |
        clinMat$tumor_stage == "stage iv"] <- "late"

sample.early <- colnames(rna.tcga)[stage == "early"]
sample.late <- colnames(rna.tcga)[stage == "late"]


####################################
# GENE-WISE CORRELATION

outDir <- file.path("../output")
dir.create(outDir)
setwd(outDir)

mygenes <- rownames(rna.tcga)
cor.early <- sapply(mygenes, function(i){
  return(cor(as.numeric(rna.tcga[i, sample.early]), as.numeric(prot.tcga[i, sample.early]),
             method = "spearman", use = "complete.obs"))
})
cor.late <- sapply(mygenes, function(i){
  return(cor(as.numeric(rna.tcga[i, sample.late]), as.numeric(prot.tcga[i, sample.late]),
             method = "spearman", use = "complete.obs"))
})


geneWiseMat <- data.frame(symbol = mygenes,
                          cor.early = cor.early,
                          cor.late = cor.late,
                          delta = cor.late - cor.early,
                          Group = ann.sub$Group)

# TEST NORMALITY
shapi <- data.frame(shapiro.pvalue = shapiro.test(geneWiseMat$delta)$p.value)
write.xlsx(shapi, "figure7ABCD_Shapiro_test.xlsx", row.names = FALSE)

# PVALUE (PNORM)
pvLess <- sapply(geneWiseMat$delta, getNormPV, v = geneWiseMat$delta, lower = TRUE)
pvGreater <- sapply(geneWiseMat$delta, getNormPV, v = geneWiseMat$delta, lower = FALSE)

geneWiseMat$pv.less <- pvLess
geneWiseMat$FDR.less <- p.adjust(geneWiseMat$pv.less, method = "BH")
geneWiseMat$pv.greater <- pvGreater
geneWiseMat$FDR.greater <- p.adjust(geneWiseMat$pv.greater, method = "BH")
geneWiseMat$score <- -log10(sapply(1:nrow(geneWiseMat),
                                   function(i) min(geneWiseMat[i, c("pv.less", "pv.greater")])))


#SAVE
setwd(file.path("~/Research/Sajib/P2_half_lives/TCGA/early_vs_late"))
write.xlsx(geneWiseMat, "genewise_late_vs_early.xlsx")

# SPLIT
dfList <- lapply(unique(geneWiseMat$Group), function(i){
  return(geneWiseMat[geneWiseMat$Group == i, ])
})
names(dfList) <- unique(geneWiseMat$Group)
write.xlsx(dfList, "genewise_late_vs_early_SPLIT.xlsx")

# SAVE EXPRESSION
setwd(file.path("~/Research/Sajib/P2_half_lives/data"))
rna.colon <- read.xlsx("rna_colon.xlsx", sheet = 1, rowNames = TRUE)
rna.colon <- rna.colon[match(mygenes, row.names(rna.colon)), ]
prot.colon <- read.xlsx("prot_colon.xlsx", sheet = 1, rowNames = TRUE)
prot.colon <- prot.colon[match(mygenes, row.names(prot.colon)), ]

rna.colon <- cbind(geneWiseMat, rna.colon)
prot.colon <- cbind(geneWiseMat, prot.colon)
rna.early <- cbind(geneWiseMat, rna.tcga[mygenes, sample.early])
prot.early <- cbind(geneWiseMat, prot.tcga[mygenes, sample.early])
rna.late <- cbind(geneWiseMat, rna.tcga[mygenes, sample.late])
prot.late <- cbind(geneWiseMat, prot.tcga[mygenes, sample.late])

toxlsx <- list(rna.colon = rna.colon,
               prot.colon = prot.colon,
               rna.early = rna.early,
               prot.early = prot.early,
               rna.late = rna.late,
               prot.late = prot.late)
toxlsx <- lapply(toxlsx, function(i) i[order(i$Group, i$symbol), ])

setwd(file.path("~/Research/Sajib/P2_half_lives/TCGA/early_vs_late"))
write.xlsx(toxlsx, "colon_COADearly_COADlate_expression.xlsx")

# PLOT
geneWiseMat <- geneWiseMat[geneWiseMat$Group != "grey.zone", ]# remove grey zone
p <- ggplot(geneWiseMat, aes(x=cor.early, y=cor.late))
p <- p + facet_wrap( ~ Group, ncol = 2)
p <- p + geom_hline(yintercept = 0, linetype= 2) + geom_vline(xintercept = 0, linetype= 2)
p <- p + geom_point(aes(colour = delta), alpha = 0.75, shape=16, size = geneWiseMat$score*2)
p <- p + scale_colour_gradient2(low = "navy", mid = "snow", high = "firebrick")
p <- p + theme_bw() + theme(legend.position="none") 
p <- p + geom_text_repel(data=geneWiseMat[geneWiseMat$score > -log10(0.05),], aes(label=symbol), cex = 2, segment.size = 0.1)

pdf("genewise_late_vs_early_scatterplot.pdf")
plot(p)
dev.off()

write.xlsx(geneWiseMat, "figure7ABCD_data.xlsx")


# DELTA NORMAL / TUMOR VS. EARLY / LATE
normal_tumor_list <- lapply(unique(as.character(geneWiseMat$Group)), function(i)
  read.xlsx("~/Research/Sajib/P2_half_lives/TCGA/gene_RNA_prot_ttest.xlsx", sheet = i))
names(normal_tumor_list) <- unique(as.character(geneWiseMat$Group))

deltaComp <- lapply(1:length(normal_tumor_list), function(i){
  genes <- normal_tumor_list[[i]]$GENE
  delta <- normal_tumor_list[[i]]$delta.median
  
  delta.lateVSearly <- geneWiseMat$delta[match(genes, geneWiseMat$symbol)]
  #cor(delta, delta.lateVSearly, method = "spearman")
  
  df <- cbind(normal_tumor_list[[i]], geneWiseMat[match(genes, geneWiseMat$symbol),])
  return(df)
})
names(deltaComp) <- names(normal_tumor_list) 

# SELECT TOP 25
deltaComp.sub <- lapply(deltaComp, function(i){
  nb <- round(nrow(i) * 0.25)
  return(i[c(order(i$delta.median)[1:nb], order(-i$delta.median)[1:nb]), ])
})

write.xlsx(deltaComp, "normal_tumor_vs_early_late.xlsx")
write.xlsx(deltaComp.sub, "normal_tumor_vs_early_late_TOP25.xlsx")

#############################
# DELTA VS. MEDIAN EXPRESSION
plotDir <- file.path(outDir, "plot")
dir.create(plotDir)
setwd(plotDir)

mygenes <- as.character(geneWiseMat$symbol)

prot.tcga[is.na(prot.tcga)] <- 0

rna.early.scaled <- scale(rna.tcga[mygenes, sample.early])
rna.late.scaled <- scale(rna.tcga[mygenes, sample.late])
prot.early.scaled <- scale(prot.tcga[mygenes, sample.early])
prot.late.scaled <- scale(prot.tcga[mygenes, sample.late])


# TWO-SIDED
rna.pv <- sapply(mygenes, function(j){
  t.test(as.numeric(rna.early.scaled[j, ]), as.numeric(rna.late.scaled[j, ]), alternative = "two.sided")$p.value
})
names(rna.pv) <- mygenes
rna.pv <- signif(rna.pv, digits = 3)
rna.qv <- p.adjust(rna.pv, method = "BH")

prot.pv <- sapply(mygenes, function(j){
  t.test(as.numeric(prot.early.scaled[j, ]), as.numeric(prot.late.scaled[j, ]), alternative = "two.sided")$p.value
})
names(prot.pv) <- mygenes
prot.pv <- signif(prot.pv, digits = 3)
prot.qv <- p.adjust(prot.pv, method = "BH")

# LESS
rna.pv.less <- sapply(mygenes, function(j){
  t.test(as.numeric(rna.early.scaled[j, ]), as.numeric(rna.late.scaled[j, ]), alternative = "less")$p.value
})
names(rna.pv.less) <- mygenes
rna.pv.less <- signif(rna.pv.less, digits = 3)
rna.qv.less <- p.adjust(rna.pv.less, method = "BH")

prot.pv.less <- sapply(mygenes, function(j){
  t.test(as.numeric(prot.early.scaled[j, ]), as.numeric(prot.late.scaled[j, ]), alternative = "less")$p.value
})
names(prot.pv.less) <- mygenes
prot.pv.less <- signif(prot.pv.less, digits = 3)
prot.qv.less <- p.adjust(prot.pv.less, method = "BH")

# GREATER
rna.pv.greater <- sapply(mygenes, function(j){
  t.test(as.numeric(rna.early.scaled[j, ]), as.numeric(rna.late.scaled[j, ]), alternative = "greater")$p.value
})
names(rna.pv.greater) <- mygenes
rna.pv.greater <- signif(rna.pv.greater, digits = 3)
rna.qv.greater <- p.adjust(rna.pv.greater, method = "BH")

prot.pv.greater <- sapply(mygenes, function(j){
  t.test(as.numeric(prot.early.scaled[j, ]), as.numeric(prot.late.scaled[j, ]), alternative = "greater")$p.value
})
names(prot.pv.greater) <- mygenes
prot.pv.greater <- signif(prot.pv.greater, digits = 3)
prot.qv.greater <- p.adjust(prot.pv.greater, method = "BH")

# plot
lapply(mygenes, function(j){
  
  ggmat <- data.frame(X = c(rep("RNA", length(as.numeric(rna.early.scaled[j, ]))),
                            rep("RNA", length(as.numeric(rna.late.scaled[j, ]))),
                            rep("Protein", length(as.numeric(prot.early.scaled[j, ]))),
                            rep("Protein", length(as.numeric(prot.late.scaled[j, ])))
  ),
  Y = c(as.numeric(rna.early.scaled[j, ]), as.numeric(rna.late.scaled[j, ]),
        as.numeric(prot.early.scaled[j, ]), as.numeric(prot.late.scaled[j, ])),
  Condition = c(rep("early", length(as.numeric(rna.early.scaled[j, ]))),
                rep("late", length(as.numeric(rna.late.scaled[j, ]))),
                rep("early", length(as.numeric(prot.early.scaled[j, ]))),
                rep("late", length(as.numeric(prot.late.scaled[j, ]))))
  )
  ggmat$X <- factor(ggmat$X, levels = c("RNA", "Protein"))
  
  # Jitter plot
  p <- ggplot(ggmat, aes(x = X, y = Y, fill = Condition)) +
    geom_point(pch = 21, position = position_jitterdodge(), alpha = 0.75) + 
    geom_boxplot(outlier.size = 0, alpha = 0.5)
  p <- p + theme_bw(base_size = 16)
  p <- p + scale_fill_manual(values=c(early = rgb(159, 138, 179, maxColorValue = 255),
                                      late = rgb(206, 171, 123, maxColorValue = 255)))
  p <- p + xlab("") + ylab("Scaled Intensity")
  p <- p + ggtitle(paste0(j, "; RNA pv: ", rna.pv[j], "; Protein pv: ", prot.pv[j]))
  
  ggsave(plot = p, filename = paste0(j, "_exp_boxplot.pdf"), width = 7, height = 6)
  write.xlsx(ggmat, paste0(j, "_exp_boxplot.xlsx"))
})

df <- data.frame(GENE = mygenes,
                 rna.pv = rna.pv, rna.qv = rna.qv,
                 prot.pv = prot.pv, prot.qv = prot.qv,
                 rna.pv.less = rna.pv.less, rna.qv.less = rna.qv.less,
                 prot.pv.less = prot.pv.less, prot.qv.less = prot.qv.less,
                 rna.pv.greater = rna.pv.greater, rna.qv.greater = rna.qv.greater,
                 prot.pv.greater = prot.pv.greater, prot.qv.greater = prot.qv.greater)

setwd(outDir)
write.xlsx(df, "gene_stages_RNA_prot_ttest.xlsx")

# SPLIT PER STABILITY GROUP
dfList <- lapply(unique(geneWiseMat$Group), function(i){
  return(df[geneWiseMat$Group == i, ])
})
names(dfList) <- unique(geneWiseMat$Group)
write.xlsx(dfList, "gene_stages_RNA_prot_ttest_SPLIT.xlsx")

# FISHER (DELTA)

stabGroup <- unique(geneWiseMat$Group)

setwd(file.path("~/Research/Sajib/P2_half_lives/TCGA/early_vs_late/Fisher"))
lapply(1:length(dbList), function(j){
  mydb <- dbList[[j]]
  mydb.name <- names(dbList)[j]
  dir.create(mydb.name)
  
  lapply(stabGroup, function(i){
    m.current <- geneWiseMat[geneWiseMat$Group == i, ]
    univ.current <- symbol2entrez(m.current$symbol)
    
    nb10 <- round(length(univ.current) * 0.1)
    nb25 <- round(length(univ.current) * 0.25)
    symbolList <- list(up10pc = m.current$symbol[order(-m.current$delta)][1:nb10],
                       up25pc = m.current$symbol[order(-m.current$delta)][1:nb25],
                       down10pc = m.current$symbol[order(m.current$delta)][1:nb10],
                       down25pc = m.current$symbol[order(m.current$delta)][1:nb25])
    entrezList <- lapply(symbolList, symbol2entrez)
    
    fhList <- lapply(entrezList, hyperG, geneSets = mydb, universe = univ.current,
                     org.library = "org.Hs.eg.db", cutoff = 1, mincount = 2)
    
    write.xlsx(fhList, file.path(mydb.name, paste0(i, "_", mydb.name, "_Fisher.xlsx")))
  })
  
})






####################################
# SAMPLE-WISE CORRELATION

mygroup <- unique(ann.sub$Group)
corList <- lapply(mygroup, function(i){
  mygenes <- ann.sub$Symbol[ann.sub$Group == i]
  
  cor.current <- sapply(1:ncol(rna.tcga), function(j){
    cor(rna.tcga[mygenes, j], prot.tcga[mygenes, j], method = "spearman", use = "complete.obs")
  })
  
  return(cor.current)
})
corMat <- do.call(rbind, corList)
rownames(corMat) <- mygroup
colnames(corMat) <- colnames(rna.tcga)

corMat <- corMat[rownames(corMat) != "grey.zone", ]

corMat.early <- corMat[, sample.early]
corMat.late <- corMat[, sample.late]


# PLOT
ggmat.early <- tidyr::gather(data.frame(t(corMat.early)))
ggmat.early$stage <- "early"

ggmat.late <- tidyr::gather(data.frame(t(corMat.late)))
ggmat.late$stage <- "late"

ggmat <- rbind(ggmat.early, ggmat.late)
ggmat$key <- factor(ggmat$key, levels = c("rS.pS", "rU.pS", "rS.pU", "rU.pU"))

# V1
p <- ggplot(ggmat, aes(x=key, y=value, fill=stage)) +
  geom_boxplot()
p <- p + scale_fill_manual(values = c(early = rgb(159, 138, 179, maxColorValue = 255),
                                      late = rgb(206, 171, 123, maxColorValue = 255)))
p <- p + theme_bw(base_size = 16)
p <- p + xlab("") + ylab("Spearman correlation")

#ggsave(plot = p, filename = "samplewise_late_vs_early_boxplot.pdf", width = 6, height = 5)
ggsave(plot = p, filename = "SupFigure6A.pdf", width = 6, height = 5)
write.xlsx(ggmat, "SupFigure6A_data.xlsx")

# TEST NORMALITY
shapi <- lapply(unique(ggmat$stage), function(i){
  sapply(unique(ggmat$key), function(j){
    shapiro.test(ggmat$value[ggmat$stage == i & ggmat$key == j])$p.value
  })
})
shapi <- do.call(cbind, shapi)
rownames(shapi) <- unique(ggmat$key)
colnames(shapi) <- unique(ggmat$stage)

shapi.qv <- apply(shapi, 2, p.adjust, method = "BH")
colnames(shapi.qv) <- paste0(colnames(shapi.qv), ".FDR")
shapi <- cbind(shapi, shapi.qv)
write.xlsx(shapi, "SupFigure6A_Shapiro_test.xlsx", row.names = TRUE)

# Wilcoxon test
wilcox.test <- sapply(unique(ggmat$key), function(i){
  x = ggmat$value[ggmat$key == i & ggmat$stage == "early"]
  y = ggmat$value[ggmat$key == i & ggmat$stage == "late"]
  
  wilcox.test(x, y)$p.value
})
wilcox.test <- data.frame(wilcox.pvalue = wilcox.test)
rownames(wilcox.test) <- unique(ggmat$key)
wilcox.test <- wilcox.test[order(wilcox.test$wilcox.pvalue), , drop = FALSE]
wilcox.test$FDR <- p.adjust(wilcox.test$"wilcox.pvalue", method = "BH")

write.xlsx(wilcox.test, "SupFigure6A_wilcox_test.xlsx", row.names = TRUE)

#t-test
ttest.res <- sapply(unique(ggmat$key), function(i){
  x = ggmat$value[ggmat$key == i & ggmat$stage == "early"]
  y = ggmat$value[ggmat$key == i & ggmat$stage == "late"]
  
  t.test(x, y, alternative = "two.sided")$p.value
})
ttest.res <- data.frame(ttest.pvalue = ttest.res)
rownames(ttest.res) <- unique(ggmat$key)
ttest.res <- ttest.res[order(ttest.res$ttest.pvalue), , drop = FALSE]
ttest.res$FDR <- p.adjust(ttest.res$"ttest.pvalue", method = "BH")

write.xlsx(ttest.res, "SupFigure6A_ttest.xlsx", row.names = TRUE)

