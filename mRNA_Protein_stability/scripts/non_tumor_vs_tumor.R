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
library(ungeviz)
library(ggrepel)
library(msigdbr)
library(circlize)
library(measurements)
library(ComplexHeatmap)

source("tools.r")
source("hyperG_parallel.R")

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


####################################
# LOAD HEALTHY COLON DATA
rna <- read.xlsx("rna_colon.xlsx", sheet = 1, rowNames = TRUE)
prot <- read.xlsx("prot_colon.xlsx", sheet = 1, rowNames = TRUE)
prot[prot == 0] <- NA

# common genes
commonGenes <- intersect(rownames(rna), rownames(prot))
rna.ctl <- rna[commonGenes, ]
prot.ctl <- prot[commonGenes, ]

colnames(rna.ctl) <- substring(gsub("TCGA-COAD.", "", colnames(rna.ctl)), 1, 12)


####################################
# INTERSECT TCGA AND CTL GENES

commonGenes <- intersect(rownames(rna.tcga), rownames(rna.ctl))

rna.tcga <- rna.tcga[commonGenes, ]
prot.tcga <- prot.tcga[commonGenes, ]

rna.ctl <- rna.ctl[commonGenes, ]
prot.ctl <- prot.ctl[commonGenes, ]

ann.sub <- ann[match(commonGenes, ann$Symbol), ]

rna.ctl[is.na(rna.ctl)] <- 0
prot.ctl[is.na(prot.ctl)] <- 0

rna.tcga[is.na(rna.tcga)] <- 0
prot.tcga[is.na(prot.tcga)] <- 0


# TCGA CLINICAL DATA
clinMat <- read.xlsx("TCGA_COAD_clinical.xlsx", sheet = 1)

####################################################################
# CHECK INFLUENCE OF AGE, GENDER, CMS.. IN TCGA SAMPLES

outDir <- file.path("../output")
dir.create(outDir)
setwd(outDir)

# Calculate mRNA / protein correaltion in TCGA
samples <- colnames(rna.tcga)

tissueNames <- colnames(rna.tcga)
corMat <- matrix(NA, nrow = length(tissueNames), ncol = 5)
rownames(corMat) <- tissueNames
colnames(corMat) <- c("rS.pS", "rU.pS", "rS.pU", "rU.pU", "grey.zone")

for(i in tissueNames){
  rna.current <- as.numeric(rna.tcga[,i])
  prot.current <- as.numeric(prot.tcga[,i])
  genes.current <- ann.sub$Symbol
  
  ann.current <- ann[match(genes.current, ann$Symbol), ]
  
  groups <- unique(ann.current$Group)
  
  # Fill in corMat
  corList <- sapply(groups, function(g){
    cor(rna.current[ann.current$Group == g],
        prot.current[ann.current$Group == g],
        method = "spearman", use = "complete.obs")
  })
  names(corList) <- groups
  
  corMat[i, names(corList)] <- corList
}

ggmat <- data.frame(SAMPLE = rownames(corMat), corMat)
ggmat <- pivot_longer(ggmat, cols = c("rS.pS", "rU.pS", "rS.pU", "rU.pU", "grey.zone"))
ggmat <- ggmat[ggmat$name != "grey.zone", ]


# Gender
ggmat$gender <- clinMat$gender[match(ggmat$SAMPLE, clinMat$submitter_id)]

# Jitter plot
p <- ggplot(ggmat, aes(x = name, y = value, fill = gender)) +
  geom_point(pch = 21, position = position_jitterdodge(jitter.width = 0.1), alpha = 0.5) + 
  geom_boxplot(outlier.size = 0, alpha = 0.1)
p <- p + theme_bw(base_size = 16)
#p <- p + scale_fill_manual(values=c("normal" = "grey80", "tumor" = "firebrick"))
p <- p + xlab("") + ylab("RNA / Protein spearman's correlation")

ggsave(plot = p, filename = "RNA_prot_correlation_jitter_GENDER.pdf", width = 7.5, height = 7)

write.xlsx(ggmat, "RNA_prot_correlation_jitter_GENDER.xlsx")

# primary_diagnosis
ggmat$primary_diag <- clinMat$primary_diagnosis[match(ggmat$SAMPLE, clinMat$submitter_id)]

# Jitter plot
p <- ggplot(ggmat[ggmat$primary_diag != "Adenosquamous carcinoma", ], aes(x = name, y = value, fill = primary_diag)) +
  geom_point(pch = 21, position = position_jitterdodge(jitter.width = 0.1), alpha = 0.5) + 
  geom_boxplot(outlier.size = 0, alpha = 0.1)
p <- p + theme_bw(base_size = 16)
#p <- p + scale_fill_manual(values=c("normal" = "grey80", "tumor" = "firebrick"))
p <- p + xlab("") + ylab("RNA / Protein spearman's correlation")

ggsave(plot = p, filename = "RNA_prot_correlation_jitter_PRIMARY_DIAG.pdf", width = 8.5, height = 7)

write.xlsx(ggmat[ggmat$primary_diag != "Adenosquamous carcinoma", ], "RNA_prot_correlation_jitter_PRIMARY_DIAG.xlsx")


# tissue_or_organ_of_origin
ggmat$tissue <- clinMat$tissue_or_organ_of_origin[match(ggmat$SAMPLE, clinMat$submitter_id)]

# Jitter plot
p <- ggplot(ggmat, aes(x = name, y = value, fill = tissue)) +
  geom_point(pch = 21, position = position_jitterdodge(jitter.width = 0.1), alpha = 0.5) + 
  geom_boxplot(outlier.size = 0, alpha = 0.1)
p <- p + theme_bw(base_size = 16)
#p <- p + scale_fill_manual(values=c("normal" = "grey80", "tumor" = "firebrick"))
p <- p + xlab("") + ylab("RNA / Protein spearman's correlation")

ggsave(plot = p, filename = "RNA_prot_correlation_jitter_TISSUE.pdf", width = 8.5, height = 7)

write.xlsx(ggmat, "RNA_prot_correlation_jitter_TISSUE.xlsx")


# below 50
ggmat$below50 <- clinMat$below50[match(ggmat$SAMPLE, clinMat$submitter_id)]

# Jitter plot
p <- ggplot(ggmat, aes(x = name, y = value, fill = below50)) +
  geom_point(pch = 21, position = position_jitterdodge(jitter.width = 0.1), alpha = 0.5) + 
  geom_boxplot(outlier.size = 0, alpha = 0.1)
p <- p + theme_bw(base_size = 16)
#p <- p + scale_fill_manual(values=c("normal" = "grey80", "tumor" = "firebrick"))
p <- p + xlab("") + ylab("RNA / Protein spearman's correlation")

ggsave(plot = p, filename = "RNA_prot_correlation_jitter_AGE.pdf", width = 7.5, height = 7)

write.xlsx(ggmat, "RNA_prot_correlation_jitter_AGE.xlsx")


# CMS
ggmat$CMS <- clinMat$CMS[match(ggmat$SAMPLE, clinMat$submitter_id)]

# Jitter plot
p <- ggplot(ggmat, aes(x = name, y = value, fill = CMS)) +
  geom_point(pch = 21, position = position_jitterdodge(jitter.width = 0.1), alpha = 0.5) + 
  geom_boxplot(outlier.size = 0, alpha = 0.1)
p <- p + theme_bw(base_size = 16)
#p <- p + scale_fill_manual(values=c("normal" = "grey80", "tumor" = "firebrick"))
p <- p + xlab("") + ylab("RNA / Protein spearman's correlation")

ggsave(plot = p, filename = "RNA_prot_correlation_jitter_CMS.pdf", width = 7.5, height = 7)

write.xlsx(ggmat, "RNA_prot_correlation_jitter_GENDER_AGE_TISSUE_CMS.xlsx")


# TEST NORMALITY

# Gender
shapi <- lapply(unique(ggmat$gender), function(i){
  sapply(unique(ggmat$name), function(j){
    shapiro.test(ggmat$value[ggmat$gender == i & ggmat$name == j])$p.value
  })
})
shapi <- do.call(cbind, shapi)
colnames(shapi) <- unique(ggmat$gender)

shapi.qv <- apply(shapi, 2, p.adjust, method = "BH")
colnames(shapi.qv) <- paste0(colnames(shapi.qv), ".FDR")
shapi <- cbind(shapi, shapi.qv)
write.xlsx(shapi, "GENDER_Shapiro_test.xlsx", row.names = TRUE)

# Primary diagnosis
primary_diag <- c("Adenocarcinoma, NOS", "Mucinous adenocarcinoma")
shapi <- lapply(primary_diag, function(i){
  sapply(unique(ggmat$name), function(j){
    shapiro.test(ggmat$value[ggmat$primary_diag == i & ggmat$name == j])$p.value
  })
})
shapi <- do.call(cbind, shapi)
colnames(shapi) <- primary_diag

shapi.qv <- apply(shapi, 2, p.adjust, method = "BH")
colnames(shapi.qv) <- paste0(colnames(shapi.qv), ".FDR")
shapi <- cbind(shapi, shapi.qv)
write.xlsx(shapi, "Primary_Diag_Shapiro_test.xlsx", row.names = TRUE)

# Tissue
tissues <- c("Ascending colon", "Cecum", "Colon, NOS", "Sigmoid colon")
shapi <- lapply(tissues, function(i){
  sapply(unique(ggmat$name), function(j){
    shapiro.test(ggmat$value[ggmat$tissue == i & ggmat$name == j])$p.value
  })
})
shapi <- do.call(cbind, shapi)
colnames(shapi) <- tissues

shapi.qv <- apply(shapi, 2, p.adjust, method = "BH")
colnames(shapi.qv) <- paste0(colnames(shapi.qv), ".FDR")
shapi <- cbind(shapi, shapi.qv)
write.xlsx(shapi, "tissue_Shapiro_test.xlsx", row.names = TRUE)

# below50
shapi <- lapply(unique(ggmat$below50), function(i){
  sapply(unique(ggmat$name), function(j){
    shapiro.test(ggmat$value[ggmat$below50 == i & ggmat$name == j])$p.value
  })
})
shapi <- do.call(cbind, shapi)
colnames(shapi) <- unique(ggmat$below50)

shapi.qv <- apply(shapi, 2, p.adjust, method = "BH")
colnames(shapi.qv) <- paste0(colnames(shapi.qv), ".FDR")
shapi <- cbind(shapi, shapi.qv)
write.xlsx(shapi, "AGE_Shapiro_test.xlsx", row.names = TRUE)

# CMS
shapi <- lapply(unique(ggmat$CMS), function(i){
  sapply(unique(ggmat$name), function(j){
    shapiro.test(ggmat$value[ggmat$CMS == i & ggmat$name == j])$p.value
  })
})
shapi <- do.call(cbind, shapi)
colnames(shapi) <- unique(ggmat$CMS)

shapi.qv <- apply(shapi, 2, p.adjust, method = "BH")
colnames(shapi.qv) <- paste0(colnames(shapi.qv), ".FDR")
shapi <- cbind(shapi, shapi.qv)
write.xlsx(shapi, "CMS_Shapiro_test.xlsx", row.names = TRUE)

###############
# WILCOXON TEST

# Gender
combos <- combn(unique(ggmat$gender),2)
wilcox.res <- adply(combos, 2, function(x) {
  
  wList <- lapply(unique(ggmat$name), function(i){
    test <- wilcox.test(ggmat$value[ggmat$gender == x[1] & ggmat$name == i],
                        ggmat$value[ggmat$gender == x[2] & ggmat$name == i], exact = TRUE)
    
    out <- data.frame("var1" = x[1]
                      , "var2" = x[2]
                      , "var1.mean" = mean(ggmat$value[ggmat$gender == x[1] & ggmat$name == i])
                      , "var2.mean" = mean(ggmat$value[ggmat$gender == x[2] & ggmat$name == i])
                      , "p.value" = sprintf("%.3f", test$p.value)
    )
    return(out)
  })
  return(do.call(rbind, wList))
})
wilcox.res$X1 <- unique(ggmat$name)
wilcox.res$p.value <- toNum(wilcox.res$p.value)
wilcox.res <- wilcox.res[order(wilcox.res$p.value), ]
wilcox.res$FDR <- p.adjust(wilcox.res$"p.value", method = "BH")
write.xlsx(wilcox.res, "GENDER_wilcox_test.xlsx", row.names = FALSE)

# Primary diagnosis
combos <- combn(primary_diag,2)
wilcox.res <- adply(combos, 2, function(x) {
  
  wList <- lapply(unique(ggmat$name), function(i){
    test <- wilcox.test(ggmat$value[ggmat$primary_diag == x[1] & ggmat$name == i],
                        ggmat$value[ggmat$primary_diag == x[2] & ggmat$name == i], exact = TRUE)
    
    out <- data.frame("var1" = x[1]
                      , "var2" = x[2]
                      , "var1.mean" = mean(ggmat$value[ggmat$primary_diag == x[1] & ggmat$name == i])
                      , "var2.mean" = mean(ggmat$value[ggmat$primary_diag == x[2] & ggmat$name == i])
                      , "p.value" = sprintf("%.3f", test$p.value)
                      
    )
    return(out)
  })
  return(do.call(rbind, wList))
})
wilcox.res$X1 <- unique(ggmat$name)
wilcox.res$p.value <- toNum(wilcox.res$p.value)
wilcox.res <- wilcox.res[order(wilcox.res$p.value), ]
wilcox.res$FDR <- p.adjust(wilcox.res$"p.value", method = "BH")
write.xlsx(wilcox.res, "PRIMARY_DIAG_wilcox_test.xlsx", row.names = FALSE)

# tissue
combos <- combn(tissues,2)
wilcox.res <- adply(combos, 2, function(x) {
  
  wList <- lapply(unique(ggmat$name), function(i){
    test <- wilcox.test(ggmat$value[ggmat$tissue == x[1] & ggmat$name == i],
                        ggmat$value[ggmat$tissue == x[2] & ggmat$name == i], exact = TRUE)
    
    out <- data.frame("var1" = x[1]
                      , "var2" = x[2]
                      , "var1.mean" = mean(ggmat$value[ggmat$tissue == x[1] & ggmat$name == i])
                      , "var2.mean" = mean(ggmat$value[ggmat$tissue == x[2] & ggmat$name == i])
                      , "p.value" = sprintf("%.3f", test$p.value)
                      
    )
    return(out)
  })
  return(do.call(rbind, wList))
})
wilcox.res$X1 <- unique(ggmat$name)
wilcox.res$p.value <- toNum(wilcox.res$p.value)
wilcox.res <- wilcox.res[order(wilcox.res$p.value), ]
wilcox.res$FDR <- p.adjust(wilcox.res$"p.value", method = "BH")
write.xlsx(wilcox.res, "TISSUE_wilcox_test.xlsx", row.names = FALSE)

#
ttest.res <- adply(combos, 2, function(x) {
  
  wList <- lapply(unique(ggmat$name), function(i){
    test <- t.test(ggmat$value[ggmat$tissue == x[1] & ggmat$name == i],
                   ggmat$value[ggmat$tissue == x[2] & ggmat$name == i], alternative = "two.sided")
    
    out <- data.frame("var1" = x[1]
                      , "var2" = x[2]
                      , "var1.mean" = mean(ggmat$value[ggmat$tissue == x[1] & ggmat$name == i])
                      , "var2.mean" = mean(ggmat$value[ggmat$tissue == x[2] & ggmat$name == i])
                      , "p.value" = sprintf("%.3f", test$p.value)
                      
    )
    return(out)
  })
  return(do.call(rbind, wList))
})
ttest.res$X1 <- unique(ggmat$name)
ttest.res$p.value <- toNum(ttest.res$p.value)
ttest.res <- ttest.res[order(ttest.res$p.value), ]
ttest.res$FDR <- p.adjust(ttest.res$"p.value", method = "BH")
write.xlsx(ttest.res, "TISSUE_ttest.xlsx", row.names = FALSE)



# BELOW 50
combos <- combn(unique(ggmat$below50),2)
wilcox.res <- adply(combos, 2, function(x) {
  
  wList <- lapply(unique(ggmat$name), function(i){
    test <- wilcox.test(ggmat$value[ggmat$below50 == x[1] & ggmat$name == i],
                        ggmat$value[ggmat$below50 == x[2] & ggmat$name == i], exact = TRUE)
    
    out <- data.frame("var1" = x[1]
                      , "var2" = x[2]
                      , "var1.mean" = mean(ggmat$value[ggmat$below50 == x[1] & ggmat$name == i])
                      , "var2.mean" = mean(ggmat$value[ggmat$below50 == x[2] & ggmat$name == i])
                      , "p.value" = sprintf("%.3f", test$p.value)
    )
    return(out)
  })
  return(do.call(rbind, wList))
})
wilcox.res$X1 <- unique(ggmat$name)
wilcox.res$p.value <- toNum(wilcox.res$p.value)
wilcox.res <- wilcox.res[order(wilcox.res$p.value), ]
wilcox.res$FDR <- p.adjust(wilcox.res$"p.value", method = "BH")
write.xlsx(wilcox.res, "AGE_wilcox_test.xlsx", row.names = FALSE)

# CMS
combos <- combn(unique(ggmat$CMS),2)
wilcox.res <- adply(combos, 2, function(x) {
  
  wList <- lapply(unique(ggmat$name), function(i){
    test <- wilcox.test(ggmat$value[ggmat$CMS == x[1] & ggmat$name == i],
                        ggmat$value[ggmat$CMS == x[2] & ggmat$name == i], exact = TRUE)
    
    out <- data.frame("var1" = x[1]
                      , "var2" = x[2]
                      , "var1.mean" = mean(ggmat$value[ggmat$CMS == x[1] & ggmat$name == i])
                      , "var2.mean" = mean(ggmat$value[ggmat$CMS == x[2] & ggmat$name == i])
                      , "p.value" = sprintf("%.3f", test$p.value)
    )
    return(out)
  })
  return(do.call(rbind, wList))
})
wilcox.res$X1 <- unique(ggmat$name)
wilcox.res$p.value <- toNum(wilcox.res$p.value)
wilcox.res <- wilcox.res[order(wilcox.res$p.value), ]
wilcox.res$FDR <- p.adjust(wilcox.res$"p.value", method = "BH")
write.xlsx(wilcox.res, "CMS_wilcox_test.xlsx", row.names = FALSE)


####################################################################
# CHECK INFLUENCE OF MOLECULAR FEATURES (PROLIFERATION, HYPOXIA, EMT...)

scMat <- read.xlsx(file.path("../data/", "H_NES_matrix.xlsx"), sheet = 1, rowNames = TRUE)
scMat <- scMat[, match(samples, colnames(scMat))]
# Define high / low groups per gene-sets
scMat.group <- lapply(1:nrow(scMat), function(i){
  sc <- as.numeric(scMat[i, ])
  sc.median <- median(sc)
  group <- ifelse(sc >= sc.median, "high", "low")
  return(group)
})

scMat.group <- lapply(1:nrow(scMat), function(i){
  sc <- as.numeric(scMat[i, ])
  sc.quantile <- quantile(sc, probs = c(1/3, 2/3))
  group <- rep("middle", length(sc))
  group[sc < sc.quantile[1]] <- "low"
  group[sc > sc.quantile[2]] <- "high"
  return(group)
})


scMat.group <- do.call(cbind, scMat.group)
rownames(scMat.group) <- samples
colnames(scMat.group) <- rownames(scMat)

# Jitter plot
lapply(1:ncol(scMat.group), function(i){
  ggmat$NES <- scMat.group[match(ggmat$SAMPLE, rownames(scMat.group)), i]
  ggmat$NES <- factor(ggmat$NES, levels = c("low", "middle", "high"))
  
  p <- ggplot(ggmat, aes(x = name, y = value, fill = NES)) +
    geom_point(pch = 21, position = position_jitterdodge(jitter.width = 0.1), alpha = 0.5) + 
    geom_boxplot(outlier.size = 0, alpha = 0.1)
  p <- p + theme_bw(base_size = 16)
  #p <- p + scale_fill_manual(values=c("normal" = "grey80", "tumor" = "firebrick"))
  p <- p + xlab("") + ylab("RNA / Protein spearman's correlation")
  

  ggsave(plot = p, filename = paste0(colnames(scMat.group)[i], "_RNA_prot_correlation_jitter.pdf"), width = 7.5, height = 7)
  
  write.xlsx(ggmat, paste0(colnames(scMat.group)[i], "_RNA_prot_correlation_jitter.xlsx"))
  
  # SHAPIRO
  shapi <- lapply(unique(ggmat$NES), function(i){
    sapply(unique(ggmat$name), function(j){
      shapiro.test(ggmat$value[ggmat$NES == i & ggmat$name == j])$p.value
    })
  })
  shapi <- do.call(cbind, shapi)
  colnames(shapi) <- unique(ggmat$NES)
  
  shapi.qv <- apply(shapi, 2, p.adjust, method = "BH")
  colnames(shapi.qv) <- paste0(colnames(shapi.qv), ".FDR")
  shapi <- cbind(shapi, shapi.qv)
  write.xlsx(shapi, paste0(colnames(scMat.group)[i], "_Shapiro_test.xlsx"), row.names = TRUE)
  
  # WILCOXON
  combos <- combn(unique(ggmat$NES),2)
  wilcox.res <- adply(combos, 2, function(x) {
    x <- as.character(x)
    wList <- lapply(unique(ggmat$name), function(i){
      test <- wilcox.test(ggmat$value[ggmat$NES == x[1] & ggmat$name == i],
                          ggmat$value[ggmat$NES == x[2] & ggmat$name == i], exact = TRUE)
      
      out <- data.frame("var1" = x[1]
                        , "var2" = x[2]
                        , "var1.mean" = mean(ggmat$value[ggmat$NES == x[1] & ggmat$name == i])
                        , "var2.mean" = mean(ggmat$value[ggmat$NES == x[2] & ggmat$name == i])
                        , "p.value" = sprintf("%.3f", test$p.value)
      )
      return(out)
    })
    return(do.call(rbind, wList))
  })
  wilcox.res$X1 <- rep(unique(ggmat$name), ncol(combos))
  wilcox.res$p.value <- as.numeric(wilcox.res$p.value)
  wilcox.res <- wilcox.res[order(wilcox.res$p.value), ]
  wilcox.res$FDR <- p.adjust(wilcox.res$"p.value", method = "BH")
  write.xlsx(wilcox.res, paste0(colnames(scMat.group)[i], "_wilcox_test.xlsx"), row.names = FALSE)
  
  
  
  
})

####################################################################
# CHECK INFLUENCE OF TUMOR PURITY

tpMat <- read.xlsx(file.path("../data", "Tumor purity.xlsx"), sheet = 1, rowNames = TRUE)
rownames(tpMat) <- sapply(strsplit(rownames(tpMat), split = "-"),
                          function(i) paste(i[1:3], collapse = "-"))
tpMat <- t(tpMat[, c("IHC", "CPE")])
tpMat <- tpMat[, match(samples, colnames(tpMat))]
# Define high / low groups per gene-sets
tpMat.group <- lapply(1:nrow(tpMat), function(i){
  sc <- as.numeric(tpMat[i, ])
  sc.median <- median(sc)
  group <- ifelse(sc >= sc.median, "high", "low")
  return(group)
})

tpMat.group <- lapply(1:nrow(tpMat), function(i){
  sc <- as.numeric(tpMat[i, ])
  sc.quantile <- quantile(sc, probs = c(1/3, 2/3))
  group <- rep("middle", length(sc))
  group[sc < sc.quantile[1]] <- "low"
  group[sc > sc.quantile[2]] <- "high"
  return(group)
})


tpMat.group <- do.call(cbind, tpMat.group)
rownames(tpMat.group) <- samples
colnames(tpMat.group) <- rownames(tpMat)

# Jitter plot
lapply(1:ncol(tpMat.group), function(i){
  ggmat$NES <- tpMat.group[match(ggmat$SAMPLE, rownames(tpMat.group)), i]
  ggmat$NES <- factor(ggmat$NES, levels = c("low", "middle", "high"))
  
  p <- ggplot(ggmat, aes(x = name, y = value, fill = NES)) +
    geom_point(pch = 21, position = position_jitterdodge(jitter.width = 0.1), alpha = 0.5) + 
    geom_boxplot(outlier.size = 0, alpha = 0.1)
  p <- p + theme_bw(base_size = 16)
  #p <- p + scale_fill_manual(values=c("normal" = "grey80", "tumor" = "firebrick"))
  p <- p + xlab("") + ylab("RNA / Protein spearman's correlation")
  
  
  ggsave(plot = p, filename = paste0(colnames(tpMat.group)[i], "_RNA_prot_correlation_jitter.pdf"), width = 7.5, height = 7)
  
  write.xlsx(ggmat, paste0(colnames(tpMat.group)[i], "_RNA_prot_correlation_jitter.xlsx"))
  
  # SHAPIRO
  shapi <- lapply(unique(ggmat$NES), function(i){
    sapply(unique(ggmat$name), function(j){
      shapiro.test(ggmat$value[ggmat$NES == i & ggmat$name == j])$p.value
    })
  })
  shapi <- do.call(cbind, shapi)
  colnames(shapi) <- unique(ggmat$NES)
  
  shapi.qv <- apply(shapi, 2, p.adjust, method = "BH")
  colnames(shapi.qv) <- paste0(colnames(shapi.qv), ".FDR")
  shapi <- cbind(shapi, shapi.qv)
  write.xlsx(shapi, paste0(colnames(tpMat.group)[i], "_Shapiro_test.xlsx"), row.names = TRUE)
  
  # WILCOXON
  combos <- combn(unique(ggmat$NES),2)
  wilcox.res <- adply(combos, 2, function(x) {
    x <- as.character(x)
    wList <- lapply(unique(ggmat$name), function(i){
      test <- wilcox.test(ggmat$value[ggmat$NES == x[1] & ggmat$name == i],
                          ggmat$value[ggmat$NES == x[2] & ggmat$name == i], exact = TRUE)
      
      out <- data.frame("var1" = x[1]
                        , "var2" = x[2]
                        , "var1.mean" = mean(ggmat$value[ggmat$NES == x[1] & ggmat$name == i])
                        , "var2.mean" = mean(ggmat$value[ggmat$NES == x[2] & ggmat$name == i])
                        , "p.value" = sprintf("%.3f", test$p.value)
      )
      return(out)
    })
    return(do.call(rbind, wList))
  })
  wilcox.res$X1 <- rep(unique(ggmat$name), ncol(combos))
  wilcox.res$p.value <- as.numeric(wilcox.res$p.value)
  wilcox.res <- wilcox.res[order(wilcox.res$p.value), ]
  wilcox.res$FDR <- p.adjust(wilcox.res$"p.value", method = "BH")
  write.xlsx(wilcox.res, paste0(colnames(tpMat.group)[i], "_wilcox_test.xlsx"), row.names = FALSE)
  
  
  
  
})






######################################################################
# UNPAIRED ANALYSIS


# GENERATE RANDOM UNPAIRED SAMPLES
nb <- 20
comb.ctl <- read.xlsx(file.path("../data", "comb20.xlsx"), sheet = "ctl")
comb.tcga <- read.xlsx(file.path("../data", "comb20.xlsx"), sheet = "tcga")

rna.ctl.sub <- rna.ctl[, unique(comb.ctl$Var1)]
prot.ctl.sub <- prot.ctl[, unique(comb.ctl$Var2)]
rna.tcga.sub <- rna.tcga[, unique(comb.tcga$Var1)]
prot.tcga.sub <- prot.tcga[, unique(comb.tcga$Var2)]


# CALCULATE CORRELATION OF REAL ADN RANDOM GROUPS
stabGroup <- unique(ann.sub$Group)

ggmat <- c()

ggmat.delta <- c()

geneRankList <- list()

for(i in stabGroup){
  mygenes <- ann.sub$Symbol[ann.sub$Group == i]
  
  # get gene ranke delta
  geneRankMat.ctl <- matrix(NA, nrow = length(mygenes), ncol = nrow(comb.ctl))
  geneRankMat.tcga <- matrix(NA, nrow = length(mygenes), ncol = nrow(comb.ctl))
  for(j in 1:nrow(comb.ctl)){
    rank.ctl <- abs(rank(rna.ctl.sub[mygenes, comb.ctl[j, 1]]) - rank(prot.ctl.sub[mygenes, comb.ctl[j, 2]]))
    rank.tcga <- abs(rank(rna.tcga.sub[mygenes, comb.tcga[j, 1]]) - rank(prot.tcga.sub[mygenes, comb.tcga[j, 2]]))                 
    geneRankMat.ctl[, j] <- rank.ctl
    geneRankMat.tcga[, j] <- rank.tcga
  }
  rownames(geneRankMat.ctl) <- mygenes
  rownames(geneRankMat.tcga) <- mygenes
  
  geneRankList[[i]] <- list(ctl = geneRankMat.ctl, tcga = geneRankMat.tcga)
  
  # get coorelation in "NORMAL"
  corList.ctl <- sapply(1:nrow(comb.ctl), function(j){
    cor(rna.ctl.sub[mygenes, comb.ctl[j, 1]], prot.ctl.sub[mygenes, comb.ctl[j, 2]], method = "spearman", use = "complete.obs")
  })
  
  # get coorelation in "TUMOR"
  corList.tcga <- sapply(1:nrow(comb.tcga), function(j){
    cor(rna.tcga.sub[mygenes, comb.tcga[j, 1]], prot.tcga.sub[mygenes, comb.tcga[j, 2]], method = "spearman", use = "complete.obs")
  })
  
  pv <- t.test(corList.ctl, corList.tcga)$"p.value"
  delta <- median(corList.ctl) - median(corList.tcga)
  print(delta)
  
  # fill in ggmat
  ggmat <- rbind(ggmat, data.frame(cor = c(corList.ctl, corList.tcga),
                                   stabGroup = i, 
                                   Condition = c(rep("normal", length(corList.ctl)), rep("tumor", length(corList.tcga)))
  ))
  
  # get correlation in "NORMAL" and "TUMOR" (random genes)
  nbRd <- 1000
  delta.list <- mclapply(1:nbRd, function(k){
    idx <- sample(rownames(rna.ctl.sub), length(mygenes))
    
    cor.ctl.current <- sapply(1:nrow(comb.ctl), function(j){
      cor(rna.ctl.sub[idx, comb.ctl[j, 1]], prot.ctl.sub[idx, comb.ctl[j, 2]], method = "spearman", use = "complete.obs")
    })
    
    cor.tcga.current <- sapply(1:nrow(comb.tcga), function(j){
      cor(rna.tcga.sub[idx, comb.tcga[j, 1]], prot.tcga.sub[idx, comb.tcga[j, 2]], method = "spearman", use = "complete.obs")
    })
    
    #return(t.test(cor.ctl.current, cor.tcga.current))
    return(median(cor.ctl.current) - median(cor.tcga.current))# delta of correlation with random genes
  }, mc.cores = 8)
  #pv.random <- sapply(t.test.list, function(k) k$p.value)
  delta.random <- unlist(delta.list)
  print(median(delta.random))
  
  delta.pv <- getEmpiricalPV(delta.random, delta, "less")
  
  ggmat.delta <- rbind(ggmat.delta, data.frame(delta = c(delta, delta.random),
                                               stabGroup = i, 
                                               Condition = c(rep("observed", length(delta)), rep("random", length(delta.random)))
  ))
  print(delta.pv)# is the observed difference between ctl and tcga more significnat than expected by chance
  
  # delta vs. avg expression
}

# Jitter plot
p <- ggplot(ggmat[ggmat$stabGroup != "grey.zone", ], aes(x = stabGroup, y = cor, fill = Condition)) +
  geom_point(pch = 21, position = position_jitterdodge(), alpha = 0.5) + 
  geom_boxplot(outlier.size = 0, alpha = 0.1)
p <- p + theme_bw(base_size = 16)
p <- p + scale_fill_manual(values=c("normal" = "grey80", "tumor" = "firebrick"))
p <- p + xlab("") + ylab("RNA / Protein spearman's correlation")

setwd(outDir)
#ggsave(plot = p, filename = "RNA_prot_correlation_jitter.pdf", width = 7.5, height = 7)
ggsave(plot = p, filename = "figure4A.pdf", width = 7.5, height = 7)
#write.xlsx(ggmat, "RNA_prot_correlation_jitter.xlsx")
write.xlsx(ggmat, "figure4A_data.xlsx")


# TEST NORMALITY
shapi <- lapply(unique(stabGroup), function(i){
  x <- ggmat$cor[ggmat$stabGroup == i & ggmat$Condition == "normal"]
  shapiro.test(x)
})
shapi <- sapply(shapi, function(i) i$p.value)
shapi.normal <- data.frame(normal.shapiro.pvalue = shapi)

shapi <- lapply(unique(stabGroup), function(i){
  x <- ggmat$cor[ggmat$stabGroup == i & ggmat$Condition == "tumor"]
  shapiro.test(x)
})
shapi <- sapply(shapi, function(i) i$p.value)
shapi.coad <- data.frame(tumor.shapiro.pvalue = shapi)

shapi <- cbind(shapi.normal, shapi.coad)
rownames(shapi) <- unique(stabGroup)
write.xlsx(shapi, "figure4A_Shapiro_test.xlsx", row.names = TRUE)

# t.test normal vs. tumor per stability group
wilcox.res <- sapply(unique(stabGroup), function(i){
  x <- ggmat$cor[ggmat$stabGroup == i & ggmat$Condition == "normal"]
  y <- ggmat$cor[ggmat$stabGroup == i & ggmat$Condition == "tumor"]
  wilcox.test(x, y)$p.value
})
wilcox.res <- data.frame(wilcox.res)
write.xlsx(wilcox.res, "figure4A_wilcox_test.xlsx", row.names = TRUE)


# delta density plot (Sup Figure 3D)
ggmat.delta <- ggmat.delta[ggmat.delta$stabGroup != "grey.zone", ]

p <- ggplot(ggmat.delta, aes(x=stabGroup, y = delta))
p <- p + geom_violin(data = subset(ggmat.delta, Condition == "random"), fill = "grey80")
p <- p + geom_hpline(data = subset(ggmat.delta, Condition != "random"), col = "firebrick")
p <- p + theme_bw(base_size = 18)
p <- p + xlab("") + ylab("Normal - Tumor correlation (RNA / Protein)")
#ggsave(plot = p, filename = "delta_observed_vs_random.pdf", width = 7, height = 6)
ggsave(plot = p, filename = "SupFigure5D.pdf", width = 7, height = 6)
write.xlsx(ggmat.delta, "SupFigure5D_data.xlsx")

# TEST NORMALITY
shapi <- lapply(unique(ggmat.delta$stabGroup), function(i){
  x <- ggmat.delta$delta[ggmat.delta$stabGroup == i & ggmat.delta$Condition == "random"]
  shapiro.test(x)
})
shapi <- sapply(shapi, function(i) i$p.value)
shapi <- data.frame(shapiro.pvalue = shapi)
write.xlsx(shapi, "SupFigure5D_Shapiro_test.xlsx")

ttest.res <- lapply(unique(ggmat.delta$stabGroup), function(i){
  x <- ggmat.delta$delta[ggmat.delta$stabGroup == i & ggmat.delta$Condition == "random"]
  y <- ggmat.delta$delta[ggmat.delta$stabGroup == i & ggmat.delta$Condition == "observed"]
  pnorm(y, mean = mean(x), sd = sd(x), lower.tail = FALSE)# UPPER TAIL
})
ttest.res <- do.call(rbind, ttest.res)
rownames(ttest.res) <- unique(ggmat.delta$stabGroup)
colnames(ttest.res) <- "ttest.pvalue"  
write.xlsx(ttest.res, "SupFigure5D_ttest.xlsx")

##############################
# GENE RANK CTL vs. TCGA DELTA

rankDeltaPV <- lapply(geneRankList, function(i){
  rank.ctl <- i$ctl
  rank.tcga <- i$tcga
  
  rank.ctl.median <- apply(rank.ctl, 1, median)
  rank.tcga.median <- apply(rank.tcga, 1, median)
  
  delta.median <- rank.ctl.median - rank.tcga.median
  
  pv <- sapply(1:nrow(rank.ctl), function(j) wilcox.test(rank.ctl[j, ], rank.tcga[j, ], alternative = "two.sided")$"p.value")
  qv <- p.adjust(pv, method = "BH")
  
  pv.less <- sapply(1:nrow(rank.ctl), function(j) wilcox.test(rank.ctl[j, ], rank.tcga[j, ], alternative = "less")$"p.value")
  qv.less <- p.adjust(pv.less, method = "BH")
  
  pv.greater <- sapply(1:nrow(rank.ctl), function(j) wilcox.test(rank.ctl[j, ], rank.tcga[j, ], alternative = "greater")$"p.value")
  qv.greater <- p.adjust(pv.greater, method = "BH")
  
  df <- data.frame(GENE = rownames(rank.ctl),
                   rank.ctl.median = rank.ctl.median,
                   rank.tcga.median = rank.tcga.median,
                   delta.median = delta.median,
                   pv = pv,
                   qv = qv,
                   pv.less = pv.less,
                   qv.less = qv.less,
                   pv.greater = pv.greater,
                   qv.greater = qv.greater)
  return(df)
})


# PLOT
nbToShow <- 10
ggmatList <- lapply(1:length(rankDeltaPV), function(i){
  ggmat <- data.frame(rankDeltaPV[[i]])
  ggmat$stabGroup <- names(rankDeltaPV)[i]
  ggmat$Score <- ifelse(ggmat$delta.median > 0, -log10(ggmat$qv), log10(ggmat$qv))
  ggmat$toShow <- "no"
  ggmat$toShow[order(ggmat$delta.median)[1:nbToShow]] <- "yes"
  ggmat$toShow[order(-ggmat$delta.median)[1:nbToShow]] <- "yes"
  return(ggmat)
})
ggmat <- do.call(rbind, ggmatList)
ggmat <- ggmat[ggmat$stabGroup != "grey.zone", ]
ggmat$GENE <- factor(ggmat$GENE, levels = ggmat$GENE[order(-ggmat$delta.median)])

p <- ggplot(data=ggmat, aes(x=GENE, y=delta.median, fill = Score, label = GENE)) +
  geom_bar(stat="identity")
p <- p + scale_fill_gradient2(midpoint=0,  low="navy", mid="snow",
                              high="firebrick",  space = "Lab")
p <- p + geom_text_repel(
  force = 3,
  nudge_y      = -0.05,
  direction    = "x",
  angle        = 45,
  vjust        = 2,
  segment.size = 0.2,
  box.padding = 1,
  max.iter = 100000,
  data = subset(ggmat, toShow == "yes" & delta.median < 0)
) 
p <- p + geom_text_repel(
  force = 5,
  nudge_y      = -0.05,
  direction    = "x",
  angle        = 45,
  vjust        = -2,
  segment.size = 0.2,
  box.padding = 1,
  max.iter = 100000,
  data = subset(ggmat, toShow == "yes" & delta.median > 0)
) 
p <- p + theme(axis.title.x=element_blank(),
               axis.text.x=element_blank(),
               axis.ticks.x=element_blank())
p <- p + ylab("absolute delta rank normal - tumor")
p <- p + facet_wrap(~stabGroup, ncol = 2, scales='free')

setwd(outDir)
#ggsave(plot = p, filename = "gene_delta_rank_barplot.pdf", width = 7, height = 7)
ggsave(plot = p, filename = "figure4B.pdf", width = 7, height = 7)

#write.xlsx(rankDeltaPV, "gene_delta_rank_Wilcox.xlsx")
write.xlsx(rankDeltaPV, "figure4B_data.xlsx")


# FISHER (DELTA RANK)

lapply(1:length(dbList), function(j){
  mydb <- dbList[[j]]
  mydb.name <- names(dbList)[j]
  dir.create(mydb.name)
  
  lapply(unique(ggmat$stabGroup), function(i){
    m.current <- ggmat[ggmat$stabGroup == i, ]
    univ.current <- symbol2entrez(m.current$GENE)
    
    nb10 <- round(length(univ.current) * 0.1)
    nb25 <- round(length(univ.current) * 0.25)
    symbolList <- list(up10pc = m.current$GENE[order(-m.current$delta.median)][1:nb10],
                       up25pc = m.current$GENE[order(-m.current$delta.median)][1:nb25],
                       down10pc = m.current$GENE[order(m.current$delta.median)][1:nb10],
                       down25pc = m.current$GENE[order(m.current$delta.median)][1:nb25])
    print(sapply(symbolList, length))
    print(length(univ.current))
    entrezList <- lapply(symbolList, symbol2entrez)
    
    fhList <- lapply(entrezList, hyperG, geneSets = mydb, universe = univ.current,
                     org.library = "org.Hs.eg.db", cutoff = 1, mincount = 2)
    
    write.xlsx(fhList, file.path(mydb.name, paste0(i, "_", mydb.name, "_Fisher.xlsx")))
  })
  
})

#############################
# DELTA VS. MEDIAN EXPRESSION


pvcutoff <- 0.05

rna.ctl.sub.save <- rna.ctl.sub
rna.tcga.sub.save <- rna.tcga.sub
prot.ctl.sub.save <- prot.ctl.sub
prot.tcga.sub.save <- prot.tcga.sub

rna.ctl.sub <- scale(rna.ctl.sub)
rna.tcga.sub <- scale(rna.tcga.sub)
prot.ctl.sub <- scale(prot.ctl.sub)
prot.tcga.sub <- scale(prot.tcga.sub)

rna.ctl.median <- scale(apply(rna.ctl.sub, 1, median))[,1]
rna.tcga.median <- scale(apply(rna.tcga.sub, 1, median))[,1]
prot.ctl.median <- scale(apply(prot.ctl.sub, 1, median))[,1]
prot.tcga.median <- scale(apply(prot.tcga.sub, 1, median))[,1]

rna.ctl.median <- apply(rna.ctl.sub, 1, median)
rna.tcga.median <- apply(rna.tcga.sub, 1, median)
prot.ctl.median <- apply(prot.ctl.sub, 1, median)
prot.tcga.median <- apply(prot.tcga.sub, 1, median)


rankDeltaPV.cor <- lapply(rankDeltaPV, function(i){
  mygenes <- as.character(i$GENE)
  delta <- i $delta.median
  
  nb25 <- round(length(mygenes) * 0.25)
  mygenes.up <- as.character(mygenes[order(i$pv.greater)[1:nb25]])
  mygenes.down <- as.character(mygenes[order(i$pv.less)[1:nb25]])
  
  tt.rna.up.greater <- t.test(rna.ctl.median[mygenes.up], rna.tcga.median[mygenes.up], alternative = "greater")$p.value
  tt.prot.up.greater <- t.test(prot.ctl.median[mygenes.up], prot.tcga.median[mygenes.up], alternative = "greater")$p.value
  
  tt.rna.down.greater <- t.test(rna.ctl.median[mygenes.down], rna.tcga.median[mygenes.down], alternative = "greater")$p.value
  tt.prot.down.greater <- t.test(prot.ctl.median[mygenes.down], prot.tcga.median[mygenes.down], alternative = "greater")$p.value
  
  tt.rna.up.less <- t.test(rna.ctl.median[mygenes.up], rna.tcga.median[mygenes.up], alternative = "less")$p.value
  tt.prot.up.less <- t.test(prot.ctl.median[mygenes.up], prot.tcga.median[mygenes.up], alternative = "less")$p.value
  
  tt.rna.down.less <- t.test(rna.ctl.median[mygenes.down], rna.tcga.median[mygenes.down], alternative = "less")$p.value
  tt.prot.down.less <- t.test(prot.ctl.median[mygenes.down], prot.tcga.median[mygenes.down], alternative = "less")$p.value
  
  
  return(c(tt.rna.up.greater, tt.rna.up.less,
           tt.prot.up.greater, tt.prot.up.less,
           tt.rna.down.greater, tt.rna.down.less,
           tt.prot.down.greater, tt.prot.down.less))
  
  #return(mycor)
})

rankDeltaPV.cor <- do.call(rbind, rankDeltaPV.cor)# Normal greater
colnames(rankDeltaPV.cor) <- c("tt.rna.up.greater", "tt.rna.up.less",
                               "tt.prot.up.greater", "tt.prot.up.less",
                               "tt.rna.down.greater", "tt.rna.down.less",
                               "tt.prot.down.greater", "tt.prot.down.less")


write.xlsx(rankDeltaPV.cor, "Normal_Tumor_delta_abundance_TOP0.25.xlsx", row.names = TRUE)


plotDir <- file.path(outDir, "plot")
dir.create(plotDir)
setwd(plotDir)

geneExpDelta <- lapply(rankDeltaPV, function(i){
  mygenes <- as.character(i$GENE)
  delta <- i $delta.median
  
  # TWO-SIDED
  rna.pv <- sapply(mygenes, function(j){
    t.test(as.numeric(rna.ctl.sub[j, ]), as.numeric(rna.tcga.sub[j, ]), alternative = "two.sided")$p.value
  })
  names(rna.pv) <- mygenes
  rna.pv <- signif(rna.pv, digits = 3)
  rna.qv <- p.adjust(rna.pv, method = "BH")
  
  prot.pv <- sapply(mygenes, function(j){
    t.test(as.numeric(prot.ctl.sub[j, ]), as.numeric(prot.tcga.sub[j, ]), alternative = "two.sided")$p.value
  })
  names(prot.pv) <- mygenes
  prot.pv <- signif(prot.pv, digits = 3)
  prot.qv <- p.adjust(prot.pv, method = "BH")
  
  # LESS
  rna.pv.less <- sapply(mygenes, function(j){
    t.test(as.numeric(rna.ctl.sub[j, ]), as.numeric(rna.tcga.sub[j, ]), alternative = "less")$p.value
  })
  names(rna.pv.less) <- mygenes
  rna.pv.less <- signif(rna.pv.less, digits = 3)
  rna.qv.less <- p.adjust(rna.pv.less, method = "BH")
  
  prot.pv.less <- sapply(mygenes, function(j){
    t.test(as.numeric(prot.ctl.sub[j, ]), as.numeric(prot.tcga.sub[j, ]), alternative = "less")$p.value
  })
  names(prot.pv.less) <- mygenes
  prot.pv.less <- signif(prot.pv.less, digits = 3)
  prot.qv.less <- p.adjust(prot.pv.less, method = "BH")
  
  # GREATER
  rna.pv.greater <- sapply(mygenes, function(j){
    t.test(as.numeric(rna.ctl.sub[j, ]), as.numeric(rna.tcga.sub[j, ]), alternative = "greater")$p.value
  })
  names(rna.pv.greater) <- mygenes
  rna.pv.greater <- signif(rna.pv.greater, digits = 3)
  rna.qv.greater <- p.adjust(rna.pv.greater, method = "BH")
  
  prot.pv.greater <- sapply(mygenes, function(j){
    t.test(as.numeric(prot.ctl.sub[j, ]), as.numeric(prot.tcga.sub[j, ]), alternative = "greater")$p.value
  })
  names(prot.pv.greater) <- mygenes
  prot.pv.greater <- signif(prot.pv.greater, digits = 3)
  prot.qv.greater <- p.adjust(prot.pv.greater, method = "BH")
  
  # plot
  if(TRUE){
    lapply(mygenes, function(j){
      
      ggmat <- data.frame(X = c(rep("RNA", length(as.numeric(rna.ctl.sub[j, ]))*2),
                                rep("Protein", length(as.numeric(prot.ctl.sub[j, ]))*2)),
                          Y = c(as.numeric(rna.ctl.sub[j, ]), as.numeric(rna.tcga.sub[j, ]),
                                as.numeric(prot.ctl.sub[j, ]), as.numeric(prot.tcga.sub[j, ])),
                          Condition = c(rep("normal", length(as.numeric(rna.ctl.sub[j, ]))),
                                        rep("tumor", length(as.numeric(rna.tcga.sub[j, ]))),
                                        rep("normal", length(as.numeric(prot.ctl.sub[j, ]))),
                                        rep("tumor", length(as.numeric(prot.tcga.sub[j, ]))))
      )
      ggmat$X <- factor(ggmat$X, levels = c("RNA", "Protein"))
      
      # Jitter plot
      p <- ggplot(ggmat, aes(x = X, y = Y, fill = Condition)) +
        geom_point(pch = 21, position = position_jitterdodge(), alpha = 0.75) + 
        geom_boxplot(outlier.size = 0, alpha = 0.5)
      p <- p + theme_bw(base_size = 16)
      p <- p + scale_fill_manual(values=c("normal" = rgb(137,167,134, maxColorValue = 255),
                                          "tumor" = rgb(196, 146, 149, maxColorValue = 255)))
      p <- p + xlab("") + ylab("Scaled Intensity")
      p <- p + ggtitle(paste0(j, "; RNA pv: ", rna.pv.less[j], "; Protein pv: ", prot.pv.less[j]))
      
      ggsave(plot = p, filename = paste0(j, "_exp_boxplot.pdf"), width = 7, height = 6)
      write.xlsx(ggmat, paste0(j, "_exp_boxplot.xlsx"))
    })
  }
  
  df <- cbind(i,
              rna.pv = rna.pv, rna.qv = rna.qv,
              prot.pv = prot.pv, prot.qv = prot.qv,
              rna.pv.less = rna.pv.less, rna.qv.less = rna.qv.less,
              prot.pv.less = prot.pv.less, prot.qv.less = prot.qv.less,
              rna.pv.greater = rna.pv.greater, rna.qv.greater = rna.qv.greater,
              prot.pv.greater = prot.pv.greater, prot.qv.greater = prot.qv.greater
  )
  return(df)
})
setwd(ourDir)
write.xlsx(geneExpDelta, "gene_RNA_prot_ttest.xlsx")

#####################################################################################
# CORRELATION WITH ALL GENES
ggmat <- c()
mygenes <- ann.sub$Symbol

# get coorelation in "NORMAL"
corList.ctl <- sapply(1:nrow(comb.ctl), function(j){
  cor(rna.ctl.sub[mygenes, comb.ctl[j, 1]], prot.ctl.sub[mygenes, comb.ctl[j, 2]], method = "spearman")
})

# get coorelation in "TUMOR"
corList.tcga <- sapply(1:nrow(comb.tcga), function(j){
  cor(rna.tcga.sub[mygenes, comb.tcga[j, 1]], prot.tcga.sub[mygenes, comb.tcga[j, 2]], method = "spearman")
})

pv <- t.test(corList.ctl, corList.tcga)$"p.value"
delta <- median(corList.ctl) - median(corList.tcga)
print(delta)

# fill in ggmat
ggmat <- rbind(ggmat, data.frame(cor = c(corList.ctl, corList.tcga),
                                 stabGroup = "ALL", 
                                 Condition = c(rep("normal", length(corList.ctl)), rep("tumor", length(corList.tcga)))
))

p <- ggplot(ggmat[ggmat$stabGroup != "grey.zone", ], aes(x = stabGroup, y = cor, fill = Condition)) +
  geom_point(pch = 21, position = position_jitterdodge(), alpha = 0.5) + 
  geom_boxplot(outlier.size = 0, alpha = 0.1)
p <- p + theme_bw(base_size = 16)
p <- p +scale_fill_manual(values=c("normal" = "grey80", "tumor" = "firebrick"))
p <- p + xlab("") + ylab("RNA / Protein spearman's correlation")

#ggsave(plot = p, filename = "ALL_RNA_prot_correlation_jitter.pdf", width = 4.5, height = 7)
ggsave(plot = p, filename = "SupFigure5C.pdf", width = 4.5, height = 7)


toxlsx <- ggmat[ggmat$stabGroup != "grey.zone", ]
write.xlsx(toxlsx, "SupFigure5C_data.xlsx")

x <- ggmat$cor[ggmat$Condition == "normal"]
y <- ggmat$cor[ggmat$Condition == "tumor"]

shapi <- data.frame(shapiro.pvalue = c(shapiro.test(x)$p.value, shapiro.test(y)$p.value))
rownames(shapi) <- c("normal", "tumor")
write.xlsx(shapi, "SupFigure5C_Shapiro_test.xlsx", row.names = TRUE)

#t.test(x, y, alternative = "less")
wilcox.res <- data.frame(wilcox.pvalue = wilcox.test(x, y, alternative = "two.sided")$p.value)
rownames(wilcox.res) <- "ALL"

write.xlsx(wilcox.res, "SupFigure5C_wilcox_test.xlsx", row.names = TRUE)


#############################################

stabCor <- lapply(stabGroup, function(i){
  mygenes <- ann.sub$Symbol[ann.sub$Group == i]
  
  # get coorelation in "NORMAL"
  corList.ctl <- sapply(1:nrow(comb.ctl), function(j){
    cor(rna.ctl.sub[mygenes, comb.ctl[j, 1]], prot.ctl.sub[mygenes, comb.ctl[j, 2]], method = "spearman")
  })
  
  # get coorelation in "TUMOR"
  corList.tcga <- sapply(1:nrow(comb.tcga), function(j){
    cor(rna.tcga.sub[mygenes, comb.tcga[j, 1]], prot.tcga.sub[mygenes, comb.tcga[j, 2]], method = "spearman")
  })
  
  return(c(median(corList.ctl), median(corList.tcga)))
})

stabCor <- matrix(unlist(stabCor), ncol = 2, byrow = TRUE)
rownames(stabCor) <- stabGroup
colnames(stabCor) <- c("normal", "tumor")
stabCorList <- lapply(1:nrow(stabCor), function(i){
  stabCor.diff <- (stabCor[,1] - stabCor[i,1]) - (stabCor[,2] - stabCor[i,2])
  return(stabCor.diff)
})
stabCorMat <- do.call(rbind, stabCorList)
rownames(stabCorMat) <- stabGroup

nbRd <- 1000
corList <- mclapply(1:nbRd, function(k){
  
  stabCor.rd <- lapply(stabGroup, function(i){
    
    mygenes <- ann.sub$Symbol[ann.sub$Group == i]
    idx <- sample(rownames(rna.ctl.sub), length(mygenes))
    
    cor.ctl.current <- sapply(1:nrow(comb.ctl), function(j){
      cor(rna.ctl.sub[idx, comb.ctl[j, 1]], prot.ctl.sub[idx, comb.ctl[j, 2]], method = "spearman")
    })
    
    cor.tcga.current <- sapply(1:nrow(comb.tcga), function(j){
      cor(rna.tcga.sub[idx, comb.tcga[j, 1]], prot.tcga.sub[idx, comb.tcga[j, 2]], method = "spearman")
    })
    
    #return(t.test(cor.ctl.current, cor.tcga.current))
    return(c(median(cor.ctl.current), median(cor.tcga.current)))
  })
  stabCor.rd <- matrix(unlist(stabCor.rd), ncol = 2, byrow = TRUE)
  return(stabCor.rd)
  
}, mc.cores = 8)

# vs. rS.pS 

pvList <- lapply(1:nrow(stabCorMat), function(k){
  
  corList.diff <- lapply(corList, function(x){
    d.ctl <- x[,1] - x[k,1]
    d.tcga <- x[,2] - x[k,2]
    return(d.ctl - d.tcga)
  })
  corList.diff <- matrix(unlist(corList.diff), ncol = 5, byrow = TRUE)
  
  return(sapply(1:ncol(stabCorMat), function(j){
    getEmpiricalPV(corList.diff[, j], stabCorMat[k, j], "greater")
  }))
  
})
pvMat <- do.call(rbind, pvList)
rownames(pvMat) <- stabGroup
colnames(pvMat) <- stabGroup
pvMat[pvMat==0] <- 1/nbRd

for(i in 1:ncol(pvMat)){
  pvMat[i,i] <- NA
}

setwd(file.path("~/Research/Sajib/P2_half_lives/TCGA/"))
write.xlsx(list(delta = stabCorMat, pv = pvMat), "stabGroup_normal_tumor_delta_cor.xlsx", row.names = TRUE)

# "CORPLOT"
setwd(file.path("~/Research/Sajib/P2_half_lives/TCGA/save_median/"))
stabCorMat <- read.xlsx("stabGroup_normal_tumor_delta_cor.xlsx", sheet = "delta", rowNames = TRUE)
pvMat <- read.xlsx("stabGroup_normal_tumor_delta_cor.xlsx", sheet = "pv", rowNames = TRUE)


keep <- c("rS.pS", "rU.pS", "rS.pU", "rU.pU")
stabCorMat <- stabCorMat[keep, keep]
pvMat <- pvMat[keep, keep]

stabCorMat <- as.matrix(stabCorMat)
pvMat <- as.matrix(pvMat)

col<- colorRampPalette(c("navy", "snow", "firebrick"))(10)

pdf("SupFigure5E.pdf", width = 4, height = 4)
corrplot(stabCorMat,
         #type="upper",
         col=col,
         diag=FALSE ,
         addCoef.col = "black",
         tl.col="black", tl.srt=45,
         cl.lim = c(-0.11, 0.11),
         is.corr = FALSE)
dev.off()

write.xlsx(list(delta= stabCorMat, pv = pvMat), "SupFigure5E_data.xlsx", row.names = FALSE)

