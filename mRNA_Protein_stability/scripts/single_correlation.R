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
library(tidyverse)

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
    cor.random <- getRandomCor(length(idx), 1000, toNum(x), toNum(y))	
    
    corre <- c(corre, cor.true)
    pv.greater <- c(pv.greater, getEmpiricalPV(cor.random, cor.true, "greater"))
    pv.less <- c(pv.less, getEmpiricalPV(cor.random, cor.true, "less"))	
  }
  return(data.frame(Sample = rep(name, length(groups)),
                    Group = names(groups),
                    Correlation = corre,
                    PV.greater = pv.greater,
                    PV.less = pv.less))
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
setwd(file.path("../data/"))
rna <- read.xlsx("rna_DS1.xlsx", sheet = 1, rowNames = TRUE)
prot <- read.xlsx("protein_DS1.xlsx", sheet = 1, rowNames = TRUE)

ann.sub <- ann[match(rownames(rna), ann$Symbol), ]

# SINGLE TISSUE CORRELATION ANALYSIS
outDir <- file.path("../output")
dir.create(outDir)
setwd(outDir)

tissueNames <- colnames(rna)
corMat <- matrix(NA, nrow = length(tissueNames), ncol = 5)
rownames(corMat) <- tissueNames
colnames(corMat) <- c("rS.pS", "rU.pS", "rS.pU", "rU.pU", "grey.zone")

for(i in tissueNames){
  rna.current <- as.numeric(rna[,i])
  prot.current <- as.numeric(prot[,i])
  genes.current <- ann.sub$Symbol
  
  # remove NA
  rna.NA <- is.na(rna.current)
  prot.NA <- is.na(prot.current)
  tokeep <- !(rna.NA | prot.NA)
  
  rna.current <- rna.current[tokeep]
  prot.current <- prot.current[tokeep]
  genes.current <- genes.current[tokeep]
  
  ann.current <- ann[match(genes.current, ann$Symbol), ]
  
  groups <- unique(ann.current$Group)
  
  # Fill in corMat
  corList <- sapply(groups, function(g){
    cor(rna.current[ann.current$Group == g],
        prot.current[ann.current$Group == g],
        method = "spearman")
  })
  names(corList) <- groups
  
  corMat[i, names(corList)] <- corList
  
  
  # Bootstrap PV
  groupList <- lapply(groups, function(g) ann.current$Symbol[ann.current$Group == g])
  names(groupList) <- groups
  corPV <- correlationPipeline(x = rna.current, y = prot.current,
                               genes = genes.current,
                               groups = groupList, name = i)
  write.xlsx(corPV, paste0(i, "_correlation_bootstrap_PV.xlsx"), row.names = FALSE)  
  
  
  # Plot
  ggList <- lapply(groups, function(g){
    ggmat <- data.frame(X = rna.current[ann.current$Group == g],
                        Y = prot.current[ann.current$Group == g],
                        Group = g
    )
    ggmat$Density <- getDensity(ggmat$X, ggmat$Y)
    return(ggmat)
  })
  ggmat <- do.call(rbind, ggList)
  ggmat$Group <- factor(ggmat$Group, levels = c("rS.pS", "rU.pS", "rS.pU", "rU.pU", "grey.zone"))
  
  densityScatterplot(ggmat, outFile = paste0(i, "_correlation.pdf"))
  
  # without grey zone
  ggmat <- do.call(rbind, ggList)
  ggmat <- ggmat[ggmat$Group != "grey.zone",]
  
  densityScatterplot(ggmat, outFile = paste0(i, "_correlation_woGZ.pdf"))
  
}
write.xlsx(corMat, "DS1_correlation.xlsx", row.names = TRUE)

pheatmap(corMat, cluster_cols = FALSE, color = colorRampPalette(c("navy", "snow", "firebrick"))(16), border_color = NA,
         cellwidth = 10, cellheight = 10,
         filename = "DS1_correlation_heatmap.pdf")

corMat <- corMat[, -ncol(corMat)]
pheatmap(corMat, cluster_cols = FALSE, color = magma(16),
         border_color = NA,
         cellwidth = 10, cellheight = 10,
         filename = "DS1_correlation_heatmap_woGZ.pdf")

write.xlsx(corMat, "figure2A_data.xlsx", row.names = TRUE)

# TEST NORMALITY
shapi <- apply(corMat, 2, shapiro.test)
shapi <- sapply(shapi, function(i) i$p.value)
shapi <- data.frame(shapiro.pvalue = shapi)
shapi$FDR <- p.adjust(shapi$shapiro.pvalue, method = "BH")
write.xlsx(shapi, "figure2A_Shapiro_test.xlsx", row.names = TRUE)


# ANOVA
anova.input <- tidyr::gather(data.frame(corMat))
anova.res <- aov(value ~ key, anova.input)# analyse de variance

sink("DS1_correlation_anova_woGZ.txt")
summary(anova.res)
sink()

# T.TEST
combos <- combn(ncol(corMat),2)
ttest.res <- adply(combos, 2, function(x) {
  test <- t.test(corMat[, x[1]], corMat[, x[2]])
  
  out <- data.frame("var1" = colnames(corMat)[x[1]]
                    , "var2" = colnames(corMat)[x[2]]
                    , "var1.mean" = mean(corMat[, x[1]])
                    , "var2.mean" = mean(corMat[, x[2]])
                    , "t.value" = sprintf("%.3f", test$statistic)
                    ,  "df"= test$parameter
                    ,  "p.value" = sprintf("%.3f", test$p.value)
  )
  return(out)
})
ttest.res$p.value <- toNum(ttest.res$p.value)
ttest.res$FDR <- p.adjust(ttest.res$p.value, method = "BH")
ttest.res <- ttest.res[, -1]
#write.xlsx(ttest.res, "DS1_correlation_ttest_woGZ.xlsx", row.names = FALSE)
write.xlsx(ttest.res, "figure2A_ttest.xlsx", row.names = FALSE)

####################################
# LOAD HUMAN TISSUE DATA (DS2)
setwd(file.path("~/Research/Sajib/P2_half_lives/data"))
rna <- read.xlsx("rna_DS2.xlsx", sheet = 1, rowNames = TRUE)
prot <- read.xlsx("protein_DS2.xlsx", sheet = 1, rowNames = TRUE)

ann.sub <- ann[match(rownames(rna), ann$Symbol), ]

# SINGLE TISSUE CORRELATION ANALYSIS
setwd(outDir)
tissueNames <- colnames(rna)
corMat <- matrix(NA, nrow = length(tissueNames), ncol = 5)
rownames(corMat) <- tissueNames
colnames(corMat) <- c("rS.pS", "rU.pS", "rS.pU", "rU.pU", "grey.zone")

for(i in tissueNames){
  rna.current <- as.numeric(rna[,i])
  prot.current <- as.numeric(prot[,i])
  genes.current <- ann.sub$Symbol
  
  # remove NA
  rna.NA <- is.na(rna.current)
  prot.NA <- is.na(prot.current)
  tokeep <- !(rna.NA | prot.NA)
  
  rna.current <- rna.current[tokeep]
  prot.current <- prot.current[tokeep]
  genes.current <- genes.current[tokeep]
  
  ann.current <- ann[match(genes.current, ann$Symbol), ]
  
  groups <- unique(ann.current$Group)
  
  # Fill in corMat
  corList <- sapply(groups, function(g){
    cor(rna.current[ann.current$Group == g],
        prot.current[ann.current$Group == g],
        method = "spearman")
  })
  names(corList) <- groups
  
  corMat[i, names(corList)] <- corList
  
  
  # Bootstrap PV
  groupList <- lapply(groups, function(g) ann.current$Symbol[ann.current$Group == g])
  names(groupList) <- groups
  corPV <- correlationPipeline(x = rna.current, y = prot.current,
                               genes = genes.current,
                               groups = groupList, name = i)
  write.xlsx(corPV, paste0(i, "_correlation_bootstrap_PV.xlsx"), row.names = FALSE)  
  
  # Plot
  ggList <- lapply(groups, function(g){
    ggmat <- data.frame(X = rna.current[ann.current$Group == g],
                        Y = prot.current[ann.current$Group == g],
                        Group = g
    )
    ggmat$Density <- getDensity(ggmat$X, ggmat$Y)
    return(ggmat)
  })
  ggmat <- do.call(rbind, ggList)
  ggmat$Group <- factor(ggmat$Group, levels = c("rS.pS", "rU.pS", "rS.pU", "rU.pU", "grey.zone"))
  
  densityScatterplot(ggmat, outFile = paste0(i, "_correlation.pdf"))
  
  # without grey zone
  ggmat <- do.call(rbind, ggList)
  ggmat <- ggmat[ggmat$Group != "grey.zone",]
  
  densityScatterplot(ggmat, outFile = paste0(i, "_correlation_woGZ.pdf"))
  
}

write.xlsx(corMat, "DS2_correlation.xlsx", row.names = TRUE)

pheatmap(corMat, cluster_cols = FALSE, color = colorRampPalette(c("navy", "snow", "firebrick"))(16), border_color = NA,
         cellwidth = 10, cellheight = 10,
         filename = "DS2_correlation_heatmap.pdf")

corMat <- corMat[, -ncol(corMat)]
pheatmap(corMat, cluster_cols = FALSE, color = magma(16), border_color = NA,
         cellwidth = 10, cellheight = 10,
         filename = "DS2_correlation_heatmap_woGZ.pdf")

write.xlsx(corMat, "figure2B_data.xlsx", row.names = TRUE)

# TEST NORMALITY
shapi <- apply(corMat, 2, shapiro.test)
shapi <- sapply(shapi, function(i) i$p.value)
shapi <- data.frame(shapiro.pvalue = shapi)
shapi$FDR <- p.adjust(shapi$shapiro.pvalue, method = "BH")
write.xlsx(shapi, "figure2B_Shapiro_test.xlsx", row.names = TRUE)

# ANOVA
anova.input <- tidyr::gather(data.frame(corMat))
anova.res <- aov(value ~ key, anova.input)# analyse de variance

sink("DS2_correlation_anova_woGZ.txt")
summary(anova.res)
sink()

# T.TEST
combos <- combn(ncol(corMat),2)
ttest.res <- adply(combos, 2, function(x) {
  test <- t.test(corMat[, x[1]], corMat[, x[2]])
  
  out <- data.frame("var1" = colnames(corMat)[x[1]]
                    , "var2" = colnames(corMat)[x[2]]
                    , "var1.mean" = mean(corMat[, x[1]])
                    , "var2.mean" = mean(corMat[, x[2]])
                    , "t.value" = sprintf("%.3f", test$statistic)
                    ,  "df"= test$parameter
                    ,  "p.value" = sprintf("%.3f", test$p.value)
  )
  return(out)
})
ttest.res$p.value <- toNum(ttest.res$p.value)
ttest.res$FDR <- p.adjust(ttest.res$p.value, method = "BH")
ttest.res <- ttest.res[, -1]
#write.xlsx(ttest.res, "DS2_correlation_ttest_woGZ.xlsx", row.names = FALSE)
write.xlsx(ttest.res, "figure2B_ttest.xlsx", row.names = FALSE)


####################################
# LOAD HUMAN TISSUE DATA (DS3)
setwd(file.path("~/Research/Sajib/P2_half_lives/data"))
rna <- read.xlsx("rna_DS3.xlsx", sheet = 1, rowNames = TRUE)
prot <- read.xlsx("protein_DS3.xlsx", sheet = 1, rowNames = TRUE)

ann.sub <- ann[match(rownames(rna), ann$Symbol), ]

# SINGLE TISSUE CORRELATION ANALYSIS
setwd(outDir)
tissueNames <- colnames(rna)
corMat <- matrix(NA, nrow = length(tissueNames), ncol = 5)
rownames(corMat) <- tissueNames
colnames(corMat) <- c("rS.pS", "rU.pS", "rS.pU", "rU.pU", "grey.zone")

for(i in tissueNames){
  rna.current <- as.numeric(rna[,i])
  prot.current <- as.numeric(prot[,i])
  genes.current <- ann.sub$Symbol
  
  # remove NA
  rna.NA <- is.na(rna.current)
  prot.NA <- is.na(prot.current)
  tokeep <- !(rna.NA | prot.NA)
  
  rna.current <- rna.current[tokeep]
  prot.current <- prot.current[tokeep]
  genes.current <- genes.current[tokeep]
  
  ann.current <- ann[match(genes.current, ann$Symbol), ]
  
  groups <- unique(ann.current$Group)
  
  # Fill in corMat
  corList <- sapply(groups, function(g){
    cor(rna.current[ann.current$Group == g],
        prot.current[ann.current$Group == g],
        method = "spearman")
  })
  names(corList) <- groups
  
  corMat[i, names(corList)] <- corList
  
  
  # Bootstrap PV
  groupList <- lapply(groups, function(g) ann.current$Symbol[ann.current$Group == g])
  names(groupList) <- groups
  corPV <- correlationPipeline(x = rna.current, y = prot.current,
                               genes = genes.current,
                               groups = groupList, name = i)
  write.xlsx(corPV, paste0(i, "_correlation_bootstrap_PV.xlsx"), row.names = FALSE)  
  
  
  
  # Plot
  ggList <- lapply(groups, function(g){
    ggmat <- data.frame(X = rna.current[ann.current$Group == g],
                        Y = prot.current[ann.current$Group == g],
                        Group = g
    )
    ggmat$Density <- getDensity(ggmat$X, ggmat$Y)
    return(ggmat)
  })
  ggmat <- do.call(rbind, ggList)
  ggmat$Group <- factor(ggmat$Group, levels = c("rS.pS", "rU.pS", "rS.pU", "rU.pU", "grey.zone"))
  
  densityScatterplot(ggmat, outFile = paste0(i, "_correlation.pdf"))
  
  # without grey zone
  ggmat <- do.call(rbind, ggList)
  ggmat <- ggmat[ggmat$Group != "grey.zone",]
  
  densityScatterplot(ggmat, outFile = paste0(i, "_correlation_woGZ.pdf"))
  
}

write.xlsx(corMat, "DS3_correlation.xlsx", row.names = TRUE)

pheatmap(corMat, cluster_cols = FALSE, color = colorRampPalette(c("navy", "snow", "firebrick"))(16), border_color = NA,
         cellwidth = 10, cellheight = 10,
         filename = "DS3_correlation_heatmap.pdf")

corMat <- corMat[, -ncol(corMat)]
pheatmap(corMat, cluster_cols = FALSE, color = magma(16), border_color = NA,
         cellwidth = 10, cellheight = 10,
         filename = "DS3_correlation_heatmap_woGZ.pdf")

write.xlsx(corMat, "figure2C_data.xlsx", row.names = TRUE)

# TEST NORMALITY
shapi <- apply(corMat, 2, shapiro.test)
shapi <- sapply(shapi, function(i) i$p.value)
shapi <- data.frame(shapiro.pvalue = shapi)
shapi$FDR <- p.adjust(shapi$shapiro.pvalue, method = "BH")
write.xlsx(shapi, "figure2C_Shapiro_test.xlsx", row.names = TRUE)

# ANOVA
anova.input <- tidyr::gather(data.frame(corMat))
anova.res <- aov(value ~ key, anova.input)# analyse de variance

sink("DS3_correlation_anova_woGZ.txt")
summary(anova.res)
sink()

# T.TEST
combos <- combn(ncol(corMat),2)
ttest.res <- adply(combos, 2, function(x) {
  test <- t.test(corMat[, x[1]], corMat[, x[2]])
  
  out <- data.frame("var1" = colnames(corMat)[x[1]]
                    , "var2" = colnames(corMat)[x[2]]
                    , "var1.mean" = mean(corMat[, x[1]])
                    , "var2.mean" = mean(corMat[, x[2]])
                    , "t.value" = sprintf("%.3f", test$statistic)
                    ,  "df"= test$parameter
                    ,  "p.value" = sprintf("%.3f", test$p.value)
  )
  return(out)
})
ttest.res$p.value <- toNum(ttest.res$p.value)
ttest.res$FDR <- p.adjust(ttest.res$p.value, method = "BH")
ttest.res <- ttest.res[, -1]
#write.xlsx(ttest.res, "DS3_correlation_ttest_woGZ.xlsx", row.names = FALSE)
write.xlsx(ttest.res, "figure2C_ttest.xlsx", row.names = FALSE)




####################################
# LOAD CELL LINE DATA
setwd(file.path("~/Research/Sajib/P2_half_lives/data"))
rna <- read.xlsx("rna_cell_lines.xlsx", sheet = 1, rowNames = TRUE)
prot <- read.xlsx("protein_cell_lines.xlsx", sheet = 1, rowNames = TRUE)

commonSamples <- intersect(colnames(rna), colnames(prot))

rna <- rna[, commonSamples]
prot <- prot[, commonSamples]

ann.sub <- ann[match(rownames(rna), ann$Symbol), ]

# SINGLE TISSUE CORRELATION ANALYSIS
setwd(outDir)

tissueNames <- colnames(rna)
corMat <- matrix(NA, nrow = length(tissueNames), ncol = 5)
rownames(corMat) <- tissueNames
colnames(corMat) <- c("rS.pS", "rU.pS", "rS.pU", "rU.pU", "grey.zone")

for(i in tissueNames){
  rna.current <- as.numeric(rna[,i])
  prot.current <- as.numeric(prot[,i])
  genes.current <- ann.sub$Symbol
  
  # remove NA
  rna.NA <- is.na(rna.current)
  prot.NA <- is.na(prot.current)
  tokeep <- !(rna.NA | prot.NA)
  
  rna.current <- rna.current[tokeep]
  prot.current <- prot.current[tokeep]
  genes.current <- genes.current[tokeep]
  
  ann.current <- ann[match(genes.current, ann$Symbol), ]
  
  groups <- unique(ann.current$Group)
  
  # Fill in corMat
  corList <- sapply(groups, function(g){
    cor(rna.current[ann.current$Group == g],
        prot.current[ann.current$Group == g],
        method = "spearman")
  })
  names(corList) <- groups
  
  corMat[i, names(corList)] <- corList
  
  
  # Bootstrap PV
  groupList <- lapply(groups, function(g) ann.current$Symbol[ann.current$Group == g])
  names(groupList) <- groups
  corPV <- correlationPipeline(x = rna.current, y = prot.current,
                               genes = genes.current,
                               groups = groupList, name = i)
  write.xlsx(corPV, paste0(i, "_correlation_bootstrap_PV.xlsx"), row.names = FALSE)  
  
  
  
  # Plot
  ggList <- lapply(groups, function(g){
    ggmat <- data.frame(X = rna.current[ann.current$Group == g],
                        Y = prot.current[ann.current$Group == g],
                        Group = g
    )
    ggmat$Density <- getDensity(ggmat$X, ggmat$Y)
    return(ggmat)
  })
  ggmat <- do.call(rbind, ggList)
  ggmat$Group <- factor(ggmat$Group, levels = c("rS.pS", "rU.pS", "rS.pU", "rU.pU", "grey.zone"))
  
  densityScatterplot(ggmat, outFile = paste0(i, "_correlation.pdf"))
  
  # without grey zone
  ggmat <- do.call(rbind, ggList)
  ggmat <- ggmat[ggmat$Group != "grey.zone",]
  
  densityScatterplot(ggmat, outFile = paste0(i, "_correlation_woGZ.pdf"))
  
}

write.xlsx(corMat, "cell_lines_correlation.xlsx", row.names = TRUE)


pheatmap(corMat, cluster_cols = FALSE, color = colorRampPalette(c("navy", "snow", "firebrick"))(16), border_color = NA,
         cellwidth = 10, cellheight = 10,
         filename = "cell_lines_correlation_heatmap.pdf")

corMat <- corMat[, -ncol(corMat)]
pheatmap(corMat, cluster_cols = FALSE, color = magma(16), border_color = NA,
         cellwidth = 10, cellheight = 10,
         filename = "cell_lines_correlation_heatmap_woGZ.pdf")

write.xlsx(corMat, "figure3A_data.xlsx", row.names = TRUE)

# TEST NORMALITY
shapi <- apply(corMat, 2, shapiro.test)
shapi <- sapply(shapi, function(i) i$p.value)
shapi <- data.frame(shapiro.pvalue = shapi)
shapi$FDR <- p.adjust(shapi$shapiro.pvalue, method = "BH")
write.xlsx(shapi, "figure3A_Shapiro_test.xlsx", row.names = TRUE)


# ANOVA
anova.input <- tidyr::gather(data.frame(corMat))
anova.res <- aov(value ~ key, anova.input)# analyse de variance

sink("cell_lines_correlation_anova_woGZ.txt")
summary(anova.res)
sink()

# T.TEST
combos <- combn(ncol(corMat),2)
ttest.res <- adply(combos, 2, function(x) {
  test <- t.test(corMat[, x[1]], corMat[, x[2]])
  
  out <- data.frame("var1" = colnames(corMat)[x[1]]
                    , "var2" = colnames(corMat)[x[2]]
                    , "var1.mean" = mean(corMat[, x[1]])
                    , "var2.mean" = mean(corMat[, x[2]])
                    , "t.value" = sprintf("%.3f", test$statistic)
                    ,  "df"= test$parameter
                    ,  "p.value" = sprintf("%.3f", test$p.value)
  )
  return(out)
})
ttest.res$p.value <- toNum(ttest.res$p.value)
ttest.res$FDR <- p.adjust(ttest.res$p.value, method = "BH")
ttest.res <- ttest.res[, -1]
#write.xlsx(ttest.res, "cell_lines_correlation_ttest_woGZ.xlsx", row.names = FALSE)
write.xlsx(ttest.res, "figure3A_ttest.xlsx", row.names = FALSE)





####################################
# LOAD TCGA COAD DATA
setwd(file.path("~/Research/Sajib/P2_half_lives/data"))
rna <- read.xlsx("rna_COAD.xlsx", sheet = 1, rowNames = TRUE)
prot <- read.xlsx("prot_COAD.xlsx", sheet = 1, rowNames = TRUE)

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
rna <- rna[commonGenes, ]
prot <- prot[commonGenes, ]

ann.sub <- ann[match(rownames(rna), ann$Symbol), ]


# Replace 0 by NA
prot[prot == 0] <- NA

# SINGLE TISSUE CORRELATION ANALYSIS
setwd(outDir)

tissueNames <- colnames(rna)
corMat <- matrix(NA, nrow = length(tissueNames), ncol = 5)
rownames(corMat) <- tissueNames
colnames(corMat) <- c("rS.pS", "rU.pS", "rS.pU", "rU.pU", "grey.zone")

for(i in tissueNames){
  rna.current <- as.numeric(rna[,i])
  prot.current <- as.numeric(prot[,i])
  genes.current <- ann.sub$Symbol
  
  # remove NA
  rna.NA <- is.na(rna.current)
  prot.NA <- is.na(prot.current)
  tokeep <- !(rna.NA | prot.NA)
  
  rna.current <- rna.current[tokeep]
  prot.current <- prot.current[tokeep]
  genes.current <- genes.current[tokeep]
  
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
  
  
  # Bootstrap PV
  groupList <- lapply(groups, function(g) ann.current$Symbol[ann.current$Group == g])
  names(groupList) <- groups
  corPV <- correlationPipeline(x = rna.current, y = prot.current,
                               genes = genes.current,
                               groups = groupList, name = i)
  write.xlsx(corPV, paste0(i, "_correlation_bootstrap_PV.xlsx"), row.names = FALSE)  
  
  
  
  # Plot
  ggList <- lapply(groups, function(g){
    ggmat <- data.frame(X = rna.current[ann.current$Group == g],
                        Y = prot.current[ann.current$Group == g],
                        Group = g
    )
    ggmat$Density <- getDensity(ggmat$X, ggmat$Y)
    return(ggmat)
  })
  ggmat <- do.call(rbind, ggList)
  ggmat$Group <- factor(ggmat$Group, levels = c("rS.pS", "rU.pS", "rS.pU", "rU.pU", "grey.zone"))
  
  densityScatterplot(ggmat, outFile = paste0(i, "_correlation.pdf"))
  
  # without grey zone
  ggmat <- do.call(rbind, ggList)
  ggmat <- ggmat[ggmat$Group != "grey.zone",]
  
  densityScatterplot(ggmat, outFile = paste0(i, "_correlation_woGZ.pdf"))
  
}

write.xlsx(corMat, "COAD_correlation.xlsx", row.names = TRUE)


pheatmap(corMat, cluster_cols = FALSE, color = colorRampPalette(c("navy", "snow", "firebrick"))(16), border_color = NA,
         cellwidth = 10, cellheight = 10,
         filename = "COAD_correlation_heatmap.pdf")

corMat <- corMat[, -ncol(corMat)]
pheatmap(corMat, cluster_cols = FALSE, color = magma(16), border_color = NA,
         cellwidth = 10, cellheight = 10,
         filename = "COAD_correlation_heatmap_woGZ.pdf")


write.xlsx(corMat, "SupFigure3B_data.xlsx", row.names = TRUE)

# TEST NORMALITY
shapi <- apply(corMat, 2, shapiro.test)
shapi <- sapply(shapi, function(i) i$p.value)
shapi <- data.frame(shapiro.pvalue = shapi)
shapi$FDR <- p.adjust(shapi$shapiro.pvalue, method = "BH")
write.xlsx(shapi, "SupFigure3B_Shapiro_test.xlsx", row.names = TRUE)



# ANOVA
anova.input <- tidyr::gather(data.frame(corMat))
anova.res <- aov(value ~ key, anova.input)# analyse de variance

sink("COAD_correlation_anova_woGZ.txt")
summary(anova.res)
sink()

# T.TEST
combos <- combn(ncol(corMat),2)
ttest.res <- adply(combos, 2, function(x) {
  test <- t.test(corMat[, x[1]], corMat[, x[2]])
  
  out <- data.frame("var1" = colnames(corMat)[x[1]]
                    , "var2" = colnames(corMat)[x[2]]
                    , "var1.mean" = mean(corMat[, x[1]])
                    , "var2.mean" = mean(corMat[, x[2]])
                    , "t.value" = sprintf("%.3f", test$statistic)
                    ,  "df"= test$parameter
                    ,  "p.value" = sprintf("%.3f", test$p.value)
  )
  return(out)
})
ttest.res$p.value <- toNum(ttest.res$p.value)
ttest.res$FDR <- p.adjust(ttest.res$p.value, method = "BH")
ttest.res <- ttest.res[, -1]
#write.xlsx(ttest.res, "COAD_correlation_ttest_woGZ.xlsx", row.names = FALSE)
write.xlsx(ttest.res, "SupFigure5B_ttest.xlsx", row.names = FALSE)


# WILCOXON TEST
combos <- combn(ncol(corMat),2)
wilcox.res <- adply(combos, 2, function(x) {
  test <- wilcox.test(corMat[, x[1]], corMat[, x[2]], exact = TRUE)
  
  out <- data.frame("var1" = colnames(corMat)[x[1]]
                    , "var2" = colnames(corMat)[x[2]]
                    , "var1.mean" = mean(corMat[, x[1]])
                    , "var2.mean" = mean(corMat[, x[2]])
                    , "p.value" = sprintf("%.3f", test$p.value)
  )
  return(out)
})
wilcox.res$p.value <- toNum(wilcox.res$p.value)
wilcox.res$FDR <- p.adjust(wilcox.res$p.value, method = "BH")
wilcox.res <- wilcox.res[, -1]
write.xlsx(wilcox.res, "SupFigure5B_wilcox_test.xlsx", row.names = FALSE)


# SUP FIGURE 5A 
ggmat <- data.frame(SAMPLE = rownames(corMat), corMat)
ggmat <- pivot_longer(ggmat, cols = c("rS.pS", "rU.pS", "rS.pU", "rU.pU"))
ggmat$name <- factor(ggmat$name, levels = c("rS.pS", "rU.pS", "rS.pU", "rU.pU"))

# Jitter plot
mycol <- c(rgb(237, 189, 79, maxColorValue = 255),
           rgb(230, 121, 123, maxColorValue = 255),
           rgb(105, 153, 196, maxColorValue = 255),
           rgb(105, 48, 48, maxColorValue = 255))

p <- ggplot(ggmat, aes(x = name, y = value, fill = name)) +
  geom_boxplot(outlier.size = 0, alpha = 0.8)
p <- p + theme_bw(base_size = 16)
p <- p + scale_fill_manual(values=mycol)
p <- p + xlab("") + ylab("RNA / Protein spearman's correlation")


ggsave(plot = p, filename = "SupFigure5A.pdf", width = 7.5, height = 7)

write.xlsx(ggmat, "SupFigure5A_data.xlsx")

