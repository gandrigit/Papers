#options(java.parameters = "-Xmx8g")

library(plyr)
library(ggplot2)
library(openxlsx)
library(GeneAnswers)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(pheatmap)
library(reshape2)

source("tools.r")

############################################
###                                      ###
###               FUNCTION               ###
###                                      ###
############################################

toNum <- function(x){
  if(is.numeric(x)) return(x)
  return(as.numeric(levels(x))[x])
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
  return(data.frame(Sample = rep(name, length(groups)),
                    Group = names(groups),
                    Correlation = corre,
                    PV.greater = pv.greater,
                    PV.less = pv.less))
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
  
  # loop groups
  for(g in groups)
  {
    idx <- which(genes %in% g)
    cor1.true <- cor(toNum(x1[idx]), toNum(y1[idx]),
                     method = "spearman", use = "complete.obs")
    cor2.true <- cor(toNum(x2[idx]), toNum(y2[idx]),
                     method = "spearman", use = "complete.obs")
    delta.true <- cor1.true - cor2.true
    
    delta.random <- getRandomDelta(length(idx), 10000, toNum(x1), toNum(y1), toNum(x2), toNum(y2))	
    
    cor1 <- c(cor1, cor1.true)
    cor2 <- c(cor2, cor2.true)
    delta <- c(delta, delta.true)
    delta.rd <- c(delta.rd, median(delta.random))
    
    pv.greater <- c(pv.greater, getEmpiricalPV(delta.random, delta.true, "greater"))
    pv.less <- c(pv.less, getEmpiricalPV(delta.random, delta.true, "less"))	
  }
  
  qv.greater <- p.adjust(pv.greater, method = "BH")
  qv.less <- p.adjust(pv.less, method = "BH")
  return(data.frame(Sample = rep(name, length(groups)),
                    Group = names(groups),
                    Correlation1 = cor1,
                    Correlation2 = cor2,
                    Delta = delta, 
                    Delta.random = delta.rd, 
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
setwd(file.path("../data/"))
ann <- read.xlsx("stability_group_human.xlsx", sheet = 1)

groups.unique <- unique(ann$Group)
ann.list <- list()
for(i in 1:length(groups.unique))
{
  ann.list[[i]] <- ann$Symbol[ann$Group == groups.unique[i]]
}
names(ann.list) <- groups.unique


####################################
# LOAD HUMAN TISSUE DATA (M.Selbach)
rna1 <- read.xlsx("rna_DS1.xlsx", sheet = 1, rowNames = TRUE)
prot1 <- read.xlsx("protein_DS1.xlsx", sheet = 1, rowNames = TRUE)

colnames(rna1) <- toupper(colnames(rna1))
colnames(prot1) <- toupper(colnames(prot1))


# LOAD CELL LINE DATA
rna2 <- read.xlsx("rna_cell_lines.xlsx", sheet = 1, rowNames = TRUE)
prot2 <- read.xlsx("protein_cell_lines.xlsx", sheet = 1, rowNames = TRUE)

colnames(rna2) <- toupper(colnames(rna2))
colnames(prot2) <- toupper(colnames(prot2))


# INTERSECT THE TWO DATASETS
commonGenes <- intersect(rownames(rna1), rownames(rna2))

rna1 <- rna1[commonGenes, ]
prot1 <- prot1[commonGenes, ]

rna2 <- rna2[commonGenes, ]
prot2 <- prot2[commonGenes, ]

ann.sub <- ann[match(commonGenes, ann$Symbol), ]


# Pairs
pairMat <- matrix(c("A549", "HEPG2", "RKO", "LNCAP", "GAMG",
                    "LUNG", "LIVER", "COLON", "PROSTATE", "FRONTAL.CORTEX"),
                  ncol = 2, byrow= FALSE)
pairMat.names <- apply(pairMat, 1, paste, collapse = "-")

# DELTA CORRELATION PVALUE
outDir <- file.path("../output")
dir.create(outDir)
setwd(outDir)

for(i in 1:nrow(pairMat)){
  pairName <- paste(pairMat[i, ], collapse = "-")
  
  # Tissue
  rna1.current <- as.numeric(rna1[,pairMat[i, 2]])
  prot1.current <- as.numeric(prot1[,pairMat[i, 2]])
  
  # Cell line
  rna2.current <- as.numeric(rna2[,pairMat[i, 1]])
  prot2.current <- as.numeric(prot2[,pairMat[i, 1]])
  
  genes.current <- ann.sub$Symbol

  ann.current <- ann[match(genes.current, ann$Symbol), ]
  
  groups <- unique(ann.current$Group)
  
  
  # Delta Bootstrap PV
  groupList <- lapply(groups, function(g) ann.current$Symbol[ann.current$Group == g])
  names(groupList) <- groups
  
  deltaPV <- deltaPipeline(x1 = rna1.current, y1 = prot1.current,
                           x2 = rna2.current, y2 = prot2.current, 
                           genes = genes.current,
                           groups = groupList,
                           name = pairName)
  
  write.xlsx(deltaPV, paste0(pairName, "_delta_correlation_bootstrap_PV.xlsx"), row.names = FALSE)  
  
}

##################################################################
# Cell line vs. tissue correlation

# GET RNA / PROTEIN CORRELATION PER GROUP IN TISSUE AND CELL LINES

stabGroup <- unique(ann.sub$Group)
stabCor <- lapply(stabGroup, function(i){
  mygenes <- ann.sub$Symbol[ann.sub$Group == i]
  
  # tissue
  cor.tissue <- sapply(1:nrow(pairMat), function(j){
    ct <- cor.test(as.numeric(rna1[mygenes, pairMat[j, 2]]),
                   as.numeric(prot1[mygenes, pairMat[j, 2]]),
                   method = "spearman", use = "complete.obs")
    return(ct$estimate["rho"])
  })
  
  # cell lines
  cor.cell <- sapply(1:nrow(pairMat), function(j){
    ct <- cor.test(as.numeric(rna2[mygenes, pairMat[j, 1]]),
                   as.numeric(prot2[mygenes, pairMat[j, 1]]),
                   method = "spearman", use = "complete.obs")
    return(ct$estimate["rho"])
  })
  
  return(data.frame(TISSUE = cor.tissue, CELL = cor.cell, GROUP = i, PAIR = pairMat.names))
  
  ct <- cor.test(cor.tissue, cor.cell)
  return(c(ct$estimate["cor"], ct$p.value))
})

stabCor <- do.call(rbind, stabCor)

write.xlsx(stabCor, "TISSUE_vs_CELLLINES_correlation.xlsx")


mycol <- c(rgb(237, 189, 79, maxColorValue = 255),
           rgb(230, 121, 123, maxColorValue = 255),
           rgb(105, 153, 196, maxColorValue = 255),
           rgb(105, 48, 48, maxColorValue = 255))

p <- ggplot(stabCor[stabCor$GROUP != "grey.zone", ], aes(TISSUE, CELL)) +
  geom_point(aes(color = GROUP, shape = PAIR), size = 5) +
  stat_smooth(method = "lm", colour = "black")
p <- p + scale_color_manual( values = mycol)
p <- p + theme_bw(base_size = 14)
#ggsave(plot = p, filename = "TISSUE_vs_CELLLINES_overall_scatterplot.pdf", width = 7, height = 5)
ggsave(plot = p, filename = "figure3C.pdf", width = 7, height = 5)

toxlsx <- stabCor[stabCor$GROUP != "grey.zone", ]
write.xlsx(toxlsx, "figure3C_data.xlsx")

# DELTA BARPLOT UPDATE
stabCor$Delta <- stabCor$TISSUE - stabCor$CELL
stabCor$GROUP <- factor(stabCor$GROUP, levels = c("rS.pS", "rU.pS", "rS.pU", "rU.pU", "grey.zone"))

p <- ggplot(data=stabCor[stabCor$GROUP != "grey.zone", ], aes(x=PAIR, y=Delta, fill=GROUP)) +
  geom_bar(stat="identity", position=position_dodge())
p <- p + theme_bw(base_size = 16)
p <- p + theme(axis.line.x = element_line(color="black", size = 0.5),
               axis.line.y = element_line(color="black", size = 0.5))	
p <- p + xlab("") + ylab("Delta (tissue vs. cell line)")
p <- p + scale_fill_manual( values = mycol)
#ggsave(plot = p, filename = "delta_barplot.pdf", width = 8, height = 4)
ggsave(plot = p, filename = "figure3B.pdf", width = 8, height = 4)


toxlsx <- stabCor[stabCor$GROUP != "grey.zone", ]
write.xlsx(toxlsx, "figure3B_data.xlsx")


