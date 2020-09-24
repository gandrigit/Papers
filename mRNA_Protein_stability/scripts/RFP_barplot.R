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


############################################
###                                      ###
###               FUNCTION               ###
###                                      ###
############################################

sem <- function(x) return(sd(x, na.rm=TRUE)/sqrt(length(x[!is.na(x)])))

data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      #sd = sd(x[[col]], na.rm=TRUE),
      sem = sem(x[[col]])
      #median = median(x[[col]], na.rm=TRUE),
      #mad = mad(x[[col]], na.rm=TRUE)
    )
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
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


############################################
###                                      ###
###                 MAIN                 ###
###                                      ###
############################################


outDir <- file.path("../output")
dir.create(outDir)
setwd(outDir)

##########################
# MOUSE
setwd(file.path("../data"))
rfpMat <- read.xlsx("RFP_data.xlsx", sheet = "RFP Mouse cells")

setwd(outDir)
stabGroup <- rfpMat$Group

ggmat <- data.frame(GENE = rfpMat$Symbol,
                    GROUP = stabGroup,
                    RFP = rfpMat$"Avg.ribosomal.foot.printing")
ggmat$GROUP <- factor(ggmat$GROUP, levels = c("rS.pS", "rU.pS", "rS.pU", "rU.pU"))

# TEST NORMALITY
shapi <- lapply(unique(stabGroup), function(i){
  x <- ggmat$RFP[ggmat$GROUP == i]
  shapiro.test(x)
})
shapi <- sapply(shapi, function(i) i$p.value)
shapi <- data.frame(Sahpiro.pvalue = shapi)
rownames(shapi) <- unique(stabGroup)

shapi$FDR <- p.adjust(shapi$Sahpiro.pvalue, method = "BH")
write.xlsx(shapi, "figure1B_left_Shapiro_test.xlsx", row.names = TRUE)

# WILCOXON TEST
combos <- combn(unique(ggmat$GROUP),2)
wilcox.res <- adply(combos, 2, function(x) {
  test <- wilcox.test(ggmat$RFP[ggmat$GROUP == as.character(x[1])],
                      ggmat$RFP[ggmat$GROUP == as.character(x[2])])
  
  out <- data.frame("var1" = x[1]
                    , "var2" = x[2]
                    , "var1.mean" = mean(ggmat$RFP[ggmat$GROUP == as.character(x[1])])
                    , "var2.mean" = mean(ggmat$RFP[ggmat$GROUP == as.character(x[2])])
                    , "p.value" = sprintf("%.3f", test$p.value)
  )
  return(out)
  
})
wilcox.res$X1 <- unique(ggmat$name)
wilcox.res$"p.value" <- as.numeric(levels(wilcox.res$"p.value"))[wilcox.res$"p.value"]
wilcox.res$FDR <- p.adjust(wilcox.res$"p.value", method = "BH")
write.xlsx(wilcox.res, "figure1B_left_wilcox_test.xlsx", row.names = FALSE)



# PLOT
mycol <- c(rgb(237, 189, 79, maxColorValue = 255),
           rgb(230, 121, 123, maxColorValue = 255),
           rgb(105, 153, 196, maxColorValue = 255),
           rgb(105, 48, 48, maxColorValue = 255))



ggMat.FP <- data_summary(ggmat, varname="RFP", 
                         groupnames=c("GROUP"))
ggMat.FP$GROUP <- factor(ggMat.FP$GROUP, levels = c("rS.pS", "rU.pS", "rS.pU", "rU.pU"))


p <- ggplot(ggMat.FP, aes(x=GROUP, y=RFP, fill=GROUP))
p <- p + theme_bw(base_size = 16) + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                          panel.grid.minor = element_blank())
p <- p + geom_bar(stat="identity", color="black", position=position_dodge())
p <- p + geom_errorbar(aes(ymin=RFP-sem, ymax=RFP+sem), width=.2, position=position_dodge(.9))                     
p <- p + theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
p <- p + scale_fill_manual(values=mycol) 
p <- p + theme(legend.position="none")  
#p <- p + theme(axis.text.x = element_text(size=14),axis.text.y = element_text(size=14))
p <- p + theme(axis.line.x = element_line(color="black", size = 0.5),
               axis.line.y = element_line(color="black", size = 0.5))	       

setwd(outDir)
pdf("RFP_Murine_barplot.pdf", width = 7, height = 8)
plot(p)
dev.off()

write.xlsx(ggMat.FP, "figure1B_left_data.xlsx")


############################################################################
# HUMAN
setwd(file.path("../data"))
rfpMat <- read.xlsx("RFP_data.xlsx", sheet = "RFP HEK-293 cell")
stabGroup <- rfpMat$Group

ggmat <- data.frame(GENE = rfpMat$Gene,
                    GROUP = stabGroup,
                    RFP = rfpMat$"RFP.Avg")
ggmat$GROUP <- factor(ggmat$GROUP, levels = c("rS.pS", "rU.pS", "rS.pU", "rU.pU"))

# TEST NORMALITY
shapi <- lapply(unique(stabGroup), function(i){
  x <- ggmat$RFP[ggmat$GROUP == i]
  shapiro.test(x)
})
shapi <- sapply(shapi, function(i) i$p.value)
shapi <- data.frame(Sahpiro.pvalue = shapi)
rownames(shapi) <- unique(stabGroup)
shapi$FDR <- p.adjust(shapi$Sahpiro.pvalue, method = "BH")
setwd(outDir)
write.xlsx(shapi, "figure1B_right_Shapiro_test.xlsx", row.names = TRUE)

# WILCOXON TEST
combos <- combn(unique(ggmat$GROUP),2)
wilcox.res <- adply(combos, 2, function(x) {
  test <- wilcox.test(ggmat$RFP[ggmat$GROUP == as.character(x[1])],
                      ggmat$RFP[ggmat$GROUP == as.character(x[2])])
  
  out <- data.frame("var1" = x[1]
                    , "var2" = x[2]
                    , "var1.mean" = mean(ggmat$RFP[ggmat$GROUP == as.character(x[1])])
                    , "var2.mean" = mean(ggmat$RFP[ggmat$GROUP == as.character(x[2])])
                    , "p.value" = sprintf("%.3f", test$p.value)
  )
  return(out)
  
})
wilcox.res$X1 <- unique(ggmat$name)
wilcox.res$"p.value" <- as.numeric(levels(wilcox.res$"p.value"))[wilcox.res$"p.value"]
wilcox.res$FDR <- p.adjust(wilcox.res$"p.value", method = "BH")
setwd(outDir)
write.xlsx(wilcox.res, "figure1B_right_wilcox_test.xlsx", row.names = FALSE)



# PLOT
ggMat.FP <- data_summary(ggmat, varname="RFP", 
                         groupnames=c("GROUP"))
ggMat.FP$GROUP <- factor(ggMat.FP$GROUP, levels = c("rS.pS", "rU.pS", "rS.pU", "rU.pU"))


p <- ggplot(ggMat.FP, aes(x=GROUP, y=RFP, fill=GROUP))
p <- p + theme_bw(base_size = 16) + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                          panel.grid.minor = element_blank())
p <- p + geom_bar(stat="identity", color="black", position=position_dodge())
p <- p + geom_errorbar(aes(ymin=RFP-sem, ymax=RFP+sem), width=.2, position=position_dodge(.9))                     
p <- p + theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
p <- p + scale_fill_manual(values=mycol) 
p <- p + theme(legend.position="none")  
#p <- p + theme(axis.text.x = element_text(size=14),axis.text.y = element_text(size=14))
p <- p + theme(axis.line.x = element_line(color="black", size = 0.5),
               axis.line.y = element_line(color="black", size = 0.5))	       

setwd(outDir)
pdf("RFP_Human_barplot.pdf", width = 7, height = 8)
plot(p)
dev.off()

write.xlsx(ggMat.FP, "figure1B_right_data.xlsx")




