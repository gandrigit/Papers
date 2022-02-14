library(org.Mm.eg.db)
library(limma)
library(openxlsx)
library(ggplot2)
library(ggrepel)
library(data.table)
library(preprocessCore)
library(UpSetR)
library(ggplot2)
library(ggrepel)

source("~/Scripts/work/tools/tools.r")
source("~/Scripts/work/tools/venn_script.r")

############################################
###                                      ###
###               FUNCTION               ###
###                                      ###
############################################

my.read_xlsx <- function(inFile)
{
  mysheets <- getSheetNames(inFile)
  mList <- lapply(mysheets, read.xlsx, xlsxFile = inFile, rowNames = TRUE)
  names(mList) <- mysheets
  return(mList)
}

symbol2entrez <- function(symbol)
{
  entrez <- mget(as.character(symbol), org.Mm.egSYMBOL2EG, ifnotfound=NA)
  entrez <- unlist(lapply(entrez, function(i) return(i[1])))
  return(entrez)
}

entrez2symbol <- function(entrez)
{
  symbol <- mget(as.character(entrez), org.Mm.egSYMBOL, ifnotfound=NA)
  symbol <- unlist(lapply(symbol, function(i) return(i[1])))
  return(symbol)
}

ggBarPoint <- function(ggmat, output)
{
  p <- ggplot(ggmat)
  p <- p + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                              panel.grid.minor = element_blank())
  p <- p + theme(axis.line.x = element_line(color="black", size = 0.5),
                 axis.line.y = element_line(color="black", size = 0.5))	
  
  p <- p + geom_point(aes(x = Var1, y = value, color = Group, group = Group),
                      size=3, position = position_dodge(-0.9))# + coord_flip()
  p <- p + scale_color_manual(values = c(WT = "blue", KO = "red"))
  
  p <- p + stat_summary(aes(x = Var1, y = value, fill = Group, group = Group),
                        fun.y=mean, geom = "bar", position = position_dodge(-0.9), alpha = 0.5)
  p <- p + scale_fill_manual(values = c(WT = "blue", KO = "red"))	
  
  p <- p + coord_cartesian(ylim=c(25,38))
  #p <- p + scale_x_discrete(limits = sort(unique(ggmat$Var1), FALSE))
  #p <- p + theme(axis.text.x = element_text(size=20),
  #	axis.text.y = element_text(face = "italic", size=20))	
  
  
  pdf(output, width = 7.5 , height = 6.5)
  plot(p)
  dev.off()			
}


log10_log2 <- function(x) return(x/log10(2))

toNum <- function(x) as.numeric(levels(x))[x]


Intersect <- function (x) {  
  # Multiple set version of intersect
  # x is a list
  if (length(x) == 1) {
    unlist(x)
  } else if (length(x) == 2) {
    intersect(x[[1]], x[[2]])
  } else if (length(x) > 2){
    intersect(x[[1]], Intersect(x[-1]))
  }
}

Union <- function (x) {  
  # Multiple set version of union
  # x is a list
  if (length(x) == 1) {
    unlist(x)
  } else if (length(x) == 2) {
    union(x[[1]], x[[2]])
  } else if (length(x) > 2) {
    union(x[[1]], Union(x[-1]))
  }
}

Setdiff <- function (x, y) {
  # Remove the union of the y's from the common x's. 
  # x and y are lists of characters.
  xx <- Intersect(x)
  yy <- Union(y)
  setdiff(xx, yy)
}

limma2Venn <- function(inputList, outName, fc, pv, adjusted = TRUE, uad = FALSE)
{
  # "uad" means "UP and DOWN"
  
  # Get genes (based on entrez ID)
  if(uad){
    geneList.up <- lapply(inputList, getGenesLimma, fcCt = fc, pvCt = pv, dsign = "up", adjusted = adjusted)
    names(geneList.up) <- paste0(names(geneList.up), ".UP")
    geneList.down <- lapply(inputList, getGenesLimma, fcCt = fc, pvCt = pv, dsign = "down", adjusted = adjusted)
    names(geneList.down) <- paste0(names(geneList.down), ".DOWN")
    geneList <- c(geneList.up, geneList.down)
  }else{
    geneList <- lapply(inputList, getGenesLimma, fcCt = fc, pvCt = pv, dsign = "both", adjusted = adjusted)
  }
  
  # UpSetPlot
  #geneList <- geneList[rev(names(geneList))]
  #pdf(paste0(outName, "_pv", pv, "_fc", fc, "_UpSetR.pdf"), width = 10, onefile = FALSE)
  #print(upset(fromList(geneList), order.by = "freq", nsets = length(geneList),
  #            sets = names(geneList), keep.order = TRUE,
  #            mainbar.y.label = "Intersection", sets.x.label = "significant loci per sample"))
  #dev.off()	
  
  
  # Do Venn
  gene.venn <- doVenn(geneList)
  
  # Info sheet
  nb <- lapply(gene.venn, length)
  nb <- unlist(nb)
  gene.venn <- gene.venn[nb!=0]
  nb <- nb[nb!=0]
  
  nbG <- lapply(names(gene.venn), gregexpr, pattern = ";")
  nbG <- lapply(nbG, function(i) return(unlist(i)))
  nbG <- unlist(lapply(nbG, function(i) return(sum(i != -1))))
  nbG <- nbG + 1
  
  venn_info <- data.frame(Group = names(gene.venn),
                          Group.number = nbG,
                          Variable.number = nb)
  
  # Details sheet
  venn.details <- lapply(1:length(gene.venn), function(i)
    cbind(Group = names(gene.venn)[i], ensembl = gene.venn[[i]])
  )
  venn.details <- as.data.frame(do.call(rbind, venn.details))
  
  genes <- venn.details$ensembl
  idx <- match(genes, inputList[[1]]$ensembl)
  venn.details <- cbind(venn.details,
                        symbol = inputList[[1]]$symbol[idx],
                        entrez = inputList[[1]]$entrez[idx]
  )
  
  for(i in 1:length(inputList))
  {
    m <- inputList[[i]][, c("logFC", "P.Value", "adj.P.Val")]
    colnames(m) <- paste(names(inputList)[i], colnames(m), sep = ".")
    idx <- match(genes, inputList[[i]]$ensembl)
    venn.details <- cbind(venn.details, m[idx, ])		
  }				
  
  # RE-ORDER
  newOrder <- order(-venn_info$Variable.number)
  venn_info <- venn_info[newOrder,]
  venn_info$ID <- paste("venn_", seq(1, nrow(venn_info)), sep = "")
  
  # SPLIT
  vennList <- lapply(venn_info$Group, function(i){
    return(venn.details[venn.details$Group == i, ])
  })
  names(vennList) <- venn_info$ID
  
  venn.details <- do.call(rbind, vennList)
  
  # SAVE
  write.xlsx(c(list(Info = venn_info, details = venn.details), vennList),
             paste0(outName, "_pv", pv, "_fc", fc, "_Venn.xlsx"),
             row.names = FALSE, overwrite = TRUE
  )
  
  
  # UpSetPlot
  keep <- FALSE
  if(uad){
    geneList.name <- names(geneList)
    geneList.name <- rev(geneList.name)
    #geneList.name <- sort(geneList.name, decreasing = TRUE)
    geneList.name <- c(grep("DOWN", geneList.name, value = TRUE), grep("UP", geneList.name, value = TRUE))
    keep <- TRUE
  }
  
  mymeta <- data.frame(sets = geneList.name,
                       reg = ifelse(grepl("UP", geneList.name), "UP", "DOWN"))
  
  print(geneList.name)
  
  pdf(paste0(outName, "_pv", pv, "_fc", fc, "_UpSetR.pdf"), width = 10, onefile = FALSE)
  print(upset(fromList(geneList), order.by = "freq", nsets = length(geneList),
              sets = geneList.name, keep.order = keep,
              text.scale = 2.5,
              point.size = 3, line.size = 1.5,
              mainbar.y.label = "Intersection", sets.x.label = "DEG per sample",
              set.metadata = list(data = mymeta, plots = list(list(type = "matrix_rows", column = "reg", assign = 5,
                                                                   colors = c("DOWN" = rgb(44,201,254, maxColorValue = 255), "UP" = rgb(255, 101, 44, maxColorValue = 255)))))))
  dev.off()		
  
}

############################################
###                                      ###
###                 MAIN                 ###
###                                      ###
############################################

# PARAMETERS
dataDir <- file.path("../data/")
limmaDir <- file.path("../protein_limma/")
dir.create(limmaDir)


# READ DATA
setwd(dataDir)
expMat.cfu <- read.xlsx("EPO_proteins_log2.xlsx", sheet = "CFU", rowNames = TRUE)
expMat.placenta <- read.xlsx("EPO_proteins_log2.xlsx", sheet = "placenta", rowNames = TRUE)
expMat.liver <- read.xlsx("EPO_proteins_log2.xlsx", sheet = "liver", rowNames = TRUE)


################
# LIMMA ANALYSIS

setwd(limmaDir)

# CFU-E
mymat <- expMat.cfu
design <- model.matrix(~ 0+factor(c(1,1,1,2,2,2)))
colnames(design) <- c("CFU.WT", "CFU.KO")
fit <- lmFit(mymat, design)

contrast.matrix <- makeContrasts(CFU.WT-CFU.KO, levels=design)

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

res <- topTable(fit2, adjust="BH", number=nrow(mymat), sort.by="P")
res$entrez <- symbol2entrez(rownames(res))
res$symbol <- rownames(res)

write.xlsx(res, "CFU_WT-KO_limma.xlsx", row.names = FALSE)


# Placenta
mymat <- expMat.placenta
design <- model.matrix(~ 0+factor(c(1,1,1,2,2)))
colnames(design) <- c("placenta.WT", "placenta.KO")
fit <- lmFit(mymat, design)


contrast.matrix <- makeContrasts(placenta.WT-placenta.KO, levels=design)

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

res <- topTable(fit2, adjust="BH", number=nrow(mymat), sort.by="P")
res$entrez <- symbol2entrez(rownames(res))
res$symbol <- rownames(res)

write.xlsx(res, "placenta_WT-KO_limma.xlsx", row.names = FALSE)


# Liver
mymat <- expMat.liver
mymat[mymat==0] <- NA
design <- model.matrix(~ 0+factor(c(1,1,1,1,2,2,2)))
colnames(design) <- c("liver.WT", "liver.KO")
fit <- lmFit(mymat, design)


contrast.matrix <- makeContrasts(liver.WT-liver.KO, levels=design)

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

res <- topTable(fit2, adjust="BH", number=nrow(mymat), sort.by="P")
res$entrez <- symbol2entrez(rownames(res))
res$symbol <- rownames(res)

write.xlsx(res, "liver_WT-KO_limma.xlsx", row.names = FALSE)


##########
# DEG VENN

setwd(limmaDir)
limmaFiles <- list.files(pattern = "limma.xlsx$")
names(limmaFiles) <- gsub("_limma.xlsx", "", limmaFiles)
limmaList <- lapply(limmaFiles, read.xlsx, sheet = 1)

pvCutoff <- 0.05
fcCutoff <- 0

limma2Venn(limmaList,
           "WT_vs_KO",
           fc = fcCutoff, pv = pvCutoff,
           adjusted = FALSE, uad = TRUE
)


#############################
# VOLCANO PLOT V2 (LFQ vs PV)

###########################
# LOAD NORMALIZED INTENSITY
setwd(dataDir)
expList <- my.read_xlsx("EPO_proteins_log2.xlsx")
expMat <- expList$liver

# GET DEGs
setwd(limmaDir)
deseqFiles <- list.files(pattern = "_limma.xlsx")
names(deseqFiles) <- gsub("_limma.xlsx", "", deseqFiles)
deseqList <- lapply(deseqFiles, read.xlsx, sheet = 1)
deseqMat <- deseqList$`liver_WT-KO`

deseqMat <- deseqMat[order(-abs(deseqMat$logFC)),]

genes <- deseqMat$symbol

ko.avg <- rowMeans(expMat[, grep("KO", colnames(expMat))])
wt.avg <- rowMeans(expMat[, grep("WT", colnames(expMat))])
all.avg <- rowMeans(data.frame(ko = ko.avg, wt = wt.avg))
all.avg <- all.avg[match(genes, names(all.avg))]

fc <- deseqMat$logFC

pv <- deseqMat$P.Value
names(pv) <- genes

isSignif <- rep("none", length(pv))
isSignif[fc > 0 & pv <= 0.05] <- "WT"
isSignif[fc < 0 & pv <= 0.05] <- "KO"

toShow <- c("Abcb10", "Abcb7", "Aco1", "Akt1", "Alas2", "Arg1", "Dna2", "Dpyd", "Epha4",
            "Htra2", "Nfs1", "Nubp2", "Pdia2", "Ppox", "Prim2", "Rfesd", "Rgn", "Slc39a8",
            "Snca", "Steap3", "Tfrc", "Ubqln1")
toShow <- c(toShow, deseqMat$symbol[1:20])

ggmat <- data.frame(Gene = genes, LFQ = all.avg, Pvalue = pv, FoldChange = -fc, Signif = isSignif)


# filter
ggmat$Gene2 <- toupper(ggmat$Gene)

p <- ggplot(ggmat, aes(x=-FoldChange, y=LFQ, color = Signif))# - FOLDCHANGE !!!
p <- p + geom_vline(xintercept = 0, colour="black", linetype = "longdash", alpha = 1)
p <- p + geom_point(data = ggmat[which(ggmat$Signif == "none"),], size=2, alpha = 1)
p <- p + geom_point(data = ggmat[which(ggmat$Signif != "none"),], size=2, alpha = 1)	
p <- p + scale_color_manual(values=c("KO" = rgb(89, 162, 248, maxColorValue = 255),
                                     "WT" = rgb(241, 182, 128, maxColorValue = 255), "none" = "grey"))
p <- p + xlab("Protein fold-change EpoR+/+ vs. EpoR-/- (Log2)")
p <- p + ylab("Averaged LFQ intensity (Log2)")
p <- p + geom_text_repel(data=ggmat[ggmat$Gene %in% toShow,],
                         aes(label=Gene2), cex = 4, segment.size = 0.25, col = "black", 
                         min.segment.length = unit(0, 'lines'))
p <- p + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                            panel.grid.minor = element_blank())
p <- p + theme(axis.text.x = element_text(size=14),axis.text.y = element_text(size=14))
p <- p + theme(axis.line.x = element_line(color="black", size = 0.5),
               axis.line.y = element_line(color="black", size = 0.5))		


setwd(limmaDir)
pdf("liver_KO_vs_WT_FC_vs_LFQ_V2.pdf", width = 9, height = 4.5, useDingbats=FALSE)
plot(p)
dev.off()





