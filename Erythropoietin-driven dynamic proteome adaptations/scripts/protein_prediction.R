# Load the libraries for annotation
library(org.Mm.eg.db)
library("GO.db")
library("annotate")
library("genefilter")
library("GOstats")
library(sva)
library(labdsv)
library(gage)
library(limma)
library(GeneAnswers)
library(mouse4302.db)
library(Cairo)
library(RColorBrewer)
library(org.Hs.eg.db)

source("/Users/gandrieux/work/epo/scripts/proper/prot_data_reader.r")

############################################
###                                      ###
###               FUNCTION               ###
###                                      ###
############################################

merge_data <- function(m1, m2)
{
  common <- intersect(rownames(m1), rownames(m2))
  common <- sort(common)
  return(cbind(m1[match(common, rownames(m1)),], m2[match(common, rownames(m2)),]))	
}

getAmplitude <- function(v)
{
	return(max(v) - min(v))	
}

getStandardError <- function(obs, pred)
{
	df <- length(obs) - 2
	dist <- obs - pred
	dist2 <- dist^2
	return(sqrt(sum(dist2)/df))
}

getStandardError2 <- function(obs, pred)
{
	n <- length(obs)
	C1=(sum((obs-mean(obs))*(pred-mean(pred))))^2
	C2=sum((obs-mean(obs))^2)
	C3=sum((pred-mean(pred))^2)
	C4=1/(n-2)
	C5=sqrt(C4*(C3-(C1/C2)))
	return(C5)
}

getR2 <- function(obs, pred)
{
	sstot <- sum((obs - mean(obs))^2)
	ssres <- sum((obs - pred)^2)
	r2 <- 1 - (ssres / sstot)
	return(r2)
}

#getR <- function(obs, pred) return(((1 / length(obs)) * sum((obs - mean(obs)) * (pred - mean(pred))) / (sd(obs) * sd(pred)))^2)


############################################
###                                      ###
###                 MAIN                 ###
###                                      ###
############################################

# PARAMETERS
dataDir <- file.path("../data")
outDir <- fil.epath("../prediction")
dir.create(outDir)

# READ DATA
setwd(dataDir)
mymat <- read.xlsx("prediction_inputs.xlsx", sheet = "RNA", rowNames = TRUE)
transTS <- mymat[,c(1:4, 6)]# early

protTS <- read.xlsx("prediction_inputs.xlsx", sheet = "PROTEIN", rowNames = TRUE)

# COMBINE RNA AND PROTEIN

tsMat <- merge_data(protTS, transTS)

testM <- tsMat[,c(1:5)]/tsMat[,c(6:10)]
testSD <- apply(testM, 1, sd)
testMedian <- apply(testM, 1, median)

setwd(outDir)
pdf(file="ratio_standard_deviation_density.pdf")
plot(density(testSD), xlab = "Standard deviation", ylab = "Density")
title("", cex.sub = 0.75)
dev.off()

trans <- tsMat[,c(6:10)]
predP <- trans * testMedian
prot <- tsMat[,c(1:5)]

# SAVE PREDICTIONS
colnames(predP) <- sub("trans", "predicted_protein", colnames(predP))
pred_mat <- cbind(trans, predP)

setwd(outDir)
write.table(pred_mat, "predicted_protein_smooth_24h.txt", quote = FALSE, sep = "\t", row.names = TRUE)

# CORRELATION
setwd(outDir)
pdf(file="predProt_Vs_prot_0h.pdf")
plot(predP[,1], prot[,1],
xlab = "predicted protein intensity", ylab = "protein intensity", pch=20, col=rgb(0,0,0,alpha=0.5))
title(round(cor(predP[,1], prot[,1], method = "spearman"), digits = 3), cex.sub = 0.75)
dev.off()

pdf(file="predProt_Vs_prot_1h.pdf")
plot(predP[,2], prot[,2],
xlab = "predicted protein intensity", ylab = "protein intensity", pch=20, col=rgb(0,0,0,alpha=0.5))
title(round(cor(predP[,2], prot[,2], method = "spearman"), digits = 3), cex.sub = 0.75)
dev.off()

pdf(file="predProt_Vs_prot_2h.pdf")
plot(predP[,3], prot[,3],
xlab = "predicted protein intensity", ylab = "protein intensity", pch=20, col=rgb(0,0,0,alpha=0.5))
title(round(cor(predP[,3], prot[,3], method = "spearman"), digits = 3), cex.sub = 0.75)
dev.off()

pdf(file="predProt_Vs_prot_3h.pdf")
plot(predP[,4], prot[,4],
xlab = "predicted protein intensity", ylab = "protein intensity", pch=20, col=rgb(0,0,0,alpha=0.5))
title(round(cor(predP[,4], prot[,4], method = "spearman"), digits = 3), cex.sub = 0.75)
dev.off()

pdf(file="predProt_Vs_prot_5h.pdf")
plot(predP[,5], prot[,5],
xlab = "predicted protein intensity", ylab = "protein intensity", pch=20, col=rgb(0,0,0,alpha=0.5))
title(round(cor(predP[,5], prot[,5], method = "spearman"), digits = 3), cex.sub = 0.75)
dev.off()


###############
# Predict protein ratio up to 24h

transTS2 <- mymat[rownames(tsMat),c(9:20)]
predP <- transTS2 * testMedian
setwd(outDir)
write.table(predP, "predicted_protein_intensity.txt", sep = "\t", quote = FALSE)


