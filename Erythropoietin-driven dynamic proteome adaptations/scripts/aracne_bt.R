############################################
###                                      ###
###               FUNCTION               ###
###                                      ###
############################################

getRowValues <- function(r)
{
  gene_ref <- r[1]
  r <- r[-1]
  genes <- r[seq(1, length(r), 2)]
  values <- as.numeric(r[seq(2, length(r), 2)])
  
  genes <- genes[!is.na(genes)]
  values <- values[!is.na(values)]
  
  names(values) <- paste(gene_ref, genes, sep = "~")
  return(values)
}

getValue <- function(a, b)
{
  if (a == b) return(0)
  return(mi[paste(a, b, sep = "~")])
}

fillAdjMatrix <- function(mi, mat)
{
  for (i in 1:nrow(mat))
  {
    for (j in 1:ncol(mat))
    {
      row_gene <- rownames(mat)[i]
      col_gene <- colnames(mat)[j]
      if (row_gene != col_gene)
      {
        miName <- paste(row_gene, col_gene, sep = "~")
        mat[i,j] <- mi[miName]
      }
    }	
  }
  return(mat)	
}

toAdj <- function(mat)
{
  genes <- as.character(mat[,1])
  
  mi <- vector("list", nrow(mat))
  for (i in 1:nrow(mat))
  {
    mi[[i]] <- getRowValues(t(mat[i,]))
  }
  mi <- unlist(mi)
  
  adjM <- matrix(0, nrow = nrow(mat), ncol = nrow(mat))
  rownames(adjM) <- genes
  colnames(adjM) <- genes
  adjM <- fillAdjMatrix(mi, adjM)
  return(adjM)
}

############################################
###                                      ###
###                 MAIN                 ###
###                                      ###
############################################

# PARAMETERS
dataDir <- file.path("../data/")
aracneIN <- file.path("../aracne_bt/inputs")
dir.create(aracneIN, recursive = TRUE)
aracneOUT <- file.path("../aracne_bt/outputs")
dir.create(aracneOUT, recursive = TRUE)

#######################
# PREPARE ARACNE INPUTS

# Load intensity matrix
setwd(dataDir)
mat <- read.delim("epo_players_predIntensity.txt")
mat <- na.omit(mat)

toRemove <- c("Epo", "Epor", "Akt1", "Dnmt1", "Jak2", "Mapk3", "Ptpn6")
mat <- mat[!(rownames(mat) %in% toRemove),]

# Create and save bootstrap matrices
setwd(aracneIN)
nb <- 10000
for(i in 1:nb)
{
  btMat <- mat[,sample(1:ncol(mat), ncol(mat), replace = TRUE)]
  write.table(btMat, paste("epo_", i, ".txt", sep = ""), sep = "\t", quote = FALSE)
}


####################
# RUN "aracne_bt.sh"
# sh arcne_bt.sh

#####################
# READ ARACNE OUTPUTS

setwd(aracneOUT)
inputFile <- "_bs10000_p1e-06.adj"
outputFile <- paste(substring(inputFile, 1, nchar(inputFile)-4), ".txt", sep ="")

rawAdjM <- read.delim(inputFile, header = FALSE, skip = 18,
                      col.names = paste0("V",seq_len(100)), fill = TRUE)
rawAdjM[rawAdjM==""] <- NA	

adjM <- toAdj(rawAdjM)

adjM[is.na(adjM)] <- 0

write.table(adjM, outputFile, quote = FALSE, sep ="\t")


