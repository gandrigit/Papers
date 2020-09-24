library(data.table)

myScale <- function(x, newMin, newMax){ (x - min(x))/(max(x)-min(x)) * (newMax - newMin) + newMin }

list2xlsx <- function(inputList, outputFile)
{	
	require(xlsx)
	wb <- createWorkbook()
	for(i in 1:length(inputList))
		{
		addDataFrame(x = data.frame(inputList[[i]]), sheet = createSheet(wb = wb, sheetName = names(inputList)[i]),
			row.names = FALSE, colnames = FALSE)
		}
	saveWorkbook(wb, outFile)	
}

saveXLSX <- function(inputList, outputFile)
{
	require(xlsx)
	wb <- createWorkbook()
	for(i in 1:length(inputList))
		{
		addDataFrame(x = inputList[[i]], sheet = createSheet(wb = wb, sheetName = names(inputList)[i]), row.names = FALSE)
		}
	saveWorkbook(wb, outputFile)	
}





log2_log10 <- function(x) return(x * log10(2))
log10_log2 <- function(x) return(x/log10(2))

is_increasing <- function(x) return(all(x == cummax(x)))
is_decreasing <- function(x) return(all(x == cummin(x)))

is_increasing_soft <- function(x) return(x[1] < x[length(x)])
is_decreasing_soft <- function(x) return(x[1] > x[length(x)])

is_increasing_lm <- function(x)
{
	lmodel <- lm(as.numeric(x)~seq(1,length(x)))
	return(lmodel$coefficients[2]>0)	
}

is_decreasing_lm <- function(x)
{
	lmodel <- lm(as.numeric(x)~seq(1,length(x)))
	return(lmodel$coefficients[2]<0)	
}

normalized <- function(vect, x, y)
{
     # Normalize to [0, 1]:
     m = min(vect)
     range = max(vect) - m
     nvect = (vect - m) / range

     # Then scale to [x,y]:
     range2 = y - x
     nvect = (nvect*range2) + x
     
     return(nvect)
}

intersect_mat <- function(m1, m2)
{
    m3 <- rbind(m1, m2)
    return(m3[duplicated(m3), , drop = FALSE])
}

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

multi.mixedorder <- function(..., na.last = TRUE, decreasing = FALSE){
    do.call(order, c(
        lapply(list(...), function(l){
            if(is.character(l)){
                factor(l, levels=mixedsort(unique(l)))
            } else {
                l
            }
        }),
        list(na.last = na.last, decreasing = decreasing)
    ))
}
    
    
my.t.test.p.value <- function(...) {
    obj<-try(t.test(...), silent=TRUE)
    if (is(obj, "try-error")) return(NA) else return(obj$p.value)
}    


mergeMatrices <- function(listOfMat)
{
	rNames <- lapply(listOfMat, rownames)
	rNames <- sort(unique(unlist(rNames)))
	
	cNames <- lapply(listOfMat, colnames)
	cNames <- sort(unique(unlist(cNames)))
	
	newMat <- matrix(0, nrow = length(rNames), ncol = length(cNames))
	rownames(newMat) <- rNames
	colnames(newMat) <- cNames
	
	for(mat in listOfMat)
		{
		newMat[rownames(mat), colnames(mat)] <- mat		
		}
	return(newMat)
}

my.read.table <- function(inFile, row.names = NULL, colClasses = NULL)
{
	header <- read.table(inFile, header = TRUE, nrow = 1)
	indata <- fread(inFile, skip = 1, header = FALSE, colClasses = colClasses)
	indata <- data.frame(indata, row.names = row.names)
	setnames(indata, names(header))
	return(indata)
}

getR2 <- function(x, y)
{
	m <- data.frame(x = x, y = y)
	model <- lm(y ~x, m)
	return(summary(model)$r.squared)
}

getR2Matrix <- function(dataMat)
{
	r2m <- matrix(NA, nrow = ncol(dataMat), ncol = ncol(dataMat))
	rownames(r2m) <- colnames(dataMat)
	colnames(r2m) <- colnames(dataMat)
	
	cbMat <- t(combn(rownames(r2m), 2))
	for(i in 1:nrow(cbMat))
		{
		cbMat.x <- cbMat[i,1]
		cbMat.y <- cbMat[i,2]
		r2 <- getR2(dataMat[,cbMat.x], dataMat[,cbMat.y])
		r2m[cbMat.x, cbMat.y] <- r2
		r2m[cbMat.y, cbMat.x] <- r2
		}
	return(r2m)
}

getJaccard <- function(x, y) return(length(intersect(x,y)) / length(union(x,y)))

my.head <- function(m)
{
	if(nrow(m) > 5 & ncol(m) > 5) return(m[1:5,1:5])
	else if(nrow(m) > 5) return(m[1:5,])
	else(return(m))
}  
  
# split 1 vector into vector of size n  
require(plyr)
plyrChunks <- function(d, n){
     is <- seq(from = 1, to = length(d), by = n)
     if(tail(is, 1) != length(d)) {
          is <- c(is, length(d)) 
     } 
     chunks <- llply(head(seq_along(is), -1), 
                     function(i){
                         start <-  is[i];
                         end <- is[i+1]-1;
                         d[start:end]})
    lc <- length(chunks)
    td <- tail(d, 1)
    chunks[[lc]] <- c(chunks[[lc]], td)
    return(chunks)
 }
 
sem <- function(x) return(sd(x)/sqrt(length(x)))
  
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
 return(data_sum)
} 
 
  
bigcor <- function(x, y = NULL, fun = c("cor", "cov"), size = 2000, verbose = TRUE, ...)
{
  fun <- match.arg(fun)
  if (fun == "cor") FUN <- cor else FUN <- cov
  if (fun == "cor") STR <- "Correlation" else STR <- "Covariance" 
  if (!is.null(y) & NROW(x) != NROW(y)) stop("'x' and 'y' must have compatible dimensions!")
   
  NCOL <- ncol(x)
  if (!is.null(y)) YCOL <- NCOL(y)
    
  ## calculate remainder, largest 'size'-divisible integer and block size
  REST <- NCOL %% size
  LARGE <- NCOL - REST  
  NBLOCKS <- NCOL %/% size
    
  ## preallocate square matrix of dimension
  ## ncol(x) in 'ff' single format
  if (is.null(y)) resMAT <- ff(vmode = "double", dim = c(NCOL, NCOL))  
  else resMAT <- ff(vmode = "double", dim = c(NCOL, YCOL))
  
  ## split column numbers into 'nblocks' groups + remaining block
  GROUP <- rep(1:NBLOCKS, each = size)
  if (REST > 0) GROUP <- c(GROUP, rep(NBLOCKS + 1, REST))
  SPLIT <- split(1:NCOL, GROUP)
  
  ## create all unique combinations of blocks
  COMBS <- expand.grid(1:length(SPLIT), 1:length(SPLIT))
  COMBS <- t(apply(COMBS, 1, sort))
  COMBS <- unique(COMBS)  
  if (!is.null(y)) COMBS <- cbind(1:length(SPLIT), rep(1, length(SPLIT)))
  
  ## initiate time counter
  timeINIT <- proc.time() 
  
  ## iterate through each block combination, calculate correlation matrix
  ## between blocks and store them in the preallocated matrix on both
  ## symmetric sides of the diagonal
  for (i in 1:nrow(COMBS)) {
    COMB <- COMBS[i, ]    
    G1 <- SPLIT[[COMB[1]]]
    G2 <- SPLIT[[COMB[2]]]    
    
    ## if y = NULL
    if (is.null(y)) {
      if (verbose) cat(sprintf("#%d: %s of Block %s and Block %s (%s x %s) ... ", i, STR,  COMB[1],
                               COMB[2], length(G1),  length(G2)))      
      RES <- FUN(x[, G1], x[, G2], ...)
      resMAT[G1, G2] <- RES
      resMAT[G2, G1] <- t(RES) 
    } else ## if y = smaller matrix or vector  
    {
      if (verbose) cat(sprintf("#%d: %s of Block %s and 'y' (%s x %s) ... ", i, STR,  COMB[1],
                               length(G1),  YCOL))    
      RES <- FUN(x[, G1], y, ...)
      resMAT[G1, ] <- RES             
    }
    
    if (verbose) {
      timeNOW <- proc.time() - timeINIT
      cat(timeNOW[3], "s\n")
    }
    
    gc()
  } 
  
  return(resMAT)
}
  
