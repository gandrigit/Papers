hyperG <- function(geneSets,DEgenes,universe, org.library, cutoff=0.1, mincount=2, parallel=T,adj.P.Val = F){
  library(org.library, character.only = T)
  library(foreach)
  library(doMC)
  if(parallel){
    registerDoMC(cores=detectCores())
    cores=detectCores()
  }else{
    cores=1
  }
  results <- mclapply(1:length(geneSets), function(i){
    results <- matrix(data=NA,ncol=9,nrow = 1)
    colnames(results) <- c('Term','Count','Size.univ', 'Size.tot','p-value','adj.P.Val','odds ratio','Entrez','Symbol')
    geneSet <- intersect(universe, geneSets[[i]])
    a <- length(intersect(DEgenes,geneSet))
    b <- length(setdiff(DEgenes,intersect(DEgenes,geneSet)))
    c <- length(setdiff(geneSet,intersect(DEgenes,geneSet)))
    d <- length(setdiff(universe,DEgenes)) - c
    contigency.matrix <- cbind(c(a,b),c(c,d))
    res <- fisher.test(contigency.matrix,alternative = 'greater')
    results[1,'Term'] <- names(geneSets)[i]
    results[1,'Count'] <- a
    results[1,'Size.univ'] <- length(intersect(geneSets[[i]], universe))
    results[1,'Size.tot'] <- length(geneSets[[i]])
    results[1,'p-value'] <- res$p.value
    results[1,'odds ratio'] <- res$estimate[[1]]
    # find genes annotated in the consensus term
    if(a > 0){
      genes <- intersect(DEgenes,geneSet)
      eid <- genes
      eid <- eid[order(eid)]
      results[1,'Entrez'] <- paste(eid,collapse="|")
    }
    return(results)
  }, mc.cores=cores)
    
  results <- as.data.frame(do.call(rbind, results))
  for(i in c(2, 3, 4, 5)){
    results[, i] <- as.numeric(as.character(results[, i]))
  }
  
  if(nrow(results) != 1){
    results <- results[order(results[,'p-value'],decreasing = FALSE),]
	results[,'adj.P.Val'] <- p.adjust(results[,'p-value'], 'BH')
	if(adj.P.Val){
		results <- as.data.frame(subset(results,results[,'adj.P.Val']<=cutoff))
	}else{
		results <- as.data.frame(subset(results,results[,'p-value']<=cutoff))
	}
    results <- as.data.frame(subset(results,results[,'Count']>=mincount))
  }else results <- as.data.frame(results)
  
  org.symb <- gsub(".db", "SYMBOL", org.library)
  # find genes 
  results$Symbol <- sapply(results$Entrez, function(x){
    y <- unlist(strsplit(as.character(x), "|", fixed=T))
    syms <- paste(unlist(mget(y, eval(parse(text=org.symb)),ifnotfound = NA)), collapse="|")
  })
  return(results)
}