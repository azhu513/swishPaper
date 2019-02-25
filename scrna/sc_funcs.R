getWilcoxStat <- function(cts, condit, paired = FALSE){
	dims <- dim(cts)
	cts.notie <- cts + 0.1 * runif(dims[1]*dims[2])
	idx <- 1:dims[1]
	wilcox.p <- sapply(idx, function(i) {
      wilcox.test(cts.notie[i, ] ~ condit)$p.value
	})
	wilcox.padj <- p.adjust(wilcox.p, method = "BH")
	names(wilcox.padj) <- rownames(cts)
	wilcox.padj
}

zeroproportion <- function(x, row = TRUE){
    if (row) {
        apply(x, 1, function(i) sum(i == 0)/length(i))
    } else {
        apply(x, 2, function(i) sum(i == 0)/length(i))
    }
}

translate <- function(ENSList = NULL, SMBList = NULL){
	if (is.null(ENSList) & is.null(SMBList)){
		stop("No input")}
	if (!is.null(ENSList) & !is.null(SMBList)){
		stop("Too many arguments")}
	if (!is.null(ENSList) & is.null(SMBList)){
		k <- gsub( "[.].*$", "", ENSList)
		ls <- mapIds(org.Mm.eg.db, column="SYMBOL", keys=k, keytype="ENSEMBL")
	}
	if (is.null(ENSList) & !is.null(SMBList)){
		ls <- mapIds(org.Mm.eg.db, column="ENSEMBL", keys=SMBList, keytype="SYMBOL")
	}
	ls
}

'%!in%' <- function(x,y)!('%in%'(x,y))

my.saveScaleFactors <- function(y, lengthCorrect=TRUE, meanDepth=NULL, sfFun=NULL, minCount=10, minN=3,
							savesf = TRUE) {
  infRepIdx <- grep("infRep",assayNames(y))
  infReps <- assays(y)[infRepIdx]
  counts <- assays(y)[["counts"]]
  length <- assays(y)[["length"]]
  nreps <- length(infReps)
  
  dims <- dim(counts)
  
  if (savesf) {
	sfs <- matrix(NA, nrow = nreps, ncol = dims[2])
  }
  
  if (is.null(meanDepth)) {
    meanDepth <- exp(mean(log(colSums(counts))))
  }
  for (k in seq_len(nreps)) {
    cat(k,"")
    # we don't technically create TPM, but something proportion to
    if (lengthCorrect) {
      tpm <- infReps[[k]] / length
    } else {
      # for 3' tagged scRNA-seq for example, don't length correct
      tpm <- infReps[[k]]
    }
    # divide out the column sum, then set all to the meanDepth
    tpm <- t(t(tpm) / colSums(tpm)) * meanDepth
    # filtering for calculting median ratio size factors
    use <- rowSums(infReps[[k]] >= minCount) >= minN
    if (is.null(sfFun)) {
      loggeomeans <- rowMeans(log(tpm[use,]))
      sf <- apply(tpm[use,], 2, function(s) {
        exp(median((log(s) - loggeomeans)[is.finite(loggeomeans)]))
      })
    } else {
      sf <- sfFun(tpm)
    }
	if (savesf) {
		sfs[k, ] <- t(as.vector(sf))
	}
    infReps[[k]] <- t( t(tpm)/sf )
  }
  cat("\n")
  assays(y)[grep("infRep",assayNames(y))] <- infReps
  
	if (savesf) {
		return(sfs)
	} else {
		return(y)
	}
}