edgerPrep <- function(cts, n.sub, condition = NULL, cov = NULL, mincount, mintotalcount) {
  # this uses countsFromAbundance (as used with limma-voom)
  # the results for edgeR were nearly indistinguishable
  # with raw counts + offset and countsFromAbundnace,
  # so for code simplicity I use countsFromAbundance 
  y <- DGEList(cts)
  if (is.null(mincount) & is.null(mintotalcount)){
	  y <- y[filterByExpr(y), ]
  } else {
	  y <- y[filterByExpr(y, min.count = mincount, min.total.count = mintotalcount), ]
  }
  y <- calcNormFactors(y)

  if (is.null(condition)){
      condition <- factor(rep(1:2,each=n.sub))
  } 

  if (is.null(cov)){  
      design <- model.matrix(~condition)
  } else {
      design <- model.matrix(~condition + cov)
  }

  list(y, design)
}

limmavoom <- function(cts, n.sub, condition = NULL, cov = NULL, mincount = NULL, mintotalcount = NULL) {
  out <- edgerPrep(cts, n.sub, condition = condition, cov = cov, mincount, mintotalcount)
  y <- out[[1]]; design <- out[[2]]
  v <- voom(y, design)
  fit <- lmFit(v,design)
  fit <- eBayes(fit)
  topTable(fit, coef=2, n=nrow(y), sort="none")
}

deseq2 <- function(cts, n.sub, condition = NULL, cov = NULL, cooksCutoff = TRUE) {

  if (is.null(condition)){
      condition <- factor(rep(1:2,each=n.sub))
  }

  if (is.null(cov)){  
      samps <- data.frame(condition = condition)
  } else {
      samps <- data.frame(condition= condition,
                          covariates = cov)
  }

  cts <- round(cts)
    suppressMessages({
        if (is.null(cov)) {
            dds <- DESeqDataSetFromMatrix(cts, samps, ~condition)
        } else {
            dds <- DESeqDataSetFromMatrix(cts, samps, ~condition + covariates)
        }
 
   # dds <- dds[,idx]
   # keep <- rowSums(counts(dds) >= 10) >= n.sub
   # table(keep)
   # dds <- dds[keep,]
    dds <- DESeq(dds, minReplicates=Inf, quiet=TRUE)
  })
  results(dds, name = "condition_2_vs_1", cooksCutoff = cooksCutoff)
}


ebseq <- function(cts, n.sub, condition, gene=TRUE) {
  #keep <- rowSums(cts >= 10) >= n.sub
  #x <- round(cts[keep,])
  x <- round(cts)
  if (is.null(condition)){
    condition <- factor(rep(1:2,each=n.sub))
  }
 
  sizes <- MedianNorm(x)
  if (gene) {
    suppressMessages({
      res <- EBTest(Data=x, Conditions=condition,
                    sizeFactors=sizes, maxround=5, Print=FALSE)
    })
  } else {
    suppressMessages({
      Ng <- GetNg(txdf[,2], txdf[,1])$IsoformNgTrun
      res <- EBTest(Data=x, NgVector=Ng, Conditions=condition,
                    sizeFactors=sizes, maxround=5, Print=FALSE)
    })
  }
  padj <- rep(1, nrow(x))
  names(padj) <- rownames(x)
  padj[match(rownames(res$PPMat), rownames(x))] <- res$PPMat[,"PPEE"]
  padj
}





