samseq <- function(cts, y) {
  x <- round(cts)
  out <- capture.output({
    samfit <- SAMseq(x, y, resp.type="Two class unpaired", fdr.output=1)
  })
  sam.padj <- rep(1,nrow(x))
  names(sam.padj) <- rownames(x)
  idx <- as.numeric(samfit$siggenes.table$genes.up[,"Gene Name"])
  sam.padj[idx] <- 1/100 * as.numeric(samfit$siggenes.table$genes.up[,"q-value(%)"])
  idx <- as.numeric(samfit$siggenes.table$genes.lo[,"Gene Name"])
  sam.padj[idx] <- 1/100 * as.numeric(samfit$siggenes.table$genes.lo[,"q-value(%)"])
  sam.padj
}

samseqRBE.round <- function(cts, y, cov) {
  x <- round(cts)
  x <- exp(limma::removeBatchEffect(log(x+.5), batch=cov))
  x <- round(x)
  out <- capture.output({
    samfit <- SAMseq(x, y, resp.type="Two class unpaired", fdr.output=1)
  })
  sam.padj <- rep(1,nrow(x))
  names(sam.padj) <- rownames(x)
  idx <- as.numeric(samfit$siggenes.table$genes.up[,"Gene Name"])
  sam.padj[idx] <- 1/100 * as.numeric(samfit$siggenes.table$genes.up[,"q-value(%)"])
  idx <- as.numeric(samfit$siggenes.table$genes.lo[,"Gene Name"])
  sam.padj[idx] <- 1/100 * as.numeric(samfit$siggenes.table$genes.lo[,"q-value(%)"])
  sam.padj
}