samseq <- function(cts, condition) {
  x <- round(cts)
  out <- capture.output({
    samfit <- SAMseq(x, condition, resp.type="Two class unpaired", fdr.output=1)
  })
  sam.padj <- rep(1,nrow(x))
  names(sam.padj) <- rownames(x)
  if (samfit$siggenes.table$ngenes.up > 0) {
    up <- samfit$siggenes.table$genes.up
    if (!is.matrix(up)) {
      idx <- as.numeric(up["Gene Name"])
      sam.padj[idx] <- 1/100 * as.numeric(up["q-value(%)"])
    } else {
      idx <- as.numeric(up[,"Gene Name"])
      sam.padj[idx] <- 1/100 * as.numeric(up[,"q-value(%)"])
    }
  }
  if (samfit$siggenes.table$ngenes.lo > 0) {
    lo <- samfit$siggenes.table$genes.lo
    if (!is.matrix(lo)) {
      idx <- as.numeric(lo["Gene Name"])
      sam.padj[idx] <- 1/100 * as.numeric(lo["q-value(%)"])
    } else {
      idx <- as.numeric(lo[,"Gene Name"])
      sam.padj[idx] <- 1/100 * as.numeric(lo[,"q-value(%)"])
    }
  }
  sam.padj <- pmin(sam.padj, 1)
}