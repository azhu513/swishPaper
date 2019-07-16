eval.highrep <- function(input, n.sub, batch = FALSE){
	y <- scaleInfReps(input$se)
	y <- labelKeep(y)

	infVar <- assays(input$vse)[["variance"]]
	mu <- assays(input$vse)[["counts"]]
	infRV <- pmax(infVar - mu, 0)/(mu + 5) + .01
	meanInfRV <- rowMeans(infRV)
	mcols(y)$InfRV <- meanInfRV

	cts <- round(input$txi[mcols(y)$keep, ])
	y <- y[mcols(y)$keep, ]
	
	if (batch) {
		cts.rbe <- exp(limma::removeBatchEffect(log(cts+.5), batch=colData(y)[["batch"]])) 
		
		set.seed(1)
		y <- swish(y, x="condition", cov="batch")
		set.seed(1)
		samseq <- samseq(cts.rbe, condition = colData(y)[["condition"]])

		tt.gene <- limmavoom(cts, n.sub, condition = colData(y)[["condition"]], cov=colData(y)[["batch"]])
	} else {
	set.seed(1)
	y <- swish(y, x="condition")
	
	set.seed(1)
	SAMseq.padj <- samseq(cts, condition = colData(y)[["condition"]])

	tt.gene <- limmavoom(cts, n.sub, condition = colData(y)[["condition"]])
	}
	padj <- data.frame(row.names=rownames(input$se))
	padj$swish <- rowData(y)$qvalue[match(rownames(padj), rownames(y))]
	padj$SAMseq <- SAMseq.padj[match(rownames(padj), names(SAMseq.padj))]
	padj$limma <- tt.gene$adj.P.Val[match(rownames(padj), rownames(tt.gene))]
	padj$infRV <- meanInfRV
	
	return(padj)	
}

falsePos <- function(x){
	swish <- sum(x$swish < 0.05, na.rm = T)/length(x$swish[!x$swish])
	SAMseq <- sum(x$SAMseq < 0.05, na.rm = T)/length(x$SAMseq[!x$SAMseq])
	limma <- sum(x$limma < 0.05, na.rm = T)/length(x$limma[!x$limma])
	as.vector(c(swish, SAMseq, limma))
}



falsePosByPanel <- function(x, group.idx){
	swish <- sum(x$swish[group.idx] < 0.05, na.rm = T)
	SAMseq <- sum(x$SAMseq[group.idx] < 0.05, na.rm = T)
	limma <- sum(x$limma[group.idx] < 0.05, na.rm = T)
	as.vector(c(swish, SAMseq, limma))
}