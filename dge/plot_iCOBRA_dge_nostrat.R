n.sub <- 6

## true values ##
library(fishpond)
library(devtools)
suppressMessages(library(SummarizedExperiment))
suppressMessages(library(iCOBRA))
library(ggplot2)
library(rafalib)

source("../plot_funcs.R")

load("data/simulate.rda")
dge <- union(dge.genes, dte.genes)

load(file="data/sim_gene_nostrat.rda") 
load(paste0("data/simres_gene_nostrat_", n.sub, ".rda"))
load(paste0("data/sleuth_gene_nostrat_", n.sub, ".rda"))


## Plot ##

## Every method is stratification adjusted
padj <- data.frame(row.names=rownames(y.swish))
padj$swish <- rowData(y.swish)$qvalue[match(rownames(padj), rownames(y.swish))]

padj$SAMseq <- y.samseq[match(rownames(padj), names(y.samseq))]
padj$SAMseq <- pmin(padj$SAMseq, 1)

padj$limma <- tt.gene$adj.P.Val[match(rownames(padj), rownames(tt.gene))]
padj$DESeq2 <- res.gene$padj[match(rownames(padj), rownames(res.gene))]

padj$EBSeq <- ebseq.padj[match(rownames(padj), names(ebseq.padj))]
padj$sleuth <- sleuth.res$qval[match(rownames(padj), sleuth.res$target_id)]

status <- as.numeric(rownames(padj) %in% dge)

brks <- c(quantile(mcols(y.swish)[["InfRV"]], probs = seq(0, 1, by = 1/3), na.rm = TRUE))
InfRV_cat = cut(mcols(y.swish)[["InfRV"]], brks)
cd <- COBRAData(padj=padj, truth=data.frame(status = status, InfRV = InfRV_cat, row.names=rownames(padj)))

pdf(file=paste0("plots/dge_nostrat_swish_q25_newpanel_", n.sub, ".pdf"))
	bigpar()
    print(myicobra(cd, lvl="Gene", overall = TRUE, ylim = c(0.45,1), xlim.right = 0.2))
dev.off()        

sessionInfo()
