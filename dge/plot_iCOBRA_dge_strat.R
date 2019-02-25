n.sub <- 12

## true values ##
library(fishpond)
library(devtools)
suppressMessages(library(SummarizedExperiment))
suppressMessages(library(iCOBRA))
library(ggplot2)
source("../plot_funcs.R")
library(rafalib)

sessionInfo()

load("data/simulate.rda")
dge <- union(dge.genes, dte.genes)

load(file="data/sim_gene_strat.rda") 
load(paste0("data/simres_gene_strat_", n.sub, ".rda"))
load(paste0("data/sleuth_gene_strat_", n.sub, ".rda"))

## Plot ##

## Every method is stratification adjusted
padj <- data.frame(row.names=rownames(y.swishstrat))
padj$swish <- rowData(y.swishstrat)$qvalue[match(rownames(padj), names(y.swishstrat))]

padj$SAMseq <- y.samseqRBE[match(rownames(padj), names(y.samseqRBE))]
padj$SAMseq <- pmin(padj$SAMseq, 1)

padj$limma <- tt.gene$adj.P.Val[match(rownames(padj), rownames(tt.gene))]
padj$DESeq2 <- res.gene$padj[match(rownames(padj), rownames(res.gene))]

padj$EBSeq <- ebseq.padjRBE[match(rownames(padj), names(ebseq.padjRBE))]

status <- as.numeric(rownames(padj) %in% dge)

brks <- c(quantile(mcols(y.swishstrat)[["InfRV"]], probs = seq(0, 1, by = 1/3), na.rm = TRUE))
InfRV_cat = cut(mcols(y.swishstrat)[["InfRV"]], brks)
cd <- COBRAData(padj=padj, truth=data.frame(status = status, InfRV = InfRV_cat, row.names=rownames(padj)))

pdf(file=paste0("plots/dge_strat_swish_q25_newpanel_nosleuth_", n.sub, ".pdf"))
    print(myicobra(cd, lvl="Gene", overall = TRUE, ylim = c(0.45,1), xlim.right = 0.2))
dev.off()   

padj$sleuth <- sleuth.res$qval[match(rownames(padj), sleuth.res$target_id)]
cd <- COBRAData(padj=padj, truth=data.frame(status = status, InfRV = InfRV_cat, row.names=rownames(padj)))
pdf(file=paste0("plots/dge_strat_swish_q25_newpanel_", n.sub, ".pdf"))
    print(myicobra(cd, lvl="Gene", overall = TRUE, ylim = c(0.45,1), xlim.right = 0.6))
dev.off()   