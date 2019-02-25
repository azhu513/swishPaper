library(Seurat)
library(fishpond)
library(tximeta)
library(SummarizedExperiment)
library(org.Mm.eg.db)
library(samr)
library(DESeq2)

source("../samseq_func.R")
source("../scrna/sc_funcs.R")

load(file = "data/mouse_900_sample_lower_filter_y_se.rda")

seed <- 131
set.seed(seed)

dds <- DESeqDataSetFromMatrix(round(assays(y)[["counts"]]), data.frame(x=1:129), ~1)
dds <- estimateSizeFactors(dds, type="poscounts")
n.cts <- counts(dds, normalized=TRUE)

y.swish <- swish(y, x = "condition")
res <- mcols(y.swish)[mcols(y.swish)[["keep"]], ]
res.wilcox <- getWilcoxStat(n.cts, condit = colData(y)[["condition"]])
res$wilcox <- res.wilcox 

save(res, y.swish, res.wilcox, file = "data/mouse_900_sample_lower_filter_DE_res.rda")
