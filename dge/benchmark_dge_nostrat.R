set.seed(1)
n.sub <- 6

library(fishpond)
library(samr)
library(SummarizedExperiment)
library(EBSeq)
library(tximeta)
library(limma)
library(DESeq2)
library(edgeR)

sessionInfo()
stopifnot(packageVersion("samr") == "2.0")

load("data/sim_gene_nostrat.rda")
load("data/simulate.rda")
source("../dte_funcs.R")
source("../samseq_func.R")

## Analysis ###
y <- y[mcols(y)$keep, ]
cts <- round(cts[mcols(y)$keep, ])

y.swish <- swish(y, x="condition")
y.samseq <- samseq(cts, y = colData(y)[["condition"]])
tt.gene <- limmavoom(cts, n.sub, condition = colData(y)[["condition"]])
res.gene <- deseq2(cts, n.sub, condition = colData(y)[["condition"]])
ebseq.padj <- ebseq(cts, n.sub, condition = colData(y)[["condition"]], gene=TRUE)

save(y.swish, y.samseq, tt.gene, res.gene, ebseq.padj, file =  paste0("data/simres_gene_nostrat_", n.sub, ".rda"))

sessionInfo()
