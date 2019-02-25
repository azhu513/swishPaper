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

load("data/sim_txp_nostrat.rda")
load("data/simulate.rda")
source("../dte_funcs.R")
source("../samseq_func.R")

## Analysis ###
y <- y[mcols(y)$keep, ]
cts <- round(cts[mcols(y)$keep, ])

y.swish <- swish(y, x="condition")
y.samseq <- samseq(cts, y = colData(y)[["condition"]])
tt.txp <- limmavoom(cts, n.sub, condition = colData(y)[["condition"]])
res.txp <- deseq2(cts, n.sub, condition = colData(y)[["condition"]])
ebseq.padj <- ebseq(cts, n.sub, condition = colData(y)[["condition"]], gene=FALSE)

save(y.swish, y.samseq, tt.txp, res.txp, ebseq.padj, file =  paste0("data/simres_txp_nostrat_", n.sub, ".rda"))

sessionInfo()
