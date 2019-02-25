set.seed(1)
n.sub <- 12

library(fishpond)
library(samr)
library(SummarizedExperiment)
library(EBSeq)
library(tximeta)
library(limma)
library(DESeq2)
library(edgeR)

stopifnot(packageVersion("samr") == "2.0")

load("data/sim_gene_strat.rda")
load("data/simulate.rda")
source("../dte_funcs.R")
source("../samseq_func.R")

## Analysis ###
y <- y[mcols(y)$keep, ]
cts <- round(cts[mcols(y)$keep, ])
cts.rbe  <- exp(limma::removeBatchEffect(log(cts+.5), batch=colData(y)[["batch"]]))

y.swishstrat <- swish(y, x="condition", cov="batch")
y.samseqRBE <- samseqRBE.round(cts, y = colData(y)[["condition"]], cov = colData(y)[["batch"]])
tt.gene <- limmavoom(cts, n.sub, condition = colData(y)[["condition"]], cov=colData(y)[["batch"]])
res.gene <- deseq2(cts, n.sub, condition = colData(y)[["condition"]], cov=colData(y)[["batch"]])
ebseq.padjRBE <- ebseq(cts.rbe, n.sub, condition = colData(y)[["condition"]], gene=TRUE)

save(y.swishstrat, y.samseqRBE, tt.gene, res.gene, ebseq.padjRBE, file =  paste0("data/simres_gene_strat_", n.sub, ".rda"))

sessionInfo()
