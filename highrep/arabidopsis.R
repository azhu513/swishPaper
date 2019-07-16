n.sub <- 8
rep <- 100 						#repitition of the analysis

library(fishpond)
library(samr)
library(SummarizedExperiment)
library(tximeta)
library(limma)
library(edgeR)
library(iCOBRA)
library(caret)

wd <- "/proj/milovelab/zhu/projects/inf-uncertainty/"
setwd(wd)

load("data/athaliana.rda")
source("inf-uncertainty-2/RCodes/samseq_func.R")
source("/proj/milovelab/zhu/projects/inf-uncertainty/inf-uncertainty-2/RCodes/dte_funcs.R")
source("/proj/milovelab/zhu/projects/inf-uncertainty/inf-uncertainty-2/RCodes/yeast/eval.highrep.R")

batch1 <- factor(ifelse(se$batch == 1, "yes", "no"))

set.seed(1)
part <- createDataPartition(batch1, times=rep, p=7/16)

test.gene <- NULL
test.txp <- NULL

for (i in 1:rep) {
	idx <- c(part[[i]], setdiff(seq_len(16),part[[i]]))
	
	se.t <- se[, idx]
	gse.t <- gse[, idx]
	colData(se.t)$condition <- as.factor(c(rep("A", n.sub), rep("B", n.sub)))
	colData(gse.t)$condition <- as.factor(c(rep("A", n.sub), rep("B", n.sub)))
	
	colData(se.t)$batch <- as.factor(ifelse(colData(se.t)$batch == 1, "yes", "no"))
	colData(gse.t)$batch <- as.factor(ifelse(colData(gse.t)$batch == 1, "yes", "no"))

	
	test.txp[[i]] <- list(idx = idx, se = se.t, vse = vse[, idx], txi = txi$counts[, idx])
	test.gene[[i]] <- list(idx = idx, se = gse.t, vse = vgse[, idx], txi = gtxi$counts[, idx])
}

save(test.gene, test.txp, file = paste0("data/arabidopsis_", n.sub, ".rda"))

res.txp <- lapply(test.txp, FUN = eval.highrep, n.sub = n.sub, batch = FALSE)
res.gene <- lapply(test.gene, FUN = eval.highrep, n.sub = n.sub, batch = FALSE)

save(res.txp, res.gene, file = paste0("data/arabidopsis_res_", n.sub, ".rda"))
