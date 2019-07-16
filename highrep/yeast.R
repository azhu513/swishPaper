args <- commandArgs(trailingOnly = TRUE)
n.sub <- as.numeric(args[1])
rep <- 100 						#repitition of the analysis

library(fishpond)
library(samr)
library(SummarizedExperiment)
library(tximeta)
library(limma)
library(edgeR)
library(iCOBRA)

wd <- "/proj/milovelab/zhu/projects/inf-uncertainty/"
setwd(wd)

load("data/yeast.rda")
source("inf-uncertainty-2/RCodes/samseq_func.R")
source("/proj/milovelab/zhu/projects/inf-uncertainty/inf-uncertainty-2/RCodes/dte_funcs.R")
source("/proj/milovelab/zhu/projects/inf-uncertainty/inf-uncertainty-2/RCodes/yeast/eval.highrep.R")

test.list <- NULL

for (i in 1:rep) {
	set.seed(i)
	idx <- sample(seq_len(ncol(se)), n.sub * 2, replace = F)
	se.t <- se[, idx]
	colData(se.t)$condition <- as.factor(c(rep("A", n.sub), rep("B", n.sub)))
	test.list[[i]] <- list(idx = idx, se = se.t, vse = vse[, idx], txi = txi$counts[, idx])
}

save(test.list, file = paste0("data/yeast_", n.sub, ".rda"))

res <- lapply(test.list, FUN = eval.highrep, n.sub = n.sub)

save(res, file = paste0("data/yeast_res_", n.sub, ".rda"))
