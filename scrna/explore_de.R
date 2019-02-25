alphas <- c(0.01, 0.05, 0.1)
alpha <- 0.01

library(Seurat)
library(fishpond)
library(tximeta)
library(SummarizedExperiment)
library(org.Mm.eg.db)
library(samr)
library(DESeq2)
library(rafalib)
library(goseq)

source("../samseq_func.R")
source("../scrna/sc_funcs.R")

seed <- 131
set.seed(seed)

load("data/mouse_900_sample_lower_filter_y_se.rda")
load("data/mouse_900_sample_lower_filter_scalingfactors.rda")
load("data/mouse_900_sample_lower_filter_DE_res.rda")

hist(sfs)
hist(res$stat, breaks=20, col="grey")
with(res, plot(stat, -log10(qvalue + .001)))

with(res, table(wilcox=wilcox < alpha, swish=qvalue < alpha)

res$wsig <- ifelse(res$wilcox < .01, "W+", "W-")
res$ssig <- ifelse(res$qvalue < .01, "S+", "S-")
res$sig <- factor(paste(res$wsig, res$ssig, sep="/"))

png(file = "plots/mouse_900_scRNA_Boxplot_InfRV.png", width = 400, height = 400, unit = "px")			
infRVList <- with(res, split(InfRV, sig))[c(3,2,4)]
names(infRVList)
with(res, boxplot(infRVList, ylim=c(0,.3), col="skyblue", xlab="method agreement", ylab="InfRV"))
points(1:3, c(.165, .125, .125), pch=20, col="white", cex=5)
text(1:3, c(.165, .125, .125), paste0("#gene=",lengths(infRVList)))
points(1:3, c(.3, .3, .3), pch=2, lwd=2)
dev.off()

sapply(infRVList, function(x) sum(x > .3))

geneList <- sub("\\..*", "", res$gene_id)
geneListI <- ifelse(res$qvalue < .01, 1, 0)
names(geneListI) <- geneList
pwf <- nullp(geneListI, "hg38", "ensGene")
fit <- goseq(pwf, "hg38", "ensGene")
head(fit[fit$ontology == "BP" & fit$numInCat > 10,])
