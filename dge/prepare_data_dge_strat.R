set.seed(1)

suppressPackageStartupMessages(library(samr))
library(tximeta)
library(fishpond)
suppressPackageStartupMessages(library(SummarizedExperiment))

dir <- "data/quants"
subdir <- rep(c("bias","unif"), each=12)
nms_bias <- as.vector(t(outer(1:6,1:2,function(x,y) paste0(x,"_",y))))
nms_unif <- as.vector(t(outer(7:12,1:2,function(x,y) paste0(x,"_",y))))
nms <- c(nms_bias, nms_unif)
files <- file.path(dir, subdir, nms, "quant.sf")
all(file.exists(files))
condition <- factor(sub(".*_(.)","\\1",nms))
batch <- factor(rep(1:2,each=12))
coldata <- data.frame(files, nms, condition, batch)
coldata$names <- paste(rep(c("bias", "unif"), each = 12), coldata$nms)

se <- tximeta(coldata)

gse <- tximeta::summarizeToGene(se, countsFromAbundance="lengthScaledTPM")
gse2 <- tximeta::summarizeToGene(se, varReduce = TRUE)
save(se, gse, gse2, file="data/sim_gene_strat_se.rda")

cts <- assays(gse)[["counts"]]

y <- scaleInfReps(gse)
y <- labelKeep(y, minN=3)

infVar <- assays(gse2)[["variance"]]
mu <- assays(gse2)[["counts"]]
infRV <- pmax(infVar - mu, 0)/(mu + 5) + .01
meanInfRV <- rowMeans(infRV)
mcols(y)$InfRV <- meanInfRV

save(y, cts, file="data/sim_gene_strat.rda") 
sessionInfo()