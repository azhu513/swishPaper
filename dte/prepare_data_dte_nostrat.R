set.seed(1)

suppressPackageStartupMessages(library(samr))
library(tximeta)
library(fishpond)
suppressPackageStartupMessages(library(SummarizedExperiment))

dir <- "data/quants"
subdir <- rep("bias", 12)
names <- as.vector(t(outer(1:6,1:2,function(x,y) paste0(x,"_",y))))
files <- file.path(dir, subdir, names, "quant.sf")
all(file.exists(files))
condition <- factor(sub(".*_(.)","\\1",names))
coldata <- data.frame(files, names, condition)
se <- tximeta(coldata)
se2 <- tximeta(coldata, varReduce=TRUE)
save(se, se2, file="data/sim_txp_nostrat_se.rda")

se_cfa <- tximeta(coldata, countsFromAbundance="lengthScaledTPM", skipMeta=TRUE)
cts <- assays(se_cfa)[["counts"]]

y <- scaleInfReps(se)
y <- labelKeep(y, minN=3)

infVar <- assays(se2)[["variance"]]
mu <- assays(se2)[["counts"]]
infRV <- pmax(infVar - mu, 0)/(mu + 5) + .01
meanInfRV <- rowMeans(infRV)
mcols(y)$InfRV <- meanInfRV

save(y, cts,  file="data/sim_txp_nostrat.rda")

sessionInfo()
