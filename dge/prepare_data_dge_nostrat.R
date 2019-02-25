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

#### Different from here with the other version #####

gse <- tximeta::summarizeToGene(se, countsFromAbundance="lengthScaledTPM")
gse2 <- tximeta::summarizeToGene(se, varReduce = TRUE)

save(se, gse, gse2, file="data/sim_gene_nostrat_se.rda")

#gse_cfa <- tximeta::summarizeToGene(se, countsFromAbundance="lengthScaledTPM")
cts <- assays(gse)[["counts"]]

y <- scaleInfReps(gse)
y <- labelKeep(y, minN=3)

infVar <- assays(gse2)[["variance"]]
mu <- assays(gse2)[["counts"]]
infRV <- pmax(infVar - mu, 0)/(mu + 5) + .01
meanInfRV <- rowMeans(infRV)
mcols(y)$InfRV <- meanInfRV

save(y, cts, file="data/sim_gene_nostrat.rda") 

sessionInfo()
