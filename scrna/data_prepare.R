###################################################################################
id.1 <- "6"
id.2 <- "7"
###################################################################################
library(Seurat)
library(tximeta)
library(fishpond)
library(DESeq2)
library(SummarizedExperiment)

seed <- 131
set.seed(seed)

source("../scrna/sc_funcs.R")

load("data/mouse_900_se.rda")
load("data/mouse_900_quant.rda")
mnh <- readRDS("data/mnh.rds")

cidx.1 <- WhichCells(mnh, ident=id.1)
cidx.2 <- WhichCells(mnh, ident=id.2)
save(cidx.1, cidx.2, file = "data/mouse_900_sample_cellIds.rda")

counts <- mat.d[, c(cidx.1, cidx.2)]
ids <- match(c(cidx.1, cidx.2), colnames(mat))
infReps <- list()
for (k in seq_along(bmat.list.d)){
    temp <- bmat.list.d[[k]]
    infReps[[k]] <- temp[, ids]
}

names(infReps) <-  paste0(rep("infRep", length(infReps)), seq(1:length(infReps)))
infReps[["counts"]] <- counts

save(infReps, file = "data/mouse_900_sample_infReps.rda")

se2 <- se[, c(cidx.1, cidx.2)]
condition <-  as.factor(c(rep(1, length(cidx.1)), rep(2, length(cidx.2))))
colData(se2)$condition <- condition

variance <- assays(se2)[["variance"]]
mu <- assays(se2)[["counts"]]

assays(se2) <- infReps

sfFun <- function(m) estimateSizeFactorsForMatrix(m, geoMeans=exp(rowSums(log(m) * as.numeric(m > 0))/ncol(m)))

y <- scaleInfReps(se2, lengthCorrect=FALSE, sfFun=sfFun, minCount=5, minN=3)

y <- labelKeep(y, minCount=5, minN=3)

infCV <- sqrt(variance)/(mu+0.5)
meanInfCV <- rowMeans(infCV)
mcols(y)$InfCV <- meanInfCV

infRV <- pmax(variance - mu, 0)/(mu + 5) + .01
meanInfRV <- rowMeans(infRV)
mcols(y)$InfRV <- meanInfRV

y <- y[mcols(y)$keep,] 

sfs <- my.saveScaleFactors(se2, lengthCorrect=FALSE, sfFun=sfFun, minCount=5, minN=3, savesf = TRUE)

save(y, se2, file = "data/mouse_900_sample_y_se.rda")

save(sfs, file = "data/mouse_900_sample_scalingfactors.rda")

rm(y, sfs, infCV, infRV)

y <- scaleInfReps(se2, lengthCorrect=FALSE, sfFun=sfFun, minCount=3, minN=5)

y <- labelKeep(y, minCount=3, minN=5)

infCV <- sqrt(variance)/(mu+0.5)
meanInfCV <- rowMeans(infCV)
mcols(y)$InfCV <- meanInfCV

infRV <- pmax(variance - mu, 0)/(mu + 5) + .01
meanInfRV <- rowMeans(infRV)
mcols(y)$InfRV <- meanInfRV

y <- y[mcols(y)$keep,] 

sfs <- my.saveScaleFactors(se2, lengthCorrect=FALSE, sfFun=sfFun, minCount=3, minN=5, savesf = TRUE)

save(y, se2, file = "data/mouse_900_sample_lower_filter_y_se.rda")

save(sfs, file = "data/mouse_900_sample_lower_filter_scalingfactors.rda")
