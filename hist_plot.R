library(tximeta)
suppressPackageStartupMessages(library(SummarizedExperiment))
dir <- "quants"
subdir <- rep(c("bias","unif"), each=12)
nms_bias <- as.vector(t(outer(1:6,1:2,function(x,y) paste0(x,"_",y))))
nms_unif <- as.vector(t(outer(7:12,1:2,function(x,y) paste0(x,"_",y))))
nms <- c(nms_bias, nms_unif)
files <- file.path(dir, subdir, nms, "quant.sf")
all(file.exists(files))
condition <- factor(sub(".*_(.)","\\1",nms))
batch <- factor(rep(1:2,each=12))
names <- paste0(nms, "_b", batch)
coldata <- data.frame(files, names, condition, batch, stringsAsFactors=FALSE)
coldata <- coldata[1:12,]

system.time({
  se <- tximeta(coldata, varReduce=TRUE)
})
system.time({
  se0 <- tximeta(coldata)
})
system.time({
  gse <- summarizeToGene(se0, varReduce=TRUE)
})

gkeep <- rowSums(assays(gse)[["counts"]] >= 10) >= 3
ginfVar <- assays(gse)[["variance"]]
gmu <- assays(gse)[["counts"]]
ginfRV <- pmax(ginfVar - gmu, 0)/(gmu + 5) + .01
gmeanInfRV <- rowMeans(ginfRV)
hist(log10(gmeanInfRV)[gkeep & gmeanInfRV > .01],breaks=100)

keep <- rowSums(assays(se)[["counts"]] >= 10) >= 3
infVar <- assays(se)[["variance"]]
mu <- assays(se)[["counts"]]
infRV <- pmax(infVar - mu, 0)/(mu + 5) + .01
meanInfRV <- rowMeans(infRV)
hist(log10(meanInfRV)[keep & meanInfRV > .01],breaks=100)

gdat <- data.frame(InfRV=log10(gmeanInfRV)[gkeep & gmeanInfRV > .01])
dat <- data.frame(InfRV=log10(meanInfRV)[keep & meanInfRV > .01 & meanInfRV < 100])
library(ggplot2)
library(cowplot)
g1 <- ggplot(gdat, aes(x=InfRV)) + geom_histogram(binwidth=.1,color="white") + ylim(0,3000) +
  ggtitle("Gene level") + xlab("log10 of mean InfRV")
g2 <- ggplot(dat, aes(x=InfRV)) + geom_histogram(binwidth=.1,color="white") + ylim(0,3000) +
  ggtitle("Transcript level") + xlab("log10 of mean InfRV")
pdf(file="hist_mean_infrv.pdf", height=3)
plot_grid(g1, g2, labels=c("A","B"))
dev.off()
