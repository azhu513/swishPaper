# cluster 7 vs cluster 5
load("mouse_900_sample_lower_filter_DE_res.rda")

library(SummarizedExperiment)
y <- y.swish
mcols(y)$wilcox <- res.wilcox

library(DESeq2)
library(fishpond)
library(ggplot2)

library(abind)
infReps <- abind(as.list(assays(y)[1:20]), along=3)
infVar <- apply(infReps, 1:2, var)
mu <- apply(infReps, 1:2, mean)
infRV <- pmax(infVar - mu, 0)/(mu + 5) + .01 

mcols(y)$InfRV <- rowMeans(infRV)

getWilcoxStat <- function(cts, condit) {
  dims <- dim(cts)
  cts.notie <- cts + 0.1 * runif(dims[1]*dims[2])
  wilcox.p <- sapply(seq_len(dims[1]), function(i) {
    wilcox.test(cts.notie[i, ] ~ condit)$p.value
  })
  wilcox.padj <- p.adjust(wilcox.p, method = "BH")
  names(wilcox.padj) <- rownames(cts)
  wilcox.padj
}

dds.cts <- DESeqDataSetFromMatrix(round(assays(y)[["counts"]]),
                                  data.frame(x=seq_len(ncol(y))), ~1)
dds.cts <- estimateSizeFactors(dds.cts, type="poscounts")
n.cts <- counts(dds.cts, normalized=TRUE)

set.seed(1)
n.sub <- 20
smp.idx <- c(sample(which(y$condition==1),n.sub),
             sample(which(y$condition==2),n.sub))
set.seed(1)
wilcox <- getWilcoxStat(n.cts[,smp.idx], condit=y$condition[smp.idx])

#plot(-log10(mcols(y)$wilcox), -log10(wilcox))

set.seed(1)
y2 <- swish(y[,smp.idx], x="condition")
#plot(-log10(mcols(y)$qvalue), -log10(mcols(y2)$qvalue))

#plot(-log10(wilcox), -log10(mcols(y2)$qvalue)); abline(0,1,col="red")

####

mcols(y)$wilcox <- wilcox
mcols(y)$qvalue <- mcols(y2)$qvalue

alpha <- .05
with(mcols(y), table(wilcox=wilcox < alpha, swish=qvalue < alpha))

mcols(y)$wsig <- ifelse(mcols(y)$wilcox < alpha, "W+", "W-")
mcols(y)$ssig <- ifelse(mcols(y)$qvalue < alpha, "S+", "S-")
mcols(y)$sig <- factor(paste(mcols(y)$wsig, mcols(y)$ssig, sep="/"))

# load in the bulk RNA-seq of VZ and CP, taken as pseudo-gold-standard

if (!"txi" %in% ls()) {
  load("txi.rda")
}

common <- intersect(rownames(y), rownames(txi$abundance))
dat <- data.frame(log2FC=mcols(y)[common,"log2FC"],
                  CP=log2(rowMeans(txi$abundance[common,c(7:11)])+5),
                  VZ=log2(rowMeans(txi$abundance[common,c(12:16)])+5),
                  category=mcols(y)[common,"sig"])
dat$size <- ifelse(dat$category %in% c("W+/S-","W-/S+","W+/S+"), 2, 1)

# dat[grep("ENSMUSG00000032446",rownames(dat)),] # Eomes
# dat[grep("ENSMUSG00000038255",rownames(dat)),] # Neurod2

col <- scale_colour_gradient(low="black",high="green")
rgb2 <- function(x,y,z) rgb(x,y,z,maxColorValue=255)
dat <- dat[c(which(dat$category == "W-/S-"),
             which(dat$category == "W+/S+"),
             which(dat$category == "W-/S+"),
             which(dat$category == "W+/S-")),]
with(dat, cor(log2FC, VZ-CP)) # cor = .61

pdf(file="scrna_lfc_scatter_20.pdf", width=5, height=5)
ggplot(dat, aes(VZ - CP, log2FC, col=category, size=size)) +
  geom_point(alpha=.75) +
  geom_hline(yintercept=0) + geom_vline(xintercept=0) +
  ylim(-2,2) + xlim(-4,4) +
  ylab("Cluster 7 (Eomes) - Cluster 5 (Neurod2)") +
  scale_size(range=c(.5, 1.5), guide="none") +
  scale_color_manual(values=c("grey",
                              rgb2(230,159,0),
                              rgb2(0,114,178),
                              rgb2(0,158,115))) +
  theme_bw()
dev.off()

if (!"dds" %in% ls()) {
  cts <- round(txi$counts[,7:16])
  dds <- DESeqDataSetFromMatrix(cts, data.frame(x=rep(c("CP","VZ"),each=5)), ~x)
  keep <- rowSums(counts(dds) >= 10) >= 5
  dds <- DESeq(dds)[keep,]
}

#
common <- intersect(rownames(y), rownames(dds))
res <- results(dds)[common,]

alpha2 <- .05

getNumberFDR <- function(alpha, m) {
  tab <- table(sig=factor(mcols(y)[common,m] < alpha, levels=c("FALSE","TRUE")),
               de=res$padj < alpha2 &
                 (sign(res$log2FoldChange) == sign(mcols(y)[common,"log2FC"])))
  c(number=sum(tab[2,]), FDR=prop.table(tab,1)[2,1])
}

alphas <- c(.001,.0025,.005,.01,.025,.05,.075,.1,.15,.2)
dat <- rbind(data.frame(t(sapply(alphas, getNumberFDR, "qvalue")),method="Swish"),
             data.frame(t(sapply(alphas, getNumberFDR, "wilcox")),method="Wilcoxon"))
alphas2 <- c(.01, .05, .1)
dat2 <- rbind(data.frame(t(sapply(alphas2, getNumberFDR, "qvalue")),method="Swish"),
              data.frame(t(sapply(alphas2, getNumberFDR, "wilcox")),method="Wilcoxon"))
dat2$shape <- ifelse(dat2$FDR < alphas2, 19, 21)
ggplot(dat, aes(FDR, number, col=method)) + geom_line() +
  ylim(0,max(dat$number)) + xlim(0,.2) +
  geom_point(data=dat2, aes(shape=shape), size=5, fill="white") +
  scale_shape_identity() +
  theme_bw()

###

alphas <- c(.001,.0025,.005,.01,.025,.05,.075,.1,.15,.2)
niter <- 20
swish.res <- array(dim=c(10,2,niter))
wilcox.res <- array(dim=c(10,2,niter))
colnames(wilcox.res) <- colnames(swish.res) <- c("number","FDR")
for (i in 1:niter) {
  cat(i)
  set.seed(i)
  n.sub <- 30
  smp.idx <- c(sample(which(y$condition==1),n.sub),
               sample(which(y$condition==2),n.sub))
  set.seed(i)
  wilcox <- getWilcoxStat(n.cts[,smp.idx], condit=y$condition[smp.idx])
  set.seed(i)
  y2 <- swish(y[,smp.idx], x="condition", quiet=TRUE)
  mcols(y)$wilcox <- wilcox
  mcols(y)$qvalue <- mcols(y2)$qvalue
  swish.res[,,i] <- t(sapply(alphas, getNumberFDR, "qvalue"))
  wilcox.res[,,i] <- t(sapply(alphas, getNumberFDR, "wilcox"))
}
dat <- rbind(data.frame(apply(swish.res, 1:2, mean), method="swish"),
             data.frame(apply(wilcox.res, 1:2, mean), method="Wilcoxon"))
dat$FDR[is.nan(dat$FDR)] <- 0
dat2 <- dat[alphas %in% alphas2,]
dat2$shape <- ifelse(dat2$FDR < alphas2, 19, 21)

pdf(file=paste0("scrna_number_fdr_",n.sub,".pdf"), width=5, height=5)
ggplot(dat, aes(FDR, number, col=method)) + geom_line() +
  ylim(0,max(dat$number)) + xlim(0,max(dat$FDR)) +
  geom_point(data=dat2, aes(shape=shape), size=5, fill="white") +
  scale_shape_identity() +
  theme_bw() + ggtitle(paste0("n=",n.sub," vs ",n.sub,", averaged over 20 iterations") )
dev.off()

