library(ensembldb)
#ensDbFromGtf("Homo_sapiens.GRCh38.95.gtf.gz")
txdb <- EnsDb("Homo_sapiens.GRCh38.95.sqlite")
tx2gene0 <- transcripts(txdb, return="DataFrame")[,c("tx_name","gene_id")]

library(splatter)

params <- newSplatParams()
# LFC centered at 2
params <- setParam(params, "de.facLoc", log(2) * 3)
params <- setParam(params, "de.facScale", log(2) * 1)
params <- setParam(params, "de.prob", 0.1)
sim <- splatSimulate(params,
                     group.prob=c(0.5,0.5),
                     method="groups",
                     nGenes=35583,
                     seed=1)
#counts(sim) <- t(t(counts(sim))/colSums(counts(sim)))*mean(colSums(counts(sim)))
sim <- sim[,order(sim$Group)]
mcols(sim)$de <- with(mcols(sim), DEFacGroup2/DEFacGroup1 != 1)
table(mcols(sim)$de)
table(sim$Group)
sim <- sim[,c(1:20,81:100)]
#log.cts <- log(counts(sim)[mcols(sim)$de,]+1)
#pheatmap::pheatmap(log.cts[order(rowVars(log.cts),decreasing=TRUE)[1:100],],
#                   annotation_col=as.data.frame(colData(sim)["Group"]),
#                   show_colnames=FALSE, show_rownames=FALSE)
#with(mcols(sim), plot(DEFacGroup1, DEFacGroup2, log="xy"))
#with(mcols(sim), hist(log2(DEFacGroup2/DEFacGroup1), col="grey"))
#with(mcols(sim)[mcols(sim)$de,], hist(log2(DEFacGroup2/DEFacGroup1), col="grey", breaks=40))
save(sim, file="sim.rda")

library(Biostrings)
txps <- Biostrings::readDNAStringSet("Homo_sapiens.GRCh38.cdna.all.fa.gz")
tx2gene <- tx2gene0
names(txps) <- sub("\\..*","",names(txps))
tx2gene <- tx2gene[tx2gene$tx_name %in% names(txps),]
tx2gene <- tx2gene[!duplicated(tx2gene$gene_id),]
txps <- txps[names(txps) %in% tx2gene$tx_name]
txps <- subseq(txps, ifelse(width(txps) >= 400, width(txps)-400+1, 1), width(txps))
# 35,583 txps
#Biostrings::writeXStringSet(txps, file="Homo_sapiens.GRCh38.cdna_400bp.fa")

txps <- Biostrings::readDNAStringSet("Homo_sapiens.GRCh38.cdna_400bp.fa")

load("sim.rda")
library(splatter)
countmat <- counts(sim)

library(polyester)
simulate_experiment_countmat(fasta="Homo_sapiens.GRCh38.cdna_400bp.fa",
                             readmat=countmat, 
                             outdir="out",
                             paired=FALSE,
                             strand_specific=TRUE,
                             seed=1)

for (i in 1:40) {
  cat(i)
  ii <- if (i < 10) paste0("0",i) else i
  reads1 <- Biostrings::readDNAStringSet(paste0("out/sample_",ii,".fasta"))
  set.seed(1)
  idx <- sample(length(reads1))
  reads1 <- reads1[idx]
  file1 <- paste0("out/shuffle_",ii,".fa")
  Biostrings::writeXStringSet(reads1, file=file1)
  R.utils::gzip(file1, overwrite=TRUE)
}

########

library(ensembldb)
txdb <- EnsDb("Homo_sapiens.GRCh38.95.sqlite")
tx2gene0 <- transcripts(txdb, return="DataFrame")[,c("tx_name","gene_id")]
txps <- Biostrings::readDNAStringSet("Homo_sapiens.GRCh38.cdna_400bp.fa")
load("sim.rda")
library(splatter)
countmat <- counts(sim)




library(tximport)
files <- file.path("quants",list.files("quants"),"quant.sf")
txi <- tximport(files, type="salmon", tx2gene=tx2gene0, ignoreTxVersion=TRUE, varReduce=TRUE)
est <- txi$counts
rownames(est) <- sub("\\..*","",rownames(est))
gene.ids <- tx2gene0[match(names(txps), tx2gene0$tx_name),"gene_id"]
idx <- match(gene.ids, rownames(est))
est <- est[idx,,drop=FALSE]
InfRV <- pmax(txi$variance[idx,,drop=FALSE] - est, 0) / (est + 5) + .01

library(genefilter)
dat <- data.frame(
  estimated=rowSums(est),
  count=rowSums(countmat),
  tstat=rowttests(est, factor(sim$Group))$statistic,
  InfRV=pmax(rowMeans(InfRV[,1:20]),rowMeans(InfRV[,21:40])),
  status=ifelse(rowData(sim)$de,"DE","null"),
  label=rownames(countmat), stringsAsFactors=FALSE)
dat$tstat <- ifelse(is.na(dat$tstat), 0, dat$tstat)

library(ggplot2)
sub.dat <- dat[dat$count > 5 & abs(dat$tstat) > 5 & dat$status == "null" & dat$InfRV > .2,]

library(cowplot)
g1 <- ggplot(dat, aes(log10(estimated+5), tstat, col=status)) +
  geom_point(size=2) +
  geom_abline(slope=0, intercept=0) +
  geom_point(data=sub.dat, size=5, shape=1, show.legend=FALSE) +
  ylab("t-statistic")
g2 <- ggplot(dat, aes(log10(count+5), tstat, col=status)) +
  geom_point(size=2) +
  geom_abline(slope=0, intercept=0) +
  geom_point(data=sub.dat, size=5, shape=1, show.legend=FALSE) +
  ylab("t-statistic")

png(width=1200, height=600, res=120, file="splatter_ma.png")
plot_grid(g1, g2)
dev.off()

genes <- as.numeric(sub("Gene","",sub.dat$label))

## gene <- genes[1]
## rbind(true=counts(sim)[gene,],
##       est=round(txi$counts[idx,][gene,]),
##       InfRV=round(InfRV[gene,],1))
## boxplot(txi$counts[idx,][gene,] ~ sim$Group)
## wilcox.test(txi$counts[idx,][gene,] ~ sim$Group)

###

library(tximeta)
coldata <- data.frame(files, names=1:40, condition=factor(rep(1:2,each=20)))
se <- tximeta(coldata)
y0 <- summarizeToGene(se)
y <- y0
mcols(y)$tx_name <- tx2gene0$tx_name[match(rownames(y), tx2gene0$gene_id)]
# re-order
y <- y[ match(names(txps), mcols(y)$tx_name), ]
mcols(y)$InfRV <- pmax(rowMeans(InfRV[,1:20]),rowMeans(InfRV[,21:40]))
mcols(y)$de <- mcols(sim)$de

getWilcoxStat <- function(cts, condit, paired = FALSE){
  dims <- dim(cts)
  cts.notie <- cts + 0.1 * runif(dims[1]*dims[2])
  idx <- 1:dims[1]
  wilcox.p <- sapply(idx, function(i) {
    wilcox.test(cts.notie[i, ] ~ condit)$p.value
  })
  wilcox.padj <- p.adjust(wilcox.p, method = "BH")
  names(wilcox.padj) <- rownames(cts)
  wilcox.padj
}

library(fishpond)
sfFun <- function(m) {
  DESeq2::estimateSizeFactorsForMatrix(m, geoMeans=exp(rowSums(log(m) * as.numeric(m > 0))/ncol(m)))
}
y <- scaleInfReps(y, lengthCorrect=FALSE, sfFun=sfFun, minCount=5, minN=3)
y <- labelKeep(y, minCount=5, minN=3)
y <- swish(y, x="condition")
res <- mcols(y)[mcols(y)$keep,]

prop.table(with(res, table(qvalue < .05, de)), 1)
prop.table(with(res, table(qvalue < .05, de)), 2)

library(DESeq2)
dds <- DESeqDataSetFromMatrix(round(assay(y)[mcols(y)$keep,]),
                              data.frame(x=seq_len(ncol(y))), ~1)
dds <- estimateSizeFactors(dds, type="poscounts")
n.cts <- counts(dds, normalized=TRUE)
res.wilcox <- getWilcoxStat(n.cts, condit=y$condition)
res$wilcox <- res.wilcox 

prop.table(with(res, table(wilcox < .05, de)), 1)
prop.table(with(res, table(wilcox < .05, de)), 2)

sub.res <- res[match(names(txps[genes]), res$tx_name),]
sub.res <- sub.res[order(sub.res$wilcox),]

pdf(file="splatter_gene_padj.pdf", width=8)
library(rafalib)
bigpar()
plot(-log10(sub.res$wilcox), ylim=c(.01,12), cex=2, pch=20,
     xlab="14 selected genes, ordered by Wilcoxon",
     ylab="-log10 adjusted p-value")
points(-log10(sub.res$qvalue), col="blue", cex=2, pch=20)
abline(h=-log10(.05), lty=2)
text(4, -log10(.05), "FDR cutoff of 0.05", pos=3, cex=1.5)
legend("topright", c("Wilcoxon","swish"), col=c("black","blue"), pch=20, inset=.05, cex=1.5)
dev.off()

pdf(file="splatter_inf_reps.pdf", width=8, height=5)
par(mfrow=c(3,3), mar=c(1,4,3,1))
for (i in c(3,4,7, 8,10,11, 12,13,14)) {
  plotInfReps(y, sub.res$gene_id[i], x="condition", xaxis=FALSE)
}
dev.off()
