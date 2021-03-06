args <- commandArgs(trailingOnly = TRUE)
gene <- args[1] == "gene"
batch <- args[2] == "2batch"

suppressPackageStartupMessages(library(SummarizedExperiment))
dir <- "quants"
subdir <- rep(c("bias","unif"), each=12)
nms_bias <- as.vector(t(outer(1:6,1:2,function(x,y) paste0(x,"_",y))))
nms_unif <- as.vector(t(outer(7:12,1:2,function(x,y) paste0(x,"_",y))))
nms <- c(nms_bias, nms_unif)
files <- file.path(dir, subdir, nms, "quant.sf")
all(file.exists(files))

## condition <- factor(sub(".*_(.)","\\1",nms))
## batch.var <- factor(rep(1:2,each=12))
## names <- paste0(nms, "_b", batch.var)
## coldata <- data.frame(files=files, names=names,
##                       condition=condition,
##                       batch=batch.var,
##                       stringsAsFactors=FALSE)
## se <- tximeta::tximeta(coldata)
## save(se, file="se_full.rda")

library(tximeta)
library(tximport)
load("se_full.rda")
load("quants/simulate.rda")
if (gene) {
  stopifnot(ncol(se) == 24)
  gse <- summarizeToGene(se)
  v <- tximport(files, type="salmon", tx2gene=txdf[,2:1], varReduce=TRUE)
  assays(gse)[["variance"]] <- v$variance
  idx <- if (batch) 1:24 else 1:12
  txi <- tximport(files[idx], type="salmon", tx2gene=txdf[,2:1],
                  countsFromAbundance="lengthScaledTPM")
  gse <- gse[,idx]
} else {
  stopifnot(ncol(se) == 24)
  v <- tximport(files, type="salmon", txOut=TRUE, varReduce=TRUE)
  assays(se)[["variance"]] <- v$variance
  idx <- if (batch) 1:24 else 1:12
  txi <- tximport(files[idx], type="salmon", txOut=TRUE,
                countsFromAbundance="lengthScaledTPM")
  se <- se[,idx]
}

###

myicobra <- function(cd, lvl="Gene", n.sub) {
  cp <- calculate_performance(cd,
                              binary_truth="status",
                              aspect=c("fdrtpr","fdrtprcurve"),
                              thrs=c(.01,.05,.1),
                              splv="InfRV", maxsplit=Inf)
  cobraplot <- prepare_data_for_plot(cp,
                                     colorscheme=cols,
                                     facetted = TRUE,
                                     incloverall = TRUE)
  yrng <- c(.4, 1)
  xrng <- c(0,max(.2, max(fdrtpr(cp)$FDR)))
  plot_fdrtprcurve(cobraplot, plottype="points",
                   pointsize=2.5,
                   xaxisrange=xrng,
                   yaxisrange=yrng,
                   title=paste0(lvl,"-level, n=",n.sub," vs ",n.sub))
}

###

if (gene) {
  infVar <- assays(gse)[["variance"]]
  mu <- assays(gse)[["counts"]]
} else {
  infVar <- assays(se)[["variance"]]
  mu <- assays(se)[["counts"]]
}
infRV <- pmax(infVar - mu, 0)/(mu + 5) + .01
meanInfRV <- rowMeans(infRV)
#hist(log10(meanInfRV)[meanInfRV > .01],breaks=100)
#plot(rowMeans(mu) + 5, meanInfRV, log="xy", cex=.1)

if (gene) {
  de <- rownames(gse) %in% union(dge.genes, dte.genes)
} else {
  iso.any <- iso.dtu | iso.dte | iso.dge
  de <- rownames(se) %in% names(iso.dtu)[iso.any]
}

library(fishpond)

if (gene) {
  y <- gse
} else {
  y <- se
}
mcols(y)$de <- de
mcols(y)$InfRV <- meanInfRV
y <- scaleInfReps(y)
y <- labelKeep(y)
y <- y[mcols(y)$keep,]

qs <- quantile(mcols(y)$InfRV, c(.33,.66,1))
mcols(y)$cut <- cut(mcols(y)$InfRV, c(0,qs))
#with(mcols(y), hist(log10(InfRV)[InfRV > .01],breaks=100))

###########
## swish ##
###########

set.seed(1)
if (batch) {
  y <- swish(y, x="condition", cov="batch")
} else {
  y <- swish(y, x="condition")
}
padj <- data.frame(row.names=rownames(y))
padj$swish <- mcols(y)$qvalue

############
## SAMseq ##
############

samseq <- function(cts, condition) {
  x <- round(cts)
  out <- capture.output({
    samfit <- SAMseq(x, condition, resp.type="Two class unpaired", fdr.output=1)
  })
  sam.padj <- rep(1,nrow(x))
  names(sam.padj) <- seq_len(nrow(x))
  if (samfit$siggenes.table$ngenes.up > 0) {
    up <- samfit$siggenes.table$genes.up
    if (!is.matrix(up)) {
      idx <- as.numeric(up["Gene Name"])
      sam.padj[idx] <- 1/100 * as.numeric(up["q-value(%)"])
    } else {
      idx <- as.numeric(up[,"Gene Name"])
      sam.padj[idx] <- 1/100 * as.numeric(up[,"q-value(%)"])
    }
  }
  if (samfit$siggenes.table$ngenes.lo > 0) {
    lo <- samfit$siggenes.table$genes.lo
    if (!is.matrix(lo)) {
      idx <- as.numeric(lo["Gene Name"])
      sam.padj[idx] <- 1/100 * as.numeric(lo["q-value(%)"])
    } else {
      idx <- as.numeric(lo[,"Gene Name"])
      sam.padj[idx] <- 1/100 * as.numeric(lo[,"q-value(%)"])
    }
  }
  sam.padj <- pmin(sam.padj, 1)
}

library(samr)
cts <- txi$counts[rownames(y),]
if (batch) {
  cts.rbe  <- exp(limma::removeBatchEffect(log(cts+.5), batch=y$batch))
  set.seed(1)
  sam.padj <- samseq(cts.rbe, y$condition)
} else {
  set.seed(1)
  sam.padj <- samseq(cts, y$condition)
}
padj$SAMseq <- sam.padj

##############
## wilcoxon ##
##############

getWilcoxStat <- function(cts, condition){
  cts.notie <- cts + 0.1 * runif(nrow(cts)*ncol(cts))
  wilcox.p <- sapply(seq_len(nrow(cts)), function(i) {
    wilcox.test(cts.notie[i, ] ~ condition)$p.value
  })
  p.adjust(wilcox.p, method="BH")
}

if (batch) {
  wilcox.padj <- samseq(cts.rbe, y$condition)
} else {
  wilcox.padj <- samseq(cts, y$condition)
}
padj$Wilcoxon <- wilcox.padj

###########
## limma ##
###########

edgerPrep <- function(cts, condition = NULL, cov = NULL) {
    y <- DGEList(cts)
    y <- y[filterByExpr(y),]
    y <- calcNormFactors(y)
    if (is.null(cov)){
      design <- model.matrix(~condition)
    } else {
      design <- model.matrix(~condition + cov)
    }
    list(y, design)
}

limmavoom <- function(cts, condition = NULL, cov = NULL) {
  out <- edgerPrep(cts, condition = condition, cov = cov)
  y <- out[[1]]; design <- out[[2]]
  v <- voom(y, design)
  fit <- lmFit(v, design)
  fit <- eBayes(fit)
  topTable(fit, coef=2, n=nrow(y), sort="none")
}

library(limma)
library(edgeR)

if (batch) {
  tt <- limmavoom(cts, condition=y$condition, cov=y$batch)
} else {
  tt <- limmavoom(cts, condition=y$condition)
}
padj$limma <- 1
padj[rownames(tt),"limma"] <- tt$adj.P.Val

############
## DESeq2 ##
############

deseq2 <- function(cts, condition = NULL, cov = NULL) {
  cts <- round(cts)
  suppressMessages({
    if (is.null(cov)) {
      samps <- data.frame(condition=condition)
      dds <- DESeqDataSetFromMatrix(cts, samps, ~condition)
    } else {
      samps <- data.frame(condition=condition, cov=cov)
      dds <- DESeqDataSetFromMatrix(cts, samps, ~condition + cov)
    }
    dds <- DESeq(dds, minReplicates=Inf, quiet=TRUE)
  })
  res <- results(dds, name="condition_2_vs_1")
  res$padj[is.na(res$padj)] <- 1
  res
}

library(DESeq2)

if (batch) {
  res <- deseq2(cts, condition=y$condition, cov=y$batch)
} else {
  res <- deseq2(cts, condition=y$condition)
}
padj$DESeq2 <- res$padj

###########
## EBSeq ##
###########

ebseq <- function(cts, condition, txdf, gene=TRUE) {
  x <- round(cts)
  sizes <- MedianNorm(x)
  if (gene) {
    suppressMessages({
      res <- EBTest(Data=x, Conditions=condition,
                    sizeFactors=sizes, maxround=5, Print=FALSE)
    })
  } else {
    suppressMessages({
      Ng <- GetNg(txdf[,2], txdf[,1])$IsoformNgTrun
      res <- EBTest(Data=x, NgVector=Ng, Conditions=condition,
                    sizeFactors=sizes, maxround=5, Print=FALSE)
    })
  }
  padj <- rep(1, nrow(x))
  names(padj) <- rownames(x)
  padj[match(rownames(res$PPMat), rownames(x))] <- res$PPMat[,"PPEE"]
  padj
}

library(EBSeq)

ebseq.file <- paste0("ebseq/ebseq_",ifelse(gene,"gene","txp"),"_",ifelse(batch,"2batch","1batch"),".rda")
if (!file.exists(ebseq.file)) {
  if (batch) {
    eb.padj <- ebseq(cts.rbe, condition=y$condition, txdf=txdf, gene=gene)
  } else {
    eb.padj <- ebseq(cts, condition=y$condition, txdf=txdf, gene=gene)
  }
  save(eb.padj, file=ebseq.file)
} else {
  load(ebseq.file)
}

padj$EBSeq <- eb.padj

############
## sleuth ##
############

sleuth.file <- paste0("sleuth_",
                      ifelse(gene,"gene","txp"),"_",
                      ifelse(batch,"strat_12","nostrat_6"),
                      ".rda")
load(sleuth.file)
padj$sleuth <- sleuth.res$qval[match(rownames(padj), sleuth.res$target_id)]

############
## iCOBRA ##
############

library(iCOBRA)
rgb2 <- function(x,y,z) rgb(x,y,z,maxColorValue=255)
cols <- c("DESeq2"="black",
          "EBSeq"=rgb2(230,159,0),
          "limma"=rgb2(86,180,233),
          "SAMseq"=rgb2(0,158,115),
          "sleuth"=rgb2(0,114,178),
          "swish"=rgb2(213,94,0),
          "Wilcoxon"=rgb2(204,121,167))

if (gene) {
  status <- as.numeric(rownames(padj) %in% union(dge.genes, dte.genes))
} else {
  iso.any <- iso.dtu | iso.dte | iso.dge
  status <- as.numeric(rownames(padj) %in% names(iso.dtu)[iso.any])
}
truth <- data.frame(status, InfRV=mcols(y)$cut,
                    row.names=rownames(padj))
filename <- paste0("icobra_",
                   ifelse(gene,"gene","txp"),"_",
                   ifelse(batch,"2batch","1batch"))
#save(padj, truth, gene, batch, file=paste0(filename, ".rda"))
load(paste0(filename,".rda"))
cd <- COBRAData(padj=padj, truth=truth)
n.sub <- if (batch) 12 else 6
pdf(paste0(filename,".pdf"))

padj$sleuth <- NULL
cd <- COBRAData(padj=padj, truth=truth)
pdf(paste0(filename,"_no_sleuth.pdf"))

padj$Wilcoxon <- NULL
cd <- COBRAData(padj=padj, truth=truth)
pdf(paste0(filename,"_no_wilcox.pdf"))

print(myicobra(cd, if (gene) "Gene" else "Transcript", n.sub))
dev.off()
