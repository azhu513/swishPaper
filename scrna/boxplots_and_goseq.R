# first analysis: clusters 6 and 7
# third analysis: clusters 4 and 6

third <- FALSE
third <- TRUE

if (third) {
  load("mouse_900_sample_lower_filter_y_se_3rdAnalysis.rda")
  load("mouse_900_sample_lower_filter_DE_res_3rdAnalysis.rda")
} else {
  load("mouse_900_sample_lower_filter_y_se.rda")
  load("mouse_900_sample_lower_filter_DE_res.rda")
}

y <- y.swish
mcols(y)$wilcox <- res.wilcox

# tabulate overlap

with(mcols(y), table(wilcox=wilcox < .01, swish=qvalue < .01))
mcols(y)$wsig <- ifelse(mcols(y)$wilcox < .01, "W+", "W-")
mcols(y)$ssig <- ifelse(mcols(y)$qvalue < .01, "S+", "S-")
mcols(y)$sig <- factor(paste(mcols(y)$wsig, mcols(y)$ssig, sep="/"))
infRVList <- with(mcols(y), split(InfRV, sig))[c(3,2,4)]
names(infRVList)

sapply(infRVList, function(x) sum(x > .25))

### draw ggplot2 boxplots

boxfile <- if (third) "scrna_infrv_3.pdf" else "scrna_infrv.pdf"
main <- if (third) "DE genes: cluster 6 vs 4" else "DE genes: cluster 7 vs 6"
dat <- data.frame(InfRV=unlist(infRVList),
                  method=rep(names(infRVList),lengths(infRVList)))
dat$method <- relevel(dat$method, "W+/S-")
nums <- data.frame(method=names(infRVList),y=c(.25,.25,.25),
                   label=paste0("#gene = ",lengths(infRVList)))
library(ggplot2)
pdf(file=boxfile, width=4, height=4)
ggplot(dat, aes(x=method,y=InfRV)) +
  geom_boxplot(outlier.color=NA) + ylim(0,.25) +
  geom_jitter(width=0.3, height=0, alpha=0.2) +
  geom_label(data=nums, aes(x=method,y=y,label=label)) +
  ggtitle(main)
dev.off()

### goseq analysis

library(goseq)
geneListI <- as.integer(with(mcols(y), qvalue < .01))
names(geneListI) <- sub("\\..*","",rownames(y))
pwf <- nullp(geneListI, "mm10", "ensGene")
fit <- goseq(pwf, "mm10", "ensGene")
fit <- fit[fit$ontology == "BP" & fit$numInCat > 10,]
fit$padj <- p.adjust(fit$over_represented_pvalue)
table(fit$padj < .01)
tab <- fit[fit$padj < .01,c(1,4,5,6)]
tab$term <- sub("nucleobase-containing", "nuc. cont.", tab$term)

geneListI <- as.integer(with(mcols(y), wilcox < .01))
names(geneListI) <- sub("\\..*","",rownames(y))
pwf2 <- nullp(geneListI, "mm10", "ensGene")
fit2 <- goseq(pwf2, "mm10", "ensGene")
fit2 <- fit2[fit2$ontology == "BP" & fit2$numInCat > 10,]
fit2$padj <- p.adjust(fit2$over_represented_pvalue)
table(fit2$padj < .01)
tab2 <- fit2[fit2$padj < .01,c(1,4,5,6)]
tab2$term <- sub("microtubule", "micro.", tab2$term)
tab2$term <- sub("cytoskeleton", "cyto.", tab2$term)

tab$both <- ifelse(tab$category %in% tab2$category, "X", "")
tab <- tab[,c(1:3,5,4)]

tab2$both <- ifelse(tab2$category %in% tab$category, "X", "")
tab2 <- tab2[,c(1:3,5,4)]
