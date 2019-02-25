library(sleuth)

suppressPackageStartupMessages(library(tximeta))
suppressPackageStartupMessages(library(SummarizedExperiment))
suppressPackageStartupMessages(library(GenomicFeatures))

load("data/simulate.rda")
tx2gene <- txdf[,2:1]
colnames(tx2gene) <- c("target_id","gene_id")

n.sub <- 12
dir <- "data/quants"
subdir <- rep(c("kallisto_bias","kallisto_unif"), each=12)pwd

nms_bias <- as.vector(t(outer(1:6,1:2,function(x,y) paste0(x,"_",y))))
nms_unif <- as.vector(t(outer(7:12,1:2,function(x,y) paste0(x,"_",y))))
names <- c(nms_bias, nms_unif)
files_dir <- file.path(dir, subdir, names)
				
condition <- factor(sub(".*_(.)","\\1",names))
batch <- factor(rep(1:2,each=12))
names2 <- paste(rep(c("bias", "unif"), each = 12), names)

s2c <- data.frame(	 	 
  sample = names2,
  condition = condition,  
  batch = batch,
  path = files_dir, 
  stringsAsFactors=FALSE
)	 	 

st.gene <- system.time({
  so <- sleuth_prep(s2c, ~ condition + batch, target_mapping = tx2gene, aggregation_column = "gene_id", num_cores = 1, gene_mode = TRUE)
  so <- sleuth_fit(so)
  so <- sleuth_fit(so, ~ batch, 'reduced')
  so <- sleuth_lrt(so, 'reduced', 'full')
  sleuth.res <- sleuth_results(so, 'reduced:full', test_type = 'lrt', pval_aggregate = FALSE)
})

save(sleuth.res, st.gene, file=paste0("data/sleuth_gene_strat_",n.sub,".rda"))

rm(list=ls())


n.sub <- 6
dir <- "data/quants"
subdir <- rep(c("kallisto_bias"), each=12)
names <- as.vector(t(outer(1:6,1:2,function(x,y) paste0(x,"_",y))))
files_dir <- file.path(dir, subdir, names)

load("data/simulate.rda")
tx2gene <- txdf[,2:1]
colnames(tx2gene) <- c("target_id","gene_id")

condition <- factor(sub(".*_(.)","\\1",names))

s2c <- data.frame(	 	 
  sample = names,	 	 
  path = files_dir, 
  condition, 
  stringsAsFactors=FALSE
)

st2.gene <- system.time({
  so <- sleuth_prep(s2c, ~ condition, target_mapping = tx2gene, aggregation_column = "gene_id", num_cores = 1, gene_mode = TRUE)
  so <- sleuth_fit(so)
  so <- sleuth_fit(so, ~ 1, 'reduced')
  so <- sleuth_lrt(so, 'reduced', 'full')
  sleuth.res <- sleuth_results(so, 'reduced:full', test_type = 'lrt', pval_aggregate = FALSE)
})

save(sleuth.res, st2.gene, file=paste0("data/sleuth_gene_nostrat_",n.sub,".rda"))

