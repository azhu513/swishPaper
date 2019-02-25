library(sleuth)
suppressPackageStartupMessages(library(tximeta))
suppressPackageStartupMessages(library(SummarizedExperiment))

n.sub <- 12
dir <- "data/quants"
subdir <- rep(c("kallisto_bias","kallisto_unif"), each=12)
nms_bias <- as.vector(t(outer(1:6,1:2,function(x,y) paste0(x,"_",y))))
nms_unif <- as.vector(t(outer(7:12,1:2,function(x,y) paste0(x,"_",y))))
names <- c(nms_bias, nms_unif)
files_dir <- file.path(dir, subdir, names)

condition <- factor(sub(".*_(.)","\\1",names))
batch <- factor(rep(1:2,each=12))
names2 <- paste(rep(c("bias", "unif"), each = 12), names)

s2c <- data.frame(	 	 
  sample = names2,	 	 
  path = files_dir, 
  condition, 
  batch, 
  stringsAsFactors=FALSE
)	 	 

st <- system.time({
  so <- sleuth_prep(s2c, ~ condition + batch, num_cores=1)
  so <- sleuth_fit(so)
  so <- sleuth_fit(so, ~ batch, 'reduced')
  so <- sleuth_lrt(so, 'reduced', 'full')
  sleuth.res <- sleuth_results(so, 'reduced:full', test_type = 'lrt')
})

save(sleuth.res, st, file=paste0("data/sleuth_txp_strat_",n.sub,".rda"))

rm(list=ls())

n.sub <- 6
dir <- "data/quants"
subdir <- rep("kallisto_bias", 12)
names <- as.vector(t(outer(1:6,1:2,function(x,y) paste0(x,"_",y))))
files_dir <- file.path(dir, subdir, names)

condition <- factor(sub(".*_(.)","\\1",names))
s2c <- data.frame(	 	 
  sample = names,	 	 
  path = files_dir, 
  condition, 
  stringsAsFactors=FALSE
)

st2 <- system.time({
  so <- sleuth_prep(s2c, ~ condition, num_cores=1)
  so <- sleuth_fit(so)
  so <- sleuth_fit(so, ~ 1, 'reduced')
  so <- sleuth_lrt(so, 'reduced', 'full')
  sleuth.res <- sleuth_results(so, 'reduced:full', test_type = 'lrt')
})

save(sleuth.res, st2, file=paste0("data/sleuth_txp_nostrat_",n.sub,".rda"))
