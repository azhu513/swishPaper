
wd <- "/proj/milovelab/zhu/projects/inf-uncertainty/"
setwd(wd)

source("/proj/milovelab/zhu/projects/inf-uncertainty/inf-uncertainty-2/RCodes/yeast/eval.highrep.R")
rep <- 100

fpr <- NULL

sss <- c(5, 10, 15, 20)
for (n.sub in sss){
  load(paste0("data/yeast_res_", n.sub, ".rda"))
  fpr[[n.sub]] <- sapply(res, FUN = falsePos) 

  fprmat <- t(fpr[[n.sub]])
   
  pdf(file = paste0("/proj/milovelab/zhu/projects/inf-uncertainty/plots/yeast_overall_", n.sub, ".pdf" ))
	boxplot(fprmat, names = c("swish", "SAMseq", "limma"), ylim = c(0, 1), main = "Overall")
  dev.off()
  
} 

# Arabidopsis #

plotarab <- function(lvl = "Gene", n.sub){
	if (lvl == "Gene") {
		fpr <- t(sapply(res.gene, FUN = falsePos))
		pdf(file = paste0("/proj/milovelab/zhu/projects/inf-uncertainty/plots/Arabidopsis_Gene_Overall_", n.sub, ".pdf" ))	
			boxplot(fpr, names = c("swish", "SAMseq", "limma"), ylim = c(0, 1))				
		dev.off()
	} else {
		fpr <- t(sapply(res.txp, FUN = falsePos))
		
		pdf(file = paste0("/proj/milovelab/zhu/projects/inf-uncertainty/plots/Arabidopsis_Txp_Overall_", n.sub, ".pdf" ))
			boxplot(fpr, names = c("swish", "SAMseq", "limma"), ylim = c(0, 1))
		dev.off()
	}
}


n.sub <- 8
load(paste0("data/arabidopsis_res_", n.sub, ".rda"))	
plotarab(lvl="Gene", n.sub)
plotarab(lvl="Txp", n.sub)

library(rlist)
n.sub <- 5
load(paste0("data/arabidopsis_res_", n.sub, ".rda"))	
res.gene <- list.remove(res.gene, 15)
plotarab(lvl="Gene", n.sub)
plotarab(lvl="Txp", n.sub)


