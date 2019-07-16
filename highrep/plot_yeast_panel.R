
wd <- "/proj/milovelab/zhu/projects/inf-uncertainty/"
setwd(wd)

source("/proj/milovelab/zhu/projects/inf-uncertainty/inf-uncertainty-2/RCodes/yeast/eval.highrep.R")
rep <- 100
load("data/yeast.rda")

fpr <- NULL
rvs <- NULL
InfRVcat <- NULL

sss <- c(5, 10, 15, 20)
for (n.sub in sss){
  load(paste0("data/yeast_res_", n.sub, ".rda"))
  length <- length(res[[1]]$SAMseq)
  fpr[[n.sub]] <- sapply(res, FUN = falsePos)/length 
  rvs[[n.sub]] <- rowMeans(sapply(res, FUN = function(x) as.vector(x$infRV)))
  brks <- c(quantile(rvs[[n.sub]], probs = seq(0, 1, by = 1/3), na.rm = TRUE))
  InfRV_cat <- cut(rvs[[n.sub]], brks)
 

  fprmat <- t(fpr[[n.sub]])
  res.bypanel <- NULL
  
  pdf(file = paste0("/proj/milovelab/zhu/projects/inf-uncertainty/plots/yeast_bypanel_zoom_", n.sub, ".pdf" ))
	par(mfrow=c(2,2))
		boxplot(fprmat, names = c("swish", "SAMseq", "limma"), ylim = c(0, 0.01), main = "Overall")
  
	  for (grp in levels(InfRV_cat)){
		idx <- which(InfRV_cat == grp)
		length.sub <- length(idx)
		res.bypanel[[grp]] <- sapply(res, FUN = falsePosByPanel, group.idx = idx)/length.sub 
		res.bypanel.t <-t(res.bypanel[[grp]])
		boxplot(res.bypanel.t, names = c("swish", "SAMseq", "limma"), ylim = c(0, 0.01), main = grp)
	  }
  
  dev.off()
  
  
  pdf(file = paste0("/proj/milovelab/zhu/projects/inf-uncertainty/plots/yeast_bypanel_", n.sub, ".pdf" ))
	par(mfrow=c(2,2))
		boxplot(fprmat, names = c("swish", "SAMseq", "limma"), ylim = c(0, 1), main = "Overall")
  
	  for (grp in levels(InfRV_cat)){
		idx <- which(InfRV_cat == grp)
		length.sub <- length(idx)
		res.bypanel[[grp]] <- sapply(res, FUN = falsePosByPanel, group.idx = idx)/length.sub  
		res.bypanel.t <-t(res.bypanel[[grp]])
		boxplot(res.bypanel.t, names = c("swish", "SAMseq", "limma"), ylim = c(0, 1), main = grp)
	  }
  
  dev.off()
  
  
  InfRVcat[[n.sub]] <- InfRVcat 

} 



# Arabidopsis #

plotarab <- function(lvl = "Gene", n.sub){
	if (lvl == "Gene") {
		length <- length(res.gene[[1]]$SAMseq)
		fpr <- t(sapply(res.gene, FUN = falsePos)/length)
		
		rvs <- rowMeans(sapply(res.gene, FUN = function(x) as.vector(x$infRV)))
		brks <- c(quantile(rvs, probs = seq(0, 1, by = 1/3), na.rm = TRUE))
		InfRV_cat <- cut(rvs, brks)

		pdf(file = paste0("/proj/milovelab/zhu/projects/inf-uncertainty/plots/Arabidopsis_Gene_Panel_zoom_", n.sub, ".pdf" ))
			par(mfrow=c(2,2))
			boxplot(fpr, names = c("swish", "SAMseq", "limma"), ylim = c(0, 0.01))
			
			  for (grp in levels(InfRV_cat)){
				idx <- which(InfRV_cat == grp)
				length.sub <- length(idx)
				res.bypanel <- t(sapply(res.gene, FUN = falsePosByPanel, group.idx = idx)/length.sub)
				boxplot(res.bypanel, names = c("swish", "SAMseq", "limma"), ylim = c(0, 0.01), main = grp)
			  }			
		
		dev.off()
		
		pdf(file = paste0("/proj/milovelab/zhu/projects/inf-uncertainty/plots/Arabidopsis_Gene_Panel_", n.sub, ".pdf" ))
			par(mfrow=c(2,2))	
			boxplot(fpr, names = c("swish", "SAMseq", "limma"), ylim = c(0, 1))
			
			  for (grp in levels(InfRV_cat)){
				idx <- which(InfRV_cat == grp)
				length.sub <- length(idx)
				res.bypanel <- t(sapply(res.gene, FUN = falsePosByPanel, group.idx = idx)/length.sub)
				boxplot(res.bypanel, names = c("swish", "SAMseq", "limma"), ylim = c(0, 1), main = grp)
			  }			
		
		dev.off()
	} else {
		length <- length(res.txp[[1]]$SAMseq)
		fpr <- t(sapply(res.txp, FUN = falsePos)/length)
		rvs <- rowMeans(sapply(res.txp, FUN = function(x) as.vector(x$infRV)))
		brks <- c(quantile(rvs, probs = seq(0, 1, by = 1/3), na.rm = TRUE))
		InfRV_cat <- cut(rvs, brks)

		pdf(file = paste0("/proj/milovelab/zhu/projects/inf-uncertainty/plots/Arabidopsis_Txp_Panel_zoom_", n.sub, ".pdf" ))
			par(mfrow=c(2,2))
			boxplot(fpr, names = c("swish", "SAMseq", "limma"), ylim = c(0, 0.01))
			
			  for (grp in levels(InfRV_cat)){
				idx <- which(InfRV_cat == grp)
				length.sub <- length(idx)
				res.bypanel <- t(sapply(res.txp, FUN = falsePosByPanel, group.idx = idx)/length.sub)
				boxplot(res.bypanel, names = c("swish", "SAMseq", "limma"), ylim = c(0, 0.01), main = grp)
			  }					
			
		dev.off()
		
		pdf(file = paste0("/proj/milovelab/zhu/projects/inf-uncertainty/plots/Arabidopsis_Txp_Panel_", n.sub, ".pdf" ))
			par(mfrow=c(2,2))
			boxplot(fpr, names = c("swish", "SAMseq", "limma"), ylim = c(0, 1))
			
			  for (grp in levels(InfRV_cat)){
				idx <- which(InfRV_cat == grp)
				length.sub <- length(idx)
				res.bypanel <- t(sapply(res.txp, FUN = falsePosByPanel, group.idx = idx)/length.sub)
				boxplot(res.bypanel, names = c("swish", "SAMseq", "limma"), ylim = c(0, 1), main = grp)
			  }					
			
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


