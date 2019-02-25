rgb2 <- function(x,y,z) rgb(x,y,z,maxColorValue=255)


myicobra <- function(cd, lvl, label=FALSE, split=TRUE, overall = FALSE, ylim, xlim.right) {
  if (split) {
    cp <- calculate_performance(cd,
                                binary_truth="status",
                                aspect=c("fdrtpr","fdrtprcurve"),
                                thrs=c(.01,.05,.1),
                                splv = "InfRV",
                                maxsplit = Inf)  
    cobraplot <- prepare_data_for_plot(cp, colorscheme=cols, facetted = TRUE,
                                       incloverall = overall)
  } else {
    cp <- calculate_performance(cd,
                                binary_truth="status",
                                aspect=c("fdrtpr","fdrtprcurve"),
                                thrs=c(.01,.05,.1))   
    cobraplot <- prepare_data_for_plot(cp, colorscheme=cols)
  }
  yrng <- ylim
  xrng <- c(0, xlim.right)
  #xrng <- c(0,max(.21,max(fdrtpr(cp)$FDR)))
  plt <- plot_fdrtprcurve(cobraplot, plottype="points",
                          xaxisrange=xrng,
                          yaxisrange=yrng,
                          #yaxisrange=c(0,max(fdrtpr(cp)$TPR)),
                          stripsize=8,
                          title=paste0(lvl,"-level, n=",n.sub," vs ",n.sub))
  if (label)
    plt <- plt + geom_text(aes(label=method,color=method,vjust=-1))
  plt
}

catPlot <- function(res, keeprule = NULL, error = TRUE) {
  tops <- c(100,150,200,250,300,350,400)
  df <- do.call(rbind, lapply(res, function(x) {
    tab <- do.call(cbind, lapply(2:length(x$lfcs), function(i) {
      sapply(tops, function(top) {
	    if (is.null(keeprule)) {
		length(intersect(order(-abs(x$lfcs$true))[1:top],
                         order(-abs(x$lfcs[[i]]))[1:top]))/top
		} else {
        length(intersect(order(-abs(x$lfcs$true[keeprule]))[1:top],
                         order(-abs(x$lfcs[[i]][keeprule]))[1:top]))/top
		}
      })
    }))
    rownames(tab) <- tops
    colnames(tab) <- names(x$lfcs)[-1]
    df <- melt(tab)
    names(df) <- c("top", "method", "concordance")
    df
  }))
  df <- fixdf(df)
  fun.y <- mean
  fun.ymax <- function(x) mean(x) + sd(x)/sqrt(length(x))
  fun.ymin <- function(x) mean(x) - sd(x)/sqrt(length(x))
  g <- ggplot(df, aes(x=top, y=concordance, color=method, group=method, linetype = method)) +
    stat_summary(fun.y=fun.y,geom="point", size=pt.size) +
    stat_summary(fun.y=fun.y,geom="line", size=line.size) +
    xlab("top genes by LFC") + ylab("concordance") +
	scale_linetype_manual(values=lt)
  g <- g +
    stat_summary(data=subset(df, method=="apeglm"),
                 aes(x=top, y=concordance, color=method),
                 fun.y=fun.y, geom="point", size=pt.size) +
    stat_summary(data=subset(df, method=="apeglm"),
                 aes(x=top, y=concordance, color=method, linetype = method),
                 fun.y=fun.y, geom="line", size=line.size) 
  if (error) {
    g <- g + stat_summary(fun.ymax=fun.ymax,fun.ymin=fun.ymin,geom="errorbar")
  }
  g
}

cols <- c("DESeq2"="black",
          "EBSeq"=rgb2(230,159,0),
          "limma"=rgb2(86,180,233),
          "SAMseq"=rgb2(0,158,115),
          "sleuth"=rgb2(0,114,178),
          "swish"=rgb2(213,94,0))
		  
