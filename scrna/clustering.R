seed <- 131

library(Seurat)
library(loomR)
library(tximeta)
library(SummarizedExperiment)
library(org.Mm.eg.db)

load("data/mouse_900_se.rda")
load("data/mouse_900_quant.rda")

mnh <- CreateSeuratObject(mat, min.cells = 1, min.features = 200, project = "mouse_900")

idx.mito <- seqnames(se) == "chrM"

mito.percent <- apply(mat[as.vector(idx.mito, mode = "any"), ], 2, sum)/apply(mat, 2, sum)

mnh <- AddMetaData(mnh, metadata = mito.percent, col.name = "percent.mito")

png(file = "plots/vlnplt.png", width = 1500, height = 500)
VlnPlot(object = mnh, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
dev.off()

quantile(mnh@meta.data$nGene, probs = seq(0,1,0.05))

quantile(mnh@meta.data$nUMI, probs = seq(0,1,0.05))

quantile(mnh@meta.data$percent.mito, probs = seq(0,1,0.05))


png(file = "plots/qc_scatter.png", width = 1000, height = 500)
par(mfrow = c(1, 2))
GenePlot(object = mnh, gene1 = "nUMI", gene2 = "percent.mito")
GenePlot(object = mnh, gene1 = "nUMI", gene2 = "nGene")
dev.off()

mnh <- FilterCells(object = mnh, subset.names = c("nGene", "percent.mito"), 
                    low.thresholds = c(1500, -Inf), high.thresholds = c(6500, 0.04))

# 835 cells left

mnh <- NormalizeData(mnh, normalization.method = "LogNormalize", 
                     scale.factor = 10000)

mnh <- FindVariableGenes(object = mnh, mean.function = ExpMean, dispersion.function = LogVMR, 
                         x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5, num.bin = 20)

# 2732 features selected
mnh <- ScaleData(object = mnh, vars.to.regress = c("nUMI", "percent.mito"))

mnh <- RunPCA(object = mnh, pc.genes = mnh@var.genes, pcs.print = 1:10, pcs.compute = 30,
              genes.print = 10)

PrintPCA(object = mnh, pcs.print = 1:30, genes.print = 5, use.full = FALSE)
VizPCA(object = mnh, pcs.use = 1:20)
PCAPlot(object = mnh, dim.1 = 1, dim.2 = 2)

mnh <- JackStraw(object = mnh, num.replicate = 100, display.progress = FALSE)

png(file = "plots/pca_jstraw.png", width = 1000, height = 1000)
JackStrawPlot(object = mnh)
dev.off()

png(file = "plots/pca_elbow.png", width = 500, height = 500)
PCElbowPlot(object = mnh, num.pc = 30)
dev.off()

mnh  <- FindClusters(object = mnh, reduction.type = "pca", dims.use = 1:10, 
                     resolution = 0.6, print.output = 0, save.SNN = TRUE, random.seed = seed)

mnh <- RunTSNE(object = mnh, dims.use = 1:10, do.fast = TRUE, k.seed = seed)

png(file = "plots/tsne_0.6.png", width = 600, height = 600)
TSNEPlot(object = mnh, do.label = TRUE, pt.size = 1.5, label.size= 6)
dev.off()

mnh <- StashIdent(object = mnh, save.name = "ClusterNames_0.6")

mnh <- FindClusters(object = mnh, resolution = 1.2)

png(file = "plots/tsne_1.2.png", width = 600, height = 600)
TSNEPlot(object = mnh,  do.label = TRUE, pt.size = 2, label.size = 10)
dev.off()

## plot tsne with known markers ##
k <- c("Neurod6","Gad2","Eomes","Mki67","Olig2","Trem2","Igfbp7")
k.ensm.nv <- mapIds(org.Mm.eg.db, column="ENSEMBL", keys=k, keytype="SYMBOL")

feature.ensm <- rownames(mnh@scale.data)
feature.ensm.nv <- gsub( "[.].*$", "", feature.ensm)

k.ensm <- feature.ensm[match( k.ensm.nv, feature.ensm.nv)]

png(file = "plots/feature_tsne_eomes_mki67.png", width = 1000, height = 500)
FeaturePlot(object = mnh , features.plot = k.ensm[3:4], cols.use = c("grey", "red"), 
    reduction.use = "tsne", pt.size = 2, no.legend = TRUE)
dev.off()

saveRDS(mnh, file = "data/mnh.rds")