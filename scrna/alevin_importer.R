boot.cell.names <- read.table("prediction/alevin/quants_boot_rows.txt", strings=FALSE)[,1]
gene.names <- read.table("prediction/alevin/quants_mat_cols.txt", strings=FALSE)[,1]
cell.names <- read.table("prediction/alevin/quants_mat_rows.txt", strings=FALSE)[,1]

library(jsonlite)
jsonPath <- file.path("prediction/cmd_info.json")
cmd_info <- jsonlite::fromJSON(jsonPath)

num.cells <- length(cell.names)
num.genes <- length(gene.names)
num.bootstraps <- as.numeric(cmd_info$numCellBootstraps)

mat.file <- "prediction/alevin/quants_mat.gz"
mean.mat.file <- "prediction/alevin/quants_mean_mat.gz"
boot.mat.file <- "prediction/alevin/quants_boot_mat.gz"

mat <- matrix(nrow=num.genes, ncol=num.cells, dimnames=list(gene.names, cell.names))
con <- gzcon(file(mat.file, "rb"))
for (j in seq_len(num.cells)) {
  mat[,j] <- readBin(con, double(), endian = "little", n=num.genes)
}
close(con)

library(Matrix)
mat <- Matrix(mat)

num.non.zero <- matrix(nrow=num.cells, ncol=num.bootstraps)
con <- gzcon(file(boot.mat.file, "rb"))
for (j in seq_len(num.cells)) {
  if (j %% 50 == 0) cat(j,"")
  for (k in seq_len(num.bootstraps)) {
    vec <- readBin(con, double(), endian = "little", n=num.genes)
    num.non.zero[j,k] <- sum(vec > 0)
  }
}
close(con)

is <- lapply(seq_len(num.bootstraps), function(k) numeric(sum(num.non.zero[,k])))
xs <- is

# easy to build this one
js <- lapply(seq_len(num.bootstraps), function(k) rep(seq_len(num.cells), num.non.zero[,k]))

con <- gzcon(file(boot.mat.file, "rb"))
ptr <- rep(1, num.bootstraps)
for (j in seq_len(num.cells)) {
  if (j %% 50 == 0) cat(j,"")
  for (k in seq_len(num.bootstraps)) {
    vec <- readBin(con, double(), endian = "little", n=num.genes)
    n <- num.non.zero[j,k]
    stopifnot(n == sum(vec > 0))
    is[[k]][ptr[k]:(ptr[k]+n-1)] <- which(vec > 0)
    xs[[k]][ptr[k]:(ptr[k]+n-1)] <- vec[vec > 0]
    ptr[k] <- ptr[k] + n
  }
}
close(con)

bmat.list <- list()
for (k in seq_len(num.bootstraps)) {
  bmat <- sparseMatrix(i=is[[k]], j=js[[k]], x=xs[[k]])
  idx <- match(cell.names, boot.cell.names)
  bmat <- bmat[,idx]
  stopifnot(all(abs(colSums(mat) - colSums(bmat)) < 1))
  bmat.list[[k]] <- bmat
}

mat.d <- as.matrix(mat)

bmat.list.d <- list()
for (k in seq_along(bmat.list)){
    tempmat.d <- as.matrix(bmat.list[[k]])
    bmat.list.d[[k]] <- tempmat.d
}

save(mat, bmat.list, mat.d, bmat.list.d, is, xs, ptr, js, file = "data/mouse_900_quant.rda")

## tximeta ##
library("tximeta")
files <- file.path("prediction/alevin/quants_mat.gz")
file.exists(files)
coldata <- data.frame(files, names="neurons")
se <- tximeta(coldata, type="alevin")

save(se, file = "data/data.mouse_900_se.rda")

