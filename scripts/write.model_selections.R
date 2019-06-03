#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=3) {
  stop("Three arguments must be supplied (<loomfile> <loo_results> <attr_name)", call.=FALSE)
}

loomfile <- args[1]
loo_results <- args[2]
attr_name <- args[3]

library(loomR)

ds <- connect(loomfile, mode='r+')
results <- readRDS(loo_results)
results <- results[!is.na(rowSums(results)), ]

best1 <- results[,1]<0 & results[,3]<0
best2 <- results[,1]>0 & results[,5]<0
best3 <- results[,3]>0 & results[,5]>0
cat(sprintf('%5d Poisson genes\n', sum(best1)))
cat(sprintf('%5d Negative Binomial genes\n', sum(best2)))
cat(sprintf('%5d Zero-Inflated Neg. Binomial genes\n', sum(best3)))
cat(sprintf('%5d Total\n', sum(best1, best2, best3)))

gname <- ds$row.attrs$GeneID[]
best1idx <- match(rownames(results[best1,]), gname)
best2idx <- match(rownames(results[best2,]), gname)
best3idx <- match(rownames(results[best3,]), gname)
bestmodel <- rep('', length(gname))
bestmodel[best1idx] <- 'P'
bestmodel[best2idx] <- 'NB'
bestmodel[best3idx] <- 'ZINB'
ra <- vector(mode="list", length=1)
names(ra) <- attr_name
ra[[attr_name]] <- bestmodel
ds$add.row.attribute(ra, overwrite=TRUE)
ds$close_all()
