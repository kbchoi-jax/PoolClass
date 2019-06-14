#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=3) {
  stop("Three arguments must be supplied (<loomfile> <loo_results> <attr_name>)", call.=FALSE)
}

loomfile <- args[1]
loo_results <- args[2]
attr_name <- args[3]

library(loomR)

results <- readRDS(loo_results)
#results <- results[!is.na(rowSums(results)), ]
gsurv <- names(results)

bestmodel <- c()
for (g in gsurv) {
  res <- results[[g]]
  bestmodel <- c(bestmodel, as.integer(gsub('model', '', rownames(res)[1])))
}
cat(sprintf('%5d Poisson genes\n', sum(bestmodel==1)))
cat(sprintf('%5d Negative Binomial genes\n', sum(bestmodel==2)))
cat(sprintf('%5d Zero-Inflated Poisson genes\n', sum(bestmodel==3)))
cat(sprintf('%5d Zero-Inflated Neg. Binomial genes\n', sum(bestmodel==4)))
cat(sprintf('%5d Total\n', length(results)))

ds <- connect(loomfile, mode='r+')
gname <- ds$row.attrs$GeneID[]
best_idx <- match(gsurv, gname)
bestmodel_fullsize <- rep(0, length(gname))
bestmodel_fullsize[best_idx] <- bestmodel
ra <- vector(mode="list", length=1)
names(ra) <- attr_name
ra[[attr_name]] <- bestmodel_fullsize
ds$add.row.attribute(ra, overwrite=TRUE)
ds$close_all()
