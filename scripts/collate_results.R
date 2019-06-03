#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=2) {
  stop("Two argument must be supplied (<indir> <outfile>)", call.=FALSE)
}

indir <- args[1]
outfile <- args[2]

flist <- Sys.glob(sprintf('%s/_poolclass*', indir))
flist <- sort(flist, decreasing = FALSE)
results <- data.frame()
for (k in 1:length(flist)) {
  newres <- readRDS(flist[k])
  results <- rbind(results, newres)
}
saveRDS(results, file = outfile)
results <- results[!is.na(rowSums(results)), ]
best1 <- sum(results[,1]<0 & results[,3]<0)
best2 <- sum(results[,1]>0 & results[,5]<0)
best3 <- sum(results[,3]>0 & results[,5]>0)
cat(sprintf('%5d Poisson genes\n', sum(best1)))
cat(sprintf('%5d Negative Binomial genes\n', sum(best2)))
cat(sprintf('%5d Zero-Inflated Neg. Binomial genes\n', sum(best3)))
cat(sprintf('%5d Total\n', sum(best1, best2, best3)))
