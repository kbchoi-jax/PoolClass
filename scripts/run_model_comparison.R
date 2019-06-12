#!/usr/bin/env Rscript
source('compare_count_models.R')

nCores <- min(4, parallel::detectCores())

args = commandArgs(trailingOnly=TRUE)
if (length(args)!=3) {
  stop("Three arguments must be supplied (<cntfile> <outfile> <seed>).n", call.=FALSE)
}

cntfile <- args[1]
outfile <- args[2]
seed <- as.integer(args[3])

load(cntfile)  # This will load 'cntmat', 'gsurv', and 'csize'
gname <- rownames(cntmat)
num_genes <- length(gsurv)
exposure <- log(csize)

m_zip  <- readRDS('zip.rds')
m_zinb <- readRDS('zinb.rds')
results <- list()
processed <- c()

for (gg in c(1:num_genes)) {
  if(gsurv[gg]) {
    y <- round(unlist(cntmat[gg,]))
    cat(sprintf("\nTesting %s\n", gname[gg]))
    tryCatch({
      res <- compare_count_models(y, exposure, m_zip, m_zinb, nCores, seed)
      results[[gname[gg]]] <- res
    }, error = function(err) {
      cat(sprintf("Error while fitting %s\n", gname[gg]))
    })
  }
}

saveRDS(results, file = outfile)
