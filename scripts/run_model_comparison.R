#!/usr/bin/env Rscript
source('compare_count_models.R')

#nCores <- parallel::detectCores()
nCores <- 4

args = commandArgs(trailingOnly=TRUE)
if (length(args)!=2) {
  stop("Two arguments must be supplied (<cntfile> <outfile>).n", call.=FALSE)
}

cntfile <- args[1]
outfile <- args[2]

load(cntfile)  # This will load 'cntmat', 'gsurv', and 'csize'
gname <- rownames(cntmat)
num_genes <- length(gsurv)
exposure <- log(csize)

m <- readRDS('zinb.rds')
results <- data.frame()
processed <- c()

for (gg in c(1:num_genes)) {
  if(gsurv[gg]) {
    y <- round(unlist(cntmat[gg,]))
    cat(sprintf("\nTesting %s\n", gname[gg]))
    tryCatch({
      res <- compare_count_models(y, exposure, m, nCores)
      results <- rbind(results, res)
      processed <- c(processed, gname[gg])
    }, error = function(err) {
      cat(sprintf("Error while fitting %s\n", gname[gg]))
    })
  }
}

rownames(results) <- processed
colnames(results) <- rep(c("elpd_diff", "se"), 3)
saveRDS(results, file = outfile)
