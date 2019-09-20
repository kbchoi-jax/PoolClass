#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=4) {
  stop("Four arguments must be supplied (<cntfile> <outfile> <cores> <seed>).n", call.=FALSE)
}

library(scRATE)

cntfile <- args[1]
outfile <- args[2]
nCores <- as.integer(args[3])
seed <- as.integer(args[4])

load(cntfile)  # This will load 'cntmat', 'gsurv', and 'csize'
gname <- rownames(cntmat)
num_genes <- length(gsurv)
exposure <- log(csize)

results.elpd_loo <- list()
results.mean_par <- list()
for (gg in c(1:num_genes)) {
  if(gsurv[gg]) {
    y <- round(unlist(cntmat[gg,]))
    cat(sprintf("\nTesting %s\n", gname[gg]))
    tryCatch({
      model_fit <- fit_count_models(y, exposure, nCores, seed, adapt_delta = 0.95, brms4zi=FALSE)
      results.elpd_loo[[gname[gg]]] <- compare_count_models(model_fit)
      results.mean_par[[gname[gg]]] <- get_model_params(model_fit)
    }, error = function(err) {
      cat(sprintf("Error while fitting %s\n", gname[gg]))
    })
  }
}

saveRDS(results.elpd_loo, file = outfile)
saveRDS(results.mean_par, file = gsub('elpd_loo', 'mean_par', outfile))
