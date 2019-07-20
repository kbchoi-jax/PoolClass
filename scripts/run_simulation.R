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

results <- list()

for (gg in c(1:num_genes)) {
  if(gsurv[gg]) {
    y <- round(unlist(cntmat[gg,]))
    cat(sprintf("\nTesting %s\n", gname[gg]))
    tryCatch({
      model_fit <- fit_count_models(y, exposure, nCores, seed, brms4zi = TRUE)
      elpd_loo <- compare_count_models(model_fit)
      # Generate simulated counts for each model
      yrep1 <- posterior_predict(model_fit$P)
      yrep2 <- posterior_predict(model_fit$NB)
      yrep3 <- posterior_predict(model_fit$ZIP)
      yrep4 <- posterior_predict(model_fit$ZINB)
      num_posterior_samp <- dim(yrep_1)[1]
      sim1 <- yrep1[sample(num_posterior_samp, 1),]
      sim2 <- yrep2[sample(num_posterior_samp, 1),]
      sim3 <- yrep3[sample(num_posterior_samp, 1),]
      sim4 <- yrep4[sample(num_posterior_samp, 1),]
      # Fit simulated counts with all four models
      model_fit_sim1 <- fit_count_models(sim1, exposure, nCores, seed)
      model_fit_sim2 <- fit_count_models(sim2, exposure, nCores, seed)
      model_fit_sim3 <- fit_count_models(sim3, exposure, nCores, seed)
      model_fit_sim4 <- fit_count_models(sim4, exposure, nCores, seed)
      elpd_loo_sim1 <- compare_count_models(model_fit_sim1)
      elpd_loo_sim2 <- compare_count_models(model_fit_sim2)
      elpd_loo_sim3 <- compare_count_models(model_fit_sim3)
      elpd_loo_sim4 <- compare_count_models(model_fit_sim4)
      res <- list(orig=elpd_loo,
                  sim1=elpd_loo_sim1,
                  sim2=elpd_loo_sim2,
                  sim3=elpd_loo_sim3,
                  sim4=elpd_loo_sim4)
      results[[gname[gg]]] <- res
    }, error = function(err) {
      cat(sprintf("Error while fitting %s\n", gname[gg]))
    })
  }
}

saveRDS(results, file = outfile)
