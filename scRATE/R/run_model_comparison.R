#!/usr/bin/env Rscript
#source('./compare_count_models.R')

run_model_comparison <- function(cntfile, outfile, zip_model, zinb_model, nCores, seed) {
  nCores <- min(4, parallel::detectCores())
  load(cntfile)  # This will load 'cntmat', 'gsurv', and 'csize'
  gname <- rownames(cntmat)
  num_genes <- length(gsurv)
  exposure <- log(csize)

  #zip_model  <- readRDS('zip.rds')
  #zinb_model <- readRDS('zinb.rds')
  results <- list()
  processed <- c()

  for (gg in c(1:num_genes)) {
    if(gsurv[gg]) {
      y <- round(unlist(cntmat[gg,]))
      cat(sprintf("\nTesting %s\n", gname[gg]))
      tryCatch({
        res <- compare_count_models(y, exposure, zip_model, zinb_model, nCores, seed)
        results[[gname[gg]]] <- res
      }, error = function(err) {
        cat(sprintf("Error while fitting %s\n", gname[gg]))
      })
    }
  }

  saveRDS(results, file = outfile)
}

collate_model_selection <- function(loomfile, loo_results, attr_name) {
  results <- readRDS(loo_results)
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
  best_idx <- match(rownames(results), gname)
  bestmodel_fullsize <- rep(0, length(gname))
  bestmodel_fullsize[best_idx] <- bestmodel
  ra <- vector(mode="list", length=1)
  names(ra) <- attr_name
  ra[[attr_name]] <- bestmodel_fullsize
  ds$add.row.attribute(ra, overwrite=TRUE)
  ds$close_all()
}
