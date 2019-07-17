#' Bayesian model selection using loo
#'
#' @export
#' @param cntfile Expression quantity file (RData format: Use 'save' and 'load').
#' @param outfile Output file name to store ELPD_loo results (RDS format).
#' @param nCores Number of cores.
#' @param seed Seed number.
#' @return A list of ELPD_loo results returned by loo::compare
#'
run_model_comparison <- function(cntfile, nCores=NULL, seed=NULL, outfile=NULL) {
  if(is.null(nCores)) {
    nCores <- parallel::detectCores()
  }
  if(is.null(seed)) {
    seed <- 1004
  }

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
        fit <- fit_count_models(y, exposure, nCores, seed)
        res <- compare_count_models(fit)
        results[[gname[gg]]] <- res
      }, error = function(err) {
        cat(sprintf("Error while fitting %s\n", gname[gg]))
      })
    }
  }

  if(!is.null(outfile)) {
    saveRDS(results, file = outfile)
  }
  return(results)
}
