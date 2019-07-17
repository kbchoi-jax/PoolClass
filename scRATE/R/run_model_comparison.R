#' Bayesian linear regression with Stan
#'
#' @export
#' @param cntfile Numeric vector of UMI counts.
#' @param outfile Numeric vector of cell sizes (total UMI counts per cell).
#' @param nCores Number of cores.
#' @param seed Seed number.
#' @param ... Arguments passed to `rstan::sampling` (e.g. iter, chains).
#' @return A list of `stanfit` returned by four models: Poisson, Negative-Binomial, Zero-Inflated Poisson, & Zero-Inflated Negative-Binomial
#'
run_model_comparison <- function(cntfile, outfile, nCores, seed) {
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
        res <- compare_count_models(y, exposure, nCores, seed)
        results[[gname[gg]]] <- res
      }, error = function(err) {
        cat(sprintf("Error while fitting %s\n", gname[gg]))
      })
    }
  }

  saveRDS(results, file = outfile)
}
