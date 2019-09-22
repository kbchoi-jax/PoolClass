#' Bayesian model selection for scRNA-seq count data
#'
#' @export
#' @param loo_dir A folder name in which leave-one-out ELPD result files reside
#' @param globstr Search string (wildcard supported) for loo result files (in RDS format)
#' @param loo_outfile Name of the file to save collated leave-one-out ELPD results (optional)
#' @return A list of collated results if 'loo_outfile' is not given
#'
collate_results <- function(loo_dir, globstr='_scrate*', loo_outfile=NULL) {
                                     
  flist <- Sys.glob(file.path(loo_dir, globstr))
  flist <- sort(flist, decreasing = FALSE)
  results <- c()
  for (i in 1:length(flist)) {
    results <- c(results, readRDS(flist[i]))
  }

  if (!is.null(loo_outfile)) {
    saveRDS(results, file=loo_outfile)
  } else {
    return(results)
  }

}
