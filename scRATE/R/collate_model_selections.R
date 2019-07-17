#' Bayesian linear regression with Stan
#'
#' @export
#' @param loomfile Numeric vector of UMI counts.
#' @param loo_results Numeric vector of cell sizes (total UMI counts per cell).
#' @param attr_name Name of the .
#' @return A list of `stanfit` returned by four models: Poisson, Negative-Binomial, Zero-Inflated Poisson, & Zero-Inflated Negative-Binomial
#'
collate_model_selections <- function(loomfile, loo_results, attr_name) {
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
