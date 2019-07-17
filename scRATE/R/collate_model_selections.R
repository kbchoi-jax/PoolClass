#' Bayesian linear regression with Stan
#'
#' @export
#' @param loomfile A expression quantity file (loom format).
#' @param loo_results Leave-one-out ELPD result file.
#' @param attr_name Name of the row attribute in loomfile for storing best model calls.
#' @return best_model_calls Best model calls are also stored in the input loomfile.
#'
collate_model_selections <- function(loomfile, loo_results, attr_name) {
  results <- readRDS(loo_results)
  gsurv <- names(results)

  bestmodel <- c()
  for (g in gsurv) {
    res <- results[[g]]
    if (rownames(res)[1] == 'model1') {
      bestmodel <- c(bestmodel, 1)
    } else if (rownames(res)[1] == 'model2') {
      if (abs(res['model1',][['elpd_diff']]) < 2 * res['model1',][['se_diff']]) {
        bestmodel <- c(bestmodel, 1)
      } else {
        bestmodel <- c(bestmodel, 2)
      }
    } else if (rownames(res)[1] == 'model3') {
      if (abs(res['model1',][['elpd_diff']]) < 2 * res['model1',][['se_diff']]) {
        bestmodel <- c(bestmodel, 1)
      } else if (abs(res['model2',][['elpd_diff']]) < 2 * res['model2',][['se_diff']]) {
        bestmodel <- c(bestmodel, 2)
      } else {
        bestmodel <- c(bestmodel, 3)
      }
    } else if (rownames(res)[1] == 'model4') {
      if (abs(res['model1',][['elpd_diff']]) < 2 * res['model1',][['se_diff']]) {
        bestmodel <- c(bestmodel, 1)
      } else if (abs(res['model2',][['elpd_diff']]) < 2 * res['model2',][['se_diff']]) {
        bestmodel <- c(bestmodel, 2)
      } else if (abs(res['model3',][['elpd_diff']]) < 2 * res['model3',][['se_diff']]) {
        bestmodel <- c(bestmodel, 3)
      } else {
        bestmodel <- c(bestmodel, 4)
      }
    }
  }
  cat(sprintf('%5d Poisson genes\n', sum(bestmodel==1)))
  cat(sprintf('%5d Negative Binomial genes\n', sum(bestmodel==2)))
  cat(sprintf('%5d Zero-Inflated Poisson genes\n', sum(bestmodel==3)))
  cat(sprintf('%5d Zero-Inflated Neg. Binomial genes\n', sum(bestmodel==4)))
  cat(sprintf('%5d Total\n', length(results)))

  ds <- connect(loomfile, mode='r+')
  gname <- ds$row.attrs$GeneID[]
  best_idx <- match(gsurv, gname)
  bestmodel_fullsize <- rep(0, length(gname))
  bestmodel_fullsize[best_idx] <- bestmodel
  ra <- vector(mode="list", length=1)
  names(ra) <- attr_name
  ra[[attr_name]] <- bestmodel_fullsize
  ds$add.row.attribute(ra, overwrite=TRUE)
  ds$close_all()
  return(bestmodel_fullsize)
}
