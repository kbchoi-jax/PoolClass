#' Bayesian model selection for scRNA-seq count data
#'
#' @export
#' @param fit_list A list of model fitting results
#' @param margin A multiplier for standard deviation (SD) in leave-one-out ELPD for calling models
#' @param loomfile A expression quantity file (loom format)
#' @param attr_name Name of the row attribute in loomfile for storing best model calls
#' @param verbose Whether to print out overall model calls
#' @return best_model_calls Best model calls are also stored in the input loomfile
#'
collate_model_selections <- function(fit_list, margin=2, loomfile=NULL, attr_name=NULL, verbose=FALSE) {
                                     
  gsurv <- names(fit_list)
  bestmodel <- c()
  for (g in gsurv) {
    bestmodel <- c(bestmodel, select_model(fit_list[[g]][['elpd_loo']], margin=margin))
  }

  if(verbose==TRUE) {
    cat(sprintf('%5d Poisson genes\n', sum(bestmodel==1)))
    cat(sprintf('%5d Negative Binomial genes\n', sum(bestmodel==2)))
    cat(sprintf('%5d Zero-Inflated Poisson genes\n', sum(bestmodel==3)))
    cat(sprintf('%5d Zero-Inflated Neg. Binomial genes\n', sum(bestmodel==4)))
    cat(sprintf('%5d Total\n', length(fit_list)))
  }

  if(is.null(loomfile)) {
    return(bestmodel)
  } else {
    ds <- connect(loomfile, mode='r+')
    gname <- ds$row.attrs$GeneID[]
    best_idx <- match(gsurv, gname)
    bestmodel_fullsize <- rep(0, length(gname))
    bestmodel_fullsize[best_idx] <- bestmodel
    ra <- vector(mode="list", length=1)
    if(is.null(attr_name)) {
      attr_name <- 'BestModel'
    }
    names(ra) <- attr_name
    ra[[attr_name]] <- bestmodel_fullsize
    ds$add.row.attribute(ra, overwrite=TRUE)
    ds$close_all()
    return(bestmodel_fullsize)
  }

}
