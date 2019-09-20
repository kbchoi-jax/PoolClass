#' Bayesian model selection for scRNA-seq count data
#'
#' @export
#' @param loo_dir A folder name in which leave-one-out ELPD result files reside
#' @param globstr Search string (wildcard supported) for loo result files (in RDS format)
#' @param margin A multiplier for standard deviation (SD) in leave-one-out ELPD for calling models.
#' @param loo_outfile Name of the file to save collated leave-one-out ELPD results (optional)
#' @param loomfile A expression quantity file (loom format)
#' @param attr_name Name of the row attribute in loomfile for storing best model calls
#' @param verbose Whether to print out overall model calls
#' @return best_model_calls Best model calls are also stored in the input loomfile
#'
collate_model_selections <- function(loo_dir, globstr='_scrate*', margin=2, loo_outfile=NULL,
                                     loomfile=NULL, attr_name=NULL, verbose=FALSE) {
  flist <- Sys.glob(file.path(loo_dir, globstr))
  flist <- sort(flist, decreasing = FALSE)
  results <- c()
  for (i in 1:length(flist)) {
    results <- c(results, readRDS(flist[i]))
  }
  gsurv <- names(results)

  bestmodel <- c()
  for (g in gsurv) {
    res <- results[[g]]
    if (rownames(res)[1] == 'model1') {
      bestmodel <- c(bestmodel, 1)
    } else if (rownames(res)[1] == 'model2') {
      if (abs(res['model1',][['elpd_diff']]) < margin * res['model1',][['se_diff']]) {
        bestmodel <- c(bestmodel, 1)
      } else {
        bestmodel <- c(bestmodel, 2)
      }
    } else if (rownames(res)[1] == 'model3') {
      if (abs(res['model1',][['elpd_diff']]) < margin * res['model1',][['se_diff']]) {
        bestmodel <- c(bestmodel, 1)
      } else if (abs(res['model2',][['elpd_diff']]) < margin * res['model2',][['se_diff']]) {
        bestmodel <- c(bestmodel, 2)
      } else {
        bestmodel <- c(bestmodel, 3)
      }
    } else if (rownames(res)[1] == 'model4') {
      if (abs(res['model1',][['elpd_diff']]) < margin * res['model1',][['se_diff']]) {
        bestmodel <- c(bestmodel, 1)
      } else if (abs(res['model2',][['elpd_diff']]) < margin * res['model2',][['se_diff']]) {
        bestmodel <- c(bestmodel, 2)
      } else if (abs(res['model3',][['elpd_diff']]) < margin * res['model3',][['se_diff']]) {
        bestmodel <- c(bestmodel, 3)
      } else {
        bestmodel <- c(bestmodel, 4)
      }
    }
  }

  if(verbose==TRUE) {
    cat(sprintf('%5d Poisson genes\n', sum(bestmodel==1)))
    cat(sprintf('%5d Negative Binomial genes\n', sum(bestmodel==2)))
    cat(sprintf('%5d Zero-Inflated Poisson genes\n', sum(bestmodel==3)))
    cat(sprintf('%5d Zero-Inflated Neg. Binomial genes\n', sum(bestmodel==4)))
    cat(sprintf('%5d Total\n', length(results)))
  }

  if (!is.null(loo_outfile)) {
    saveRDS(results, file=loo_outfile)
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
