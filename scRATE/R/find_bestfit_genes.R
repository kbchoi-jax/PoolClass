#' Bayesian model selection for scRNA-seq count data
#'
#' @export
#' @param fit_file A file name that contains model fitting result reside
#' @param margin A multiplier for standard deviation (SD) in leave-one-out ELPD for calling models
#' @param loomfile A expression quantity file (loom format)
#' @param attr_name Name of the row attribute in loomfile for storing best model calls
#' @param verbose Whether to print out overall model calls
#' @return best_model_calls Genes that Best fits to a given model are also stored in the input loomfile
#'
find_bestfit_genes <- function(fit_file, model2check, margin=2, loomfile=NULL, attr_name=NULL, verbose=FALSE) {
                                     
  results <- readRDS(fit_file)
  gsurv <- names(results)

  best_genes <- c()
  for (g in gsurv) {
    if (is_bestfit(results[[g]][['elpd_loo']], model2check=model2check, margin=margin)) {
      best_genes <- c(best_genes, g)
    }
  }

  if(verbose==TRUE) {
    if (model2check == 1) {
      cat(sprintf('%5d Poisson genes found\n', length(best_genes)))
    } else if (model2check == 2) {
      cat(sprintf('%5d Negative Binomial genes found\n', length(best_genes)))
    } else if (model2check == 3) {
      cat(sprintf('%5d Zero-Inflated Poisson genes found\n', length(best_genes)))
    } else if (model2check == 4) {
      cat(sprintf('%5d Zero-Inflated Neg. Binomial genes found\n', length(best_genes)))
    } else {
      cat(('Error: model2check should be either of 1, 2, 3, or 4.'))
    }
  }

  if(is.null(loomfile)) {
    return(best_genes)
  } else {
    ds <- connect(loomfile, mode='r+')
    gname <- ds$row.attrs$GeneID[]
    gene_idx <- match(best_genes, gname)
    best_genes_fullsize <- rep(0, length(gname))
    best_genes_fullsize[gene_idx] <- model2check
    ra <- vector(mode="list", length=1)
    if(is.null(attr_name)) {
      attr_name <- sprintf('FitBest_%d', model2check)
    }
    names(ra) <- attr_name
    ra[[attr_name]] <- best_genes_fullsize
    ds$add.row.attribute(ra, overwrite=TRUE)
    ds$close_all()
    return(best_genes_fullsize)
  }

}
