#' Bayesian model selection for scRNA-seq count data
#'
#' @export
#' @param loo_results ELPD_loo for models
#' @param model2check Model to check whether a gene bestfits to it (Either of 1, 2, 3, or 4)
#' @param margin A multiplier for confidence interval
#' @return Whether a given gene best fits to the specified model (boolean)
#
is_bestfit <- function(loo_results, model2check, margin=2) {

  m3idx <- which(rownames(loo_results) == 'model_fit$ZIP')
  if(length(m3idx)) { rownames(loo_results)[m3idx] <- 'model3' }
  m4idx <- which(rownames(loo_results) == 'model_fit$ZINB')
  if(length(m4idx)) { rownames(loo_results)[m4idx] <- 'model4' }

  if (model2check == 1) {
    if (rownames(loo_results)[1] == 'model1') {
      if (abs(loo_results['model2',][['elpd_diff']]) > margin * loo_results['model2',][['se_diff']] &&
          abs(loo_results['model3',][['elpd_diff']]) > margin * loo_results['model3',][['se_diff']] &&
          abs(loo_results['model4',][['elpd_diff']]) > margin * loo_results['model4',][['se_diff']]) 
      {
        return(TRUE)
      } else {
        return(FALSE)
      }
    } else {
      return(FALSE)
    }
  } else if (model2check == 2) {
    if (rownames(loo_results)[1] == 'model2') {
      if (abs(loo_results['model3',][['elpd_diff']]) > margin * loo_results['model3',][['se_diff']] &&
          abs(loo_results['model4',][['elpd_diff']]) > margin * loo_results['model4',][['se_diff']] &&
          abs(loo_results['model1',][['elpd_diff']]) > margin * loo_results['model1',][['se_diff']]) 
      {
        return(TRUE)
      } else {
        return(FALSE)
      }
    } else {
      return(FALSE)
    }
  } else if (model2check == 3) {
    if (rownames(loo_results)[1] == 'model3') {
      if (abs(loo_results['model4',][['elpd_diff']]) > margin * loo_results['model4',][['se_diff']] &&
          abs(loo_results['model1',][['elpd_diff']]) > margin * loo_results['model1',][['se_diff']] &&
          abs(loo_results['model2',][['elpd_diff']]) > margin * loo_results['model2',][['se_diff']]) 
      {
        return(TRUE)
      } else {
        return(FALSE)
      }
    } else {
      return(FALSE)
    }
  } else if (model2check == 4) {
    if (rownames(loo_results)[1] == 'model4') {
      if (abs(loo_results['model1',][['elpd_diff']]) > margin * loo_results['model1',][['se_diff']] &&
          abs(loo_results['model2',][['elpd_diff']]) > margin * loo_results['model2',][['se_diff']] &&
          abs(loo_results['model3',][['elpd_diff']]) > margin * loo_results['model3',][['se_diff']]) 
      {
        return(TRUE)
      } else {
        return(FALSE)
      }
    } else {
      return(FALSE)
    }
  }
        
}
