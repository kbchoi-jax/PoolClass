#' Bayesian model selection for scRNA-seq count data
#'
#' @export
#' @param model_fit A list of four model fits
#' @return A list of mean parameter values for each model
#
get_model_params <- function(model_fit) {
  result <- list()
  result[['P']]    <- model_fit[['P']]$coefficients
  result[['NB']]   <- colMeans(as.data.frame(model_fit[['NB']]))
  result[['ZIP']]  <- colMeans(as.matrix(model_fit[['ZIP']], pars = c("b_Intercept", "b_zi_Intercept")))
  result[['ZINB']] <- colMeans(as.matrix(model_fit[['ZINB']], pars = c("b_Intercept", "b_zi_Intercept", "shape")))    
  return(result)
}
