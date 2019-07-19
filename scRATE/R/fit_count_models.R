#' Bayesian model selection for scRNA-seq count data
#'
#' @export
#' @param y Numeric vector of UMI counts.
#' @param exposure Numeric vector of cell sizes (total UMI counts per cell).
#' @param nCores Number of cores.
#' @param seed Seed number.
#' @return A list of `stanfit` returned by four models: Poisson, Negative-Binomial, Zero-Inflated Poisson, & Zero-Inflated Negative-Binomial
#'
fit_count_models <- function(y, exposure, nCores=NULL, seed=NULL, brms4zi=FALSE) {
  if(is.null(nCores)) {
    nCores <- parallel::detectCores()
  }
  if(is.null(seed)) {
    seed <- 1004
  }

  gexpr <- data.frame(y, exposure)

  fit_1 <- stan_glm(y ~ 1,
                    family=poisson,
                    offset=exposure,
                    data=gexpr,
                    cores = nCores,
                    seed=seed,
                    refresh=0)

  fit_2 <- stan_glm(y ~ 1,
                    family=neg_binomial_2,
                    offset=exposure,
                    data=gexpr,
                    cores = nCores,
                    seed=seed,
                    refresh=0)

  myprior_3 <- get_prior(bf(y ~ 1 + offset(exposure), zi ~ 1 + offset(exposure)),
                         data = gexpr,
                         family = zero_inflated_poisson())
  myprior_3_values <- eval(parse(text=gsub("student_t", "c", myprior_3$prior[1])))
  if(brms4zi) {
    fit_3 <- brm(bf(y ~ 1 + offset(exposure), zi ~ 1 + offset(exposure)),
                 family = zero_inflated_poisson(),
                 data = gexpr,
                 prior = myprior_3,
                 cores = nCores,
                 seed = seed,
                 refresh=500)
  } else {
    fit_3 <- rstan::sampling(stanmodels$zip,
                            data=list(N=length(y),
                                      Y=y,
                                      offset=exposure,
                                      offset_zi=exposure,
                                      prior_only=0,
                                      df=myprior_3_values[1],
                                      loc=myprior_3_values[2],
                                      scale=myprior_3_values[3]),
                            cores = nCores,
                            seed=seed)
  }

  myprior_4 <- get_prior(bf(y ~ 1 + offset(exposure), zi ~ 1 + offset(exposure)),
                         data = gexpr,
                         family = zero_inflated_negbinomial())
  myprior_4_values <- eval(parse(text=gsub("student_t", "c", myprior_4$prior[1])))
  if (brms4zi) {
    fit_4 <- brm(bf(y ~ 1 + offset(exposure), zi ~ 1 + offset(exposure)),
                 family = zero_inflated_negbinomial(),
                 data = gexpr,
                 prior = myprior_4,
                 cores = nCores,
                 seed = seed,
                 refresh=500)
  } else {
    fit_4 <- rstan::sampling(stanmodels$zinb,
                            data=list(N=length(y),
                                      Y=y,
                                      offset=exposure,
                                      offset_zi=exposure,
                                      prior_only=0,
                                      df=myprior_4_values[1],
                                      loc=myprior_4_values[2],
                                      scale=myprior_4_values[3]),
                            cores = nCores,
                            seed=seed)
  }

  return(list("P"=fit_1, "NB"=fit_2, "ZIP"=fit_3, "ZINB"=fit_4))

}
