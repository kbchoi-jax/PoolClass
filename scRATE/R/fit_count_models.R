#' Bayesian model selection for scRNA-seq count data
#'
#' @export
#' @param y Numeric vector of UMI counts.
#' @param exposure Numeric vector of cell sizes (total UMI counts per cell).
#' @param nCores Number of cores.
#' @param seed Seed number.
#' @param model2fit A specific model to fit (1:P, 2:NB, 3:ZIP, 4:ZINB, NULL:All)
#' @return A list of `stanfit` returned by model(s)
#'
fit_count_models <- function(y, exposure, nCores=NULL, seed=NULL, model2fit=NULL, brms4zi=FALSE) {

  if(is.null(nCores)) {
    nCores <- min(4, parallel::detectCores())
  }
  if(is.null(seed)) {
    seed <- 1004
  }

  gexpr <- data.frame(y, exposure)
  fitting <- list()

  if(is.null(model2fit) | model2fit==1) {
    fitting[["P"]] <- stan_glm(y ~ 1,
                               family = poisson,
                               offset = exposure,
                               data = gexpr,
                               cores = nCores,
                               seed = seed,
                               refresh = 0)
  }

  if(is.null(model2fit) | model2fit==2) {
    fitting[["NB"]] <- stan_glm(y ~ 1,
                                family = neg_binomial_2,
                                offset = exposure,
                                data = gexpr,
                                cores = nCores,
                                seed = seed,
                                refresh = 0)
  }

  if(is.null(model2fit) | model2fit==3) {
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
                  refresh = 500)
    } else {
      fit_3 <- rstan::sampling(stanmodels$zip,
                              data=list(N = length(y),
                                        Y = y,
                                        offset = exposure,
                                        offset_zi = exposure,
                                        prior_only = 0,
                                        df = myprior_3_values[1],
                                        loc = myprior_3_values[2],
                                        scale = myprior_3_values[3]),
                              cores = nCores,
                              seed = seed)
    }
    fitting[["ZIP"]] <- fit_3
  }

  if(is.null(model2fit) | model2fit==4) {
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
                  refresh = 500)
    } else {
      fit_4 <- rstan::sampling(stanmodels$zinb,
                              data=list(N = length(y),
                                        Y = y,
                                        offset = exposure,
                                        offset_zi = exposure,
                                        prior_only = 0,
                                        df = myprior_4_values[1],
                                        loc = myprior_4_values[2],
                                        scale = myprior_4_values[3]),
                              cores = nCores,
                              seed = seed)
    }
    fitting[["ZINB"]] <- fit_4
  }

  return(fitting)

}
