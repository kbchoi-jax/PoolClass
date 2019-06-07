library(rstan)
library(rstanarm)
library(rstantools)
library(brms)
library(loo)

compare_count_models <- function(y, exposure, zinb_model, nCores) {
  
  gexpr <- data.frame(y, exposure)
  
  fit_1 <- stan_glm(y ~ 1,
                    family=poisson,
                    offset=exposure,
                    data=gexpr,
                    cores = nCores,
                    #seed=SEED,
                    refresh=0)
  loo_1 <- loo(fit_1)
  
  fit_2 <- stan_glm(y ~ 1,
                    family=neg_binomial_2,
                    offset=exposure, 
                    data=gexpr,
                    cores = nCores,
                    #seed=SEED, 
                    refresh=0)
  loo_2 <- loo(fit_2)
  
  myprior <- get_prior(bf(y ~ 1 + offset(exposure), zi ~ 1 + offset(exposure)),
                       data = gexpr, 
                       family = zero_inflated_negbinomial())
  myprior_values <- eval(parse(text=gsub("student_t", "c", myprior$prior[1])))
  fit_3 <- sampling(zinb_model, data=list(N=length(y), 
                                          Y=y, 
                                          offset=exposure,
                                          offset_zi=exposure,
                                          prior_only=0,
                                          df=myprior_values[1],
                                          loc=myprior_values[2],
                                          scale=myprior_values[3]),
                    cores = nCores)
  loo_3 <- loo(fit_3)
  
  cmp12 <- compare(loo_1, loo_2)
  cmp13 <- compare(loo_1, loo_3)
  cmp23 <- compare(loo_2, loo_3)
  
  res <- c(cmp12, cmp13, cmp23)
  return(res)
  
}


fit_count_models <- function(y, exposure, zinb_model, nCores) {
  
  gexpr <- data.frame(y, exposure)
  
  fit_1 <- stan_glm(y ~ 1,
                    family=poisson,
                    offset=exposure,
                    data=gexpr,
                    cores = nCores,
                    #seed=SEED,
                    refresh=0)

  fit_2 <- stan_glm(y ~ 1,
                    family=neg_binomial_2,
                    offset=exposure, 
                    data=gexpr,
                    cores = nCores,
                    #seed=SEED, 
                    refresh=0)

  myprior <- get_prior(bf(y ~ 1 + offset(exposure), zi ~ 1 + offset(exposure)),
                       data = gexpr, 
                       family = zero_inflated_negbinomial())
  myprior_values <- eval(parse(text=gsub("student_t", "c", myprior$prior[1])))
  fit_3 <- sampling(zinb_model, data=list(N=length(y), 
                                          Y=y, 
                                          offset=exposure,
                                          offset_zi=exposure,
                                          prior_only=0,
                                          df=myprior_values[1],
                                          loc=myprior_values[2],
                                          scale=myprior_values[3]),
                    cores = nCores)

  return(list("P"=fit_1, "NB"=fit_2, "ZINB"=fit_3))
  
}

