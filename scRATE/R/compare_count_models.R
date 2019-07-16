library(rstan)
library(rstanarm)
library(rstantools)
library(brms)
library(loo)

compare_count_models <- function(y, exposure, zip_model, zinb_model, nCores, seed) {
  
  gexpr <- data.frame(y, exposure)
  
  fit_1 <- stan_glm(y ~ 1,
                    family=poisson,
                    offset=exposure,
                    data=gexpr,
                    cores = nCores,
                    seed=seed,
                    refresh=0)
  loo_1 <- loo(fit_1)
  
  fit_2 <- stan_glm(y ~ 1,
                    family=neg_binomial_2,
                    offset=exposure, 
                    data=gexpr,
                    cores = nCores,
                    seed=seed, 
                    refresh=0)
  loo_2 <- loo(fit_2)

  myprior_3 <- get_prior(bf(y ~ 1 + offset(exposure), zi ~ 1 + offset(exposure)),
                         data = gexpr, 
                         family = zero_inflated_poisson())
  myprior_3_values <- eval(parse(text=gsub("student_t", "c", myprior_3$prior[1])))
  fit_3 <- sampling(zip_model, 
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
  loo_3 <- loo(fit_3)
  
  myprior_4 <- get_prior(bf(y ~ 1 + offset(exposure), zi ~ 1 + offset(exposure)),
                       data = gexpr, 
                       family = zero_inflated_negbinomial())
  myprior_4_values <- eval(parse(text=gsub("student_t", "c", myprior_4$prior[1])))
  fit_4 <- sampling(zinb_model,
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
  loo_4 <- loo(fit_4)
  
  res <- loo_compare(loo_1, loo_2, loo_3, loo_4)
  return(res)
  
}


fit_count_models <- function(y, exposure, zip_model, zinb_model, nCores, seed) {
  
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
  fit_3 <- sampling(zip_model, 
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

  myprior_4 <- get_prior(bf(y ~ 1 + offset(exposure), zi ~ 1 + offset(exposure)),
                       data = gexpr, 
                       family = zero_inflated_negbinomial())
  myprior_4_values <- eval(parse(text=gsub("student_t", "c", myprior_4$prior[1])))
  fit_4 <- sampling(zinb_model,
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

  return(list("P"=fit_1, "NB"=fit_2, "ZIP"=fit_3, "ZINB"=fit_4))
  
}

