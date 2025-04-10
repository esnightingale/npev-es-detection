fit_npafp_model <- function(data, W, formula){

  # Fit brms model
  fit <- brms::brm(formula = formula,
                       data = data, 
                       data2 = list(W = W),
                       family = binomial(link = "logit"),
                       prior = c(prior(normal(-1,1), class = "Intercept"), # Prior for baseline prevalence centred on ~25%
                                 prior(normal(0, 2), class = "b")),  # Priors for fixed effects
                       iter = 2000, chains = 4, cores = 4,
                       control = list(adapt_delta = 0.95)  # To avoid divergences
                       )
  
  return(fit)
}
