priors <- list(
  
  baseline = c(
    prior(normal(-1,1), class = "Intercept"), # Prior for baseline prevalence centred on ~25%
    prior(student_t(3, 0, 2.5), class = "sds") # Prior for spline std dev.
  ),
  
  kroiss = c(
    prior(normal(-1,1), class = "Intercept"), # Prior for baseline prevalence centred on ~25%
    prior(normal(0,0.5), class = "ar"), # Prior for AR(1) term
    # prior(beta(0.5,0.5),class="rhocar"), # BYM2 mixing parameter
    # prior(normal(0, 1), class = "sdcar"), # BYM2 std. dev.
    prior(student_t(3, 0, 2.5), class = "sds") # Prior for spline std dev.
  ),
  
  alt1 = c(
    prior(normal(-1,1), class = "Intercept"), # Prior for baseline prevalence centred on ~25%
    prior(normal(0, 2), class = "b"), # Priors for fixed effects
    prior(normal(0,0.5), class = "ar"), # Prior for AR(1) term
    # prior(beta(0.5,0.5),class="rhocar"), # BYM2 mixing parameter
    # prior(normal(0, 1), class = "sdcar"), # BYM2 std. dev.
    prior(student_t(3, 0, 2.5), class = "sds") # Prior for spline std dev.
  ),
  
  alt2 = c(
    prior(normal(-1,1), class = "Intercept"), # Prior for baseline prevalence centred on ~25%
    prior(normal(0, 2), class = "b"), # Priors for fixed effects
    prior(normal(0,0.5), class = "ar"), # Prior for AR(1) term
    # prior(beta(0.5,0.5),class="rhocar"), # BYM2 mixing parameter
    # prior(normal(0, 1), class = "sdcar"), # BYM2 std. dev.
    prior(student_t(3, 0, 2.5), class = "sds") # Prior for spline std dev. 
  ),
  
  alt3 = c(
    prior(normal(-1,1), class = "Intercept"), # Prior for baseline prevalence centred on ~25%
    prior(normal(0, 2), class = "b"), # Priors for fixed effects
    prior(normal(0,0.5), class = "ar"), # Prior for AR(1) term
    # prior(beta(0.5,0.5),class="rhocar"), # BYM2 mixing parameter
    # prior(normal(0, 1), class = "sdcar"), # BYM2 std. dev.
    prior(student_t(3, 0, 2.5), class = "sds") # Prior for spline std dev. 
  )
)
