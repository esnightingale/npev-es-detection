################################################################################
################################################################################
# Fit NPEV model
################################################################################

pacman::p_load(here, rstan, brms)
rm(list = ls(all = TRUE))  # clear all data
gc()

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

dir <- "data/Pakistan/analysis/prevalence"
outdir <- "output/prevalence/final"

# Read data ---------------------------------------------------------------
# - Set up in npev_prev_model_exploration.Rmd

fitdata <- readRDS(here(dir, "fitdata.rds"))
W <- readRDS(here(dir, "../W.rds"))

# Model specification -----------------------------------------------------

f = bf(
  n_npev | trials(n_npafp) ~     
    # Overall average prevalence
    1 + 
    (1|adm1_name) + (1|t) +
    # Correlated/uncorrelated random deviations by district 
    car(W, gr = guid, type = "bym2") + 
    # Long-term spline over time
    s(t, k = 10) +
    # Seasonal temporal splines overall and by province
    s(month_of_year, bs = "cc", k = 12) +
    # Long-term deviation by guid (penalise wiggliness now we're fitting many)
    s(t, guid, bs = "fs", k = 6)
)

# Prior specification -----------------------------------------------------

# Check default priors given this formula
priors <- c(prior(normal(-1,1), class = Intercept), # More informative on logit scale
            # Province and time deviations - few province levels and other temporal effects so keep these tight
            prior(exponential(2), class = sd,  group = "adm1_name"),
            prior(exponential(2), class = sd,  group = "t"),
            # BYM2 parameters
            # More density on non-correlated deviation
            prior(beta(1,2), class = rhocar), 
            # Moderately regularising prior on CAR sd (many levels, potential heterogeneity)
            prior(exponential(1), class = sdcar), 
            # Overall time effects
            # Spline coefficients
            prior(normal(0,1), class = "b", coef = "st_1"), 
            # These splines should be well-informed by the data so moderately regularising prior 
            prior(exponential(1), class = sds, coef = 's(t, k = 10)'), 
            prior(exponential(1), class = sds, coef = 's(month_of_year, bs = "cc", k = 12)'), 
            # District-specific time effects
            # Tighter prior on spline sd for district deviations
            prior(exponential(2.5), class = sds, coef = 's(t, guid, bs = "fs", k = 6)'))

# Prior check -------------------------------------------------------------

fit_prior <- brm(formula = f,
                 prior = priors,
                 data = fitdata,
                 data2 = list(W = W),
                 family = "binomial",
                 sample_prior = "only",
                 cores = parallel::detectCores(),
                 control = list(adapt_delta = 0.98),
                 save_pars = save_pars(all = T),
                 file = here(outdir,"fit_prior.rds"))

# Fitting -----------------------------------------------------------------

fit <- brm(formula = f,
           prior = priors,
           data = fitdata,
           data2 = list(W = W),
           family = "binomial",
           refresh=250,
           iter = 5000, warmup = 2000, 
           cores = parallel::detectCores(),
           control = list(adapt_delta = 0.98),
           save_pars = save_pars(all = T),
           file = here(outdir,"fit_posterior.rds"))

################################################################################
################################################################################
