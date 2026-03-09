################################################################################
################################################################################
# Fit ES model - primary and sensitivity
################################################################################

pacman::p_load(tidyverse, sf, janitor, units, here, brms)
rm(list = ls(all = TRUE))  # clear all data

options(mc.cores = parallel::detectCores())

source(here("R/config.R"))
# ES samples linked to district NPEV prevalence, by site coordinates
dir <- file.path(dir, "analysis/detection")
es_npev_xy <- readRDS(here(dir, "es_npev_xy.rds"))

# Data set up -------------------------------------------------------------

fitdata <- es_npev_xy |> 
  mutate(y = as.numeric(!is.na(npev)),
         t = as.numeric(as.factor(month)),
         moy = as.factor(as.numeric(month_of_year)),
         year_f = factor(year),
         season = as.factor(lubridate::quarter(month)),
         dist_t_id = as.factor(paste(guid,t, sep = ":")),
         # Add 1 to avoid 0 catchment
         log_pop_5k = log(catchment_pop_5k+1),
         log_pop_5k_s = scale(log_pop_5k),
         built_s = scale(site_built),
         precip_s = scale(site_precip21)
  )

# Model structure ---------------------------------------------------------

# Formula
f = bf(y ~ 1 + (1|site_id) + log_pop_5k_s + me(pred_mean, pred_sd, gr = dist_t_id))

# Priors
priors <- c(prior(normal(0,1), class = Intercept),
            prior(exponential(3), class = sd),
            prior(normal(0,1), class = b),
            prior(normal(0,1), class = meanme),
            prior(exponential(2), class = sdme))

# Fitting -----------------------------------------------------------------

# Function to fit with given data set, for main and sensitivity analyses

fit_es_model <- function(data, name = "primary"){
  
  save_latent = F
  save_all = F
  if (name == "primary"){save_latent = T; save_all = T}
  
  fit <- brm(f, 
             # Exclude pv+ positive samples (WPV/VDPV) from primary analysis
             data = data,
             family = "bernoulli",
             prior = priors,
             refresh=250,
             empty=FALSE, init = NULL,
             iter = 4000, thin = 2, chains = 4, 
             cores = parallel::detectCores(),
             control = list(adapt_delta = 0.96),
             save_pars = save_pars(all = save_all, latent = save_latent))
  return(fit)
}

# Sensitivity to excluded samples -----------------------------------------

# Define alternative data sets to include PV+ samples
alt_fitdata <- list(
  # Exclude PV+ samples (WPV/VDPV) from primary analysis
  primary = filter(fitdata, !pv),
  SA_pv_obs = fitdata,
  SA_pv_neg = mutate(fitdata, y = if_else(pv, 0, y)),
  SA_pv_pos = mutate(fitdata, y = if_else(pv, 1, y)))

# Quick summary of how this changes NPEV positivity
lapply(alt_fitdata, \(x){ summary(x$y) })

# Fit above model to all four datasets
fits_all <- purrr::map2(alt_fitdata, names(alt_fitdata), fit_es_model)

lapply(fits_all, \(x){ summary(x) })

# Save primary fit, and supplementary for comparison
outdir_fit <- file.path("output/detection", country)
dir.create(outdir_fit, recursive = TRUE, showWarnings = FALSE)
saveRDS(fits_all[[1]], here(outdir_fit, "fit_primary.rds"))
saveRDS(fits_all, here(outdir_fit, "fits_all.rds"))

# Also save alternative datasets
saveRDS(alt_fitdata, here(dir,"fitdata_SA.rds"))

################################################################################
################################################################################