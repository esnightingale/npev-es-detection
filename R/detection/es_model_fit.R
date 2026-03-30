################################################################################
################################################################################
# Fit ES model - primary and sensitivity
################################################################################

pacman::p_load(tidyverse, sf, janitor, units, here, brms)
rm(list = ls(all = TRUE))  # clear all data

options(mc.cores = parallel::detectCores())

source(here("R/config.R"))
# ES samples linked to district EV prevalence, by site coordinates
dir <- file.path(dir, "analysis/detection")
es_ev_xy <- readRDS(here(dir, "es_ev_xy.rds"))

# Data set up -------------------------------------------------------------
# Outcome: EV detection (any virus - PV, vaccine, or NPEV) among all ES samples
fitdata <- es_ev_xy |> 
  mutate(
    # catchment_pop_5k = coalesce(catchment_pop_5k, site_catchment_pop, 0),
    y = as.numeric(ev),
    t = as.numeric(as.factor(month)),
    moy = as.factor(as.numeric(month_of_year)),
    year_f = factor(year),
    season = as.factor(lubridate::quarter(month)),
    dist_t_id = as.factor(paste(guid,t, sep = ":")),
    # Covariates
    catchment_s = scale(site_catchment_pop)[1,],
    # Add 1 to avoid 0 catchment
    log_pop_5k = log(site_catchment_pop+1),
         log_pop_5k_s = scale(log_pop_5k)[1,],
         built_s = scale(site_built)[1,],
         precip_s = scale(site_precip21)[1,]
  )

# Raw relationship
plot(jitter(fitdata$pred_p_logit), fitdata$y, 
     xlab = "logit(prevalence)", ylab = "ES positivity")

# Model structure ---------------------------------------------------------

# Formula
f_base = bf(y ~ 1 + (1|site_id) + me(pred_p_logit, sd_p_logit))
f = bf(y ~ 1 + (1|site_id) + me(pred_p_logit, sd_p_logit) +
         moy + catchment_s + site_class + built_s + precip_s)
# Priors
priors <- c(prior(normal(0,1), class = Intercept),
            prior(exponential(3), class = sd),
            prior(normal(0,1), class = b),
            prior(normal(0,1), class = meanme),
            prior(exponential(2), class = sdme))

# Fitting -----------------------------------------------------------------

# Function to fit with given data set, for main and sensitivity analyses

fit_es_model <- function(f, data, name = "primary"){
  
  save_latent = F
  save_all = F
  if (name == "primary"){save_latent = T; save_all = T}
  
  fit <- brm(f, 
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


# Primary analysis: EV detection across all samples ----------------------

# Fit above model to all four datasets
fit <- fit_es_model(fitdata)

summary(fit)

outdir_fit <- file.path("output/detection", country)
dir.create(outdir_fit, recursive = TRUE, showWarnings = FALSE)
saveRDS(fit, here(outdir_fit, "fit_primary.rds"))
saveRDS(fitdata, here(dir, "fitdata.rds"))


# Primary analysis: Base and covariate-adjusted model ---------------------

f_list <- list(f_base, f)
fits <- lapply(f_list, fit_es_model, data = fitdata)

lapply(fits, \(x){ summary(x) })

# Quick summary:

# Conditional effects (should show flat line now, sloped after fix)
plot(conditional_effects(fits[[1]], "pred_p_logit"))

# Predicted vs observed
pred <- posterior_epred(fits[[1]])
plot(colMeans(pred), colMeans(fitdata$y), xlim = c(0,1), ylim = c(0,1))
abline(0,1, col = "red")

# Site random effects
site_effects <- ranef(fits[[2]])$site_id[,, "Intercept"]
site_summary <- site_effects %>%
  as_tibble() %>%
  mutate(site_id = rownames(site_effects)) %>%
  pivot_longer(-site_id, names_to = "draw", values_to = "logit_sens") %>%
  group_by(site_id) %>%
  summarise(
    mean_sens = plogis(mean(logit_sens)),
    low_sens = plogis(quantile(logit_sens, 0.025)),
    hi_sens = plogis(quantile(logit_sens, 0.975)),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_sens)) |> 
  left_join(fitdata |> group_by(site_id) |> summarise(pred_p_logit = mean(pred_p_logit)))
View(site_summary)

ggplot(site_summary, aes(mean_sens, xmin = low_sens, xmax = hi_sens,
                         y = site_id)) +
  geom_pointrange() +
  theme(axis.text.y = element_blank())


# Calculate site-level probability of detection at average district:month prevalence
# for comparability: 
# 1. Define AFP prevalence grid (your range: use quantiles or fixed)
p_afp_grid <- quantile(plogis(fitdata$pred_p_logit), c(0.1, 0.5, 0.9))  # 10th/50th/90th
logit_p_grid <- qlogis(p_afp_grid)
print("Grid AFP EV %:") 
round(p_afp_grid * 100, 1)

# 2. Newdata: all sites × grid (average covariates within site)
# Model expects pred_p_logit (not logit_p_grid) for me(pred_p_logit, sd_p_logit)
newfitdata_sites <- fitdata %>%
  group_by(site_id, site_class) %>%
  summarise(
    moy = factor(round(median(as.numeric(moy))), levels = levels(moy)),
    catchment_s = mean(catchment_s),
    built_s = mean(built_s), 
    precip_s = mean(precip_s),
    sd_p_logit = mean(sd_p_logit),  # site-average ME
    .groups = "drop"
  ) %>%
  crossing(logit_p_grid = logit_p_grid) %>%
  mutate(pred_p_logit = logit_p_grid)  # brms me() requires this column name

# 3. Posterior predictions (probability scale)
site_probs_raw <- posterior_epred(
  fits[[2]],                    # your fitted model
  newdata = newfitdata_sites,
  ndraws = 1000,
  transform = TRUE           # P(Y=1)
) 

# 4. Convert to tidy, summarise
# Row order from crossing: site1-g1, site1-g2, site1-g3, site2-g1, ... so grid index cycles 1,2,3,1,2,3
n_sites <- length(unique(newfitdata_sites$site_id))
colnames(site_probs_raw) <- paste0(newfitdata_sites$site_id, "_", 
  rep(seq_along(logit_p_grid), times = n_sites))

grid_lookup <- tibble(logit_p_idx = seq_along(logit_p_grid), p_afp = p_afp_grid)

site_summary <- site_probs_raw %>%
  as_tibble() %>%
  mutate(draw = row_number()) %>%
  pivot_longer(-draw, names_to = "site_grid", values_to = "prob") %>%
  mutate(
    site_id = sub("_\\d+$", "", site_grid),
    logit_p_idx = as.integer(sub(".*_", "", site_grid))
  ) %>%
  left_join(grid_lookup, by = "logit_p_idx") %>%
  group_by(site_id, p_afp) %>%
  summarise(
    mean_prob = mean(prob),
    low_prob  = quantile(prob, 0.025),
    hi_prob   = quantile(prob, 0.975),
    n_draws   = n(),
    .groups = "drop"
  ) %>%
  arrange(p_afp, desc(mean_prob))

# 5. View/export
print("Site rankings at median AFP EV:")
site_summary %>% 
  filter(near(p_afp, median(p_afp_grid))) %>%
  select(site_id, p_afp, mean_prob, low_prob, hi_prob) %>%
  print(n = 20)

# 6. Save full table
write_csv(site_summary, "site_sensitivity_by_afp_prevalence.csv")


# Sensitivity to prevalence -----------------------------------------------

fitdata2 <- filter(fitdata, pred_p_logit < -1)

fit <- fit_es_model(f, fitdata2)
summary(fit)


# Sensitivity to excluded samples -----------------------------------------

# Primary: EV detection among ALL ES samples
# Sensitivity analyses: compare handling of PV+ samples (for method comparison)
alt_fitdata <- list(
  primary = fitdata,
  SA_pv_obs = fitdata,
  SA_pv_neg = mutate(fitdata, y = if_else(pv, 0, y)),
  SA_pv_pos = mutate(fitdata, y = if_else(pv, 1, y)))

# Quick summary of how this changes EV positivity
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