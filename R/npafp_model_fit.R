################################################################################
# Fit and compare model structures for areal NPEV prevalence among NPAFP cases 
# - In particular explore specification of RW temporal structure as in:
#   https://discourse.mc-stan.org/t/random-walk-or-autoregressive-prior-in-brms/28411/11 
################################################################################

# Env setup -----------------------------------------------------------
pacman::p_load(tidyverse, here, cmdstanr, brms, gtsummary, units, loo)
theme_set(theme_minimal())

dir <- "data/Pakistan"

# Load data ------------------------------------------------------------------

# NP-AFP analysis data
data <- readRDS(here(dir, "afp_analysis.rds"))

# Set up additional variables
fitdata <- mutate(data, 
                  month_of_year = as.numeric(month_of_year),
                  year_F = factor(year),
                  guid_F = factor(guid, levels = unique(data$guid)),
                  t_F = factor(t, levels = 1:max(data$t)),
                  # standardise pop density for stability
                  log_pop_dens = log(guid_pop_dens),
                  log_pop_dens_s = (log_pop_dens - mean(log_pop_dens))/sd(log_pop_dens)) 

# Shapefiles and adjacency matrix
shape2 <- readRDS(here(dir,"shape2.rds")) 
W <- readRDS(here(dir,"W.rds"))

# Define formulae and priors--------------------------------------------------

# Set up stanvars and priors
sv <- stanvar(scode = "real rw1_lpdf(vector alpha, real sigma) {
    int N = num_elements(alpha);
    return normal_lpdf(alpha[2:N] - alpha[1:N - 1] | 0, sigma);
  }",
              block = "functions") +
  stanvar(scode = "r_1_1 ~ rw1(sd_t_rw);", block = "model", position = "end") +
  stanvar(scode = "real<lower=0> sd_t_rw;", block = "parameters", position = "end") +
  stanvar(scode = "sd_t_rw ~ normal(0,1);", block = "model", position = "end") +
  stanvar(scode = "sum(r_1_1) ~ normal(0, N_1*.001);", block = "model", position = "end")

# Formula
f1 <- bf(n_npev | trials(n_npafp) ~
               # Fixed effect on population density
               log_pop_dens_s +
               # 1st order random walk by month
               (1|t))
  
# Priors
pr1 <- prior(constant(10000), class = "sd", group = "t") +
  prior(normal(-1, 2), class=Intercept) +
  prior(normal(0, 2), class=b)

# Fit model
fit <- brm(f1, 
           family = "binomial",
           data = fitdata, 
           prior = pr1,
           stanvars = sv, 
           chains=1)

summary(fit)


# Alternative approach ----------------------------------------------------

# Define a function to set up the RW prior
prior_rw <- function(idx_inter = 1, 
                     sigma_inter, 
                     prior_scale_sigma_rw,...
                     ) {
  
  pr <- prior_string(paste0("rw1(", idx_inter, ", ", sigma_inter, ", sigma_rw)"), ...) +
    prior_string(paste0("target += normal_lpdf(sigma_rw| 0, ", prior_scale_sigma_rw, ")"), check=FALSE)
  
  sv_rw1 <- stanvar(scode = "real rw1_lpdf(vector alpha, int idx_inter, real sigma_inter, real sigma_rw) {
    int N = num_elements(alpha);
    return normal_lpdf(alpha[2:N] - alpha[1:N - 1] | 0, sigma_rw) + normal_lpdf(alpha[idx_inter]| 0, sigma_inter);
  }", 
                    block = "functions")
  
  sv_sigma_rw <- stanvar(scode="real<lower=0> sigma_rw;", 
                         block="parameters")
  
  list(prior=pr, stanvar=sv_rw1 + sv_sigma_rw)
}

# prior_rw <- function(idx_inter = 1, prior_inter = "normal(0, 5)",
#                       prior_scale_sigma_rw = "exponential(1)", ...) {
#   stanvars <- stanvar(
#     scode = paste(
#       "vector[N_t_factor] b_rw;",
#       "real sigma_rw;",
#       "b_rw[", idx_inter, "] = b_rw_Intercept;",  # Anchor first time point
#       "for (i in 1:(N_t_factor-1)) {",
#       "  b_rw[i + 1] = b_rw[i] + normal_rng(0, sigma_rw);",  # RW(1) prior
#       "}"
#     ),
#     block = "model"
#   )
#   
#   prior_list <- c(
#     prior_string(prior_inter, class = "b", nlpar = "rw", coef = "t_F1"),  # Prior for first time point
#     prior_string(prior_scale_sigma_rw, class = "sigma_rw")  # Prior for RW scale
#   )
#   
#   list(prior = prior_list, stanvars = stanvars)
# }

# # Generate RW prior components
# rw_prior <- prior_rw(
#   idx_inter = 1,                   # Anchor first time point (t=1)
#   prior_inter = "normal(-1, 2)",    # Prior for initial time effect
#   prior_scale_sigma_rw = "exponential(0.5)"  # Prior for RW volatility
# )
# 
# # Combine with other priors (e.g., fixed effects)
# full_prior <- c(
#   rw_prior$prior,
#   set_prior("normal(0, 1)", class = "b", nlpar = "lin")  # Priors for fixed effects
# )

# Define formula with linear affect of log population density and temporal RW
f_base = bf(
  n_npev | trials(n_npafp) ~ lin + rw,
  lin ~ 0 + log_pop_dens_s +
    # Periodic spline on 1-12 for seasonality
    s(month_of_year, bs = "cc", k = 12),
  rw ~ (1|guid) + t_F, nl=TRUE, center=FALSE)

# Check default priors given this formula
get_prior(f_base, fitdata, family = binomial("logit"))

# Define the RW prior
rw_prior <- prior_rw(1, 5, 0.5, class="b", nlpar="rw")

# Define the complete set of priors
pr_base = c(
  # prior(normal(0,1), class = "b", coef = "lin_logguid_pop_dens_s"), # Prior for fixed effects
  rw_prior$prior
)

# Fit the model
alt_fit <- brm(f_base, data = fitdata,
           family = "binomial",
               prior = pr_base,
               stanvars = rw_prior$stanvar,
               refresh=250,
               empty=FALSE, init=0.2,
           control = list(adapt_delta = 0.95)) ##, sample_prior="only")

alt_fit
plot_pp_checks(alt_fit)
plot(alt_fit)

# Predicted values for each time point
model_pred <- tibble(n_npafp = mean(fitdata$n_npafp),
                     t,
                     t_F=levels(fitdata2$t_F)) # x=0

rvar(posterior_epred(alt_fit, newdata=model_pred))


