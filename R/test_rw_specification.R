# we recommend running this in a fresh R session or restarting your current session
install.packages("cmdstanr", repos = c('https://stan-dev.r-universe.dev', getOption("repos")))

library(here)
library(brms)
library(dplyr)
library(tidyr)
library(posterior)

# instruct brms to use cmdstanr as backend and cache all Stan binaries
options(brms.backend="cmdstanr", cmdstanr_write_stan_file_dir=here("_brms-cache"))

# create cache directory if not yet available
dir.create(here("_brms-cache"), FALSE)

# set the R random seed of the session for reproducibility
set.seed(5886935)

## PREVIOUS APPROACH

# Simulate Data
N <- 100
d <- tibble(agecat = sample.int(8, N, replace = TRUE)) %>%
  mutate(mu = sin(pi*agecat/7), x=rnorm(N), 
         Y = mu + 2 * x + rnorm(N, sd = 1),
  )

# Set up stanvars and priors
sv <- stanvar(scode = "real rw1_lpdf(vector alpha, real sigma) {
    int N = num_elements(alpha);
    return normal_lpdf(alpha[2:N] - alpha[1:N - 1] | 0, sigma);
  }",
              block = "functions") +
  stanvar(scode = "r_1_1 ~ rw1(sd_agecat_rw);", block = "model", position = "end") +
  stanvar(scode = "real<lower=0> sd_agecat_rw;", block = "parameters", position = "end") +
  stanvar(scode = "sd_agecat_rw ~ normal(0,1);", block = "model", position = "end") +
  stanvar(scode = "sum(r_1_1) ~ normal(0, N_1*.001);", block = "model", position = "end")

pr <- prior(constant(10000), class = "sd", group = "agecat") +
  prior(normal(0, 5), class=Intercept) +
  prior(normal(0, 2), class=b)

# Fit model
fit <- brm(Y ~ 1 + x + (1|agecat), data = d, prior = pr, stanvars = sv, chains=1)


## ALTERNATIVE APPROACH

## Defines an RW prior for brms. This requires to define in brms a
## non-linear model such that the RW term is it's own linear model
## which must include the overall intercept (so set the overall
## intercept to 0 for the remaining linear model). The user can pick
## which of the categories for the random walk is used as overall
## intercept by setting the idx_inter to the respective index. The
## function returns a list with the elements prior and stanvar, which
## need to be passed to the brm call respectivley.
prior_rw <- function(idx_inter, sigma_inter, prior_scale_sigma_rw, ...) {
  pr <- prior_string(paste0("rw1(", idx_inter, ", ", sigma_inter, ", sigma_rw)"), ...) +
    prior_string(paste0("target += normal_lpdf(sigma_rw| 0, ", prior_scale_sigma_rw, ")"), check=FALSE)
  
  sv_rw1 <- stanvar(scode = "real rw1_lpdf(vector alpha, int idx_inter, real sigma_inter, real sigma_rw) {
    int N = num_elements(alpha);
    return normal_lpdf(alpha[2:N] - alpha[1:N - 1] | 0, sigma_rw) + normal_lpdf(alpha[idx_inter]| 0, sigma_inter);
  }", block = "functions")
  
  sv_sigma_rw <- stanvar(scode="real<lower=0> sigma_rw;", block="parameters")
  
  list(prior=pr, stanvar=sv_rw1 + sv_sigma_rw)
}

age_rw_prior <- prior_rw(1, 5, 0.5, class="b", nlpar="rw")

alt_model <- bf(Y ~ lin + rw, lin ~ 0 + x, rw ~ 1 + agecatF, nl=TRUE, center=FALSE)

d2 <- mutate(d, agecatF=factor(agecat))

get_prior(alt_model, d2)

alt_prior <- prior(normal(0, 2), class=b, coef=x, nlpar=lin) +
  age_rw_prior$prior

alt_fit <- brm(alt_model, data = d2,
               prior = alt_prior,
               stanvars = age_rw_prior$stanvar,
               refresh=250,
               empty=FALSE, init=0.2) ##, sample_prior="only")

alt_fit

model_pred <- tibble(agecatF=levels(d2$agecatF), x=0)

rvar(posterior_epred(alt_fit, newdata=model_pred))
