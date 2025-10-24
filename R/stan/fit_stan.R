################################################################################
################################################################################

install.packages("rstan")
install.packages("sf")         # For spatial data
install.packages("fields")     # For distance calculations
install.packages("Matrix")     # For sparse adjacency matrices (if using ICAR)
install.packages("ggplot2")    # For visualization

library(rstan)
library(sf)
library(fields)
library(Matrix)
library(ggplot2)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# Simulate data -----------------------------------------------------------

# Simulate spatial coordinates for N locations
set.seed(42)
N_sites <- 50
coords <- cbind(runif(N_sites, 0, 10), runif(N_sites, 0, 10)) # (X,Y)

# Assign each location to one of N_areas
N_districts <- 10
district <- sample(1:N_districts, N_sites, replace = TRUE)

# Generate areal centroids for Gaussian Process
district_coords <- aggregate(coords, by=list(district), FUN=mean)[,2:3]

# Temporal structure: N_months with variable observations
N_months <- 12
N_obs <- 200
site <- sample(1:N_sites, N_obs, replace = TRUE)
month <- sample(1:N_months, N_obs, replace = TRUE)
samples <- sample(5:20, N_obs, replace = TRUE)  # Number of samples per site-month
x <- rbinom(N_obs, samples, prob = 0.3)  # Simulated detection data

# Pack data into a list for Stan
stan_data <- list(
  N_obs = N_obs,
  N_sites = N_sites,
  N_districts = N_districts,
  N_months = N_months,
  x = x,
  samples = samples,
  site = site,
  month = month,
  district = district,
  coords = coords,
  district_coords = district_coords
)


# Fitting ------------------------------------------------------------------

# Compile the model
model <- stan_model("point_es_model.stan")

# Fit the model
fit <- sampling(model, data = stan_data, iter = 2000, chains = 4, warmup = 1000, 
                thin = 2, control = list(adapt_delta = 0.95))

# Print summary
print(fit, pars = c("mu", "beta", "sigma_dist", "sigma_space", "rho_space"))

# Summarise and plot -------------------------------------------------------

# Extract posterior samples
posterior_samples <- extract(fit)

# Plot posterior distributions
hist(posterior_samples$beta, 
     main = "Posterior of Beta (Effect of Z)", 
     xlab = "Beta", 
     col = "skyblue", breaks = 30)

# Extract spatial effect means
spatial_means <- colMeans(posterior_samples$spatial_effect)

ggplot(data = data.frame(x = coords[,1], 
                         y = coords[,2], 
                         effect = spatial_means),
       aes(x = x, y = y, color = effect)) +
  geom_point(size = 3) +
  scale_color_viridis_c() +
  theme_minimal() +
  labs(title = "Spatial Effect Estimates", color = "Effect")

# Check convergence --------------------------------------------------------

summary(fit)$summary[, "Rhat"]

# If ~1, then the chains have converged
# If >1.1, then the chains have not converged; try increasing adapt_delta / iter

################################################################################
################################################################################