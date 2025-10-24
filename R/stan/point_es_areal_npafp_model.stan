//

data {
    int<lower=1> N_obs;             // Number of observations (site:month)
    int<lower=1> N_sites;       // Number of unique site locations
    int<lower=1> N_districts;           // Number of districts
    int<lower=0> N_edges;     // Borders between districts (**?**)
    int<lower=1> N_months;          // Number of months
    int<lower=0> x[N_obs];          // Observed NPEV detections
    int<lower=0> samples[N_obs];     // Number of samples
    int<lower=1, upper=N_sites> site[N_obs]; // Site index
    int<lower=1, upper=N_months> month[N_obs];       // Month index
    int<lower=1, upper=N_districts> district[N_locations];   // District index for each site
    int<lower=1, upper=N_districts> node1[N_edges];     // For each district, coordinates of adjacent district centroids
    int<lower=1, upper=N_districts> node2[N_edges];
    matrix[N_sites, 2] coords;  // GPS coordinates (latitude, longitude)
}

parameters {
    real mu;                        // Overall intercept
    real beta;                      // Effect of district prevalence Z
    vector[N_districts] Z;               // Latent true values of Z_j
    vector[N_districts] eta_dist;        // Random effects for districts
    vector[N_sites] spatial_raw; // GP spatial effect (latent)
    vector[N_months] temp_raw;       // Temporal effect (latent)
    
    real<lower=0> tau_Z;  // Precision parameter
    real<lower=0> sigma_dist;        // SD of district effects
    real<lower=0> sigma_space;       // SD of spatial effect
    real<lower=0> rho_space;         // Length scale for spatial GP
    real<lower=0> sigma_time;        // SD of temporal effect
    
}

transformed parameters {
    vector[N_sites] spatial_effect;
    vector[N_months] temporal_effect;

    // Spatial effect via Gaussian Process
    {
        matrix[N_sites, N_sites] K = cov_exp_quad(coords, sigma_space, rho_space);
        for (i in 1:N_sites)
            K[i, i] += 1e-6; // Jitter for numerical stability
        spatial_effect = cholesky_decompose(K) * spatial_raw;
    }

    // Temporal effect modeled as a simple random walk
    temporal_effect[1] = temp_raw[1];
    for (t in 2:N_months)
        temporal_effect[t] = temporal_effect[t-1] + sigma_time * temp_raw[t];
}

model {
    // Priors
    mu ~ normal(0, 1);
    beta ~ normal(0, 1);
    sigma_dist ~ normal(0, 1);
    sigma_space ~ normal(0, 1);
    rho_space ~ normal(0, 1);
    sigma_time ~ normal(0, 1);
    
    tau_Z ~ normal(0, 1);
    Z ~ normal(0, 1); // Weak prior to center

    eta_dist ~ normal(0, sigma_dist);
    spatial_raw ~ normal(0, 1);
    temp_raw ~ normal(0, 1);

    // Latent ICAR model (with sum-to-zero constraint) for district prevalence Z:
    for (i in 1:N_edges) {
        target += -0.5 * tau_Z * square(Z[node1[i]] - Z[node2[i]]);
    }
    sum(Z) ~ normal(0, 0.001 * N_districts);

    // Likelihood
    for (n in 1:N_obs) {
        real logit_p = mu 
                     + beta * Z[district[site[n]]] // Latent district-level predictor effect (prevalence)
                     + spatial_effect[site[n]] // Spatial effect
                     + temporal_effect[month[n]]; // Temporal effect
        
        x[n] ~ binomial_logit(samples[n], logit_p);
    }
}
