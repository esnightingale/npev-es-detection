// generated with brms 2.22.0
functions {
}
data {
  int<lower=1> N;  // Total number of observations (districts*months)
  array[N] int Y;  // Number of NPEV positives 
  array[N] int trials;  // Number of non-polio AFP cases
  
  //int<lower=1> Nprov; // Number of provinces
  int<lower=1> Nloc; // Number of districts
  //array[N] int<lower=1> Jprov; // Province index
  array[N] int<lower=1> Jloc; // District index
  
  int<lower=1> Ntime; // Number of time points
  array[N] int<lower=1> Jtime;  // Time index 
  
  // data for spline 1 (seasonal)
  int nb_1;  // number of bases
  array[nb_1] int knots_1;  // number of knots
  // basis function matrices
  matrix[N, knots_1[1]] Zs_1_1;
  
  // data for the CAR structure
  int<lower=0> Nedges; 
  array[Nedges] int<lower=1> edges1;
  array[Nedges] int<lower=1> edges2;
  real<lower=0> car_scale;  // scaling factor of the spatial CAR component
  
  int prior_only;  // should the likelihood be ignored?
}

parameters {
  real Intercept;  // temporary intercept for centered predictors
  
  // parameters for spline 1
  vector[knots_1[1]] zs_1_1;  // standardized penalized spline coefficients
  vector<lower=0>[nb_1] sds_1;  // SDs of penalized spline coefficients
  
  // parameters for the BYM2 structure
  vector[Nloc] zcar;  // spatial part
  vector[Nloc] nszcar;  // non-spatial part
  real<lower=0> sdcar;  // SD of the CAR structure
  real<lower=0,upper=1> rhocar;  // proportion of variance in the spatial part
  
  // parameters for district-level random walks
  vector[Ntime] rw_raw;  // raw increments
  vector[Ntime] rw0;  // initial values for each time point
  real<lower=0> sdrw;  // SD of random walk
}
transformed parameters {

  vector[Nloc] rcar;
  vector[Ntime] rw;
  
  real lprior = 0;  // prior contributions to the log posterior
  lprior += normal_lpdf(Intercept | 0, 1);
  
  // Compute penalized spline coefficients
  vector[knots_1[1]] s_1_1;
  s_1_1 = sds_1[1] * zs_1_1;
  // Prior for spline SD
  lprior += exponential_lpdf(sds_1 | 0.5);
  
  // Join the spatial and the non-spatial CAR component
  rcar = (sqrt(1 - rhocar) * nszcar + sqrt(rhocar * inv(car_scale)) * zcar) * sdcar;
  // Priors for BYM2 parameters
  lprior += exponential_lpdf(sdcar | 0.5);
  lprior += beta_lpdf(rhocar | 0.5, 0.5);
  
  // Compute uncentered RW1 values by accumulating increments
  rw[1] = rw0[1] + rw_raw[1];
    for (t in 2:Ntime) {
      rw[t] = rw[t - 1] + rw_raw[t];
    }
  // Prior for RW SD
  lprior += exponential_lpdf(sdrw | 1);
  
}
model {
  
 if (!prior_only) {
   
    vector[N] mu = rep_vector(0.0, N);
    mu += Intercept + Zs_1_1 * s_1_1;
    
    // add BYM2 and RW terms to the linear predictor
    for (n in 1:N) {
      mu[n] += rcar[Jloc[n]] + rw[Jtime[n]];
    }
    
    // Likelihood
    target += binomial_logit_lpmf(Y | trials, mu);
  }
  
  // priors including constants
  target += lprior;
  target += std_normal_lpdf(zs_1_1);
  
  // improper prior on the spatial BYM2 component
  target += -0.5 * dot_self(zcar[edges1] - zcar[edges2]);
  // soft sum-to-zero constraint
  target += normal_lpdf(sum(zcar) | 0, 0.001 * Nloc);
  // proper prior on the non-spatial BYM2 component
  target += std_normal_lpdf(nszcar);
  
  // RW1 increments prior
  target += std_normal_lpdf(to_vector(rw_raw));
  target += normal_lpdf(rw0 | 0, sdrw);
  
}
generated quantities {

  // Log-likelihood and posterior prediction
  array[N] real log_lik;
  vector[N] linpred = Intercept + Zs_1_1 * s_1_1;
    for (n in 1:N) {
      // add more terms to the linear predictor
      linpred[n] += rcar[Jloc[n]] + rw[Jtime[n]];
      log_lik[n] = binomial_logit_lpmf(Y[n] | trials[n], linpred[n]);
    }
  // Back-transform linpred and draw from binomial pmf for prediction
  vector[N] epred = inv_logit(linpred); 
  array[N] real y_pred = binomial_rng(trials, epred);
  
  real logit_rhocar = log(rhocar / (1.0 - rhocar));
  
  // RW deviations from overall seasonal trend
  vector[Ntime] rw_dev;
  real mean_t = mean(rw);
  for (t in 1:Ntime) {
      rw_dev[t] = rw[t] - mean_t;
  }
}
