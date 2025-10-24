// generated with brms 2.22.0
functions {
}
data {
  int<lower=1> N;  // total number of observations
  array[N] int Y;  // response variable
  array[N] int trials;  // number of trials
  // data for spline 1 (seasonal)
  int nb_1;  // number of bases
  array[nb_1] int knots_1;  // number of knots
  // basis function matrices
  matrix[N, knots_1[1]] Zs_1_1;
  // data for the CAR structure
  int<lower=1> Nloc;
  array[N] int<lower=1> Jloc;
  int<lower=0> Nedges;
  array[Nedges] int<lower=1> edges1;
  array[Nedges] int<lower=1> edges2;
  // scaling factor of the spatial CAR component
  real<lower=0> car_scale;
  // time data for random walk
  int<lower=1> Nprov; // number of provinces
  int<lower=1> Ntime;  // number of time points
  array[N] int<lower=1> Jprov;
  array[N] int<lower=1> time_id;  // time index for each observation
  int prior_only;  // should the likelihood be ignored?
}

parameters {
  real Intercept;  // temporary intercept for centered predictors
  // parameters for spline 1
  // standardized penalized spline coefficients
  vector[knots_1[1]] zs_1_1;
  vector<lower=0>[nb_1] sds_1;  // SDs of penalized spline coefficients
  real<lower=0> sdcar;  // SD of the CAR structure
  // parameters for the BYM2 structure
  vector[Nloc] zcar;  // spatial part
  vector[Nloc] nszcar;  // non-spatial part
  // proportion of variance in the spatial part
  real<lower=0,upper=1> rhocar;
  // parameters for province-level random walks
  matrix[Nprov, Ntime] rw_p_raw;  // raw increments
  vector[Ntime] rw0;  // initial values for each time point (shared across provinces)
  real<lower=0> sdrw_p;  // shared SD of random walks
}
transformed parameters {

  // penalized spline coefficients
  vector[knots_1[1]] s_1_1;
  
  // scaled parameters for the BYM2 structure
  vector[Nloc] rcar;
  real lprior = 0;  // prior contributions to the log posterior

  // join the spatial and the non-spatial CAR component
  rcar = (sqrt(1 - rhocar) * nszcar + sqrt(rhocar * inv(car_scale)) * zcar) * sdcar;

  // compute penalized spline coefficients
  s_1_1 = sds_1[1] * zs_1_1;

  // compute uncentered RW1 values by accumulating increments
  matrix[Nprov, Ntime] rw_p;
  for (p in 1:Nprov) {
    rw_p[p, 1] = rw0[1] + rw_p_raw[p, 1];
    for (t in 2:Ntime) {
      rw_p[p, t] = rw_p[p, t - 1] + rw_p_raw[p, t];
    }
  }

  lprior += student_t_lpdf(Intercept | 3, 0, 2.5);
  lprior += student_t_lpdf(sds_1 | 3, 0, 2.5) - 1 * student_t_lccdf(0 | 3, 0, 2.5);

  // Priors for BYM2
  // PC prior: P(sdcar > U) = alpha  →  λ = -log(alpha) / U
  real lambda_sdcar = -log(0.01) / 1;  // 1% chance sdcar > 1
  lprior += exponential_lpdf(sdcar | lambda_sdcar);
  lprior += beta_lpdf(rhocar | 0.5, 1); // PC prior shrinking to non-spatial correlation

  // Compute uncentered RW1 values by accumulating increments
  for (d in 1:Nloc) {
    rw_d[d, 1] = rw0[1] + rw_d_raw[d, 1];
    for (t in 2:Ntime) {
      rw_d[d, t] = rw_d[d, t - 1] + rw_d_raw[d, t];
    }
  }
  // Priors for RW
  real lambda_sdrw = -log(0.01) / 1;
  lprior += exponential_lpdf(sdrw_d | lambda_sdrw);
  
  }
model {
  // likelihood including constants
  if (!prior_only) {
    // initialize linear predictor term
    vector[N] mu = rep_vector(0.0, N);
    mu += Intercept + Zs_1_1 * s_1_1;
    for (n in 1:N) {
      // add more terms to the linear predictor
      mu[n] += rcar[Jloc[n]] + rw_p[Jprov[n], time_id[n]];
    }
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
  target += std_normal_lpdf(to_vector(rw_p_raw));
  target += normal_lpdf(rw0 | 0, sdrw_p);
}
generated quantities {

  // Log-likelihood and posterior prediction
  array[N] real log_lik;
  vector[N] linpred = Intercept + Zs_1_1 * s_1_1;
    for (n in 1:N) {
      // add more terms to the linear predictor
      linpred[n] += rcar[Jloc[n]] + rw_p[Jprov[n], time_id[n]];
      log_lik[n] = binomial_logit_lpmf(Y[n] | trials[n], linpred[n]);
    }
  // Back-transform linpred and draw from binomial pmf for prediction
  vector[N] epred = inv_logit(linpred); 
  array[N] real y_pred = binomial_rng(trials, epred);
  
  // Province-level RW deviations from overall seasonal trend
  matrix[Nprov, Ntime] rw_dev;
  for (t in 1:Ntime) {
    real mean_t = mean(to_vector(rw_p[, t]));
    for (p in 1:Nprov) {
      rw_dev[p, t] = rw_p[p, t] - mean_t;
    }
  }
}
