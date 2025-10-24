// generated with brms 2.22.0
functions {
}
data {
  int<lower=1> N;  // total number of observations
  array[N] int Y;  // response variable
  array[N] int trials;  // number of trials
  int<lower=1> K;  // number of population-level effects
  matrix[N, K] X;  // population-level design matrix
  int<lower=1> Kc;  // number of population-level effects after centering
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
  // data for group-level effects of ID 1
  int<lower=1> N_1;  // number of grouping levels
  int<lower=1> M_1;  // number of coefficients per level
  array[N] int<lower=1> J_1;  // grouping indicator per observation
  // group-level predictor values
  vector[N] Z_1_1;
  // time data for random walk
  int<lower=1> T;  // number of time points
  array[N] int<lower=1> time_id;  // time index for each observation
  int prior_only;  // should the likelihood be ignored?
}
transformed data {
  matrix[N, Kc] Xc;  // centered version of X without an intercept
  vector[Kc] means_X;  // column means of X before centering
  for (i in 2:K) {
    means_X[i - 1] = mean(X[, i]);
    Xc[, i - 1] = X[, i] - means_X[i - 1];
  }
}
parameters {
  vector[Kc] b;  // regression coefficients
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
  vector<lower=0>[M_1] sd_1;  // group-level standard deviations
  array[M_1] vector[N_1] z_1;  // standardized group-level effects
  // parameters for district-level random walks
  matrix[Nloc, T] rw_d;  // RW1 per district
  real<lower=0> sdrw_d;  // shared SD of random walks
}
transformed parameters {
  // penalized spline coefficients
  vector[knots_1[1]] s_1_1;
  // scaled parameters for the BYM2 structure
  vector[Nloc] rcar;
  vector[N_1] r_1_1;  // actual group-level effects
  real lprior = 0;  // prior contributions to the log posterior
  // join the spatial and the non-spatial CAR component
  rcar = (sqrt(1 - rhocar) * nszcar + sqrt(rhocar * inv(car_scale)) * zcar) * sdcar;
  // compute penalized spline coefficients
  s_1_1 = sds_1[1] * zs_1_1;
  r_1_1 = (sd_1[1] * (z_1[1]));
  lprior += student_t_lpdf(Intercept | 3, 0, 2.5);
  lprior += student_t_lpdf(sds_1 | 3, 0, 2.5)
    - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  lprior += student_t_lpdf(sdcar | 3, 0, 2.5)
    - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  lprior += beta_lpdf(rhocar | 1, 1);
  lprior += student_t_lpdf(sd_1 | 3, 0, 2.5)
    - 1 * student_t_lccdf(0 | 3, 0, 2.5);
}
model {
  // likelihood including constants
  if (!prior_only) {
    // initialize linear predictor term
    vector[N] mu = rep_vector(0.0, N);
    mu += Intercept + Xc * b + Zs_1_1 * s_1_1;
    for (n in 1:N) {
      // add more terms to the linear predictor
      mu[n] += rcar[Jloc[n]] + r_1_1[J_1[n]] * Z_1_1[n] + rw_d[Jloc[n], time_id[n]];
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
  target += std_normal_lpdf(z_1[1]);
  // RW1 prior per district
  for (d in 1:Nloc) {
    target += normal_lpdf(rw_d[d, 2:T] - rw_d[d, 1:(T-1)] | 0, sdrw_d);
    target += normal_lpdf(sum(rw_d[d]) | 0, 0.001 * T);
  }
  target += student_t_lpdf(sdrw_d | 3, 0, 2.5);
}
generated quantities {
  // actual population-level intercept
  real b_Intercept = Intercept - dot_product(means_X, b);
}
