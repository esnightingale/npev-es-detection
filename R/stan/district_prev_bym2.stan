// generated with brms 2.22.0
functions {
}
data {
  int<lower=1> N;  // total number of observations
  array[N] int Y;  // response variable
  array[N] int trials;  // number of trials
  // data for spline 1
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
  int prior_only;  // should the likelihood be ignored?
}
transformed data {
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
  lprior += normal_lpdf(Intercept | 0, 1);
  lprior += exponential_lpdf(sds_1 | 0.5);
  lprior += exponential_lpdf(sdcar | 0.5);
  lprior += beta_lpdf(rhocar | 0.5, 0.5);
}
model {
  // likelihood including constants
  if (!prior_only) {
    // initialize linear predictor term
    vector[N] mu = rep_vector(0.0, N);
    mu += Intercept + Zs_1_1 * s_1_1;
    for (n in 1:N) {
      // add more terms to the linear predictor
      mu[n] += rcar[Jloc[n]];
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
}
generated quantities {
  // actual population-level intercept
  real b_Intercept = Intercept;
}

