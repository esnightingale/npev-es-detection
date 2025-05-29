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
  // data for spline 1
  int nb_1;  // number of bases
  array[nb_1] int knots_1;  // number of knots
  // basis function matrices
  matrix[N, knots_1[1]] Zs_1_1;
  // data for spline 2
  int nb_2;  // number of bases
  array[nb_2] int knots_2;  // number of knots
  // basis function matrices
  matrix[N, knots_2[1]] Zs_2_1;
  // data for spline 3
  int nb_3;  // number of bases
  array[nb_3] int knots_3;  // number of knots
  // basis function matrices
  matrix[N, knots_3[1]] Zs_3_1;
  // data for spline 4
  int nb_4;  // number of bases
  array[nb_4] int knots_4;  // number of knots
  // basis function matrices
  matrix[N, knots_4[1]] Zs_4_1;
  // data for spline 5
  int nb_5;  // number of bases
  array[nb_5] int knots_5;  // number of knots
  // basis function matrices
  matrix[N, knots_5[1]] Zs_5_1;
  // data for spline 6
  int nb_6;  // number of bases
  array[nb_6] int knots_6;  // number of knots
  // basis function matrices
  matrix[N, knots_6[1]] Zs_6_1;
  // data for spline 7
  int nb_7;  // number of bases
  array[nb_7] int knots_7;  // number of knots
  // basis function matrices
  matrix[N, knots_7[1]] Zs_7_1;
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
  // parameters for spline 2
  // standardized penalized spline coefficients
  vector[knots_2[1]] zs_2_1;
  vector<lower=0>[nb_2] sds_2;  // SDs of penalized spline coefficients
  // parameters for spline 3
  // standardized penalized spline coefficients
  vector[knots_3[1]] zs_3_1;
  vector<lower=0>[nb_3] sds_3;  // SDs of penalized spline coefficients
  // parameters for spline 4
  // standardized penalized spline coefficients
  vector[knots_4[1]] zs_4_1;
  vector<lower=0>[nb_4] sds_4;  // SDs of penalized spline coefficients
  // parameters for spline 5
  // standardized penalized spline coefficients
  vector[knots_5[1]] zs_5_1;
  vector<lower=0>[nb_5] sds_5;  // SDs of penalized spline coefficients
  // parameters for spline 6
  // standardized penalized spline coefficients
  vector[knots_6[1]] zs_6_1;
  vector<lower=0>[nb_6] sds_6;  // SDs of penalized spline coefficients
  // parameters for spline 7
  // standardized penalized spline coefficients
  vector[knots_7[1]] zs_7_1;
  vector<lower=0>[nb_7] sds_7;  // SDs of penalized spline coefficients
  real<lower=0> sdcar;  // SD of the CAR structure
  // parameters for the BYM2 structure
  vector[Nloc] zcar;  // spatial part
  vector[Nloc] nszcar;  // non-spatial part
  // proportion of variance in the spatial part
  real<lower=0,upper=1> rhocar;
  real<lower=0> sigma_w;  // Random walk SD
}
transformed parameters {
  // penalized spline coefficients
  vector[knots_1[1]] s_1_1;
  // penalized spline coefficients
  vector[knots_2[1]] s_2_1;
  // penalized spline coefficients
  vector[knots_3[1]] s_3_1;
  // penalized spline coefficients
  vector[knots_4[1]] s_4_1;
  // penalized spline coefficients
  vector[knots_5[1]] s_5_1;
  // penalized spline coefficients
  vector[knots_6[1]] s_6_1;
  // penalized spline coefficients
  vector[knots_7[1]] s_7_1;
  // scaled parameters for the BYM2 structure
  vector[Nloc] rcar;
  real lprior = 0;  // prior contributions to the log posterior
  // join the spatial and the non-spatial CAR component
  rcar = (sqrt(1 - rhocar) * nszcar + sqrt(rhocar * inv(car_scale)) * zcar) * sdcar;
  // compute penalized spline coefficients
  s_1_1 = sds_1[1] * zs_1_1;
  // compute penalized spline coefficients
  s_2_1 = sds_2[1] * zs_2_1;
  // compute penalized spline coefficients
  s_3_1 = sds_3[1] * zs_3_1;
  // compute penalized spline coefficients
  s_4_1 = sds_4[1] * zs_4_1;
  // compute penalized spline coefficients
  s_5_1 = sds_5[1] * zs_5_1;
  // compute penalized spline coefficients
  s_6_1 = sds_6[1] * zs_6_1;
  // compute penalized spline coefficients
  s_7_1 = sds_7[1] * zs_7_1;
  lprior += student_t_lpdf(Intercept | 3, 0, 2.5);
  lprior += student_t_lpdf(sds_1 | 3, 0, 2.5)
    - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  lprior += student_t_lpdf(sds_2 | 3, 0, 2.5)
    - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  lprior += student_t_lpdf(sds_3 | 3, 0, 2.5)
    - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  lprior += student_t_lpdf(sds_4 | 3, 0, 2.5)
    - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  lprior += student_t_lpdf(sds_5 | 3, 0, 2.5)
    - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  lprior += student_t_lpdf(sds_6 | 3, 0, 2.5)
    - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  lprior += student_t_lpdf(sds_7 | 3, 0, 2.5)
    - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  lprior += student_t_lpdf(sdcar | 3, 0, 2.5)
    - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  lprior += beta_lpdf(rhocar | 1, 1);
}
model {
  // likelihood including constants
  if (!prior_only) {
    // initialize linear predictor term
    vector[N] mu = rep_vector(0.0, N);
    mu += Intercept + rw + Xc * b + Zs_1_1 * s_1_1 + Zs_2_1 * s_2_1 + Zs_3_1 * s_3_1 + Zs_4_1 * s_4_1 + Zs_5_1 * s_5_1 + Zs_6_1 * s_6_1 + Zs_7_1 * s_7_1;
    for (n in 1:N) {
      // add more terms to the linear predictor
      mu[n] += rcar[Jloc[n]];
    }
    target += binomial_logit_lpmf(Y | trials, mu);
  }
  // priors including constants
  target += lprior;
  target += std_normal_lpdf(zs_1_1);
  target += std_normal_lpdf(zs_2_1);
  target += std_normal_lpdf(zs_3_1);
  target += std_normal_lpdf(zs_4_1);
  target += std_normal_lpdf(zs_5_1);
  target += std_normal_lpdf(zs_6_1);
  target += std_normal_lpdf(zs_7_1);
  // improper prior on the spatial BYM2 component
  target += -0.5 * dot_self(zcar[edges1] - zcar[edges2]);
  // soft sum-to-zero constraint
  target += normal_lpdf(sum(zcar) | 0, 0.001 * Nloc);
  // proper prior on the non-spatial BYM2 component
  target += std_normal_lpdf(nszcar);
  // Random walk prior
  rw[1] ~ normal(0, 1);
  for(t in 2:N) {
    rw[t] ~ normal(rw[t-1], sigma_w);
  }
}
generated quantities {
  // actual population-level intercept
  real b_Intercept = Intercept + rw - dot_product(means_X, b);
}