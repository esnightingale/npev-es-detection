functions {
}
data {
  int<lower=1> N;  // total number of observations
  array[N] int Y;  // response variable
  array[N] int trials;  // number of trials
  int<lower=1> K;  // number of population-level effects
  matrix[N, K] X;  // population-level design matrix
  int<lower=1> Kc;  // number of population-level effects after centering
  // data for splines
  int Ks;  // number of linear effects
  matrix[N, Ks] Xs;  // population-level design matrix for splines
  // data for spline 1
  int nb_1;  // number of bases
  array[nb_1] int knots_1;  // number of knots
  matrix[N, knots_1[1]] Zs_1_1;
  // data for spline 2
  int nb_2;
  array[nb_2] int knots_2;
  matrix[N, knots_2[1]] Zs_2_1;
  // data for spline 3
  int nb_3;
  array[nb_3] int knots_3;
  matrix[N, knots_3[1]] Zs_3_1;
  // data for spline 4
  int nb_4;
  array[nb_4] int knots_4;
  matrix[N, knots_4[1]] Zs_4_1;
  // data for spline 5
  int nb_5;
  array[nb_5] int knots_5;
  matrix[N, knots_5[1]] Zs_5_1;
  // data for spline 6
  int nb_6;
  array[nb_6] int knots_6;
  matrix[N, knots_6[1]] Zs_6_1;
  // data for spline 7
  int nb_7;
  array[nb_7] int knots_7;
  matrix[N, knots_7[1]] Zs_7_1;
  // data for the CAR structure
  int<lower=1> Nloc;
  array[N] int<lower=1> Jloc;
  int<lower=0> Nedges;
  array[Nedges] int<lower=1> edges1;
  array[Nedges] int<lower=1> edges2;
  real<lower=0> car_scale;
  // ---- Random walk (RW1) indexing ----
  int<lower=1> N_rw;
  array[N] int<lower=1> J_rw;
  // ------------------------------------
  int prior_only;
}
transformed data {
  matrix[N, Kc] Xc;
  vector[Kc] means_X;
  for (i in 2:K) {
    means_X[i - 1] = mean(X[, i]);
    Xc[, i - 1] = X[, i] - means_X[i - 1];
  }
}
parameters {
  vector[Kc] b;
  real Intercept;
  vector[Ks] bs;
  vector[knots_1[1]] zs_1_1;
  vector<lower=0>[nb_1] sds_1;
  vector[knots_2[1]] zs_2_1;
  vector<lower=0>[nb_2] sds_2;
  vector[knots_3[1]] zs_3_1;
  vector<lower=0>[nb_3] sds_3;
  vector[knots_4[1]] zs_4_1;
  vector<lower=0>[nb_4] sds_4;
  vector[knots_5[1]] zs_5_1;
  vector<lower=0>[nb_5] sds_5;
  vector[knots_6[1]] zs_6_1;
  vector<lower=0>[nb_6] sds_6;
  vector[knots_7[1]] zs_7_1;
  vector<lower=0>[nb_7] sds_7;
  real<lower=0> sdcar;
  vector[Nloc] zcar;
  vector[Nloc] nszcar;
  real<lower=0,upper=1> rhocar;
  vector[N_rw] rw;
  real<lower=0> sdrw;
}
transformed parameters {
  vector[knots_1[1]] s_1_1;
  vector[knots_2[1]] s_2_1;
  vector[knots_3[1]] s_3_1;
  vector[knots_4[1]] s_4_1;
  vector[knots_5[1]] s_5_1;
  vector[knots_6[1]] s_6_1;
  vector[knots_7[1]] s_7_1;
  vector[Nloc] rcar;
  real lprior = 0;
  rcar = (sqrt(1 - rhocar) * nszcar + sqrt(rhocar * inv(car_scale)) * zcar) * sdcar;
  s_1_1 = sds_1[1] * zs_1_1;
  s_2_1 = sds_2[1] * zs_2_1;
  s_3_1 = sds_3[1] * zs_3_1;
  s_4_1 = sds_4[1] * zs_4_1;
  s_5_1 = sds_5[1] * zs_5_1;
  s_6_1 = sds_6[1] * zs_6_1;
  s_7_1 = sds_7[1] * zs_7_1;
  lprior += student_t_lpdf(sdrw | 3, 0, 2.5) - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  lprior += student_t_lpdf(Intercept | 3, 0, 2.5);
  lprior += student_t_lpdf(sds_1 | 3, 0, 2.5) - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  lprior += student_t_lpdf(sds_2 | 3, 0, 2.5) - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  lprior += student_t_lpdf(sds_3 | 3, 0, 2.5) - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  lprior += student_t_lpdf(sds_4 | 3, 0, 2.5) - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  lprior += student_t_lpdf(sds_5 | 3, 0, 2.5) - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  lprior += student_t_lpdf(sds_6 | 3, 0, 2.5) - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  lprior += student_t_lpdf(sds_7 | 3, 0, 2.5) - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  lprior += student_t_lpdf(sdcar | 3, 0, 2.5) - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  lprior += beta_lpdf(rhocar | 1, 1);
}
model {
  if (!prior_only) {
    vector[N] mu = rep_vector(0.0, N);
    mu += Intercept + Xc * b + Xs * bs +
          Zs_1_1 * s_1_1 + Zs_2_1 * s_2_1 + Zs_3_1 * s_3_1 + Zs_4_1 * s_4_1 +
          Zs_5_1 * s_5_1 + Zs_6_1 * s_6_1 + Zs_7_1 * s_7_1;
    for (n in 1:N) {
      mu[n] += rcar[Jloc[n]];
      mu[n] += rw[J_rw[n]] - mean(rw);
    }
    target += binomial_logit_lpmf(Y | trials, mu);
  }
  target += lprior;
  target += std_normal_lpdf(zs_1_1);
  target += std_normal_lpdf(zs_2_1);
  target += std_normal_lpdf(zs_3_1);
  target += std_normal_lpdf(zs_4_1);
  target += std_normal_lpdf(zs_5_1);
  target += std_normal_lpdf(zs_6_1);
  target += std_normal_lpdf(zs_7_1);
  target += -0.5 * dot_self(zcar[edges1] - zcar[edges2]);
  target += normal_lpdf(sum(zcar) | 0, 0.001 * Nloc);
  target += std_normal_lpdf(nszcar);
  target += normal_lpdf(rw[2:N_rw] - rw[1:(N_rw-1)] | 0, sdrw);
  target += normal_lpdf(sum(rw) | 0, 0.001 * N_rw);
}
generated quantities {
  real b_Intercept = Intercept - dot_product(means_X, b);
  vector[N_rw] rw_centered;
  rw_centered = rw - mean(rw);
}
