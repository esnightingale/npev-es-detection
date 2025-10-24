//

data {
    int<lower=1> N_districts;
    int<lower=0> N_edges;
    int<lower=1, upper=N_districts> node1[N_edges];
    int<lower=1, upper=N_districts> node2[N_edges];
}

parameters {
    vector[N_districts] Z;
    real<lower=0> tau_Z;  // Precision parameter
}

model {
    tau_Z ~ normal(0, 1);
    Z ~ normal(0, 1); // Weak prior to center

    // ICAR model (sum-to-zero constraint)
    for (i in 1:N_edges) {
        target += -0.5 * tau_Z * square(Z[node1[i]] - Z[node2[i]]);
    }
    sum(Z) ~ normal(0, 0.001 * N_districts);
}
