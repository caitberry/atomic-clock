data {
  int<lower=1> N;
  vector<lower=0>[N] f;      // frequency
  vector<lower=0>[N] S;      // measured PSD
}
parameters {
  vector<lower=0>[3] A;       // amplitudes for white, flicker, RW
  vector<lower=0>[3] fc;      // corner frequencies
  real<lower=0> sigma;        // log noise
}
model {
  vector[N] log_S_model;

  // Priors (tune as needed)
  A ~ lognormal(-20, 2);
  fc ~ lognormal(0, 2);
  sigma ~ normal(0, 1);

  for (i in 1:N) {
    real S_i = 0;
    S_i += (A[1]) / (1 + square(f[i] / fc[1]));                   // white
    S_i += (A[2] * inv(f[i])) / (1 + square(f[i] / fc[2]));       // flicker
    S_i += (A[3] * inv_square(f[i])) / (1 + square(f[i] / fc[3])); // random walk
    log_S_model[i] = log(S_i);
  }

  log(S) ~ normal(log_S_model, sigma);
}
