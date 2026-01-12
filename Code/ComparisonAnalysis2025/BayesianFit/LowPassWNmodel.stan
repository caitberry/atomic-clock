data {
  int<lower=1> N;            // number of frequency points
  vector<lower=0>[N] f;      // frequency values (Hz)
  vector<lower=0>[N] S;      // measured PSD values (e.g. Hz^2/Hz)
}

parameters {
  real<lower=0> S_w;         // white noise level
  real<lower=0> f_c;         // corner frequency (Hz)
  real<lower=0> sigma;       // observation noise in log domain
}

model {
  vector[N] log_S_model;

  // Prior (adjust as needed)
  S_w ~ lognormal(-20, 2);   // weak prior, adjust for scale
  f_c ~ lognormal(0, 1);     // prior on corner frequency
  sigma ~ normal(0, 1);      // log-noise SD

  // Model
  for (i in 1:N) {
    log_S_model[i] = log(S_w) - log1p(square(f[i] / f_c));
  }

  // Likelihood
  log(S) ~ normal(log_S_model, sigma);
}
