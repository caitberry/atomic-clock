# Simulate data and run Bayesian model with sum-to-zero constraint and clock-specific uncertainties in R with Stan

library(rstan)
set.seed(123)

# Settings
D <- 9       # days
K <- 3       # clocks
N <- D * choose(K, 2)  # all pairwise ratios per day

# True between-day std dev per clock
tau_true <- c(0.2, 0.3, 0.4)

# Simulate latent clock values (sum-to-zero constraint)
theta_raw <- matrix(rnorm(D * (K - 1), sd = rep(tau_true[1:(K - 1)], each = D)), D, K - 1)
theta <- matrix(0, D, K)
for (d in 1:D) {
  for (k in 1:(K - 1)) {
    theta[d, k] <- theta_raw[d, k]
  }
  theta[d, K] <- -sum(theta[d, 1:(K - 1)])  # enforce sum-to-zero
}

# Clock-specific known uncertainties (vary by day)
sigma_clock <- matrix(runif(D * K, 0.01, 0.03), D, K)  # small uncertainties per clock per day

# Observation noise (e.g., instrument noise)
sigma_obs <- 0.01

# Construct observed log ratios (all pairs, per day)
y_obs <- numeric(N)
sigma2_obs <- rep(sigma_obs^2, N)
k1 <- integer(N)
k2 <- integer(N)
day <- integer(N)
i <- 1
for (d in 1:D) {
  for (a in 1:(K - 1)) {
    for (b in (a + 1):K) {
      total_sd <- sqrt(sigma_obs^2 + sigma_clock[d, a]^2 + sigma_clock[d, b]^2)
      y_obs[i] <- theta[d, a] - theta[d, b] + rnorm(1, sd = total_sd)
      k1[i] <- a
      k2[i] <- b
      day[i] <- d
      i <- i + 1
    }
  }
}

# Stan model with clock-specific uncertainty
stan_code <- "
data {
  int<lower=1> N;
  int<lower=1> D;
  int<lower=1> K;
  vector[N] y_obs;
  vector[N] sigma2_obs;
  matrix[D, K] sigma2_clock;
  int<lower=1, upper=D> day[N];
  int<lower=1, upper=K> k1[N];
  int<lower=1, upper=K> k2[N];
}
parameters {
  matrix[D, K - 1] theta_raw;
  cholesky_factor_corr[K] L_corr;
  vector<lower=0>[K] tau;
}
transformed parameters {
  matrix[D, K] theta;
  matrix[K, K] L_Sigma;
  L_Sigma = diag_pre_multiply(tau, L_corr);

  for (d in 1:D) {
    for (k in 1:(K - 1))
      theta[d, k] = theta_raw[d, k];
    theta[d, K] = -sum(theta_raw[d]);
  }
}
model {
  tau ~ cauchy(0, 0.5);
  L_corr ~ lkj_corr_cholesky(2.0);
  for (d in 1:D)
    theta[d] ~ multi_normal_cholesky(rep_vector(0, K), L_Sigma);

  for (n in 1:N) {
    real var_total = sigma2_obs[n] + sigma2_clock[day[n], k1[n]] + sigma2_clock[day[n], k2[n]];
    y_obs[n] ~ normal(theta[day[n], k1[n]] - theta[day[n], k2[n]], sqrt(var_total));
  }
}
"

# Compile and run
stan_data <- list(N = N, D = D, K = K, y_obs = y_obs, sigma2_obs = sigma2_obs,
                  sigma2_clock = sigma_clock^2, day = day, k1 = k1, k2 = k2)
fit <- stan(model_code = stan_code, data = stan_data, seed = 42, iter = 2000, chains = 4)

print(fit, pars = c("tau"))
