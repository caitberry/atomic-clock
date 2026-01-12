data {
  int<lower=1> N;              // number of observations
  int<lower=1> D;              // number of days
  int<lower=1> K;              // number of clocks
  vector[N] y_obs;             // observed log ratios
  vector[N] sigma2_obs;        // known observation variances
  int<lower=1, upper=D> day[N];// day index
  int<lower=1, upper=K> k1[N]; // clock 1
  int<lower=1, upper=K> k2[N]; // clock 2
}
parameters {
  matrix[D, K] theta;              // latent clock offsets per day
  cholesky_factor_corr[K] L_corr;  // correlation matrix (Cholesky factor)
  vector<lower=0>[K] tau;          // standard deviations
}
transformed parameters {
  matrix[K, K] L_Sigma;
  L_Sigma = diag_pre_multiply(tau, L_corr); // Cholesky of full covariance
}
model {
  // Priors
  tau ~ cauchy(0, 0.5);            // half-Cauchy on std dev
  L_corr ~ lkj_corr_cholesky(2.0); // weakly informative prior on correlation
  
  // Latent clock offsets per day
  for (d in 1:D)
    theta[d] ~ multi_normal_cholesky(rep_vector(0, K), L_Sigma);
  
  // Observed log ratios
  for (n in 1:N)
    y_obs[n] ~ normal(theta[day[n], k1[n]] - theta[day[n], k2[n]], sqrt(sigma2_obs[n]));
}
