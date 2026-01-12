data {
  // --- 1D Data (z only) ---
  int<lower=0> N_1d;              // Number of 1D observations
  vector[N_1d] z_obs_1d;          // Observed z
  vector[N_1d] b_z_1d;            // Coefficient b_z
  vector[N_1d] c_z_1d;            // Coefficient c_z
  vector[N_1d] sigma_z_1d;        // Statistical uncertainty sigma_z

  // --- 3D Data (x, y, z) ---
  int<lower=0> N_3d;              // Number of 3D observations
  array[N_3d] vector[3] obs_3d;   // Observed vectors [x, y, z]

  // Coefficients for 3D
  vector[N_3d] a_x_3d;
  vector[N_3d] c_x_3d;
  vector[N_3d] a_y_3d;
  vector[N_3d] b_y_3d;
  vector[N_3d] b_z_3d;
  vector[N_3d] c_z_3d;

  // Covariance inputs for 3D (Statistical errors/covariances)
  vector[N_3d] sigma_x_3d;
  vector[N_3d] sigma_y_3d;
  vector[N_3d] sigma_z_3d;
  vector[N_3d] cov_xy_3d;
  vector[N_3d] cov_xz_3d;
  vector[N_3d] cov_yz_3d;
}

transformed data {
  // Pre-calculate the fixed sigmas for the eta priors
  real sigma_eta_x = sqrt(0.2^2 + 0.4^2);
  real sigma_eta_y = sqrt(0.3^2 + 0.2^2);
  real sigma_eta_z = sqrt(0.2^2 + 0.4^2);
}

parameters {
  // "True means"
  real mu_x;
  real mu_y;
  real mu_z;

  // "Dark uncertainties"
  // lower=0 handles the truncation logic automatically for the prior
  // real<lower=0> xi_x;
  // real<lower=0> xi_y;
  // real<lower=0> xi_z;
  real<lower=0> xi_al;
  real<lower=0> xi_yb;
  real<lower=0> xi_sr;

  vector[N_3d] lambdaAl;
  vector[N_1d] lambdaYb;
  vector[N_1d] lambdaSr;
  
  // Systematic uncertainty coefficients
  real alpha;
  real beta;
  real gamma;

  // Latent variables (link + geopotential)
  real eta_x;
  real eta_y;
  real eta_z;
}

model {
  // --- Priors ---
  mu_x ~ normal(0, 10000);
  mu_y ~ normal(0, 10000);
  mu_z ~ normal(0, 10000);

  // Truncated Normal priors
  // In Stan, T[0,] is implied by the parameter declaration <lower=0>
  xi_al ~ normal(500, 1000); 
  xi_yb ~ normal(500, 1000);
  xi_sr ~ normal(500, 1000);
  
  lambdaAl ~ normal(0,xi_al);
  lambdaYb ~ normal(0,xi_yb);
  lambdaSr ~ normal(0,xi_sr);
  

  alpha ~ normal(0, 1);
  beta  ~ normal(0, 1);
  gamma ~ normal(0, 1);

  // Hierarchical priors for eta
  eta_x ~ normal(mu_x, sigma_eta_x);
  eta_y ~ normal(mu_y, sigma_eta_y);
  eta_z ~ normal(mu_z, sigma_eta_z);

  // --- Likelihood 1: 1D Data Points ---
  // Vectorized operations are faster here
  {
    vector[N_1d] mu_vec_1d;
    vector[N_1d] sigma_total_1d;

    mu_vec_1d = eta_z + beta * b_z_1d - gamma * c_z_1d + lambdaYb[1:4];
    // Combine statistical sigma and dark uncertainty xi
    sigma_total_1d = sqrt(square(sigma_z_1d) + square(xi_z));

    z_obs_1d ~ normal(mu_vec_1d, sigma_total_1d);
  }

  // --- Likelihood 2: 3D Data Points ---
  // Loop is required to construct specific covariance matrices per point
  {
    vector[3] mu_vec_3d;
    matrix[3, 3] Sigma;

    for (i in 1:N_3d) {
      // 1. Construct Mean Vector
      mu_vec_3d[1] = eta_x + alpha * a_x_3d[i] - gamma * c_x_3d[i];
      mu_vec_3d[2] = eta_y + alpha * a_y_3d[i] - beta  * b_y_3d[i];
      mu_vec_3d[3] = eta_z + beta  * b_z_3d[i] - gamma * c_z_3d[i];

      // 2. Construct Covariance Matrix
      // Add 'xi' terms (dark uncertainty) to diagonal elements only
      Sigma[1, 1] = square(sigma_x_3d[i]) + square(xi_x);
      Sigma[2, 2] = square(sigma_y_3d[i]) + square(xi_y);
      Sigma[3, 3] = square(sigma_z_3d[i]) + square(xi_z);

      // Off-diagonals (stat covariance only)
      Sigma[1, 2] = cov_xy_3d[i];
      Sigma[2, 1] = cov_xy_3d[i];
      
      Sigma[1, 3] = cov_xz_3d[i];
      Sigma[3, 1] = cov_xz_3d[i];

      Sigma[2, 3] = cov_yz_3d[i];
      Sigma[3, 2] = cov_yz_3d[i];

      // Multivariate Normal Likelihood
      obs_3d[i] ~ multi_normal(mu_vec_3d, Sigma);
    }
  }
}
