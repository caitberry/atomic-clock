data {
  // --- Dimensions ---
  int<lower=0> N_1d;              // Number of 1D observations
  int<lower=0> N_3d;              // Number of 3D observations
  
  // separate counts for the random effects
  int<lower=1> N_A;               // e.g., 9
  int<lower=1> N_BC;              // e.g., 13

  // --- Indices ---
  // For 1D data (only needs B and C)
  array[N_1d] int<lower=1> idx_B_1d; 
  array[N_1d] int<lower=1> idx_C_1d; 

  // For 3D data (needs A, B, and C)
  array[N_3d] int<lower=1> idx_A_3d; 
  array[N_3d] int<lower=1> idx_B_3d; 
  array[N_3d] int<lower=1> idx_C_3d; 

  // --- 1D Data Inputs ---
  vector[N_1d] z_obs_1d;
  vector[N_1d] b_z_1d;
  vector[N_1d] c_z_1d;
  vector[N_1d] sigma_z_1d;

  // --- 3D Data Inputs ---
  array[N_3d] vector[3] obs_3d;

  vector[N_3d] a_x_3d;
  vector[N_3d] c_x_3d;
  vector[N_3d] a_y_3d;
  vector[N_3d] b_y_3d;
  vector[N_3d] b_z_3d;
  vector[N_3d] c_z_3d;

  vector[N_3d] sigma_x_3d;
  vector[N_3d] sigma_y_3d;
  vector[N_3d] sigma_z_3d;
  vector[N_3d] cov_xy_3d;
  vector[N_3d] cov_xz_3d;
  vector[N_3d] cov_yz_3d;
}

transformed data {
  real sigma_eta_x = sqrt(0.2^2 + 0.4^2);
  real sigma_eta_y = sqrt(0.3^2 + 0.2^2);
  real sigma_eta_z = sqrt(0.2^2 + 0.4^2);
}

parameters {
  real mu_x;
  real mu_y;
  real mu_z;

  real<lower=0> xi_A;
  real<lower=0> xi_B;
  real<lower=0> xi_C;

  vector[N_A] lambda_A_raw;
  vector[N_BC] lambda_B_raw;
  vector[N_BC] lambda_C_raw;

  real alpha;
  real beta;
  real gamma;

  real eta_x;
  real eta_y;
  real eta_z;
}

transformed parameters {
  // lambda = raw * scale (xi)
  vector[N_A] lambda_A = lambda_A_raw * xi_A;
  vector[N_BC] lambda_B = lambda_B_raw * xi_B;
  vector[N_BC] lambda_C = lambda_C_raw * xi_C;
}

model {
  // --- Priors ---
  mu_x ~ normal(0, 10000);
  mu_y ~ normal(0, 10000);
  mu_z ~ normal(0, 10000);

  xi_A ~ normal(500, 1000);
  xi_B ~ normal(500, 1000);
  xi_C ~ normal(500, 1000);

  lambda_A_raw ~ std_normal();
  lambda_B_raw ~ std_normal();
  lambda_C_raw ~ std_normal();

  alpha ~ normal(0, 1);
  beta  ~ normal(0, 1);
  gamma ~ normal(0, 1);

  eta_x ~ normal(mu_x, sigma_eta_x);
  eta_y ~ normal(mu_y, sigma_eta_y);
  eta_z ~ normal(mu_z, sigma_eta_z);

  // --- Likelihood 1: 1D Data Points ---
  for (i in 1:N_1d) {
    // Look up the specific index for B and C for this row
    int iB = idx_B_1d[i];
    int iC = idx_C_1d[i];

    real mu_val = eta_z + beta * b_z_1d[i] - gamma * c_z_1d[i] 
                  + lambda_B[iB] - lambda_C[iC];
    
    z_obs_1d[i] ~ normal(mu_val, sigma_z_1d[i]);
  }

  // --- Likelihood 2: 3D Data Points ---
  {
    vector[3] mu_vec;
    matrix[3, 3] Sigma;

    for (i in 1:N_3d) {
      // Look up indices for A, B, and C
      int iA = idx_A_3d[i];
      int iB = idx_B_3d[i];
      int iC = idx_C_3d[i];

      // Mean Vector using the specific indices
      mu_vec[1] = eta_x + alpha * a_x_3d[i] - gamma * c_x_3d[i] 
                  + lambda_A[iA] - lambda_C[iC];
                  
      mu_vec[2] = eta_y + alpha * a_y_3d[i] - beta  * b_y_3d[i] 
                  + lambda_A[iA] - lambda_B[iB];
                  
      mu_vec[3] = eta_z + beta  * b_z_3d[i] - gamma * c_z_3d[i] 
                  + lambda_B[iB] - lambda_C[iC];

      // Covariance Matrix
      Sigma[1, 1] = square(sigma_x_3d[i]);
      Sigma[2, 2] = square(sigma_y_3d[i]);
      Sigma[3, 3] = square(sigma_z_3d[i]);

      Sigma[1, 2] = cov_xy_3d[i];
      Sigma[2, 1] = cov_xy_3d[i];
      Sigma[1, 3] = cov_xz_3d[i];
      Sigma[3, 1] = cov_xz_3d[i];
      Sigma[2, 3] = cov_yz_3d[i];
      Sigma[3, 2] = cov_yz_3d[i];

      obs_3d[i] ~ multi_normal(mu_vec, Sigma);
    }
  }
}
