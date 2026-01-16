library(rstan)

library(readr)
library(dplyr)

# 1. Read the Data
# -------------------------------------------------------------------
# Replace with your actual filename
raw_data <- read_csv("Code/ComparisonAnalysis2025/BayesianClockComparison/bayesian_model_input_data.csv")

# 2. Split the Data
# -------------------------------------------------------------------
# 1D Data: First 4 rows (contains z, b_z, c_z, sigma_z). 
# We ignore x, y, and other columns here as they likely contain NAs.
df_1d <- raw_data[1:4, ]

# 3D Data: Row 5 onwards
df_3d <- raw_data[5:nrow(raw_data), ]

# 4. Prepare the Data List for Stan
# -------------------------------------------------------------------
# We map the dataframe columns to the exact variable names in your 'data {}' block
stan_data <- list(
  # --- 1D Data Inputs ---
  N_1d       = nrow(df_1d),
  z_obs_1d   = df_1d$z,
  b_z_1d     = df_1d$b_z,
  c_z_1d     = df_1d$c_z,
  sigma_z_1d = df_1d$sigma_z,
  
  # --- 3D Data Inputs ---
  N_3d       = nrow(df_3d),
  
  # 'array[N] vector[3]' in Stan maps to an (N x 3) Matrix in R
  obs_3d     = as.matrix(df_3d[, c("x", "y", "z")]), 
  
  # Coefficients
  a_x_3d = df_3d$a_x,
  c_x_3d = df_3d$c_x,
  a_y_3d = df_3d$a_y,
  b_y_3d = df_3d$b_y,
  b_z_3d = df_3d$b_z,
  c_z_3d = df_3d$c_z,
  
  # Statistical Uncertainties (Sigmas)
  sigma_x_3d = df_3d$sigma_x,
  sigma_y_3d = df_3d$sigma_y,
  sigma_z_3d = df_3d$sigma_z,
  
  # Covariances
  cov_xy_3d = df_3d$cov_xy,
  cov_xz_3d = df_3d$cov_xz,
  cov_yz_3d = df_3d$cov_yz
)

# 5. Run the Model
# -------------------------------------------------------------------
# Assumes you saved your Stan code in a file named "model.stan"
fit <- stan(
  file = "Code/ComparisonAnalysis2025/BayesianClockComparison/allRatioAnalysis_withCov.stan", 
  data = stan_data, 
  iter = 2000, 
  chains = 4, 
  cores = 4
)

# 6. Inspect Results
# -------------------------------------------------------------------
print(fit, pars = c("mu_x", "mu_y", "mu_z", "xi_x", "alpha", "beta", "gamma"))

# Assuming you have df_1d and df_3d already loaded in your R environment

df_1d
df_3d

# Prepare the data list for Stan
stan_data <- list(
  # --- 1D Data ---
  N_1d = nrow(df_1d),
  z_obs_1d = df_1d$z,
  b_z_1d = df_1d$b_z,
  c_z_1d = df_1d$c_z,
  sigma_z_1d = df_1d$sigma_z,
  
  # --- 3D Data ---
  N_3d = nrow(df_3d),
  # obs_3d must be a matrix or array (N_3d x 3)
  obs_3d = as.matrix(df_3d[, c("x", "y", "z")]), 
  
  # Coefficients
  a_x_3d = df_3d$a_x,
  c_x_3d = df_3d$c_x,
  a_y_3d = df_3d$a_y,
  b_y_3d = df_3d$b_y,
  b_z_3d = df_3d$b_z,
  c_z_3d = df_3d$c_z,
  
  # Variances/Covariances
  sigma_x_3d = df_3d$sigma_x,
  sigma_y_3d = df_3d$sigma_y,
  sigma_z_3d = df_3d$sigma_z,
  cov_xy_3d = df_3d$cov_xy,
  cov_xz_3d = df_3d$cov_xz,
  cov_yz_3d = df_3d$cov_yz
)

# Compile and Run the model
fit <- stan(
  file = "model.stan",   # The file name you saved above
  data = stan_data,
  iter = 2000,
  chains = 4,
  cores = 4              # Parallel processing
)

# Inspect results
print(fit)#, pars = c("mu_x", "mu_y", "mu_z", "xi_x", "alpha", "beta", "gamma",""))

posterior_samples <- as.data.frame(fit)

pairs(fit, pars = c("mu_x", "mu_y", "mu_z"))

cor(posterior_samples$mu_x,posterior_samples$mu_y)
cor(posterior_samples$mu_x,posterior_samples$mu_z)
cor(posterior_samples$mu_y,posterior_samples$mu_z)

