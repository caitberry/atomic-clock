library(readr)
library(dplyr)
library(rstan)

raw_data <- read_csv("Code/ComparisonAnalysis2025/BayesianClockComparison/bayesian_model_input_data.csv")
df_1d <- raw_data[1:4, ]
df_3d <- raw_data[5:nrow(raw_data), ]

# --- Generate Indices ---
# 1. Define the unique days for B/C (All days present in the file)
all_dates <- unique(raw_data$date) # Should be length 13
N_BC <- length(all_dates)

# 2. Define the unique days for A (Only days present in 3D data)
# Since A is required for X and Y observations, it only exists in df_3d rows
dates_A <- unique(df_3d$date)      # Should be length 9
N_A <- length(dates_A)

# 3. Create mapping functions
# match(x, table) returns the position of x in table
get_idx_BC <- function(d) match(d, all_dates)
get_idx_A  <- function(d) match(d, dates_A)

# --- Prepare Data List ---
stan_data <- list(
  # Dimensions
  N_1d = nrow(df_1d),
  N_3d = nrow(df_3d),
  N_A  = N_A,
  N_BC = N_BC,
  
  # Indices
  # For 1D: only B and C needed
  idx_B_1d = get_idx_BC(df_1d$date),
  idx_C_1d = get_idx_BC(df_1d$date),
  
  # For 3D: A, B, and C needed
  idx_A_3d = get_idx_A(df_3d$date),
  idx_B_3d = get_idx_BC(df_3d$date),
  idx_C_3d = get_idx_BC(df_3d$date),
  
  # 1D Data
  z_obs_1d   = df_1d$z,
  b_z_1d     = df_1d$b_z,
  c_z_1d     = df_1d$c_z,
  sigma_z_1d = df_1d$sigma_z,
  
  # 3D Data
  obs_3d     = as.matrix(df_3d[, c("x", "y", "z")]), 
  a_x_3d = df_3d$a_x, c_x_3d = df_3d$c_x,
  a_y_3d = df_3d$a_y, b_y_3d = df_3d$b_y,
  b_z_3d = df_3d$b_z, c_z_3d = df_3d$c_z,
  
  # Covariances (Handle NA check for rho if needed)
  sigma_x_3d = df_3d$sigma_x,
  sigma_y_3d = df_3d$sigma_y,
  sigma_z_3d = df_3d$sigma_z,
  cov_xy_3d = ifelse(is.na(df_3d$cov_xy), df_3d$rho_xy*df_3d$sigma_x*df_3d$sigma_y, df_3d$cov_xy),
  cov_xz_3d = ifelse(is.na(df_3d$cov_xz), df_3d$rho_xz*df_3d$sigma_x*df_3d$sigma_z, df_3d$cov_xz),
  cov_yz_3d = ifelse(is.na(df_3d$cov_yz), df_3d$rho_yz*df_3d$sigma_y*df_3d$sigma_z, df_3d$cov_yz)
)

# 5. Run
fit <- stan(
  file = "Code/ComparisonAnalysis2025/BayesianClockComparison/allRatioAnalysis_withCov_lambdas.stan", 
  data = stan_data, 
  warmup = 15000,
  thin = 10,
  iter = 30000, 
  chains = 4, 
  cores = 4#,control = list(adapt_delta=.99)
)

# Inspect results
print(fit)#, pars = c("mu_x", "mu_y", "mu_z", "xi_x", "alpha", "beta", "gamma",""))

posterior_samples <- as.data.frame(fit)

pairs(fit, pars = c("mu_x", "mu_y", "mu_z"))
pairs(fit, pars = c("xi_A", "xi_B", "xi_C"))
traceplot(object = fit,pars = c("lambda_C"))
pairs(fit, pars = c("lambda_C"))

cor(posterior_samples$mu_x,posterior_samples$mu_y)
cor(posterior_samples$mu_x,posterior_samples$mu_z)
cor(posterior_samples$mu_y,posterior_samples$mu_z)

