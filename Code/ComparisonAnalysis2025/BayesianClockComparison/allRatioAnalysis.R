# Load libraries
library(rstan)
library(readr)

# Set Stan options (optional, but good for reproducibility and performance)
rstan_options(auto_write = TRUE) # Saves compiled Stan models to disk
options(mc.cores = parallel::detectCores()) # Use all available cores for parallel sampling

# --- 1. Load Data ---
df_AlSr <- read_csv("Data/ClockComparison2025/BayesianAnalysisData/ErYb_AlSr_data.csv")
df_AlYb <- read_csv("Data/ClockComparison2025/BayesianAnalysisData/ErYb_AlYb_data.csv")
df_YbSr <- read_csv("Data/ClockComparison2025/BayesianAnalysisData/ErYb_YbSr_data.csv")


# --- 2. Prepare Data for Stan ---
# These names must exactly match the 'data' block in your Stan model.
stan_data <- list(
  N_AlSr = nrow(df_AlSr),
  N_AlYb = nrow(df_AlYb),
  N_YbSr = nrow(df_YbSr),
  
  x_j_data = df_AlSr$offset,
  y_j_data = df_AlYb$offset,
  z_j_data = df_YbSr$offset,
  
  sigma_x_j_data = df_AlSr$statistical_unc,
  sigma_y_j_data = df_AlYb$statistical_unc,
  sigma_z_j_data = df_YbSr$statistical_unc,
  
  a_x_j_data = df_AlSr$systematic_unc_numerator,
  a_y_j_data = df_AlYb$systematic_unc_numerator,
  a_z_j_data = df_YbSr$systematic_unc_numerator,
  
  b_x_j_data = df_AlSr$systematic_unc_denominator,
  b_y_j_data = df_AlYb$systematic_unc_denominator,
  b_z_j_data = df_YbSr$systematic_unc_denominator
)

# --- 3. Compile and Run the Stan Model ---
# This step might take a few minutes the first time as it compiles the C++ code.
# Subsequent runs will be faster due to auto_write = TRUE.
cat("Compiling Stan model...\n")
stan_model <- stan_model("Code/ComparisonAnalysis2025/BayesianClockComparison/allRatioAnalysis.stan")
cat("Starting MCMC sampling...\n")

fit <- sampling(
  stan_model,
  data = stan_data,
  chains = 4,         # Number of independent MCMC chains
  iter = 100000,        # Total iterations per chain (warmup + sampling)
  warmup = 50000,      # Number of warmup iterations per chain
  thin = 1,           # Thinning interval
  # control = list(adapt_delta = 0.95) # Adjust adapt_delta if you get 'divergent transitions' warnings
  # seed = 1234         # For reproducibility
)

# --- 4. Analyze Results ---

# Print summary of the fit
print(fit, pars = c("mu_x", "mu_y", "mu_z", "xi_x", "xi_y", "xi_z",
                    "alpha", "beta", "gamma", "eta_x", "eta_y", "eta_z"),
      probs = c(0.025, 0.5, 0.975))

# Extract samples (e.g., for further custom analysis or plotting)
posterior_samples <- as.data.frame(fit)

pairs(fit, pars = c("mu_x", "mu_y", "mu_z"))

cor(posterior_samples$mu_x,posterior_samples$mu_y)
cor(posterior_samples$mu_x,posterior_samples$mu_z)
cor(posterior_samples$mu_y,posterior_samples$mu_z)

cor(posterior_samples$eta_x,posterior_samples$eta_y)
cor(posterior_samples$eta_x,posterior_samples$eta_z)
cor(posterior_samples$eta_y,posterior_samples$eta_z)


mean(posterior_samples$alpha)
mean(posterior_samples$beta)
mean(posterior_samples$gamma)
