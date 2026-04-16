# install.packages("Rmpfr")
library(Rmpfr)
library(rstan)

# Set precision bits (128 bits gives ~38 decimal digits of precision, which is plenty)
prec <- 128

# 1. Define Nu_0 using STRINGS to prevent early truncation
nu_0 <- mpfr(c(
  "1128575290808154.4", # 1: Hg (nu3)
  "518295836590863.6",  # 2: Yb (nu8)
  "429228004229873.0"   # 3: Sr (nu12)
), precBits = prec)

# 2. Define Q using STRINGS
Q <- mpfr(c(
  "1128575290808154.62", # 1. Absolute Hg (q7)
  "518295836590863.59",  # 2. Absolute Yb (q24)
  "429228004229872.92",  # 3. Absolute Sr (q47)
  "2.17747319413456507",   # 4. Ratio Hg / Yb (q79)
  "1.20750703934333841",   # 5. Ratio Yb / Sr (q81)
  "2.62931420989890960"    # 6. Ratio Hg / Sr (q59)
), precBits = prec)

# 3. Define Uncertainties using STRINGS
U <- mpfr(c(
  "0.41",       # 1. Hg absolute unc
  "0.31",       # 2. Yb absolute unc
  "0.12",       # 3. Sr absolute unc
  "1.92e-16",   # 4. Hg/Yb unc 
  "3.40e-16",   # 5. Yb/Sr unc 
  "2.20e-16"    # 6. Hg/Sr unc 
), precBits = prec)

# 4. Calculate Reference Ratios (R_0)
# Because nu_0 is an mpfr object, R_0 automatically becomes an mpfr object
R_0 <- c(nu_0[1]/nu_0[2], nu_0[2]/nu_0[3], nu_0[1]/nu_0[3])

# 5. Calculate scaled fractional offsets (Y matrix equivalent)
scale_factor <- mpfr(1e15, precBits = prec)

Y_mpfr <- c(
  ((Q[1] - nu_0[1]) / nu_0[1]) * scale_factor,
  ((Q[2] - nu_0[2]) / nu_0[2]) * scale_factor,
  ((Q[3] - nu_0[3]) / nu_0[3]) * scale_factor,
  ((Q[4] - R_0[1]) / R_0[1]) * scale_factor,
  ((Q[5] - R_0[2]) / R_0[2]) * scale_factor,
  ((Q[6] - R_0[3]) / R_0[3]) * scale_factor
)

# 6. Calculate scaled fractional uncertainties
u_scaled_mpfr <- c(
  (U[1] / nu_0[1]) * scale_factor,
  (U[2] / nu_0[2]) * scale_factor,
  (U[3] / nu_0[3]) * scale_factor,
  (U[4] / R_0[1]) * scale_factor,
  (U[5] / R_0[2]) * scale_factor,
  (U[6] / R_0[3]) * scale_factor
)

# 7. Convert back to standard numeric for Stan
# Now that the precise subtractions are done and the numbers are scaled to O(1), 
# it is safe to cast them back to standard double-precision.
Y_stan <- as.numeric(Y_mpfr)
u_scaled_stan <- as.numeric(u_scaled_mpfr)

# 8. Construct Covariance Matrix (V) using the numeric values
V_stan <- diag(u_scaled_stan^2)

# Insert the 3 explicit off-diagonal correlations
V_stan[1, 3] <- 0.438 * u_scaled_stan[1] * u_scaled_stan[3]; V_stan[3, 1] <- V_stan[1, 3]
V_stan[2, 5] <- 0.088 * u_scaled_stan[2] * u_scaled_stan[5]; V_stan[5, 2] <- V_stan[2, 5]
V_stan[4, 6] <- 0.826 * u_scaled_stan[4] * u_scaled_stan[6]; V_stan[6, 4] <- V_stan[4, 6]

# Check the cleanly extracted data to send to rstan
print(Y_stan)
print(V_stan)

# Prepare the data list (Notice we pass V now)
stan_data <- list(
  N_clocks = 3,
  N_meas = 6,
  Y = Y,
  V = V, # Passing the raw covariance matrix
  is_ratio = c(0, 0, 0, 1, 1, 1),
  clock_1 = c(1, 2, 3, 1, 2, 1),
  clock_2 = c(0, 0, 0, 2, 3, 3)
)

# Compile and run
fit <- stan(
  file = "Code/Margolis/learningByDoing/clock_network_corr.stan", 
  data = stan_data, 
  iter = 20000, 
  warmup = 15000, 
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  chains = 4, 
  cores = 4,
  seed = 42
)

# Extract parameters
posterior_samples <- extract(fit)

# Get the fractional offsets just like before
X_hat_bayes <- colMeans(posterior_samples$x)
cov_X_bayes <- cov(posterior_samples$x)

# High-Precision Reconstruction
optimized_freqs_bayes <- nu_0 * (1 + X_hat_bayes * 1e-15)
optimized_unc_bayes <- nu_0 * (sqrt(diag(cov_X_bayes)) * 1e-15)

cat("Bayesian Optimized Frequencies:\n", sprintf("%.3f\n", optimized_freqs_bayes))
cat("Bayesian Uncertainties:\n", sprintf("%3f\n", optimized_unc_bayes))

# Extract the estimated dark uncertainty
dark_unc_estimate <- mean(posterior_samples$sigma_dark)
dark_unc_ci <- quantile(posterior_samples$sigma_dark, probs = c(0.025, 0.975))

cat("Estimated Global Dark Uncertainty (parts per 10^15):", sprintf("%.3f", dark_unc_estimate), "\n")
cat("95% Credible Interval:", sprintf("[%.3f, %.3f]", dark_unc_ci[1], dark_unc_ci[2]), "\n")

