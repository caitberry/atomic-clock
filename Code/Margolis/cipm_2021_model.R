library(rstan)
library(stringr)
library(Rmpfr)

# =====================================================================
# 1. READ AND PARSE THE RAW DATA
# =====================================================================
raw_lines <- readLines("C:/Users/aak3/Documents/atomic-clock/Data/Margolis2024DataMetrologia/ClockInputData2020_7.dat")

N <- 106
N_transitions <- 14
Q <- numeric(N)
U <- numeric(N)
num_idx <- integer(N)
den_idx <- integer(N)

# The 2017 CIPM Reference Frequencies (nu_0)
# Index 15 is Cesium (exactly 9192631770 Hz)
nu_0 <- c(
  1267402452901040, 1233030706593518, 1128575290808154, 1121015393207851,
  1064721609899145, 688358979309308,  642121496772645,  518295836590864,
  455986240494140,  444779044095485,  429228066418008,  429228004229873,
  411042129776395,  6834682610.90431, 9192631770
)

# Parse measurements (q1 to q106)
q_lines <- raw_lines[str_detect(raw_lines, "^q\\d+")]
for (line in q_lines) {
  parts <- str_split(line, "\\s+")[[1]]
  parts <- parts[parts != ""]
  q_idx <- as.integer(str_extract(parts[1], "\\d+"))
  
  trans_str <- parts[2]
  if (str_detect(trans_str, "_over_")) {
    num_idx[q_idx] <- as.integer(str_extract(str_split(trans_str, "_over_")[[1]][1], "\\d+"))
    den_idx[q_idx] <- as.integer(str_extract(str_split(trans_str, "_over_")[[1]][2], "\\d+"))
  } else {
    num_idx[q_idx] <- as.integer(str_extract(trans_str, "\\d+"))
    den_idx[q_idx] <- 15 # Implicitly relative to Cesium
  }
  
  Q[q_idx] <- as.numeric(parts[3])
  U[q_idx] <- as.numeric(parts[4])
}

# =====================================================================
# 2. TRANSFORM TO SCALED FRACTIONAL OFFSETS
# =====================================================================
Q_scaled <- numeric(N)
U_scaled <- numeric(N)

for (i in 1:N) {
  num_0 <- nu_0[num_idx[i]]
  if (den_idx[i] == 15) {
    frac_offset <- (Q[i] - num_0) / num_0
    frac_u <- U[i] / num_0
  } else {
    den_0 <- nu_0[den_idx[i]]
    R_0 <- num_0 / den_0
    frac_offset <- (Q[i] - R_0) / R_0
    frac_u <- U[i] / R_0
  }
  # Scale by 1e15 to center the data around ~1.0 for the sampler
  Q_scaled[i] <- frac_offset * 1e15
  U_scaled[i] <- frac_u * 1e15
}

# Build Scaled Covariance Matrix
V_scaled <- diag(U_scaled^2)
r_lines <- raw_lines[str_detect(raw_lines, "^r\\(q\\d+,q\\d+\\)")]

for (line in r_lines) {
  parts <- str_split(line, "\\s+")[[1]]
  indices <- str_extract_all(parts[1], "\\d+")[[1]]
  i <- as.integer(indices[1])
  j <- as.integer(indices[2])
  cor_val <- as.numeric(parts[2])
  
  # cov = r * u_i * u_j (using scaled uncertainties)
  cov_val <- cor_val * U_scaled[i] * U_scaled[j]
  V_scaled[i, j] <- cov_val
  V_scaled[j, i] <- cov_val
}

# =====================================================================
# 3. RUN THE STAN MODEL
# =====================================================================
stan_data <- list(
  N_measurements = N,
  N_transitions = N_transitions,
  Q_scaled = Q_scaled,
  V_scaled = V_scaled,
  num_idx = num_idx,
  den_idx = den_idx
)

cat("Starting Stan MCMC sampling...\n")
fit <- stan(
  file = "C:/Users/aak3/Documents/atomic-clock/Code/Margolis/cipm_2021_model.stan",
  data = stan_data,
  iter = 4000,
  warmup = 1000,
  chains = 4,
  cores = 4, # Run chains in parallel if your machine supports it
  control = list(adapt_delta = 0.99) # Helps traverse highly correlated ridges safely
)

# =====================================================================
# 4. HIGH-PRECISION RECONSTRUCTION
# =====================================================================
# Extract the posterior samples for the scaled offsets (x)
posterior_x <- rstan::extract(fit)$x

# Calculate mean and standard deviation of the scaled offsets
x_mean <- colMeans(posterior_x)
x_sd <- apply(posterior_x, 2, sd)

cat("\n======================================================\n")
cat("FINAL BAYESIAN ABSOLUTE FREQUENCIES (128-bit precision)\n")
cat("======================================================\n\n")

# Use 128-bit precision for the reconstruction to safely hold >25 decimal digits
nu_0_mpfr <- mpfr(nu_0[1:N_transitions], precBits = 128)
x_mean_mpfr <- mpfr(x_mean, precBits = 128)
x_sd_mpfr <- mpfr(x_sd, precBits = 128)
scaling_factor <- mpfr(1e-15, precBits = 128)

# Reconstruct: nu = nu_0 * (1 + x * 1e-15)
nu_estimate_mpfr <- nu_0_mpfr * (1 + x_mean_mpfr * scaling_factor)
nu_unc_mpfr <- nu_0_mpfr * (x_sd_mpfr * scaling_factor)

labels <- c(
  "115In+", "1H", "199Hg", "27Al+", "199Hg+", "171Yb+(E2)", 
  "171Yb+(E3)", "171Yb", "40Ca", "88Sr+", "88Sr", "87Sr", "40Ca+", "87Rb"
)

# Print out the results with exactly 3 decimal places for Hz
for (k in 1:N_transitions) {
  # Format using Rmpfr's format method to prevent truncation
  est_str <- format(nu_estimate_mpfr[k], digits = 18, nsmall = 3)
  unc_str <- format(nu_unc_mpfr[k], digits = 4, scientific = TRUE)
  
  cat(sprintf("%-12s: %s Hz  (Unc: %s Hz)\n", labels[k], est_str, unc_str))
}
