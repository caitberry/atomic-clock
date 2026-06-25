# install.packages("Rmpfr")
library(Rmpfr)
library(rstan)
library(ggplot2)

# Set precision bits 
prec <- 128

# 1. Define Nu_0 
nu_0 <- mpfr(c(
  "1128575290808154.4", # 1: Hg 
  "518295836590863.6",  # 2: Yb 
  "429228004229873.0"   # 3: Sr 
), precBits = prec)

# 2. Define Q (All 46 Measurements)
Q <- mpfr(c(
  # [Hg Abs]
  "1128575290808155.1", "1128575290808154.62", 
  # [Yb Abs]
  "518295836590864.0", "518295836590863.1", "518295836590865.2", "518295836590863.5", "518295836590863.59", "518295836590863.38", "518295836590863.30", "518295836590863.71", "518295836590863.61", "518295836590863.54",
  # [Sr Abs]
  "429228004229874.0", "429228004229873.65", "429228004229873.6", "429228004229874.1", "429228004229872.9", "429228004229873.9", "429228004229872.0", "429228004229873.56", "429228004229873.7", "429228004229873.13", "429228004229873.10", "429228004229872.92", "429228004229872.97", "429228004229873.04", "429228004229872.97", "429228004229872.99", "429228004229873.1", "429228004229873.00", "429228004229873.082", "429228004229873.13", "429228004229873.19",
  # [Hg/Sr]
  "2.62931420989890960", "2.62931420989890915",
  # [Yb/Sr]
  "1.2075070393433412", "1.20750703934333776", "1.207507039343337749", "1.20750703934333790", "1.20750703934333841", "1.20750703934333805", "1.20750703934333738", "1.20750703934333858", "1.20750703934333782", "1.2075070393433378482",
  # [Hg/Yb]
  "2.17747319413456507"
), precBits = prec)

# 3. Define Uncertainties (All 46 Measurements)
U <- mpfr(c(
  "6.4", "0.41", 
  "28.0", "2.0", "0.7", "8.1", "0.31", "0.57", "0.378", "0.11", "0.13", "0.259", 
  "1.1", "0.37", "1.1", "2.4", "0.5", "1.4", "1.6", "0.49", "1.4", "0.17", "0.13", "0.12", "0.16", "0.11", "0.40", "18.0", "0.429", "0.0708", "0.0773", "0.4", "1.50e-01",
  "2.2e-16", "4.6e-16",
  "1.7e-15", "2.9e-16", "5.5e-17", "7.00e-16", "3.40e-16", "3.38e-16", "3.02e-16", "4.95e-16", "7.49e-16", "8.20e-18",
  "1.92e-16"
), precBits = prec)

# 4. Map Reference Values & Scale
Ref_Vals <- c(
  rep(nu_0[1], 2),               # Abs Hg
  rep(nu_0[2], 10),              # Abs Yb
  rep(nu_0[3], 21),              # Abs Sr
  rep(nu_0[1] / nu_0[3], 2),     # Ratio Hg/Sr
  rep(nu_0[2] / nu_0[3], 10),    # Ratio Yb/Sr
  rep(nu_0[1] / nu_0[2], 1)      # Ratio Hg/Yb
)

scale_factor <- mpfr(1e15, precBits = prec)

# Calculate vectorized offsets and uncertainties
Y_mpfr <- ((Q - Ref_Vals) / Ref_Vals) * scale_factor
u_scaled_mpfr <- (U / Ref_Vals) * scale_factor

# Convert to standard numeric for Stan
Y_stan <- as.numeric(Y_mpfr)
u_scaled_stan <- as.numeric(u_scaled_mpfr)

# 5. Construct Base Covariance Matrix (Diagonal)
V_stan <- diag(u_scaled_stan^2)

# Define your exact 46 labels to match your data vectors
q_labels <- c(
  "q6", "q7",                                                                        # Abs Hg
  "q20", "q21", "q22", "q23", "q24", "q25", "q70", "q75", "q76", "q89",              # Abs Yb
  "q36", "q37", "q38", "q39", "q40", "q41", "q42", "q43", "q44", "q45", "q46",       # Abs Sr (1-11)
  "q47", "q48", "q49", "q50", "q51", "q72", "q73", "q90", "q91", "q96",              # Abs Sr (12-21)
  "q59", "q60",                                                                      # Ratio Hg/Sr
  "q64", "q65", "q66", "q80", "q81", "q82", "q83", "q93", "q94", "q102",             # Ratio Yb/Sr
  "q79"                                                                              # Ratio Hg/Yb
)

# Read in the raw correlation data
# 1. Read the file directly into a vector of lines
corr_lines <- readLines("Data/Margolis2024DataMetrologia/ClockInputData2020_7.dat")

# 2. Loop through each line
for (line in corr_lines) {
  # This automatically skips everything in the file until it hits the "r(...,...)" lines
  if (!grepl("^r\\(", trimws(line))) next 
  
  # Split the line by tabs or spaces
  parts <- strsplit(trimws(line), "\t|\\s+")[[1]]
  pair_str <- parts[1]             # e.g., "r(q7,q47)"
  r_val <- as.numeric(parts[2])    # e.g., 0.438
  
  # Extract the two q-labels using string replacement
  qs <- gsub("r\\(|\\)", "", pair_str) # turns "r(q7,q47)" into "q7,q47"
  q_pair <- strsplit(qs, ",")[[1]]     # c("q7", "q47")
  
  # Find their numeric index (1 to 46) in our subset
  idx1 <- match(q_pair[1], q_labels)
  idx2 <- match(q_pair[2], q_labels)
  
  # If BOTH measurements are in our 46-element model, add them to the matrix
  if (!is.na(idx1) && !is.na(idx2)) {
    cov_val <- r_val * u_scaled_stan[idx1] * u_scaled_stan[idx2]
    V_stan[idx1, idx2] <- cov_val
    V_stan[idx2, idx1] <- cov_val
  }
}

# Quick check: Validate that at least some off-diagonals populated successfully
cat("Number of inserted off-diagonal correlations:", sum(V_stan[lower.tri(V_stan)] != 0), "\n")


# 6. Prepare Stan Mapping Arrays 
# We use rep() to build the arrays matching the exact order of your 46 measurements
stan_data <- list(
  N_clocks = 3,
  N_meas = 46,
  Y = Y_stan,
  V = V_stan, 
  
  # 0 for the 33 Absolutes, 1 for the 13 Ratios
  is_ratio = c(rep(0, 33), rep(1, 13)),
  
  # The primary clock (Absolute clock, or Numerator for ratios)
  clock_1 = c(
    rep(1, 2),   # Hg Abs
    rep(2, 10),  # Yb Abs
    rep(3, 21),  # Sr Abs
    rep(1, 2),   # Hg/Sr (Numerator Hg)
    rep(2, 10),  # Yb/Sr (Numerator Yb)
    rep(1, 1)    # Hg/Yb (Numerator Hg)
  ),
  
  # The secondary clock (0 for Absolutes, Denominator for ratios)
  clock_2 = c(
    rep(0, 33),  # Null for Abs
    rep(3, 2),   # Hg/Sr (Denominator Sr)
    rep(3, 10),  # Yb/Sr (Denominator Sr)
    rep(2, 1)    # Hg/Yb (Denominator Yb)
  )
)

# # 7. Compile and run
# fit <- stan(
#   file = "Code/Margolis/learningByDoing/clock_network_corr.stan", 
#   data = stan_data, 
#   iter = 20000,   # Lowered from 20000 for a quicker test run
#   warmup = 15000, 
#   chains = 4, 
#   cores = 4,
#   seed = 42
# )
# 
# fit

# 7. Compile and run
fit <- stan(
  file = "Code/Margolis/learningByDoing/clock_network_corr_darkEachClock.stan", 
  data = stan_data, 
  iter = 20000,   
  warmup = 15000, 
  chains = 4, 
  cores = 4,
  control = list(adapt_delta = 0.99, max_treedepth = 12),
  seed = 42
)

###think about prior for taus, right now N(0,1), saying I strongly believe the unaccounted noise is small (likely under 1 part in 10^{15} (or adjusted for my current scale))
fit


# Extract parameters
posterior_samples <- rstan::extract(fit)

# Get the fractional offsets just like before
X_hat_bayes <- colMeans(posterior_samples$x)
cov_X_bayes <- cov(posterior_samples$x)

# Extract the individual clock dark uncertainties
dark_unc_means <- colMeans(posterior_samples$sigma_dark)
dark_unc_ci_lower <- apply(posterior_samples$sigma_dark, 2, quantile, probs = 0.025)
dark_unc_ci_upper <- apply(posterior_samples$sigma_dark, 2, quantile, probs = 0.975)

clock_names <- c("Hg", "Yb", "Sr")

cat("\n--- CLOCK-SPECIFIC DARK UNCERTAINTIES (parts per 10^15) ---\n")
for(i in 1:3) {
  cat(sprintf("%s Clock: %.3f  (95%% CI: [%.3f, %.3f])\n", 
              clock_names[i], 
              dark_unc_means[i], 
              dark_unc_ci_lower[i], 
              dark_unc_ci_upper[i]))
}


# Extract the fractional offsets (x)
X_hat_bayes <- colMeans(posterior_samples$x)
cov_X_bayes <- cov(posterior_samples$x)

# 95% Credible Intervals for the fractional offsets
x_ci_lower <- apply(posterior_samples$x, 2, quantile, probs = 0.025)
x_ci_upper <- apply(posterior_samples$x, 2, quantile, probs = 0.975)

# Calculate the new Absolute Frequencies and Uncertainties
# (Converting back out of the 1e-15 scaling)
optimized_freqs_bayes <- nu_0 * (1 + X_hat_bayes * 1e-15)
optimized_unc_bayes <- nu_0 * (sqrt(diag(cov_X_bayes)) * 1e-15)

# Print the final physics results
clock_names <- c("Hg", "Yb", "Sr")

cat("\n======================================================\n")
cat(" FINAL BAYESIAN ADJUSTED FREQUENCIES (Hz)\n")
cat("======================================================\n")

for(i in 1:3) {
  # Format the large numbers safely
  freq_str <- format(as.numeric(optimized_freqs_bayes[i]), digits=16)
  unc_str <- format(as.numeric(optimized_unc_bayes[i]), digits=3)
  
  cat(sprintf("%s Clock:\n", clock_names[i]))
  cat(sprintf("  New Freq: %s +/- %s Hz\n", freq_str, unc_str))
  cat(sprintf("  Fractional Offset (x): %.3f  (95%% CI: [%.3f, %.3f])\n\n", 
              X_hat_bayes[i], x_ci_lower[i], x_ci_upper[i]))
}



### look at results

# Load necessary plotting and data wrangling libraries
library(ggplot2)
library(tidyr)
library(dplyr)

# 1. Extract the dark uncertainty posterior samples into a data frame
# (Assuming your stanfit object is named 'fit')
posterior_samples <- rstan::extract(fit)
sigma_dark_samples <- as.data.frame(posterior_samples$sigma_dark)

# Rename columns to match the clock types
colnames(sigma_dark_samples) <- c("Hg", "Yb", "Sr")

# 2. Pivot the data to a "long" format, which ggplot loves
df_long <- pivot_longer(
  sigma_dark_samples, 
  cols = everything(), 
  names_to = "Clock", 
  values_to = "sigma_dark"
)

mean_df <- df_long %>%
  group_by(Clock) %>%
  summarize(mean_val = mean(sigma_dark))

# 3. Define the theoretical Prior function (Half-Normal, mean=0, sd=1)
prior_density <- function(x) {
  # Multiply by 2 because it's truncated at 0
  ifelse(x >= 0, 2 * dnorm(x, mean = 0, sd = 1), 0)
}

# 4. Generate the Overlay Plot
ggplot(df_long, aes(x = sigma_dark)) +
  
  # A. Plot the theoretical Prior as a thick dashed line
  stat_function(
    fun = prior_density, 
    aes(linetype = "Prior HN(0,1)"), 
    color = "black", 
    linewidth = 1.2
  ) +
  geom_vline(
    data = mean_df,
    aes(xintercept = mean_val, color = Clock),
    # linetype = "twodash", # distinct from the prior's dashed line
    linewidth = .5
  ) +
  # B. Plot the Posteriors as semi-transparent filled densities
  geom_density(
    aes(fill = Clock, color = Clock), 
    alpha = 0.5, 
    linewidth = 0.8
  ) +
  
  # C. Formatting and Aesthetics
  scale_fill_manual(values = c("Hg" = "#F8766D", "Yb" = "#00BA38", "Sr" = "#619CFF")) +
  scale_color_manual(values = c("Hg" = "#F8766D", "Yb" = "#00BA38", "Sr" = "#619CFF")) +
  scale_linetype_manual(name = "", values = "dashed") +
  
  # Zoom the x-axis in to where the action is (0 to 2)
  coord_cartesian(xlim = c(0, 2)) +
  
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.title = element_text(face = "bold"),
    legend.position = "right",
    legend.title = element_text(face = "bold")
  ) +
  labs(
    title = "Bayesian Dark Uncertainty Estimates: Posterior vs Prior",
    x = expression("Dark Uncertainty " * sigma[dark] * " (parts per " * 10^{15} * ")"),
    y = "Probability Density",
    fill = "Clock Posterior",
    color = "Clock Posterior"
  )
