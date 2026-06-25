library(Rmpfr)
library(rstan)
library(ggplot2)

# --- 1. Settings & Data Loading ---
prec <- 128
scale_factor <- mpfr(1e15, precBits = prec)
raw_data <- readLines("Data/Margolis2024DataMetrologia/ClockInputData2020_7.dat")

# --- 2. Reference Frequencies (nu1 to nu14) ---

# Define Nu_0 using STRINGS to prevent early truncation
# These reflect the 2017 CIPM recommended standard frequencies
nu_0 <- mpfr(c(
  "1267402452901050.0", # 1: 115In+ 
  "1233030706593509.0", # 2: H (1S-2S two-photon transition / 2)
  "1128575290808154.4", # 3: 199Hg
  "1121015393207851.0", # 4: 27Al+
  "1064721609899145.0", # 5: 199Hg+
  "688358979309308.3",  # 6: 171Yb+ (E2 Quadrupole)
  "642121496772645.0",  # 7: 171Yb+ (E3 Octupole)
  "518295836590863.6",  # 8: 171Yb
  "455986240494140.0",  # 9: 40Ca
  "444779044095486.5",  # 10: 88Sr+
  "429228066418007.0",  # 11: 88Sr
  "429228004229873.0",  # 12: 87Sr
  "411042129776399.8",  # 13: 40Ca+
  "6834682610.9043126"  # 14: 87Rb (Microwave transition)
), precBits = prec)


# --- 3. Parse All 106 Measurements ---
# Filters lines starting with 'q' followed by a number
q_lines <- raw_data[grep("^q\\d+", raw_data)]

# Use a list-based approach to handle inconsistent whitespace/tabs
parsed_list <- lapply(q_lines, function(line) {
  # Split by any whitespace (tabs or spaces)
  parts <- strsplit(trimws(line), "\\s+")[[1]]
  
  # Return only the first 4 columns: Label, Type, Value, Uncertainty
  # This ignores the [Reference/Comments] column which varies in length
  return(parts[1:4])
})

# Combine into a clean matrix
parsed_q <- do.call(rbind, parsed_list)

q_labels <- parsed_q[, 1]  # e.g., "q1", "q2" 
q_types  <- parsed_q[, 2]  # e.g., "nu1" or "nu3_over_nu12" 
Q_mpfr   <- mpfr(parsed_q[, 3], precBits = prec) # Measurement values 
U_mpfr   <- mpfr(parsed_q[, 4], precBits = prec) # Uncertainties

# --- 4. Logic to Determine Reference Values and Clock IDs ---
N_meas <- length(Q_mpfr)
N_clocks <- length(nu_0)

ref_vals <- mpfr(rep(0, N_meas), precBits = prec)
is_ratio <- rep(0, N_meas)
clock_1  <- rep(0, N_meas)
clock_2  <- rep(0, N_meas)

for (i in 1:N_meas) {
  type_str <- q_types[i]
  
  if (grepl("_over_", type_str)) {
    # It's a ratio: e.g., "nu3_over_nu12"
    is_ratio[i] <- 1
    parts <- strsplit(type_str, "_over_")[[1]]
    idx1 <- as.numeric(gsub("nu", "", parts[1]))
    idx2 <- as.numeric(gsub("nu", "", parts[2]))
    
    clock_1[i]  <- idx1
    clock_2[i]  <- idx2
    ref_vals[i] <- nu_0[idx1] / nu_0[idx2]
  } else {
    # It's an absolute: e.g., "nu1"
    is_ratio[i] <- 0
    idx1 <- as.numeric(gsub("nu", "", type_str))
    
    clock_1[i]  <- idx1
    clock_2[i]  <- 0 # Stan model uses 0 or N_clocks for nulls
    ref_vals[i] <- nu_0[idx1]
  }
}

# --- 5. Scale Data for Stan ---
Y_stan <- as.numeric(((Q_mpfr - ref_vals) / ref_vals) * scale_factor)
u_scaled_stan <- as.numeric((U_mpfr / ref_vals) * scale_factor)
V_stan <- diag(u_scaled_stan^2)


# --- Sanity Check: Create DataFrame for Plotting ---
measurement_df <- data.frame(
  Q_ID = q_labels,
  Y_Offset_Scaled = Y_stan,
  U_Unc_Scaled = u_scaled_stan,
  Nu_Type = q_types,
  # Create a friendly label for the facet_wrap
  Measurement_Type = ifelse(is_ratio == 1, "Frequency Ratios", "Absolute Frequencies")
)

# Force Q_ID to order numerically instead of alphabetically
q_numbers <- as.numeric(gsub("q", "", measurement_df$Q_ID))
measurement_df$Q_ID <- factor(measurement_df$Q_ID, levels = measurement_df$Q_ID[order(q_numbers)])

# Create the Uncertainty Tier column
cutoff <- 5.0

measurement_df$Unc_Tier <- ifelse(
  measurement_df$U_Unc_Scaled > cutoff, 
  "High Uncertainty", 
  "Lower Uncertainty"
)

# 2. Convert to a factor to control the order the panels appear in
measurement_df$Unc_Tier <- factor(
  measurement_df$Unc_Tier, 
  levels = c("Lower Uncertainty", "High Uncertainty")
)


# --- Plot the Data ---
pd <- position_jitter(width = 0.15, height = 0, seed = 42)

sanity_plot <- ggplot(measurement_df, aes(x = Nu_Type, y = Y_Offset_Scaled, color = Nu_Type)) +
  geom_point(size = 2,position = pd) +
  geom_errorbar(aes(ymin = Y_Offset_Scaled - U_Unc_Scaled, 
                    ymax = Y_Offset_Scaled + U_Unc_Scaled), 
                width = 0.1,position = pd) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
    legend.position = "right",
    strip.text = element_text(face = "bold", size = 11)
  ) +
  labs(
    title = "Sanity Check: Input Data Offsets",
    x = "Measurement ID",
    y = expression("Fractional Offset " ~ (Y %*% 10^{15})),
    color = "Transition"
  ) +
  # scales = "free" unlocks BOTH the X and Y axes for each panel
  facet_wrap(~ Measurement_Type+Unc_Tier, scales = "free", ncol = 1)

# Display the plot
print(sanity_plot)


# --- 6. Populate Correlations ---
corr_lines <- raw_data[grep("^r\\(", raw_data)]
for (line in corr_lines) {
  parts <- strsplit(trimws(line), "\\s+")[[1]]
  pair_str <- parts[1]
  r_val <- as.numeric(parts[2])
  
  qs <- gsub("r\\(|\\)", "", pair_str)
  q_pair <- strsplit(qs, ",")[[1]]
  
  idx1 <- match(q_pair[1], q_labels)
  idx2 <- match(q_pair[2], q_labels)
  
  if (!is.na(idx1) && !is.na(idx2)) {
    cov_val <- r_val * u_scaled_stan[idx1] * u_scaled_stan[idx2]
    V_stan[idx1, idx2] <- cov_val
    V_stan[idx2, idx1] <- cov_val
  }
}

cat("Successfully parsed", N_meas, "measurements for", N_clocks, "clocks.\n")

# --- 7. Run Stan ---
stan_data <- list(
  N_clocks = N_clocks,
  N_meas = N_meas,
  Y = Y_stan,
  V = V_stan,
  is_ratio = is_ratio,
  clock_1 = clock_1,
  clock_2 = clock_2
)

# --- 1. Compile and Run ---
# Using the updated stan_data which now contains 106 measurements and 14 clocks
fit <- stan(
  file = "Code/Margolis/learningByDoing/clock_network_corr_oneDarkUnc.stan", 
  data = stan_data, 
  iter = 20000, 
  warmup = 15000, 
  chains = 4, 
  cores = 4,
  # control = list(adapt_delta = 0.99, max_treedepth = 15),
  seed = 42
)

# --- 2. Extract Posterior Samples ---
posterior_samples <- rstan::extract(fit)

# Extract fractional offsets (x) for all 14 clocks 
X_hat_bayes <- colMeans(posterior_samples$x)
cov_X_bayes <- cov(posterior_samples$x)

# --- 3. High-Precision Frequency Reconstruction ---
# nu_0 now contains all 14 reference frequencies parsed from 's' lines 
optimized_freqs_bayes <- nu_0 * (1 + mpfr(X_hat_bayes, prec) * 1e-15)
optimized_unc_bayes   <- nu_0 * (mpfr(sqrt(diag(cov_X_bayes)), prec) * 1e-15)

# --- 4. Display Results for All 14 Clocks ---
# Convert the mpfr objects to strings using format() or formatBin() 
# so they can safely live inside a standard data.frame
results_df <- data.frame(
  Clock   = paste0("nu", 1:14),
  # Force non-scientific notation for the main frequency so you see all digits
  Freq_Hz = format(optimized_freqs_bayes, scientific = FALSE),
  # Scientific notation is usually fine for the uncertainties
  Unc_Hz  = format(optimized_unc_bayes, scientific = TRUE) 
)

# Print the resulting dataframe
print(results_df)

# --- 5. Extract Dark Uncertainty ---
# This is the global parameter sigma_dark from the Stan model 
dark_unc_estimate <- mean(posterior_samples$sigma_dark)
dark_unc_ci <- quantile(posterior_samples$sigma_dark, probs = c(0.025, 0.975))

cat("\n--- Dark Uncertainty Analysis ---\n")
cat("Estimated Global Dark Uncertainty (parts per 10^15):", 
    sprintf("%.3f", dark_unc_estimate), "\n")
cat("95% Credible Interval:", 
    sprintf("[%.3f, %.3f]", dark_unc_ci[1], dark_unc_ci[2]), "\n")



