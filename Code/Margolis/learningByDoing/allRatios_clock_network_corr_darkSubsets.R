library(Rmpfr)
library(rstan)
library(ggplot2)

# --- 1. Settings & Data Loading ---
prec <- 128
scale_factor <- mpfr(1e15, precBits = prec)
raw_data <- readLines("Data/Margolis2024DataMetrologia/ClockInputData2020_7.dat")

# --- 2. Reference Frequencies (nu1 to nu14) ---
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
q_lines <- raw_data[grep("^q\\d+", raw_data)]

parsed_list <- lapply(q_lines, function(line) {
  parts <- strsplit(trimws(line), "\\s+")[[1]]
  return(parts[1:4])
})

parsed_q <- do.call(rbind, parsed_list)
q_labels <- parsed_q[, 1]  
q_types  <- parsed_q[, 2]  
Q_mpfr   <- mpfr(parsed_q[, 3], precBits = prec) 
U_mpfr   <- mpfr(parsed_q[, 4], precBits = prec) 

# --- 4. Logic to Determine Reference Values and Clock IDs ---
N_meas <- length(Q_mpfr)
N_clocks <- length(nu_0)

ref_vals <- mpfr(rep(0, N_meas), precBits = prec)
is_ratio <- rep(0, N_meas)
clock_1  <- rep(0, N_meas)
clock_2  <- rep(0, N_meas)

# Track measurement counts for the Dark Uncertainty grouping
clock_counts <- rep(0, N_clocks)

for (i in 1:N_meas) {
  type_str <- q_types[i]
  
  if (grepl("_over_", type_str)) {
    is_ratio[i] <- 1
    parts <- strsplit(type_str, "_over_")[[1]]
    idx1 <- as.numeric(gsub("nu", "", parts[1]))
    idx2 <- as.numeric(gsub("nu", "", parts[2]))
    
    clock_1[i]  <- idx1
    clock_2[i]  <- idx2
    ref_vals[i] <- nu_0[idx1] / nu_0[idx2]
    
    # Increment counters for both clocks involved in ratio
    clock_counts[idx1] <- clock_counts[idx1] + 1
    clock_counts[idx2] <- clock_counts[idx2] + 1
  } else {
    is_ratio[i] <- 0
    idx1 <- as.numeric(gsub("nu", "", type_str))
    
    clock_1[i]  <- idx1
    clock_2[i]  <- 0 
    ref_vals[i] <- nu_0[idx1]
    
    # Increment counter
    clock_counts[idx1] <- clock_counts[idx1] + 1
  }
}

# --- NEW: Map clocks to their dark uncertainty parameter based on threshold ---
minmeas=7

clock_to_dark_idx <- rep(0, N_clocks)
has_global <- any(clock_counts < minmeas)
next_idx <- 1

if (has_global) {
  global_idx <- 1
  next_idx <- 2
  for (c in 1:N_clocks) {
    if (clock_counts[c] < minmeas) {
      clock_to_dark_idx[c] <- global_idx
    }
  }
}

for (c in 1:N_clocks) {
  if (clock_counts[c] >= minmeas) {
    clock_to_dark_idx[c] <- next_idx
    next_idx <- next_idx + 1
  }
}
N_dark <- next_idx - 1

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
  Measurement_Type = ifelse(is_ratio == 1, "Frequency Ratios", "Absolute Frequencies")
)

q_numbers <- as.numeric(gsub("q", "", measurement_df$Q_ID))
measurement_df$Q_ID <- factor(measurement_df$Q_ID, levels = measurement_df$Q_ID[order(q_numbers)])

cutoff <- 5.0
measurement_df$Unc_Tier <- ifelse(measurement_df$U_Unc_Scaled > cutoff, "High Uncertainty", "Lower Uncertainty")
measurement_df$Unc_Tier <- factor(measurement_df$Unc_Tier, levels = c("Lower Uncertainty", "High Uncertainty"))

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
  facet_wrap(~ Measurement_Type+Unc_Tier, scales = "free", ncol = 1)

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
cat("Estimating", N_dark, "unique dark uncertainty parameters.\n")

# --- 7. Run Stan ---
stan_data <- list(
  N_clocks = N_clocks,
  N_meas = N_meas,
  Y = Y_stan,
  V = V_stan,
  is_ratio = is_ratio,
  clock_1 = clock_1,
  clock_2 = clock_2,
  N_dark = N_dark,                      # NEW
  clock_to_dark_idx = clock_to_dark_idx # NEW
)

# Compile and Run
fit <- stan(
  file = "Code/Margolis/learningByDoing/clock_network_corr_darkSubsets.stan", 
  data = stan_data, 
  iter = 20000, 
  warmup = 15000, 
  chains = 4, 
  cores = 4,
  seed = 42
)

# --- Extract Posterior Samples ---
posterior_samples <- rstan::extract(fit)

# Extract fractional offsets (x)
X_hat_bayes <- colMeans(posterior_samples$x)
cov_X_bayes <- cov(posterior_samples$x)

# --- NEW: Map dynamic dark uncertainties back to the 14 individual clocks ---
# Safely handle dimensions if N_dark == 1
sigma_dark_samples <- posterior_samples$sigma_dark
if (is.null(dim(sigma_dark_samples))) {
  sigma_dark_samples <- matrix(sigma_dark_samples, ncol = 1)
}

dark_unc_estimates <- rep(0, N_clocks)
dark_unc_ci_lower  <- rep(0, N_clocks)
dark_unc_ci_upper  <- rep(0, N_clocks)

for (c in 1:N_clocks) {
  idx <- clock_to_dark_idx[c]
  dark_unc_estimates[c] <- mean(sigma_dark_samples[, idx])
  dark_unc_ci_lower[c]  <- quantile(sigma_dark_samples[, idx], probs = 0.025)
  dark_unc_ci_upper[c]  <- quantile(sigma_dark_samples[, idx], probs = 0.975)
}

# High-Precision Frequency Reconstruction
optimized_freqs_bayes <- nu_0 * (1 + mpfr(X_hat_bayes, prec) * 1e-15)
optimized_unc_bayes   <- nu_0 * (mpfr(sqrt(diag(cov_X_bayes)), prec) * 1e-15)

# Display Results
results_df <- data.frame(
  Clock   = paste0("nu", 1:14),
  Meas_Count = clock_counts,         # NEW: See exactly how many meas. a clock had
  Dark_Group = clock_to_dark_idx,    # NEW: Index 1 is the shared "Global" group
  Freq_Hz = format(optimized_freqs_bayes, scientific = FALSE),
  Unc_Hz  = format(optimized_unc_bayes, scientific = TRUE),
  Dark_Unc_1e15 = sprintf("%.3f", dark_unc_estimates),
  Dark_Unc_95CI = sprintf("[%.3f, %.3f]", dark_unc_ci_lower, dark_unc_ci_upper)
)

cat("\n### Bayesian Optimized Frequencies & Grouped Dark Uncertainties ###\n")
print(results_df)
