# ==============================================================================
# 05_simple_fit_analysis.R
# Load prepared data and run simple linear fit
# ==============================================================================

# Load Data
if(!file.exists(DATA_PREP_FILE)) stop("Run 01_prep_data.R first!")
spectralEstDF <- readRDS(DATA_PREP_FILE)

for(i in 1:length(unique(spectralEstDF$Date))){
  d=unique(spectralEstDF$Date)[i]
  
  message(paste("--- Fitting Date:", d, "---"))
  
  # Filter Data
  fit_data <- spectralEstDF %>% 
    filter(Date == d, freq > F_MIN, freq < F_MAX)
  
  # Calculate the bias
  log_bias_val <- calculate_log_bias(MY_K)
  
  F_MIN_BAND <- 0.0002
  F_MAX_BAND <- 0.008

  band_data <- fit_data %>% filter(freq >= F_MIN_BAND & freq <= F_MAX_BAND)
  fit_freq  <- lm(log(spectrum) ~ 1, data = band_data)

  # Extract the Raw (Biased) Intercept
  beta0_biased <- coef(fit_freq)[1]
  ci_biased    <- confint(fit_freq, level = 0.95)

  # Correct for the Multitaper Bias (Shift Up)
  beta0_corrected <- beta0_biased - log_bias_val
  ci_corrected    <- ci_biased - log_bias_val

  # Convert to Linear Scale
  h0_est <- exp(beta0_corrected)
  h0_ci  <- exp(ci_corrected)

  cat(sprintf("\n--- PRIMARY FIT (0 - %.4f Hz) ---\n", F_MAX_BAND))
  cat(sprintf("h0 Estimate: %.2e\n95%% CI:      [%.2e, %.2e]\n", h0_est, h0_ci[1], h0_ci[2]))

  # SENSITIVITY ANALYSIS (Stability Check) 
  # Scan upper cutoffs from 0.001 to 0.05 Hz
  cutoffs <- seq(0.001, 0.05, by = 0.0005)
  
  sensitivity_df <- do.call(rbind, lapply(cutoffs, function(fc) {
    sub_d <- fit_data %>% filter(freq >= F_MIN_BAND & freq <= fc)
    if(nrow(sub_d) < 5) return(NULL)
  
    m     <- lm(log(spectrum) ~ 1, data = sub_d)
    # Extract the Raw (Biased) Intercept
    beta0_biased <- coef(m)[1]
    ci_biased    <- confint(m, level = 0.95)
  
    # Correct for the Multitaper Bias (Shift Up)
    beta0_corrected <- beta0_biased - log_bias_val
    ci_corrected    <- ci_biased - log_bias_val
    
    # Convert to Linear Scale
    h0_est <- exp(beta0_corrected)
    h0_ci  <- exp(ci_corrected)
  
    data.frame(
      cutoff = fc,
      h0     = exp(beta0_corrected),
      lower  = exp(ci_corrected)[1],
      upper  = exp(ci_corrected)[2]
    )
}))

# # --- PLOTS ---
# # Plot A: The Primary Fit
# p1 <- ggplot() +
#   geom_point(data=fit_data, aes(x=freq, y=spectrum), color="grey80", alpha=0.5) +
#   geom_point(data=band_data, aes(x=freq, y=spectrum), color="black", alpha=0.8) +
#   geom_segment(aes(x=min(fit_data$freq), xend=F_MAX_BAND, 
#                    y=exp(beta0), yend=exp(beta0)), color="red", size=1.5) +
#   scale_x_log10() + scale_y_log10() +
#   labs(title = "Simple Fit", x="Freq (Hz)", y="PSD") + theme_bw()
# 
# # Plot B: Sensitivity / Stability
# p2 <- ggplot(sensitivity_df, aes(x=cutoff, y=h0)) +
#   geom_ribbon(aes(ymin=lower, ymax=upper), fill="blue", alpha=0.2) +
#   geom_line(color="blue") + geom_point(size=0.5) +
#   geom_vline(xintercept=F_MAX_BAND, linetype="dashed", color="red") +
#   scale_y_log10() +
#   labs(title="Sensitivity Analysis", subtitle="Effect of Upper Cutoff on h0",
#        x="Upper Cutoff (Hz)", y="Estimated h0") + theme_bw()
# 
