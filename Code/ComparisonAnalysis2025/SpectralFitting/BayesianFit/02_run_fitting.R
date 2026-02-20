# ==============================================================================
# 02_run_fitting.R
# Load prepared data and run Bayesian Fits
# ==============================================================================


# Load Data
if(!file.exists(DATA_PREP_FILE)) stop("Run 01_prep_data.R first!")
spectralEstDF <- readRDS(DATA_PREP_FILE)

# Compile Model (only needs to happen once)
stan_model_obj <- stan_model(file = STAN_MODEL)

for(i in 1:length(unique(spectralEstDF$Date))){
  d=unique(spectralEstDF$Date)[i]
  
  message(paste("--- Fitting Date:", d, "---"))
  
  # Filter Data
  fit_data <- spectralEstDF %>% 
    filter(Date == d, freq > F_MIN, freq < F_MAX)
  
  # Calculate the bias
  log_bias_val <- calculate_log_bias(MY_K)
  # New (Correct Variance for Log-Chi-Sq)
  # The SD of log(PSD) is constant: sqrt(trigamma(K))
  log_sd_const <- sqrt(trigamma(MY_K))
  
  stan_data <- list(
    N      = length(fit_data$freq),
    freq   = fit_data$freq,
    y_obs  = fit_data$spectrum,
    sigma_log = rep(log_sd_const,length(fit_data$freq)),#fit_data$sdish / fit_data$spectrum,
    tp     = TP_VAL,
    bias   = log_bias_val
  )
  
  # # Run Sampling
  fit <- sampling(
    stan_model_obj,
    data   = stan_data,
    iter   = STAN_ITER,
    warmup = STAN_WARMUP,
    chains = STAN_CORES,
    # init = init_fun
    control = list(
      adapt_delta   = 0.98  
      # max_treedepth = 15    
    )
  )
  
  # Save Fit Object
  save_name <- paste0(OUTPUT_FIT_DIR, "fit_",dataName,"_Date_", gsub("-", "_", d), ".rds")
  saveRDS(fit, file = save_name)
  message(paste("Saved:", save_name))
}