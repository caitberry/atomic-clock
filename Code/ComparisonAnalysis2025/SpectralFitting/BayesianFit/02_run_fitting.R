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
  
  # Prepare Stan Data List
  stan_data <- list(
    N      = length(fit_data$freq),
    freq   = fit_data$freq,
    y_obs  = fit_data$spectrum,
    rel_sd = fit_data$sdish / fit_data$spectrum,
    tp     = TP_VAL
  )
  
  # Run Sampling
  fit <- sampling(
    stan_model_obj,
    data   = stan_data,
    iter   = STAN_ITER,
    warmup = STAN_WARMUP,
    chains = STAN_CORES
  )
  
  # Save Fit Object
  save_name <- paste0(OUTPUT_FIT_DIR, "fit_",dataName,"_Date_", gsub("-", "_", d), ".rds")
  saveRDS(fit, file = save_name)
  message(paste("Saved:", save_name))
}