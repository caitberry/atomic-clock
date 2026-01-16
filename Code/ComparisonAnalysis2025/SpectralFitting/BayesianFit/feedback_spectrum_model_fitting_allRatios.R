# --- SCRIPT 1: GENERALIZED STAN FITTING (Updated) ---

library(rstan)
options(mc.cores = 4)
library(dplyr)
library(stringr)

# --- 1. SETUP ---
folder <- "/home/aak3/NIST/atomic-clock/Results/ClockComp2025"  
output_folder <- "/home/aak3/NIST/atomic-clock/Code/ComparisonAnalysis2025/BayesianFit/Posteriors"
stan_model_file <- "Code/ComparisonAnalysis2025/BayesianFit/feedback_spectrum_model.stan"

if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}

dataset_names <- c("YbSr","AlYb", "AlSr")
myK <- 5

# --- 2. MAIN LOOP TO PROCESS EACH DATASET ---
for (dataset_name in dataset_names) {
  
  cat("\n====================================================\n")
  cat("=== Starting Stan Fitting for:", dataset_name, "===\n")
  cat("====================================================\n")
  
  # --- a. Load data for the current dataset ---
  spectral_files <- list.files(
    path = folder,
    pattern = paste0("spectral.*", dataset_name, ".*Apr15"),
    full.names = TRUE
  )
  
  if (length(spectral_files) == 0) {
    warning(paste("No spectral files found for", dataset_name, "- skipping."))
    next
  }
  
  extract_from_list <- function(x, filename) {
    list(spectrum = data.frame(freq = x[[1]], spectrum = x[[2]], File = filename))
  }
  
  all_data <- lapply(spectral_files, function(file) {
    x <- readRDS(file)
    extract_from_list(x, basename(file))
  })
  
  spectralEstDF <- do.call(rbind, lapply(all_data, `[[`, "spectrum"))
  
  spectralEstDF <- spectralEstDF %>%
    mutate(
      Ratio = str_extract(File, "(?<=spectralEstFor)[A-Za-z]+"),
      Month = str_extract(File, "(?<=For[A-Za-z]{4})(\\d{2})"),
      Day   = str_extract(File, "(?<=For[A-Za-z]{4}\\d{2})(\\d{2})"),
      Date = as.Date(paste0("2025-", Month, "-", Day))
    ) %>%
    mutate(
      sdish = spectrum / sqrt(myK)
    )
  
  # --- b. Loop through each day and run Stan ---
  unique_dates <- unique(spectralEstDF$Date)
  
  for(i in 1:length(unique_dates)){
    d=unique_dates[i]
    cat("  -> Processing Date:", d, "\n")
    
    plot_data <- filter(spectralEstDF, Date == d)
    fit_data <- filter(plot_data, freq > 0.0001, freq < .1)
    
    if (nrow(fit_data) == 0) {
      cat("     No data in frequency range for this date. Skipping.\n")
      next
    }
    
    # --- MODIFICATION: Determine tp_value and add it to the data list ---
    if (grepl("Al", dataset_name)) {
      tp_value <- 1
    } else {
      tp_value <- 1
    }
    
    stan_data <- list(
      N = length(fit_data$freq),
      omega = 2 * pi * fit_data$freq,
      y_obs = fit_data$spectrum,
      rel_sd = fit_data$sdish / fit_data$spectrum,
      tp = tp_value  # Pass the correct tp value to Stan
    )
    
    fit <- stan(
      file = stan_model_file,
      data = stan_data,
      iter = 4000,
      warmup = 3000,
      chains = 4  
    )
    
    # --- c. Save the output ---
    output_filename <- file.path(
      output_folder,
      paste0("fit_", dataset_name, "_Date", gsub("-", "_", d), ".rds")
    )
    
    saveRDS(fit, file = output_filename)
    cat("     Fit saved to:", basename(output_filename), "\n")
  }
}

cat("\nAll datasets processed.\n")