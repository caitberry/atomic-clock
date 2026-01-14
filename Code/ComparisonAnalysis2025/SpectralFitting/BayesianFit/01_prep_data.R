# ==============================================================================
# 01_prep_data.R
# Organize data (Real or Simulated) and save to disk
# ==============================================================================
# source("Code/ComparisonAnalysis2025/SpectralFitting/BayesianFit/00_common.R")

# --- OPTION A: LOAD REAL DATA ---

dataName="AlYbApr15"
files_list <- list.files(path = DATA_RAW_DIR, pattern = "spectral.*AlYb.*Apr15", full.names = TRUE)
message(paste("Loading", length(files_list), "files..."))

all_data <- lapply(files_list, function(file) {
  x <- readRDS(file)
  extract_from_list(x, basename(file))
})

spectralEstDF <- do.call(rbind, lapply(all_data, `[[`, "spectrum"))

# Add Metadata & Clean
spectralEstDF <- spectralEstDF %>%
  mutate(
    Ratio = str_extract(File, "(?<=spectralEstFor)[A-Za-z]+"),
    Month = str_extract(File, "(?<=For[A-Za-z]{4})(\\d{2})"),
    Day   = str_extract(File, "(?<=For[A-Za-z]{4}\\d{2})(\\d{2})"),
    Date  = as.Date(paste0("2025-", Month, "-", Day)),
    # Calculate approx SD for weighting
    sdish = sqrt((1 * spectrum^2)/MY_K)
  )

# --- OPTION B: GENERATE SIMULATED DATA ---
# 
# dataName="ComplexModelSimDat"
# f_sim <- exp(seq(log(F_MIN), log(F_MAX), length.out=100))
# true_params <- list(h0=1e-31, h_m1=1.5e-33, Kp=10, Ki=5, tau=14)
# psd_true <- model_psd_r(f_sim, true_params$h0, true_params$h_m1, 
#                         true_params$Kp, true_params$Ki, true_params$tau, TP_VAL)
# # Add multiplicative noise (Chi-squared like)
# psd_noisy <- psd_true * rchisq(length(f_sim), df=2*MY_K) / (2*MY_K)
# 
# spectralEstDF <- data.frame(
#   freq = f_sim, spectrum = psd_noisy, 
#   sdish = psd_noisy/sqrt(MY_K), 
#   Date = as.Date("2025-01-01"), File = "Simulated"
# )

# --- SAVE ---
saveRDS(spectralEstDF, file = DATA_PREP_FILE)
message(paste("Data prepared and saved to:", DATA_PREP_FILE))
