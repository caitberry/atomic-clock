# source("/home/aak3/NIST/atomic-clock/Code/ComparisonAnalysis2025/SpectralFitting/BayesianFit/00_common.R")

# ==============================================================================
# 00_common.R
# Shared configuration and utility functions
# ==============================================================================
rm(list=ls())

library(rstan)
options(mc.cores = parallel::detectCores())

library(stringr)
library(dplyr)
library(ggplot2)
library(broom)
library(tidyr)
library(gridExtra)

# --- CONFIGURATION ---
# Paths
BASE_FILEPATH  <- "/home/aak3/NIST/atomic-clock/"

DATA_RAW_DIR   <- paste0(BASE_FILEPATH,"Results/ClockComp2025")
OUTPUT_FIT_DIR <- paste0(BASE_FILEPATH,"Code/ComparisonAnalysis2025/SpectralFitting/BayesianFit/Posteriors/")
OUTPUT_PLOT_DIR <- paste0(BASE_FILEPATH,"Code/ComparisonAnalysis2025/SpectralFitting/BayesianFit/Plots/")
DATA_PREP_FILE <- paste0(BASE_FILEPATH,"Code/ComparisonAnalysis2025/SpectralFitting/BayesianFit/processedData/data_prepared_for_fitting.rds")
STAN_MODEL     <- paste0(BASE_FILEPATH,"Code/ComparisonAnalysis2025/SpectralFitting/BayesianFit/complex_spectrum_model.stan")

# Processing Constants
MY_K    <- 5    # Number of tapers
TP_VAL  <- 1.0  # Probe time (s)
F_MIN   <- 0.0001
F_MAX   <- 0.1

# MCMC Settings
STAN_CORES  <- 4
STAN_ITER   <- 4000
STAN_WARMUP <- 3000

# # Ensure output dir exists
# if(!dir.exists(OUTPUT_FIT_DIR)) dir.create(OUTPUT_FIT_DIR, recursive = TRUE)

# --- HELPER FUNCTIONS ---

# 1. R implementation of the Physics Model (for posterior plotting)
model_psd_r <- function(freq, h0, h_m1, Kp, Ki, tau, tp) {
  omega <- 2 * pi * freq
  n <- length(omega)
  y_model <- numeric(n)
  
  for (i in seq_len(n)) {
    w <- omega[i]
    angle <- -w * tp / 2
    c_val <- cos(angle); s_val <- sin(angle)
    
    # Numerator G(w)
    re_num <- Kp * c_val - (-Ki/w) * s_val
    im_num <- Kp * s_val + (-Ki/w) * c_val
    
    # Denominator G(w)
    re_den <- 1; im_den <- w * tau
    den_sq <- re_den^2 + im_den^2
    
    # G = num/den
    re_G <- (re_num * re_den + im_num * im_den) / den_sq
    im_G <- (im_num * re_den - re_num * im_den) / den_sq
    
    # 1+G
    re_1pG <- 1 + re_G; im_1pG <- im_G
    one_plus_G_sq <- re_1pG^2 + im_1pG^2
    
    # Ratios
    abs_inv_sq <- 1 / one_plus_G_sq
    re_rat <- (re_G * re_1pG + im_G * im_1pG) / one_plus_G_sq
    im_rat <- (im_G * re_1pG - re_G * im_1pG) / one_plus_G_sq
    abs_ratio_sq <- re_rat^2 + im_rat^2
    
    # Result (Using cyclic freq logic: h_m1/f)
    y_model[i] <- abs_ratio_sq * h0 + abs_inv_sq * (h_m1 / freq[i])
  }
  return(y_model)
}

# 2. Extract Data Helper
extract_from_list <- function(x, filename) {
  df1 <- data.frame(freq = x[[1]], spectrum = x[[2]], File = filename)
  df2 <- data.frame(eigenval = x[[4]], K=1:MY_K, File = filename)
  list(spectrum = df1, eigenvals = df2)
}