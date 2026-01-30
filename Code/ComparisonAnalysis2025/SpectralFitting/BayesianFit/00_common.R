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
TP_VAL  <- 1.0  # Probe time (s)  ### CHANGES FOR DIFFERENT CLOCKS
F_MIN   <- 0.0001
F_MAX   <- 0.1

# MCMC Settings
STAN_CORES  <- 4
STAN_ITER   <- 10000
STAN_WARMUP <- 8000

# # Ensure output dir exists
# if(!dir.exists(OUTPUT_FIT_DIR)) dir.create(OUTPUT_FIT_DIR, recursive = TRUE)

# --- HELPER FUNCTIONS ---

# 1. R implementation of the Physics Model (for posterior plotting)
model_psd_r <- function(f, h0, h_m1, Kp, Ki, tau, tp) {
  # Angular frequency
  w <- 2 * pi * f

  # 1. Calculate Open Loop Gain G(w)
  #    Using 1i for imaginary unit
  s <- 1i * w
  G <- ( (Kp + Ki/s) / (1 + s*tau) ) * exp(-s * tp / 2)
  
  # 2. Calculate Closed Loop Shapes
  #    Mod(z)^2 gives |z|^2
  noise_transfer_h0   <- Mod(G / (1 + G))^2
  noise_transfer_hm1  <- Mod(1 / (1 + G))^2
  
  # 3. Combine
  Sy <- noise_transfer_h0 * h0 + noise_transfer_hm1 * (h_m1 / f)
  
  return(Sy)
}

# 2. Extract Data Helper
extract_from_list <- function(x, filename) {
  df1 <- data.frame(freq = x[[1]], spectrum = x[[2]], File = filename)
  df2 <- data.frame(eigenval = x[[4]], K=1:MY_K, File = filename)
  list(spectrum = df1, eigenvals = df2)
}

calculate_log_bias <- function(K) {
  # digamma(K) - log(K) is the theoretical bias for log(chi-sq/2K)
  return(digamma(K) - log(K))
}