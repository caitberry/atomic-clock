# ==============================================================================
# 06_compare_different_method_results.R
# Visualize Bayesian and Simple Fit results together
# ==============================================================================
# source("00_common.R")

# Load the original data for comparison
spectralEstDF <- readRDS(DATA_PREP_FILE)

for(i in 1:length(unique(spectralEstDF$Date))){
  target_date=unique(spectralEstDF$Date)[i]
  
  bayesian_fit_file <- paste0(OUTPUT_FIT_DIR, "fit_", dataName, "_Date_", gsub("-", "_", target_date), ".rds")
  simple_fit_file <- paste0(Output_SimpleFit_Dir, "fit_", dataName, "_Date_", gsub("-", "_", target_date), ".rds")

  bayesian_fit <- readRDS(bayesian_fit_file)
  simple_fit <- readRDS(simple_fit_file)
  
  # --- EXTRACT BAYESIAN ESTIMATES ---
  bayes_draws <- rstan::extract(bayesian_fit)$log_h0
  bayes_h0 <- exp(bayes_draws)
  
  bayes_est <- mean(bayes_h0)
  bayes_ci_lower <- quantile(bayes_h0, probs = c(0.025))
  bayes_ci_upper <- quantile(bayes_h0, probs = c(0.975))
  
  # --- EXTRACT SIMPLE FIT ESTIMATES ---
  simple_est <- simple_fit$simplefit$h0
  simple_ci_lower <- simple_fit$simplefit$CIlower
  simple_ci_upper <- simple_fit$simplefit$CIupper
  
  # --- BUILD DATA FRAME ---
  resultDF <- data.frame(
    h0     = c(bayes_est, simple_est),
    lower  = c(bayes_ci_lower, simple_ci_lower),
    upper  = c(bayes_ci_upper, simple_ci_upper),
    method = factor(c("Bayesian", "Simple Fit"), levels = c("Simple Fit", "Bayesian")) # Reorder for logic
  )
  
  # --- MAKE THE PRETTY PLOT ---
  p <- ggplot(resultDF, aes(x = method, y = h0, color = method)) +
    
    # 1. The Error Bars (Width controls the horizontal caps, size controls thickness)
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.15, size = 1) +
    
    # 2. The Main Points (Size 5 makes them pop)
    geom_point(size = 5) +
    
    # 3. Manual Colors (Matches your previous spaghetti plots)
    scale_color_manual(values = c("Bayesian" = "blue", "Simple Fit" = "red")) +
    
    # 4. Log Scale with Scientific Notation labels
    scale_y_log10(labels = scales::scientific) +
    
    # 5. Clean Theme and Labels
    theme_bw(base_size = 14) +
    labs(
      y = expression(paste(h[0])),
      x = NULL # Remove x-axis label since the ticks explain it
    ) +
    
    # 6. Minor Tweaks
    theme(
      legend.position = "none",           # Legend is redundant with x-axis labels
      panel.grid.minor.x = element_blank() # Remove vertical clutter
    )

  # Check if 'true_params' object exists and has an 'h0' component
  if (exists("true_params") && !is.null(true_params$h0)) {
    p <- p +
      geom_hline(yintercept = true_params$h0,
                 linetype = "dashed",
                 color = "black",
                 size = 1,
                 alpha = 0.7) 
  }
  # Print to screen (optional)
  print(p)
  
  # Save the plot
  save_filename <- paste0(Output_SimpleFit_Dir, "compare_h0_", gsub("-", "_", target_date), ".png")
  ggsave(save_filename, plot = p, width = 6, height = 4)
  message(paste("Saved plot to:", save_filename))
  
}




##### PLOT TOGETHER

# ==============================================================================
# 06_compare_different_method_results_combined.R
# Visualize Bayesian vs Simple Fit results for ALL dates in one plot
# ==============================================================================
library(ggplot2)
library(scales)
library(dplyr)

# source("00_common.R") # Ensure this is loaded if needed for paths/constants

# Load the original data to get the list of dates
spectralEstDF <- readRDS(DATA_PREP_FILE)
# unique_dates  <- unique(spectralEstDF$Date)
unique_dates  <- c("2025-03-07","2025-03-18", "2025-03-20")

# Initialize an empty list to store results
results_list <- list()

# --- 1. LOOP TO COLLECT DATA ---
for(i in seq_along(unique_dates)){
  target_date <- unique_dates[i]

  # Define filenames
  bayesian_fit_file <- paste0(OUTPUT_FIT_DIR, "fit_", dataName, "_Date_", gsub("-", "_", target_date), ".rds")
  simple_fit_file   <- paste0(Output_SimpleFit_Dir, "fit_", dataName, "_Date_", gsub("-", "_", target_date), ".rds")

  # Skip if files are missing
  if(!file.exists(bayesian_fit_file) | !file.exists(simple_fit_file)) {
    warning(paste("Missing fit files for date:", target_date)); next
  }

  # Load fits
  bayesian_fit <- readRDS(bayesian_fit_file)
  simple_fit   <- readRDS(simple_fit_file)

  # --- EXTRACT BAYESIAN ESTIMATES ---
  # Note: Standardize extraction. Some Stan objects use 'beta0' or 'log_h0' depending on your model version.
  # Assuming 'log_h0' represents the true unbiased physical value:
  bayes_draws    <- rstan::extract(bayesian_fit)$log_h0
  bayes_h0_chain <- exp(bayes_draws)

  bayes_est      <- median(bayes_h0_chain)
  bayes_ci_lower <- quantile(bayes_h0_chain, probs = 0.025)
  bayes_ci_upper <- quantile(bayes_h0_chain, probs = 0.975)

  # --- EXTRACT SIMPLE FIT ESTIMATES ---
  # Assuming simple_fit structure is list(simplefit = data.frame(...), ...)
  simple_est      <- simple_fit$simplefit$h0
  simple_ci_lower <- simple_fit$simplefit$CIlower
  simple_ci_upper <- simple_fit$simplefit$CIupper

  # --- STORE IN LIST ---
  # We create a small data frame for this specific date
  results_list[[i]] <- data.frame(
    Date   = as.factor(target_date),
    Method = c("Bayesian", "Simple Fit"),
    h0     = c(bayes_est, simple_est),
    lower  = c(bayes_ci_lower, simple_ci_lower),
    upper  = c(bayes_ci_upper, simple_ci_upper)
  )
}

# Combine all days into one master data frame
all_results <- do.call(rbind, results_list)

# Reorder factor levels so "Simple Fit" comes before "Bayesian" (or vice versa)
all_results$Method <- factor(all_results$Method, levels = c("Simple Fit", "Bayesian"))

# --- 2. CREATE THE COMBINED PLOT ---
p <- ggplot(all_results, aes(x = Date, y = h0, color = Method, group = Method)) +

  # A. Error Bars
  # position_dodge(width) moves the bars apart so they don't overlap
  geom_errorbar(aes(ymin = lower, ymax = upper),
                width = 0.2, size = 1,
                position = position_dodge(width = 0.5)) +

  # B. Main Points
  geom_point(size = 4, position = position_dodge(width = 0.5)) +

  # C. Styling
  scale_color_manual(values = c("Bayesian" = "blue", "Simple Fit" = "red")) +
  scale_y_log10(labels = scales::scientific) +

  theme_bw(base_size = 14) +
  labs(
    y = expression(paste(h[0])),
    x = "Date"
  ) +

  theme(
    legend.position = "bottom",
    panel.grid.minor.x = element_blank()
  )

# Print and Save
print(p)

save_filename <- paste0(Output_SimpleFit_Dir, "compare_h0_SELECT_DATES",dataName,".png")
ggsave(save_filename, plot = p, width = 8, height = 4)
message(paste("Saved combined plot to:", save_filename))



### add in adev approach comparison
library(tidyverse)
library(readxl)
library(dplyr)
library(lubridate)

# # to convert h0 to ADEV at a specific averaging time tau; this is a know function, page 74 of HFSA
# h0_to_adev <- function(h0, tau) {
#   return(sqrt(h0 / (2 * tau)))
# }

### reads in Nick's stability numbers, from 6/4/2025 email from KK
Ratio_rearranged <- read_excel("/home/aak3/NIST/atomic-clock/Data/ClockComparison2025/Ratio_rearranged.xlsx",
                               sheet = "NN")

### reads in the data sizes for each ratio measurement
alldatNvals=read.csv("/home/aak3/NIST/atomic-clock/Code/ComparisonAnalysis2025/alldatNvals.csv")

df=Ratio_rearranged[c(11:23),c(1,3,7,11)]

AvarStabilities <- df %>%
  pivot_longer(
    cols = -date,
    names_to = "variable",
    values_to = "Stability"
  ) %>%
  mutate(
    label = gsub("/", "", sub("_.*", "", variable)),  # extract before "_" and remove "/"
    VarianceFromAvar=Stability^2
  )

# Clean up 'alldatNvals' to ensure dates match
alldatNvals_clean <- alldatNvals %>%
  mutate(date = as.Date(date)) %>%
  rename(label = ratio) # Rename 'ratio' to 'label' to match AvarStabilities

# Clean up 'AvarStabilities', 'date' is <dttm> (datetime), so we convert to Date
AvarStabilities_clean <- AvarStabilities %>%
  mutate(date = as.Date(date))

# Merge and Calculate
AvarApproach_h0Est <- inner_join(alldatNvals_clean, AvarStabilities_clean, by = c("label", "date")) %>%
  mutate(AvarApproach_h0Est = VarianceFromAvar * non_na_meas_count*2)

# New plot

all_results$Date <- as.Date(all_results$Date)
AvarApproach_h0Est$date <- as.Date(AvarApproach_h0Est$date)

short_label <- substr(dataName, 1, 4)

avar_rows <- AvarApproach_h0Est %>%
  filter(label == short_label) %>%   # Ensure dataName matches 'label' in Avar table
  filter(date %in% unique(all_results$Date)) %>% # Only keep dates present in the main results
  transmute(
    Date   = date,
    Method = "Avar Derived",      # Create a new method label
    h0     = AvarApproach_h0Est,
    lower  = NA,                  # No CI for this method yet
    upper  = NA
  )

all_results$Method <- as.character(all_results$Method)
all_results_augmented <- bind_rows(all_results, avar_rows)

all_results_augmented$Method <- factor(all_results_augmented$Method,
                                       levels = c("Simple Fit", "Bayesian", "Avar Derived"))

p <- ggplot(all_results_augmented, aes(x = Date, y = h0, color = Method, group = Method)) +

  # A. Error Bars
  # position_dodge(width) moves the bars apart so they don't overlap
  geom_errorbar(aes(ymin = lower, ymax = upper),
                width = 0.2, size = 1,
                position = position_dodge(width = 0.5)) +

  # B. Main Points
  geom_point(size = 4, position = position_dodge(width = 0.5)) +

  # C. Styling
  # scale_color_manual(values = c("Bayesian" = "blue", "Simple Fit" = "red")) +
  scale_y_log10(labels = scales::scientific) +

  theme_bw(base_size = 14) +
  labs(
    y = expression(paste(h[0])),
    x = NULL
  ) +

  theme(
    legend.position = "bottom",
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank()
  ) +
  facet_wrap(~Date, scales = "free_x")

# Print and Save
print(p)

save_filename <- paste0(Output_SimpleFit_Dir, "compare_h0_all_SELECT_DATES",dataName,".png")
ggsave(save_filename, plot = p, width = 8, height = 4)
message(paste("Saved combined plot to:", save_filename))
