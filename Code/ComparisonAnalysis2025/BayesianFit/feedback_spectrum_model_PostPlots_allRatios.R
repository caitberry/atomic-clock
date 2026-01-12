# source("/home/aak3/NIST/atomic-clock/Code/ComparisonAnalysis2025/BayesianFit/feedback_spectrum_model_PostPlots_allRatios.R")

# --- GENERALIZED PLOTTING ---

library(rstan)
library(dplyr)
library(ggplot2)
library(stringr)
library(tidyr)
library(scales) 

# --- 1. MODEL FUNCTION  ---
model_psd_r <- function(omega, h0, h_m1, Kp, Ki, tau, tp) {
  # ... function code ...
  n <- length(omega); y_model <- numeric(n)
  for (i in seq_len(n)) {
    w <- omega[i]; re_num <- Kp * cos(-w * tp / 2) + (Ki / w) * sin(-w * tp / 2); im_num <- Kp * sin(-w * tp / 2) - (Ki / w) * cos(-w * tp / 2); re_den <- 1; im_den <- w * tau; denom <- re_den^2 + im_den^2; re_G <- (re_num * re_den + im_num * im_den) / denom; im_G <- (im_num * re_den - re_num * im_den) / denom; re_1pG <- 1 + re_G; im_1pG <- im_G; denom2 <- re_1pG^2 + im_1pG^2; re_ratio <- (re_G * re_1pG + im_G * im_1pG) / denom2; im_ratio <- (im_G * re_1pG - re_G * im_1pG) / denom2; abs_ratio_sq <- re_ratio^2 + im_ratio^2; abs_inv_sq <- 1 / denom2; y_model[i] <- abs_ratio_sq * h0 + abs_inv_sq * (2 * pi * h_m1 / w)
  }
  return(y_model)
}

# --- 2. SETUP  ---
results_folder <- "/home/aak3/NIST/atomic-clock/Results/ClockComp2025"
stan_fits_folder <- "/home/aak3/NIST/atomic-clock/Code/ComparisonAnalysis2025/BayesianFit/Posteriors"
output_plots_folder <- "/home/aak3/NIST/atomic-clock/Code/ComparisonAnalysis2025/BayesianFit/Plots"
if (!dir.exists(output_plots_folder)) { dir.create(output_plots_folder, recursive = TRUE) }
dataset_names <- c("YbSr","AlYb", "AlSr")
myK <- 5

# --- 3. MAIN LOOP TO PROCESS EACH DATASET ---
for (dataset_name in dataset_names) {
  
  cat("\n====================================================\n")
  cat("=== Starting Plot Generation for:", dataset_name, "===\n")
  cat("====================================================\n")
  
  if (grepl("Al", dataset_name)) { tp_value <- 1 } else { tp_value <- 1 }
  
  # --- Load original data  ---
  spectral_files <- list.files(path = results_folder, pattern = paste0("spectral.*", dataset_name, ".*Apr15"), full.names = TRUE)
  if (length(spectral_files) == 0) { warning(paste("No spectral files for", dataset_name, "- skipping.")); next }
  extract_from_list <- function(x, filename) { list(spectrum = data.frame(freq = x[[1]], spectrum = x[[2]], File = filename)) }
  all_data <- lapply(spectral_files, function(file) readRDS(file) %>% extract_from_list(basename(file)))
  spectralEstDF <- do.call(rbind, lapply(all_data, `[[`, "spectrum"))
  spectralEstDF <- spectralEstDF %>%
    mutate(Date = as.Date(paste0("2025-", str_extract(File, "(?<=For[A-Za-z]{4})(\\d{2})"), "-", str_extract(File, "(?<=For[A-Za-z]{4}\\d{2})(\\d{2})"))))
  
  stan_fit_files <- list.files(path = stan_fits_folder, pattern = paste0("fit_", dataset_name, "_Date.*\\.rds"), full.names = TRUE)
  if (length(stan_fit_files) == 0) { warning(paste("No Stan fit files for", dataset_name, "- skipping.")); next }
  
  # --- Collect plot data ---
  all_obs_data <- list(); all_pred_data <- list(); all_median_data <- list(); all_parameter_summaries <- list()
  
  for (i in 1:length(stan_fit_files)) {
    fit_file_path <- stan_fit_files[i] # Get the file path for the current iteration
    
    date_str_underscore <- str_extract(basename(fit_file_path), "(?<=Date).*(?=\\.rds)")
    if(is.na(date_str_underscore)) next
    
    current_date <- as.Date(gsub("_", "-", date_str_underscore))
    cat("  -> Collecting data for date:", as.character(current_date), "\n")
    
    stan_fit <- readRDS(fit_file_path)
    fit_input_data <- filter(spectralEstDF, Date == current_date, freq > 0.0001, freq < 0.1)
    
    freq <- fit_input_data$freq; omega <- 2 * pi * freq; post <- rstan::extract(stan_fit)
    y_post_pred <- matrix(NA, nrow = 100, ncol = length(omega))
    
    for (k in 1:100) { # Use a different index 'k' to avoid confusion with 'i'
      j <- sample(length(post$h0), 1)
      y_post_pred[k, ] <- model_psd_r(omega, post$h0[j], post$h_m1[j], post$Kp[j], post$Ki[j], post$tau[j], tp = tp_value)
    }
    
    all_obs_data[[i]] <- data.frame(freq = freq, psd = fit_input_data$spectrum, Date = current_date)
    all_pred_data[[i]] <- data.frame(freq = rep(freq, 100), psd = as.vector(t(y_post_pred)), draw = factor(rep(1:100, each=length(freq))), Date = current_date)
    all_median_data[[i]] <- data.frame(freq = freq, median_psd = apply(y_post_pred, 2, median), Date = current_date)
    param_summary <- rstan::summary(stan_fit, pars = c("h0", "h_m1", "Kp", "Ki", "tau"))$summary
    all_parameter_summaries[[i]] <- as.data.frame(param_summary) %>% tibble::rownames_to_column("Parameter") %>% mutate(Date = current_date)
  }
  
  # --- Create and save plots  ---
  if (length(all_obs_data) > 0) {
    combined_obs_df <- bind_rows(all_obs_data)
    combined_pred_df <- bind_rows(all_pred_data)
    combined_median_df <- bind_rows(all_median_data)
    
    p_faceted <- ggplot(combined_pred_df, aes(x = freq, y = psd)) +
      geom_line(aes(group = draw), alpha = 0.15, color = "dodgerblue") +
      geom_point(data = combined_obs_df, color = "black", size = 1.5, alpha = 0.8) +
      geom_line(data = combined_median_df, aes(y = median_psd), color = "red", size = 0.8) +
      facet_wrap(~ Date, scales = "free_y") +
      scale_y_log10(labels = trans_format("log10", math_format(10^.x))) + scale_x_log10(labels = trans_format("log10", math_format(10^.x))) +
      labs(title = paste("Daily Posterior Predictive Fits for", dataset_name), subtitle = "Black points: Observed data. Blue lines: Posterior draws. Red line: Posterior median.", x = "Frequency (Hz)", y = "Power Spectral Density") + theme_bw()
    
    ggsave(file.path(output_plots_folder, paste0("all_days_faceted_fit_plot_", dataset_name, ".png")), plot = p_faceted, width = 12, height = 8, dpi = 300)
    cat("  -> Faceted plot saved.\n")
  }
  
  if(length(all_parameter_summaries) > 0) {
    parameter_evolution_df <- bind_rows(all_parameter_summaries)
    param_plot <- ggplot(parameter_evolution_df, aes(x = Date, y = `50%`, group = Parameter)) +
      geom_pointrange(aes(ymin = `2.5%`, ymax = `97.5%`, color = Parameter)) +
      facet_wrap(~Parameter, scales = "free_y") +
      labs(title = paste("Evolution of Model Parameters Over Time for", dataset_name), x = "Date", y = "Parameter Value (with 95% CI)") + theme_bw() + theme(legend.position = "none")
    ggsave(file.path(output_plots_folder, paste0("parameter_evolution_", dataset_name, ".png")), plot = param_plot, width = 10, height = 7, dpi = 300)
    cat("  -> Parameter evolution plot saved.\n")
  }
}

cat("\nAll processing and plotting complete for all datasets.\n")

###################################### posterior WN part

# --- FINAL ANALYSIS: POSTERIOR OF SQRT(h0 / (2*N)) ---

cat("\n====================================================\n")
cat("=== Generating Allan Deviation Summary Plots ===\n")
cat("====================================================\n")

# --- STEP 1: DEFINE YOUR 'N' VALUES FOR EACH DAY ---

alldatNvals=read.csv("/home/aak3/NIST/atomic-clock/Code/ComparisonAnalysis2025/alldatNvals.csv")
alldatNvals$N=alldatNvals$non_na_meas_count
alldatNvals$date
# The script will now loop through each dataset (AlYb, AlSr, YbSr) and create a summary for each.
for (dataset_name in dataset_names) {
  
  cat("\n--- Processing for dataset:", dataset_name, "---\n")
  
  # Find the Stan fit files for the current dataset
  stan_fit_files <- list.files(
    path = stan_fits_folder,
    pattern = paste0("fit_", dataset_name, "_Date.*\\.rds"),
    full.names = TRUE
  )
  
  if (length(stan_fit_files) == 0) {
    warning(paste("No Stan fit files found for", dataset_name, "- skipping summary."))
    next
  }
  
  # List to store the summary for each day
  allan_dev_summaries <- list()
  
  # --- STEP 2: LOOP THROUGH FITS TO CALCULATE THE NEW QUANTITY ---
  for (i in 1:length(stan_fit_files)) {
    fit_file_path <- stan_fit_files[i]
    
    # Get the date for the current file
    date_str_underscore <- str_extract(basename(fit_file_path), "(?<=Date).*(?=\\.rds)")
    if(is.na(date_str_underscore)) next
    current_date <- as.Date(gsub("_", "-", date_str_underscore))
    
    # Look up the 'N' value for this specific date
    N_val <- alldatNvals %>% filter(date == current_date, ratio==dataset_name) %>% pull(N)
    
    if (length(N_val) == 0) {
      cat("  -> No 'N' value found for", as.character(current_date), "- skipping.\n")
      next
    }
    
    cat("  -> Calculating for", as.character(current_date), "with N =", N_val, "\n")
    
    # Load the Stan fit and extract posterior samples for h0
    stan_fit <- readRDS(fit_file_path)
    post_samples <- rstan::extract(stan_fit)
    h0_samples <- post_samples$h0
    
    # --- STEP 3: TRANSFORM EVERY POSTERIOR SAMPLE ---
    # This creates a new vector representing the posterior distribution of our target quantity.
    # This is the Allan Deviation, sigma_y(tau=N).
    adev_samples <- sqrt(h0_samples / (2 * N_val))
    
    # --- STEP 4: SUMMARIZE THE NEW DISTRIBUTION ---
    # Calculate the mean and 95% credible interval from the transformed samples.
    summary_stats <- data.frame(
      Date = current_date,
      mean_adev = mean(adev_samples),
      lower_ci_95 = quantile(adev_samples, probs = 0.025),
      upper_ci_95 = quantile(adev_samples, probs = 0.975)
    )
    
    # Add the results to our list
    allan_dev_summaries[[i]] <- summary_stats
  }
  
  # --- STEP 5: COMBINE, PRINT, AND PLOT THE RESULTS ---
  if (length(allan_dev_summaries) > 0) {
    
    # Combine the daily summaries into one data frame for this dataset
    final_adev_df <- bind_rows(allan_dev_summaries)
    
    # Print the summary table to the console
    cat("\n--- Summary for", dataset_name, "---\n")
    print(final_adev_df)
    
    final_adev_df$ratio=dataset_name
    
    adevEsts = final_adev_df %>% 
      mutate(adev=mean_adev,
             type="Complex Spectral Model",
             adev_lower = lower_ci_95, 
             adev_upper = upper_ci_95) %>%
      dplyr::select(Date,ratio,adev,adev_lower,adev_upper,type)
    
    write.csv(adevEsts,paste0("/home/aak3/NIST/atomic-clock/Code/ComparisonAnalysis2025/",dataset_name,"ComplexModel_avarEsts.csv"),row.names = F)
    
    # Create the summary plot
    p_adev <- ggplot(final_adev_df, aes(x = Date, y = mean_adev)) +
      geom_point(color = "dodgerblue4", size = 2) +
      geom_errorbar(aes(ymin = lower_ci_95, ymax = upper_ci_95), width = 0.5, color = "dodgerblue4") +
      scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
      labs(
        title = paste("Allan Deviation at τ=N from PSD fits for", dataset_name),
        subtitle = "Points are posterior means; bars are 95% credible intervals.",
        x = "Date",
        y = expression(paste("Allan Deviation ", sigma_y(tau=N)))
      ) +
      theme_minimal()
    
    # Save the plot
    output_filename_adev <- file.path(output_plots_folder, paste0("allan_deviation_summary_", dataset_name, ".png"))
    ggsave(output_filename_adev, plot = p_adev, width = 8, height = 5, dpi = 300)
    cat("  -> Summary plot saved to:", basename(output_filename_adev), "\n")
    
  } else {
    cat("\nNo results to summarize for", dataset_name, ".\n")
  }
}

cat("\n\nFinal summary analysis complete.\n")






#################### plot all together

allFiles=c(
  "Code/ComparisonAnalysis2025/YbSravarEsts.csv",
  "Code/ComparisonAnalysis2025/AlSravarEsts.csv",
  "Code/ComparisonAnalysis2025/AlYbavarEsts.csv",
  "Code/ComparisonAnalysis2025/YbSrComplexModel_avarEsts.csv",
  "Code/ComparisonAnalysis2025/AlYbComplexModel_avarEsts.csv",
  "Code/ComparisonAnalysis2025/AlSrComplexModel_avarEsts.csv"
)

allAdevRes=data.frame()
for( i in 1:6){
  temp=read.csv(allFiles[i])
  
  allAdevRes=bind_rows(temp,allAdevRes)
}

allAdevRes$Date=as.Date(allAdevRes$Date)
head(allAdevRes)
ggplot(allAdevRes,aes(Date,adev,color=type,ymin=adev_lower,ymax=adev_upper))+
  geom_point()+
  geom_errorbar()+
  facet_wrap(~ratio, scales = "free_y",ncol=1)
