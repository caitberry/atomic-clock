#########################################################################
### NEW PLOT (CORRECTED): H0 PRIOR VS POSTERIOR VS DATA COMPARISON
### This section generates a plot to diagnose if the prior on h0
### is overly restrictive by comparing its quantiles to the posterior
### and the lowest values of the observed spectral data.
###
### CORRECTED to use the actual model prior: h0 ~ lognormal(log(1e-31), 1.0)
#########################################################################

cat("\n========================================================\n")
cat("=== Generating h0 Prior/Posterior Diagnostic Plots ===\n")
cat("========================================================\n")

# --- Loop through each dataset again to create the new plot ---
for (dataset_name in dataset_names) {
  
  cat("\n--- Creating h0 diagnostic plot for:", dataset_name, "---\n")
  
  # --- Find all necessary files for this dataset ---
  spectral_files <- list.files(path = results_folder, pattern = paste0("spectral.*", dataset_name, ".*Apr15"), full.names = TRUE)
  stan_fit_files <- list.files(path = stan_fits_folder, pattern = paste0("fit_", dataset_name, "_Date.*\\.rds"), full.names = TRUE)
  
  if (length(stan_fit_files) == 0 || length(spectral_files) == 0) {
    warning(paste("Missing Stan fit or spectral files for", dataset_name, "- skipping diagnostic plot."))
    next
  }
  
  # --- Reload spectral data to ensure correct mapping ---
  extract_from_list <- function(x, filename) { list(spectrum = data.frame(freq = x[[1]], spectrum = x[[2]], File = filename)) }
  all_data <- lapply(spectral_files, function(file) readRDS(file) %>% extract_from_list(basename(file)))
  spectralEstDF <- do.call(rbind, lapply(all_data, `[[`, "spectrum"))
  spectralEstDF <- spectralEstDF %>%
    mutate(Date = as.Date(paste0("2025-", str_extract(File, "(?<=For[A-Za-z]{4})(\\d{2})"), "-", str_extract(File, "(?<=For[A-Za-z]{4}\\d{2})(\\d{2})"))))
  
  # --- List to store comparison data for each day ---
  h0_comparison_list <- list()
  
  # --- Loop through each day's fit ---
  for (fit_file_path in stan_fit_files) {
    
    # Extract date and find corresponding data
    date_str_underscore <- str_extract(basename(fit_file_path), "(?<=Date).*(?=\\.rds)")
    if(is.na(date_str_underscore)) next
    current_date <- as.Date(gsub("_", "-", date_str_underscore))
    
    fit_input_data <- filter(spectralEstDF, Date == current_date, freq > 0.0001, freq < 0.1)
    if(nrow(fit_input_data) == 0) next
    
    cat("  -> Processing date:", as.character(current_date), "\n")
    
    # --- 1. GET POSTERIOR QUANTILES ---
    stan_fit <- readRDS(fit_file_path)
    post_samples <- rstan::extract(stan_fit, pars = "h0")$h0
    posterior_summary <- data.frame(
      type = "Posterior",
      median = median(post_samples),
      lower = quantile(post_samples, 0.025),
      upper = quantile(post_samples, 0.975),
      Date = current_date
    )
    
    # --- 2. GET PRIOR QUANTILES (Corrected) ---
    # The prior from the Stan model is h0 ~ lognormal(log(1e-31), 1.0).
    # We use R's qlnorm() function with the corresponding parameters.
    prior_meanlog <- log(1e-31)
    prior_sdlog <- 1.0
    prior_summary <- data.frame(
      type = "Prior",
      median = qlnorm(0.5, meanlog = prior_meanlog, sdlog = prior_sdlog),
      lower  = qlnorm(0.025, meanlog = prior_meanlog, sdlog = prior_sdlog),
      upper  = qlnorm(0.975, meanlog = prior_meanlog, sdlog = prior_sdlog),
      Date   = current_date
    )
    
    # --- 3. GET DATA QUANTILES (from lowest 10 PSD values) ---
    # This provides an empirical estimate of the white noise floor from the data itself.
    lowest_10_psd <- sort(fit_input_data$spectrum)[1:10]
    data_summary <- data.frame(
      type = "Data (Lowest 10)",
      median = median(lowest_10_psd),
      lower = quantile(lowest_10_psd, 0.025),
      upper = quantile(lowest_10_psd, 0.975),
      Date = current_date
    )
    
    # Combine into a single data frame for this day
    h0_comparison_list[[as.character(current_date)]] <- bind_rows(prior_summary, data_summary, posterior_summary)
  }
  
  # --- Combine all daily data and create the plot ---
  if (length(h0_comparison_list) > 0) {
    
    comparison_df <- bind_rows(h0_comparison_list)
    
    # Set factor levels for a logical plot order
    comparison_df$type <- factor(comparison_df$type, levels = c("Prior", "Data (Lowest 10)", "Posterior"))
    
    p_h0_diag <- ggplot(comparison_df, aes(x = type, y = median, ymin = lower, ymax = upper, color = type)) +
      geom_pointrange(size = 0.8, fatten = 2.5) +
      facet_wrap(~ Date, scales = "free_y") +
      scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
      scale_color_manual(values = c("Prior" = "red", "Data (Lowest 10)" = "black", "Posterior" = "dodgerblue")) +
      labs(
        title = paste("Comparison of h0 (White Noise Level) for", dataset_name),
        subtitle = "Points are medians; bars represent 95% quantiles. Prior is lognormal(log(1e-31), 1.0).",
        x = "Source of Distribution",
        y = "Value of h0 (Power Spectral Density)"
      ) +
      theme_bw() +
      theme(legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1))
    
    # Save the plot
    output_filename_h0 <- file.path(output_plots_folder, paste0("h0_prior_posterior_comparison_", dataset_name, ".png"))
    ggsave(output_filename_h0, plot = p_h0_diag, width = 12, height = 9, dpi = 300)
    cat("  -> h0 diagnostic plot saved to:", basename(output_filename_h0), "\n")
    
  } else {
    cat("\nNo results to plot for h0 diagnostic for", dataset_name, ".\n")
  }
}

cat("\n\nAll h0 diagnostic plotting complete.\n")