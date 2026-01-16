# ==============================================================================
# 03_plot_results.R
# Visualize Posteriors and Fits
# ==============================================================================
# source("00_common.R")

# Load the original data for comparison
spectralEstDF <- readRDS(DATA_PREP_FILE)

for(i in 1:length(unique(spectralEstDF$Date))){
  target_date=unique(spectralEstDF$Date)[i]
  
  # Construct the expected fit filename
  fit_file <- paste0(OUTPUT_FIT_DIR, "fit_", dataName, "_Date_", gsub("-", "_", target_date), ".rds")
  
  # Skip if fit file is missing
  if(!file.exists(fit_file)) {
    warning(paste("Fit file not found for date:", target_date, "- Skipping."))
    next
  }
  
  message(paste("Processing plots for:", target_date))
  fit <- readRDS(fit_file)
  data_obs <- spectralEstDF %>% filter(Date == target_date, freq > F_MIN, freq < F_MAX)
  
  # # --- Looking at posterior results ---
  # 
  # params_of_interest <- c("h0", "h_m1", "Kp", "Ki", "tau")
  # 
  # # Extract summary matrix
  # fit_summary <- summary(fit, pars = params_of_interest)$summary
  # 
  # # Convert to data frame and select columns
  # results_df <- data.frame(
  #   Parameter = rownames(fit_summary),
  #   Mean      = fit_summary[, "mean"],
  #   SD        = fit_summary[, "sd"]
  # )
  # 
  # print(results_df)
  
  # --- PLOT 1: Posterior Predictive Check (Spaghetti Plot) ---
  post_draws <- rstan::extract(fit)
  ndraws <- 100
  idx_draws <- sample(length(post_draws$h0), ndraws)
  
  # Calculate model lines for random draws
  pred_mat <- matrix(NA, nrow=ndraws, ncol=length(data_obs$freq))
  for(k in seq_along(idx_draws)) {
    j <- idx_draws[k]
    # Note: Ensure model_psd_r is loaded from 00_common.R
    pred_mat[k,] <- model_psd_r(
      data_obs$freq, post_draws$h0[j], post_draws$h_m1[j], 
      post_draws$Kp[j], post_draws$Ki[j], post_draws$tau[j], TP_VAL
    )
  }
  
  df_lines <- data.frame(
    freq = rep(data_obs$freq, ndraws),
    psd  = as.vector(t(pred_mat)),
    draw = factor(rep(1:ndraws, each=length(data_obs$freq)))
  )
  
  df_median <- data.frame(
    freq = data_obs$freq,
    median = apply(pred_mat, 2, median)
  )
  
  p_fit <- ggplot() +
    geom_line(data=df_lines, aes(freq, psd, group=draw), alpha=0.1, col="skyblue") +
    geom_point(data=data_obs, aes(freq, spectrum), size=1, alpha=0.6) +
    geom_line(data=df_median, aes(freq, median), col="blue", size=1) +
    scale_x_log10() + scale_y_log10() +
    labs(title=paste("Posterior Fit:", target_date), y="PSD", x="Freq (Hz)") +
    theme_minimal()
  
  # SAVE PLOT 1
  plot1_filename <- paste0(OUTPUT_PLOT_DIR, "fit_plot_", dataName, "_", gsub("-", "_", target_date), ".png")
  ggsave(plot1_filename, plot = p_fit, width = 8, height = 6)
  message(paste("Saved fit plot to:", plot1_filename))
  
  
  # --- PLOT 2: Prior vs Posterior Densities ---
  # Generate Prior Samples (Hardcoded to match Stan file)
  n_p <- 4000
  priors <- data.frame(
    h0   = rlnorm(n_p, -72.78, 2.49),
    h_m1 = rlnorm(n_p, -75.79, 0.18),
    Kp   = abs(rnorm(n_p, 10, 5)),
    Ki   = rlnorm(n_p, 0, 1.17),
    tau  = rlnorm(n_p, 2.65, 1.35),
    Type = "Prior"
  )
  
  posteriors <- as.data.frame(fit) %>% 
    select(h0, h_m1, Kp, Ki, tau) %>% 
    mutate(Type = "Posterior")
  
  comp_data <- bind_rows(priors, posteriors)
  
  plot_param <- function(dat, param, log_scale=FALSE) {
    p <- ggplot(dat, aes_string(x=param, fill="Type")) +
      geom_density(alpha=0.4) +
      scale_fill_manual(values=c("Posterior"="#1f78b4", "Prior"="#b2df8a")) +
      theme_minimal() + 
      theme(legend.position="none", 
            axis.text.y=element_blank(),
            axis.title.x = element_text(size=10))
    if(log_scale) p <- p + scale_x_log10()
    return(p)
  }
  
  g1 <- plot_param(comp_data, "h0", TRUE)
  g2 <- plot_param(comp_data, "h_m1", TRUE)
  g3 <- plot_param(comp_data, "Kp", FALSE)
  g4 <- plot_param(comp_data, "Ki", FALSE)
  g5 <- plot_param(comp_data, "tau", FALSE)
  
  # SAVE PLOT 2
  # Use arrangeGrob to create an object suitable for ggsave
  g_combined <- arrangeGrob(g1, g2, g3, g4, g5, ncol=3, top=paste("Prior vs Post:", target_date))
  
  plot2_filename <- paste0(OUTPUT_PLOT_DIR, "prior_post_", dataName, "_", gsub("-", "_", target_date), ".png")
  ggsave(plot2_filename, plot = g_combined, width = 10, height = 6)
  message(paste("Saved prior/post plot to:", plot2_filename))
}
