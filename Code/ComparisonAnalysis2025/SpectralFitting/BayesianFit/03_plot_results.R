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
  
  # --- PLOT 1: Posterior Predictive Check (Spaghetti Plot) ---
  post_draws <- rstan::extract(fit)
  ndraws <- 100
  idx_draws <- sample(length(post_draws$log_h0), ndraws)
  
  # Calculate model lines for random draws
  pred_mat <- matrix(NA, nrow=ndraws, ncol=length(data_obs$freq))
  for(k in seq_along(idx_draws)) {
    j <- idx_draws[k]
    # Note: Ensure model_psd_r is loaded from 00_common.R
    pred_mat[k,] <- model_psd_r(
      data_obs$freq, exp(post_draws$log_h0[j]), exp(post_draws$log_h_m1[j]), 
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
    labs(y="PSD", x="Frequency") +
    theme_minimal()
  
  # SAVE PLOT 1
  plot1_filename <- paste0(OUTPUT_PLOT_DIR, "fit_plot_", dataName, "_", gsub("-", "_", target_date), ".png")
  ggsave(plot1_filename, plot = p_fit, width = 5, height = 4)
  message(paste("Saved fit plot to:", plot1_filename))
  
  
  # --- PLOT 2: Prior vs Posterior Densities ---
  # Generate Prior Samples (Hardcoded to match Stan file)
  n_p <- 4000
  priors <- data.frame(
    log_h0   = rnorm(n_p, -72, 2),
    log_h_m1 = rnorm(n_p,-75.63, 0.17),
    tau = runif(n_p,1,200),
    # Kp  = runif(n_p,0, 20),
    Kp   = rlnorm(n_p, 0, 4), #Needs to be truncated at 20, unsure how to do...
    Ki  = runif(n_p,0, 10),

    Type = "Prior"
  )

  posteriors <- as.data.frame(fit) %>% 
    select(log_h0, log_h_m1, Kp, Ki, tau) %>% 
    mutate(Type = "Posterior")
  
  comp_data <- bind_rows(priors, posteriors)
  
  # 2. Updated Plotting Function
  plot_param <- function(dat, param, log_scale=F) {
    
    # Filter for specific parameter
    d_sub <- dat %>% filter(!is.na(!!sym(param)))
    
    # A. Calculate Prior 95% Interval for the Label
    prior_vals <- d_sub$prior_vals <- d_sub[d_sub$Type == "Prior", param]
    prior_CI   <- quantile(prior_vals, probs = c(0.025, 0.975))
    
    # Format nicely (scientific notation for small numbers)
    fmt <- function(x) format(x, digits=2, scientific=TRUE)
    label_text <- paste0("Prior 95%: [", fmt(prior_CI[1]), ", ", fmt(prior_CI[2]), "]")
    

    p <- ggplot(d_sub, aes_string(x=param, fill="Type")) +
      geom_density(alpha=0.4, size=0.2) +
      scale_fill_manual(values=c("Posterior"="#1f78b4", "Prior"="#b2df8a")) +
      
      # Add the 95% range info in the subtitle
      labs(x = NULL, y = NULL, title = param, subtitle = label_text) +
      
      theme_minimal() + 
      theme(
        legend.position="none", 
        axis.text.y=element_blank(),
        plot.title = element_text(size=11, face="bold"),
        plot.subtitle = element_text(size=8, color="gray40")
      )
    
    if(log_scale) {
      p <- p + scale_x_log10() + annotation_logticks(sides="b")
    }
    return(p)
  }
  
  g1 <- plot_param(comp_data, "log_h0", F)
  g2 <- plot_param(comp_data, "log_h_m1", F)
  g3 <- plot_param(comp_data, "Kp", F)
  g4 <- plot_param(comp_data, "Ki", F)
  g5 <- plot_param(comp_data, "tau",F)
  
  # SAVE PLOT
  g_combined <- arrangeGrob(g1, g2, g3, g4, g5, ncol=3) 

  
  plot2_filename <- paste0(OUTPUT_PLOT_DIR, "prior_post_", dataName, "_", gsub("-", "_", target_date), ".png")
  ggsave(plot2_filename, plot = g_combined, width = 5, height = 4)
  message(paste("Saved prior/post plot to:", plot2_filename))
}
