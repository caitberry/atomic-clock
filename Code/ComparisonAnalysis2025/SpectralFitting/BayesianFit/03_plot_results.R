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
    log_h0   = rnorm(n_p, -70, 10),
    log_h_m1 = rnorm(n_p,-75.63, 0.17),
    tau = runif(n_p,1,200),
    Kp  = runif(n_p,0, 20),
    Ki  = runif(n_p,0, 10),
    # h_m1 = rlnorm(n_p, -75.79, 1),
    # Kp   = rlnorm(n_p, -1, 2),#abs(rnorm(n_p, 10, 5)),
    # Ki   = rlnorm(n_p, -1, 2),
    # tau  = rlnorm(n_p, 1, 1.5),
    # Kp   = rlnorm(n_p, 0, 5),#abs(rnorm(n_p, 10, 5)),
    # Ki   = rlnorm(n_p, 0, 5),
    # tau  = rlnorm(n_p, 0, 5),
    Type = "Prior"
  )

  posteriors <- as.data.frame(fit) %>% 
    select(log_h0, log_h_m1, Kp, Ki, tau) %>% 
    mutate(Type = "Posterior")
  
  comp_data <- bind_rows(priors, posteriors)
  
  # limits   <- range(posteriors$h0)
  # ggplot(posteriors,aes(h0))+
  #   geom_density()+
  #   geom_density(data=priors,aes(h0))+
  #   coord_cartesian(xlim=limits)
  # 
  
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
    
    # B. Determine Plot Limits (Zoom in to avoid long tails)
    # We use the 0.5% and 99.5% quantiles of the COMBINED data to set limits.
    # This ignores the extreme <1% outliers that squash the plot.
    # vals_all <- d_sub[d_sub$Type == "Posterior", param]#d_sub[[param]]
    # limits   <- range(vals_all)#quantile(vals_all, probs = c(0.005, 0.995))
    
    # Build Plot
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
    
    # Apply Log Scale if requested
    if(log_scale) {
      p <- p + scale_x_log10() + annotation_logticks(sides="b")
    }
    
    # Apply Smart Zoom (coord_cartesian keeps the data but zooms the view)
    # p <- p + coord_cartesian(xlim = limits)
    
    return(p)
  }
  
  g1 <- plot_param(comp_data, "log_h0", F)
  g2 <- plot_param(comp_data, "log_h_m1", F)
  g3 <- plot_param(comp_data, "Kp", F)
  g4 <- plot_param(comp_data, "Ki", F)
  g5 <- plot_param(comp_data, "tau",F)
  
  # SAVE PLOT
  g_combined <- arrangeGrob(g1, g2, g3, g4, g5, ncol=3) 
                            # top = textGrob(paste("Prior vs Post:", target_date), 
                                           # gp=gpar(fontsize=14, fontface="bold")))  
  # posteriors <- as.data.frame(fit) %>% 
  #   select(h0, h_m1, Kp, Ki, tau) %>% 
  #   mutate(Type = "Posterior")
  # 
  # comp_data <- bind_rows(priors, posteriors)
  # 
  # plot_param <- function(dat, param, log_scale=FALSE) {
  #   p <- ggplot(dat, aes_string(x=param, fill="Type")) +
  #     geom_density(alpha=0.4) +
  #     scale_fill_manual(values=c("Posterior"="#1f78b4", "Prior"="#b2df8a")) +
  #     theme_minimal() + 
  #     theme(legend.position="none", 
  #           axis.text.y=element_blank(),
  #           axis.title.x = element_text(size=10))
  #   if(log_scale) p <- p + scale_x_log10()
  #   return(p)
  # }
  # 
  # g1 <- plot_param(comp_data, "h0", TRUE)
  # g2 <- plot_param(comp_data, "h_m1", TRUE)
  # g3 <- plot_param(comp_data, "Kp", FALSE)
  # g4 <- plot_param(comp_data, "Ki", FALSE)
  # g5 <- plot_param(comp_data, "tau", FALSE)
  # 
  # # SAVE PLOT 2
  # # Use arrangeGrob to create an object suitable for ggsave
  # g_combined <- arrangeGrob(g1, g2, g3, g4, g5, ncol=3, top=paste("Prior vs Post:", target_date))
  # 
  
  
  plot2_filename <- paste0(OUTPUT_PLOT_DIR, "prior_post_", dataName, "_", gsub("-", "_", target_date), ".png")
  ggsave(plot2_filename, plot = g_combined, width = 10, height = 6)
  message(paste("Saved prior/post plot to:", plot2_filename))
}
