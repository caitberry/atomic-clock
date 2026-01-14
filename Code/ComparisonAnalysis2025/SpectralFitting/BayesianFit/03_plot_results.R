# ==============================================================================
# 03_plot_results.R
# Visualize Posteriors and Fits
# ==============================================================================
# source("00_common.R")

# Load the original data for comparison
spectralEstDF <- readRDS(DATA_PREP_FILE)

# Pick a date to visualize (or loop through them)
target_date <- unique(spectralEstDF$Date)[1] 

save_name <- paste0(OUTPUT_FIT_DIR, "fit_",dataName,"_Date_", gsub("-", "_", target_date), ".rds")

# readRDS("Code/ComparisonAnalysis2025/SpectralFitting/BayesianFit/Posteriors/")
fit_file <- paste0(OUTPUT_FIT_DIR, "fit_",dataName,"_Date_", gsub("-", "_", target_date), ".rds")

if(!file.exists(fit_file)) stop(paste("Fit file not found:", fit_file))

message(paste("Visualizing results for:", target_date))
fit <- readRDS(fit_file)
data_obs <- spectralEstDF %>% filter(Date == target_date, freq > F_MIN, freq < F_MAX)


# --- PLOT 1: Posterior Predictive Check (The Spaghetti Plot) ---
post_draws <- rstan::extract(fit)
ndraws <- 100
idx_draws <- sample(length(post_draws$h0), ndraws)

# Calculate model lines for random draws
pred_mat <- matrix(NA, nrow=ndraws, ncol=length(data_obs$freq))
for(i in seq_along(idx_draws)) {
  j <- idx_draws[i]
  pred_mat[i,] <- model_psd_r(
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

print(p_fit)


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
    theme_minimal() + theme(legend.position="none", axis.text.y=element_blank())
  if(log_scale) p <- p + scale_x_log10()
  return(p)
}

g1 <- plot_param(comp_data, "h0", TRUE)
g2 <- plot_param(comp_data, "h_m1", TRUE)
g3 <- plot_param(comp_data, "Kp", FALSE)
g4 <- plot_param(comp_data, "Ki", FALSE)
g5 <- plot_param(comp_data, "tau", FALSE)

grid.arrange(g1, g2, g3, g4, g5, ncol=3, top="Prior (Green) vs Posterior (Blue)")
