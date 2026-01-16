library(rstan)
library(dplyr)
library(ggplot2)
library(tidyr)

specEstFile="/home/aak3/NIST/atomic-clock/Results/ClockComp2025/spectralEstForAlSr0227_Apr15.Rds"

res=readRDS(specEstFile)

resDF=data.frame(freq=res$freq,spectrum=res$spec.hat)

# resDFfiltered=filter(resDF,freq>(8*10^-4),freq<(2*10^-2))
# > f_c_median
# [1] 0.01718324

# resDFfiltered=filter(resDF,freq>(8*10^-4),freq<(2.5*10^-2))
# > f_c_median
# [1] 0.01806874

# resDFfiltered=filter(resDF,freq>(7*10^-4),freq<(2.5*10^-2))

ggplot(resDFfiltered,aes(freq, spectrum))+
  geom_line(alpha=.4) +
  scale_y_log10() +
  scale_x_log10() +
  theme(legend.position = "bottom",legend.text = element_text(size = 10)) +
  guides(color = guide_legend(title = NULL))

data_list <- list(
  N = length(resDFfiltered$freq),
  f = resDFfiltered$freq,
  S = resDFfiltered$spectrum
)

# Compile and fit the model
fit <- stan(file = "Code/ComparisonAnalysis2025/BayesianFit/LowPassWNmodel.stan",
            data = data_list,
            iter = 10000,
            warmup = 8000,
            chains = 4,
            cores = 4,
            thin = 10,
            seed = 123)


library(dplyr)
library(ggplot2)
library(tidyr)

# Extract posterior samples
posterior_samples <- extract(fit)
n_draws <- length(posterior_samples$S_w)

# Frequencies from your data
freq_vals <- resDFfiltered$freq
n_freq <- length(freq_vals)

# Create a matrix of fitted PSD values: rows = draws, cols = frequencies
psd_matrix <- matrix(NA, nrow = n_draws, ncol = n_freq)

for (i in 1:n_draws) {
  S_w_i <- posterior_samples$S_w[i]
  f_c_i <- posterior_samples$f_c[i]
  psd_matrix[i, ] <- S_w_i / (1 + (freq_vals / f_c_i)^2)
}

# Compute pointwise credible intervals
psd_summary <- data.frame(
  freq = freq_vals,
  psd_median = apply(psd_matrix, 2, median),
  psd_lower = apply(psd_matrix, 2, quantile, probs = 0.025),
  psd_upper = apply(psd_matrix, 2, quantile, probs = 0.975)
)

# Merge with your original data if needed
plot_data <- left_join(resDFfiltered, psd_summary, by = "freq")

# Plot with credible ribbon
ggplot(plot_data, aes(freq, spectrum)) +
  geom_line(alpha = 0.4) +
  geom_ribbon(aes(ymin = psd_lower, ymax = psd_upper), fill = "blue", alpha = 0.2) +
  geom_line(aes(y = psd_median), color = "blue", size = 1) +
  scale_y_log10() +
  scale_x_log10() +
  theme(legend.position = "bottom", legend.text = element_text(size = 10)) +
  guides(color = guide_legend(title = NULL))




# Extract posterior samples
posterior_samples <- extract(fit)
n_draws <- length(posterior_samples$S_w)

# Subsample (e.g., 50 random draws)
set.seed(123)
sample_indices <- sample(1:n_draws, 50)

# Frequencies from your data
freq_vals <- resDFfiltered$freq
n_freq <- length(freq_vals)

# Create sample curves
sample_curves <- lapply(sample_indices, function(i) {
  S_w_i <- posterior_samples$S_w[i]
  f_c_i <- posterior_samples$f_c[i]
  tibble(
    freq = freq_vals,
    psd = S_w_i / (1 + (freq_vals / f_c_i)^2),
    draw = paste0("draw_", i)
  )
}) %>% bind_rows()

# Compute pointwise median curve from *all* posterior samples
psd_matrix <- sapply(1:n_draws, function(i) {
  S_w_i <- posterior_samples$S_w[i]
  f_c_i <- posterior_samples$f_c[i]
  S_w_i / (1 + (freq_vals / f_c_i)^2)
})

psd_median_df <- tibble(
  freq = freq_vals,
  psd = apply(psd_matrix, 1, median)
)

# Step 1: Compute summary stats for f_c
f_c_median <- median(posterior_samples$f_c)
f_c_lower <- quantile(posterior_samples$f_c, 0.025)
f_c_upper <- quantile(posterior_samples$f_c, 0.975)
S_w_est=mean(posterior_samples$S_w)

# Step 2: Add vertical lines to the existing plot
ggplot(resDFfiltered, aes(freq, spectrum)) +
  geom_line(alpha = 0.4) +
  geom_line(data = sample_curves, aes(x = freq, y = psd, group = draw), 
            color = "blue", alpha = 0.3) +
  geom_line(data = psd_median_df, aes(x = freq, y = psd), 
            color = "darkblue", size = 1.2) +
  geom_vline(xintercept = f_c_median, color = "red", linetype = "solid", size = 1) +
  geom_vline(xintercept = f_c_lower, color = "red", linetype = "dashed") +
  geom_vline(xintercept = f_c_upper, color = "red", linetype = "dashed") +
  scale_y_log10() +
  scale_x_log10() +
  theme(legend.position = "none")+
  geom_hline(yintercept = S_w_est, color = "green", linetype = "solid", size = 1)
  


# Step 1: Compute log-residuals
residuals_df <- resDFfiltered %>%
  left_join(psd_median_df, by = "freq") %>%
  mutate(log_residual = log(spectrum) - log(psd))

# Step 2: Residual plot
ggplot(residuals_df, aes(x = freq, y = log_residual)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  geom_point(alpha = 0.6) +
  scale_x_log10() +
  labs(y = "Log Residual", x = "Frequency (Hz)",
       title = "Residuals: log(observed PSD) - log(fitted PSD)") +
  theme_minimal()




###########trying mixed noise

resDFfiltered=filter(resDF,freq>(5*10^-4),freq<(2*10^-1))

ggplot(resDFfiltered,aes(freq, spectrum))+
  geom_line(alpha=.4) +
  scale_y_log10() +
  scale_x_log10() +
  theme(legend.position = "bottom",legend.text = element_text(size = 10)) +
  guides(color = guide_legend(title = NULL))

data_list <- list(
  N = length(resDFfiltered$freq),
  f = resDFfiltered$freq,
  S = resDFfiltered$spectrum
)

# Compile and fit the model
fit <- stan(file = "Code/ComparisonAnalysis2025/BayesianFit/LowPassMixedNoiseModel.stan",
            data = data_list,
            iter = 10000,
            warmup = 8000,
            chains = 4,
            cores = 4,
            thin = 10,
            seed = 123)

# Posterior samples assumed to come from a Stan model
posterior_samples <- extract(fit)
posterior_samples$A_white=posterior_samples$A[,1]
posterior_samples$A_flicker=posterior_samples$A[,2]
posterior_samples$A_rw=posterior_samples$A[,3]

posterior_samples$fc_white=posterior_samples$fc[,1]
posterior_samples$fc_flicker=posterior_samples$fc[,2]
posterior_samples$fc_rw=posterior_samples$fc[,3]

n_draws <- length(posterior_samples$A_white)

# Subsample for visualization
set.seed(123)
sample_indices <- sample(1:n_draws, 50)

# Frequencies (from your data)
freq_vals <- resDFfiltered$freq

# Step 1: Create posterior draw PSDs
sample_curves <- lapply(sample_indices, function(i) {
  f <- freq_vals
  S <- (
    posterior_samples$A_white[i] / (1 + (f / posterior_samples$fc_white[i])^2) +
      (posterior_samples$A_flicker[i] / f) / (1 + (f / posterior_samples$fc_flicker[i])^2) +
      (posterior_samples$A_rw[i] / f^2) / (1 + (f / posterior_samples$fc_rw[i])^2)
  )
  tibble(freq = f, psd = S, draw = paste0("draw_", i))
}) %>% bind_rows()

# Step 2: Median PSD
psd_matrix <- sapply(1:n_draws, function(i) {
  f <- freq_vals
  (
    posterior_samples$A_white[i] / (1 + (f / posterior_samples$fc_white[i])^2) +
      (posterior_samples$A_flicker[i] / f) / (1 + (f / posterior_samples$fc_flicker[i])^2) +
      (posterior_samples$A_rw[i] / f^2) / (1 + (f / posterior_samples$fc_rw[i])^2)
  )
})

psd_median_df <- tibble(
  freq = freq_vals,
  psd = apply(psd_matrix, 1, median)
)

# Step 3: Plot
ggplot(resDFfiltered, aes(freq, spectrum)) +
  geom_line(alpha = 0.4) +
  geom_line(data = sample_curves, aes(x = freq, y = psd, group = draw), 
            color = "blue", alpha = 0.3) +
  geom_line(data = psd_median_df, aes(x = freq, y = psd), 
            color = "darkblue", size = 1.2) +
  scale_y_log10() +
  scale_x_log10() +
  theme_minimal() +
  theme(legend.position = "none") +
  ggtitle("PSD with Three-Component Smoothed Model")
