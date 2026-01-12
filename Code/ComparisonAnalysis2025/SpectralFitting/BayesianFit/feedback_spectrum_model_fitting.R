library(rstan)
options(mc.cores = 4)


library(dplyr)
library(ggplot2)
library(broom)


# Set your folder path
folder <- "/home/aak3/NIST/atomic-clock/Results/ClockComp2025"  

###################################################################################
###change this stuff
myK=5
# List all files 
files_with_AlYb_spectral <- list.files(
  path = folder,
  pattern = "spectral.*AlYb.*Apr15",
  full.names = TRUE
)
####################################################################################
extract_from_list <- function(x, filename) {
  df1 <- data.frame(
    freq = x[[1]],
    spectrum = x[[2]],
    File = filename
  )
  
  df2 <- data.frame(
    eigenval = x[[4]],
    K=1:myK,
    File = filename
  )
  
  list(spectrum = df1, eigenvals = df2)
}

all_data <- lapply(files_with_AlYb_spectral, function(file) {
  x <- readRDS(file)
  extract_from_list(x, basename(file))
})

# Combine all df1s and df2s if desired
spectralEstDF <- do.call(rbind, lapply(all_data, `[[`, "spectrum"))
EVdf <- do.call(rbind, lapply(all_data, `[[`, "eigenvals"))

spectralEstDF <- spectralEstDF %>%
  mutate(
    Ratio = str_extract(File, "(?<=spectralEstFor)[A-Za-z]+"),
    Month = str_extract(File, "(?<=For[A-Za-z]{4})(\\d{2})"),
    Day   = str_extract(File, "(?<=For[A-Za-z]{4}\\d{2})(\\d{2})"),
    Date = as.Date(paste0("2025-", Month, "-", Day))
  )

EVdf <- EVdf %>%
  mutate(
    Ratio = str_extract(File, "(?<=spectralEstFor)[A-Za-z]+"),
    Month = str_extract(File, "(?<=For[A-Za-z]{4})(\\d{2})"),
    Day   = str_extract(File, "(?<=For[A-Za-z]{4}\\d{2})(\\d{2})"),
    Date = as.Date(paste0("2025-", Month, "-", Day))
  )


ggplot(EVdf,aes(K,eigenval,color=factor(Date)))+
  geom_point()+
  guides(color = guide_legend(title = NULL))

ggplot(spectralEstDF,aes(freq, spectrum, col = factor(Date)))+
  geom_line(alpha=.4) +
  scale_y_log10() +
  scale_x_log10() +
  theme(legend.position = "bottom",legend.text = element_text(size = 10)) +
  guides(color = guide_legend(title = NULL))



K <- 5

# Confidence intervals
spectralEstDF <- spectralEstDF %>%
  mutate(
    lower = spectrum * qchisq(0.025, df = 2 * K) / (2 * K),
    upper = spectrum * qchisq(0.975, df = 2 * K) / (2 * K),
    sdish=sqrt((2 * spectrum^2)/K)
  )



for(i in 1:length(unique(spectralEstDF$Date))){
  d=unique(spectralEstDF$Date)[i]
  plot_data <- filter(spectralEstDF, Date == d)
  fit_data <- filter(plot_data, freq > 0.0001, freq < .1)
  
  # fit <- lm(spectrum ~ I(freq^2), data = fit_data)
  # fit_summary <- summary(fit)
  
  # fit_data=fit_data[-1,]
  stan_data <- list(
    N = length(fit_data$freq),
    omega = 2 * pi * fit_data$freq,   # Convert to angular frequency
    y_obs = fit_data$spectrum,
    # y_sd = fit_data$sdish
    rel_sd = fit_data$sdish/fit_data$spectrum
  )
  plot(stan_data$rel_sd)
  
  # data_list <- list(
  #   N = length(spectralEstDF$freq),
  #   f = spectralEstDF$freq,
  #   y_obs = spectralEstDF$spectrum,
  #   sigma = 10^-31  # or a value you believe is realistic for log-spectrum error
  # )
  
  fit <- stan(
    file = "Code/ComparisonAnalysis2025/BayesianFit/feedback_spectrum_model.stan",
    data = stan_data,
    iter = 4000,
    warmup = 3000,
    chains = 4  
    # control = list(adapt_delta = 0.96)
  )
  
  saveRDS(fit, file = paste0("/home/aak3/NIST/atomic-clock/Code/ComparisonAnalysis2025/BayesianFit/Posteriors/fit_AlYb_Date", gsub("-", "_", d), ".rds"))
  
}



#############looking at results


print(fit, pars = c("h0", "h_m1", "Kp", "Ki", "tau"))

# Extract h0 samples and compute 95% CI
posterior_samples <- rstan::extract(fit, pars = "h0")$h0
ci <- quantile(posterior_samples, probs = c(0.025, 0.975))
print(ci)


model_psd_r <- function(omega, h0, h_m1, Kp, Ki, tau) {
  tp <- 1  # fixed value
  n <- length(omega)
  y_model <- numeric(n)
  
  for (i in seq_len(n)) {
    w <- omega[i]
    
    # Numerator G(w) components
    re_num <- Kp * cos(-w * tp / 2) + (Ki / w) * sin(-w * tp / 2)
    im_num <- Kp * sin(-w * tp / 2) - (Ki / w) * cos(-w * tp / 2)
    
    # Denominator G(w) components
    re_den <- 1
    im_den <- w * tau
    
    # Complex division num/den
    denom <- re_den^2 + im_den^2
    re_G <- (re_num * re_den + im_num * im_den) / denom
    im_G <- (im_num * re_den - re_num * im_den) / denom
    
    # 1 + G(w)
    re_1pG <- 1 + re_G
    im_1pG <- im_G
    denom2 <- re_1pG^2 + im_1pG^2
    
    # |G/(1+G)|^2 and |1/(1+G)|^2
    re_ratio <- (re_G * re_1pG + im_G * im_1pG) / denom2
    im_ratio <- (im_G * re_1pG - re_G * im_1pG) / denom2
    abs_ratio_sq <- re_ratio^2 + im_ratio^2
    abs_inv_sq <- 1 / denom2
    
    # Model PSD
    y_model[i] <- abs_ratio_sq * h0 + abs_inv_sq * (2 * pi * h_m1 / w)
  }
  
  return(y_model)
}

freq=fit_data$freq
psd=fit_data$spectrum
omega <- 2 * pi * fit_data$freq
post <- rstan::extract(fit)
n_post <- length(post$h0)
ndraws <- 100
draws_idx <- sample(n_post, ndraws)

y_post_pred <- matrix(NA, nrow = ndraws, ncol = length(omega))

for (i in seq_along(draws_idx)) {
  j <- draws_idx[i]
  y_post_pred[i, ] <- model_psd_r(
    omega,
    post$h0[j],
    post$h_m1[j],
    post$Kp[j],
    post$Ki[j],
    post$tau[j]
  )
}

df_obs <- data.frame(freq = freq, psd = psd)

df_pred_all <- data.frame(
  freq = rep(freq, times = ndraws),
  psd = as.vector(t(y_post_pred)),
  draw = factor(rep(1:ndraws, each = length(freq)))
)

library(ggplot2)
ggplot(df_pred_all, aes(x = freq, y = psd)) +
  geom_line(aes(group = draw), alpha = 0.1, color = "blue") +
  geom_point(data = df_obs, aes(x = freq, y = psd), color = "black") +
  scale_y_log10() +
  scale_x_log10() +
  labs(x = "Frequency (Hz)", y = "PSD (linear scale)",
       title = "Posterior predictive lines and observed PSD") +
  theme_minimal()

library(dplyr)
df_pred <- data.frame(
  freq = rep(freq, ndraws),
  psd = as.vector(t(y_post_pred))
)

df_summary <- df_pred %>%
  group_by(freq) %>%
  summarize(
    median = median(psd),
    lower = quantile(psd, 0.025),
    upper = quantile(psd, 0.975),
    .groups = "drop"
  )

ggplot() +
  geom_point(data = df_obs, aes(x = freq, y = psd), color = "black") +
  geom_line(data = df_summary, aes(x = freq, y = median), color = "blue") +
  geom_ribbon(data = df_summary, aes(x = freq, ymin = lower, ymax = upper),
              alpha = 0.2, fill = "blue") +
  scale_y_log10() +
  scale_x_log10() +
  labs(x = "Frequency (Hz)", y = "PSD (linear scale)",
       title = "Posterior predictive vs observed PSD") +
  theme_minimal()

################## # est tp
# print(fit, pars = c("h0", "h_m1", "Kp", "Ki", "tp", "tau"))
# 
# # Assuming your fit object is called 'fit'
# posterior_samples <- rstan::extract(fit, pars = "h0")$h0
# 
# # Compute 95% credible interval
# ci <- quantile(posterior_samples, probs = c(0.025, 0.975))
# print(ci)
# 
# ###########################################################3
# ### post curves
# ###########################################################3
# freq=fit_data$freq
# psd= fit_data$spectrum
# 
# 
# post <- rstan::extract(fit)
# n_post <- length(post$h0)  # number of posterior draws
# 
# model_psd_r <- function(omega, h0, h_m1, Kp, Ki, tp, tau) {
#   n <- length(omega)
#   y_model <- numeric(n)
#   
#   for (i in seq_len(n)) {
#     w <- omega[i]
#     # Numerator G(w) components
#     re_num <- Kp * cos(-w * tp / 2) + (Ki / w) * sin(-w * tp / 2)
#     im_num <- Kp * sin(-w * tp / 2) - (Ki / w) * cos(-w * tp / 2)
#     
#     # Denominator G(w) components
#     re_den <- 1
#     im_den <- w * tau
#     
#     # Complex division num/den
#     denom <- re_den^2 + im_den^2
#     re_G <- (re_num * re_den + im_num * im_den) / denom
#     im_G <- (im_num * re_den - re_num * im_den) / denom
#     
#     # 1 + G(w)
#     re_1pG <- 1 + re_G
#     im_1pG <- im_G
#     
#     denom2 <- re_1pG^2 + im_1pG^2
#     
#     # |G/(1+G)|^2
#     re_ratio <- (re_G * re_1pG + im_G * im_1pG) / denom2
#     im_ratio <- (im_G * re_1pG - re_G * im_1pG) / denom2
#     abs_ratio_sq <- re_ratio^2 + im_ratio^2
#     
#     # |1/(1+G)|^2
#     abs_inv_sq <- 1 / denom2
#     
#     y_model[i] <- abs_ratio_sq * h0 + abs_inv_sq * (2 * pi * h_m1 / w)
#   }
#   
#   return(y_model)
# }
# 
# 
# omega <- 2 * pi * freq  # your frequency vector in rad/s
# 
# # Sample some posterior draws (e.g., 100 draws for plotting)
# ndraws <- 100
# draws_idx <- sample(n_post, ndraws)
# 
# y_post_pred <- matrix(NA, nrow = ndraws, ncol = length(omega))
# 
# for (i in seq_along(draws_idx)) {
#   j <- draws_idx[i]
#   y_post_pred[i, ] <- model_psd_r(
#     omega,
#     post$h0[j],
#     post$h_m1[j],
#     post$Kp[j],
#     post$Ki[j],
#     post$tp[j],
#     post$tau[j]
#   )
# }
# 
# 
# df_obs <- data.frame(freq = freq, psd = psd)
# 
# df_pred <- data.frame(
#   freq = rep(freq, ndraws),
#   psd = as.vector(t(y_post_pred)),
#   draw = rep(1:ndraws, each = length(freq))
# )
# df_pred_all <- data.frame(
#   freq = rep(freq, times = ndraws),
#   psd = as.vector(t(y_post_pred)),
#   draw = factor(rep(1:ndraws, each = length(freq)))
# )
# 
# 
# ggplot(df_pred_all, aes(x = freq, y = psd)) +
#   geom_line(aes(group=draw),alpha = 0.1, color = "blue") +  # faint blue lines for all posterior draws
#   geom_point(data = df_obs, aes(x = freq, y = psd), color = "black") +
#   scale_y_log10() +
#   scale_x_log10() +
#   labs(x = "Frequency (Hz)", y = "PSD (linear scale)",
#        title = "Posterior predictive lines and observed PSD") +
#   theme_minimal()
# 
# # Compute median and 95% CI for predictions at each freq
# library(dplyr)
# df_summary <- df_pred %>%
#   group_by(freq) %>%
#   summarize(
#     median = median(psd),
#     lower = quantile(psd, 0.025),
#     upper = quantile(psd, 0.975)
#   )
# 
# ggplot() +
#   geom_point(data = df_obs, aes(x = freq, y = psd), color = "black") +
#   geom_line(data = df_summary, aes(x = freq, y = median), color = "blue") +
#   geom_ribbon(data = df_summary, aes(x = freq, ymin = lower, ymax = upper),
#               alpha = 0.2, fill = "blue") +
#   scale_y_log10() +
#   scale_x_log10() +
#   labs(x = "Frequency (Hz)", y = "PSD (linear scale)",
#        title = "Posterior predictive vs observed PSD") +
#   theme_minimal()
# 
