library(dplyr)
library(purrr)
library(broom)
library(ggplot2)

K <- 5

# Add confidence bounds
spectralEstDF <- spectralEstDF %>%
  mutate(
    lower = spectrum * qchisq(0.025, df = 2 * K) / (2 * K),
    upper = spectrum * qchisq(0.975, df = 2 * K) / (2 * K)
  )

# Filter to one date
plot_data <- filter(spectralEstDF, Date == "2025-01-16")

# Frequency cutoffs
upper_freqs <- seq(0.005, 0.1, by = 0.005)

# Fit loop with intercept, std.error, and fit stats
fit_results <- map_dfr(upper_freqs, function(f_upper) {
  fit_data <- filter(plot_data, freq > 0.0001, freq < f_upper) %>%
    mutate(weight = K / (2 * spectrum^2))
  
  if (nrow(fit_data) < 5) return(NULL)
  
  fit <- lm(spectrum ~ freq + I(freq^2), data = fit_data, weights = weight)
  
  # Extract model and summary info
  intercept_info <- tidy(fit) %>% filter(term == "(Intercept)")
  model_stats <- glance(fit)
  
  tibble(
    upper_freq = f_upper,
    intercept = intercept_info$estimate,
    std_error = intercept_info$std.error,
    r_squared = model_stats$r.squared,
    adj_r_squared = model_stats$adj.r.squared,
    sigma = model_stats$sigma  # residual standard error
  )
})


# Intercept with error band
p1 <- ggplot(fit_results, aes(x = upper_freq, y = intercept)) +
  geom_line(color = "blue") +
  geom_ribbon(aes(ymin = intercept - 2*std_error, ymax = intercept + 2*std_error), 
              fill = "blue", alpha = 0.2) +
  labs(title = "Intercept vs. Max Frequency", y = "Intercept (±2 SE)", x = "Max Frequency") +
  theme_minimal()

# R-squared
p2 <- ggplot(fit_results, aes(x = upper_freq, y = r_squared)) +
  geom_line(color = "darkgreen") +
  labs(title = "R-squared vs. Max Frequency", y = "R-squared", x = "Max Frequency") +
  theme_minimal()

p3 <- ggplot(fit_results, aes(x = upper_freq, y = sigma)) +
  geom_line(color = "darkgreen") +
  # labs(title = "R-squared vs. Max Frequency", y = "R-squared", x = "Max Frequency") +
  theme_minimal()

# Plot both side by side if desired
library(patchwork)
p1 + p2+p3

View(fit_results)


### looks like around .05 makes sense for this day,"2025-01-16", YbSr 

########################################################
########################################################
### all fitting playing around
########################################################
########################################################


K=5
spectralEstDF <- spectralEstDF %>%
  mutate(
    lower = spectrum * qchisq(0.025, df = 2 * K) / (2 * K),
    upper = spectrum * qchisq(0.975, df = 2 * K) / (2 * K)
  )

# Filter the data
plot_data <- filter(spectralEstDF, Date == "2025-01-16")
range(plot_data$freq)
fit_data=filter(plot_data,freq>.0001, freq<.05)

### I don't think it makes sense to include weights, too volitale. think about this more... 
# fit_data <- fit_data %>%
#   mutate(
#     weight_raw = K / (2 * spectrum^2),
#     weight = weight_raw / max(weight_raw)
#     # weight = weight_raw / mean(weight_raw),
#     # weight = pmax(weight, quantile(weight, 0.05))
#     )

fit <- lm(spectrum ~ freq+I(freq^2), data = fit_data)#, weights = weight)
summary(fit)

freq_seq <- seq(0, max(fit_data$freq), length.out = 200)
pred_df <- data.frame(freq = freq_seq)
pred_df$predicted_spectrum <- predict(fit, newdata = pred_df)


ggplot(plot_data, aes(x = freq, y = spectrum, col = factor(Date), ymin = lower, ymax = upper)) +
  geom_line(alpha = 0.4) +
  # geom_errorbar(width = 0, alpha = 0.4) +  # Add some width for visibility if you want
  geom_hline(yintercept = mean_val, linetype = "dashed", color = "black") +
  geom_line(data = pred_df, aes(x = freq, y = predicted_spectrum), inherit.aes = FALSE, color = "red", size = 1) +
  scale_y_log10() +
  scale_x_log10() +
  theme(legend.position = "bottom", legend.text = element_text(size = 10)) +
  guides(color = guide_legend(title = NULL))


ggplot(fit_data, aes(x = freq, y = weight)) +
  geom_point()




# Calculate Cook's distance
fit_data$cooksD <- cooks.distance(fit)

# Quick look at largest values
head(fit_data[order(-fit_data$cooksD), ], 10)

# Plot Cook's Distance vs Frequency
library(ggplot2)
ggplot(fit_data, aes(x = freq, y = cooksD)) +
  geom_point() +
  scale_x_log10() +
  labs(title = "Cook's Distance vs Frequency", y = "Cook's Distance") +
  geom_hline(yintercept = 4 / nrow(fit_data), linetype = "dashed", color = "red") +
  theme_minimal()




### testing log model for better residuals, but i don't like the fit
# fit_orig <- lm(spectrum ~ freq + I(freq^2), data = fit_data, weights = weight)
# fit_log <- lm(log(spectrum) ~ log(freq) + I(log(freq)^2), data = fit_data, weights = rep(K, nrow(fit_data)))
# 
# freq_seq <- seq(min(fit_data$freq), max(fit_data$freq), length.out = 200)
# pred_df <- data.frame(freq = freq_seq)
# 
# # Original scale prediction
# pred_df$fit_orig <- predict(fit_orig, newdata = pred_df)
# 
# # Log-scale prediction, back-transformed to original scale
# log_pred <- predict(fit_log, newdata = data.frame(
#   freq = freq_seq,
#   `log(freq)` = log(freq_seq)
# ))
# pred_df$fit_log_backtransformed <- exp(log_pred)
# 
# 
# ggplot(fit_data, aes(x = freq, y = spectrum)) +
#   geom_point(alpha = 0.3) +
#   geom_line(data = pred_df, aes(x = freq, y = fit_orig), color = "red", size = 1, linetype = "solid") +
#   geom_line(data = pred_df, aes(x = freq, y = fit_log_backtransformed), color = "blue", size = 1, linetype = "dashed") +
#   scale_x_log10() +
#   scale_y_log10() +
#   labs(title = "Comparison of Original vs Log-scale Model Fit",
#        subtitle = "Red = original scale fit | Blue dashed = log-log fit (back-transformed)",
#        y = "Spectrum",
#        x = "Frequency") +
#   theme_minimal()


library(ggplot2)

fit_data$residuals <- residuals(fit)
fit_data$fitted <- fitted(fit)

ggplot(fit_data, aes(x = freq, y = residuals)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Fitted values", y = "Residuals") +
  theme_minimal()
# scale_x_log10() +
# scale_y_log10() 

ggplot(fit_data, aes(x = fitted, y = residuals)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Fitted values", y = "Residuals") +
  theme_minimal()
# scale_x_log10() +
# scale_y_log10() 

ggplot(fit_data, aes(sample = residuals)) +
  stat_qq() +
  stat_qq_line() +
  labs(
    title = "Normal Q-Q Plot of Residuals"
  ) +
  theme_minimal()


