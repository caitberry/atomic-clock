# ==============================================================================
# 04_plot_simulation_recovery.R
# Compare Posterior Estimates to "True" Simulation Values
# ==============================================================================

# --- CONFIGURATION ---
# Define the TRUE values used in 01_prep_data.R for simulation
# (Ensure these match exactly what you used to generate the data!)
# true_params <- list(h0=5.191944e-31, h_m1=1.617518e-33, Kp=1.736019e-02, Ki=1.279804e-01, tau=2.5)
true_params <- list(log_h0=log(5.580327e-31),log_h_m1= log(1.013525e-33), Kp= 2.110357e-02, Ki= 1.212464e-01, tau= 3.112662e+00) ##TRY THESE TUESDAY

# File path for the simulated fit
# (Adjust the date/name if you saved your simulation differently)
SIM_DATE <- "2026-01-14" 
dataName="ComplexModelSimDat"

FIT_FILE <- paste0(OUTPUT_FIT_DIR, "fit_",dataName,"_Date_", gsub("-", "_", SIM_DATE), ".rds")
PLOT_FILE <- paste0(OUTPUT_FIT_DIR, "../Plots/simulation_recovery.png")

# --- LOAD DATA ---
if(!file.exists(FIT_FILE)) stop(paste("Simulation fit file not found:", FIT_FILE))

message("Loading simulation fit...")
fit_sim <- readRDS(FIT_FILE)

# Extract samples
post_samples <- as.data.frame(fit_sim) %>%
  select(log_h0, log_h_m1, Kp, Ki, tau) %>%
  pivot_longer(cols = everything(), names_to = "Parameter", values_to = "Value")

# Create a dataframe for True Values to map to the facets
true_vals_df <- data.frame(
  Parameter = names(true_params),
  TrueValue = as.numeric(true_params)
)

# --- PLOTTING ---

p1=post_samples %>%
  # filter(Parameter %in% c("Kp", "Ki", "tau")) %>%
  ggplot(aes(x = Value)) +
  # Posterior Density
  geom_density(fill = "#1f78b4", alpha = 0.5, color = NA) +
  # True Value Line
  geom_vline(data = true_vals_df,#filter(true_vals_df, Parameter %in% c("Kp", "Ki", "tau")), 
             aes(xintercept = TrueValue), color = "red", linetype = "dashed", size = 1) +
  # Formatting
  facet_wrap(~Parameter, scales = "free", ncol = 3) +
  labs(y = "Posterior Density", x = "Parameter Value") +
  theme_bw() +
  theme(strip.background = element_rect(fill = "gray90"))
plot_filename <- paste0(OUTPUT_PLOT_DIR, "postVtruth_", dataName, "_", gsub("-", "_", target_date), ".png")
ggsave(plot_filename, plot = p1, width = 5, height = 4)


# --- PRINT NUMERICAL SUMMARY ---
# This prints a table to the console checking if Truth is within the 95% CI
summary_stats <- post_samples %>%
  group_by(Parameter) %>%
  summarize(
    Mean = mean(Value),
    Lower_95 = quantile(Value, 0.025),
    Upper_95 = quantile(Value, 0.975),
    .groups = 'drop'
  ) %>%
  left_join(true_vals_df, by = "Parameter") %>%
  mutate(
    Recovered = (TrueValue >= Lower_95 & TrueValue <= Upper_95)
  )

print(summary_stats)
