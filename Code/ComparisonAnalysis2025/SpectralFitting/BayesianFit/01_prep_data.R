# ==============================================================================
# 01_prep_data.R
# Organize data (Real or Simulated) and save to disk
# ==============================================================================
# source("Code/ComparisonAnalysis2025/SpectralFitting/BayesianFit/00_common.R")

# --- OPTION A: LOAD REAL DATA ---
# 
# dataName="AlYbApr15"
# files_list <- list.files(path = DATA_RAW_DIR, pattern = "spectral.*AlYb.*Apr15", full.names = TRUE)
# message(paste("Loading", length(files_list), "files..."))
# 
# all_data <- lapply(files_list, function(file) {
#   x <- readRDS(file)
#   extract_from_list(x, basename(file))
# })
# 
# spectralEstDF <- do.call(rbind, lapply(all_data, `[[`, "spectrum"))
# 
# spectralEstDF <- spectralEstDF %>%
#   mutate(
#     Ratio = str_extract(File, "(?<=spectralEstFor)[A-Za-z]+"),
#     Month = str_extract(File, "(?<=For[A-Za-z]{4})(\\d{2})"),
#     Day   = str_extract(File, "(?<=For[A-Za-z]{4}\\d{2})(\\d{2})"),
#     Date = as.Date(paste0("2025-", Month, "-", Day))
#   )
# 
# # Confidence intervals
# spectralEstDF <- spectralEstDF %>%
#   mutate(
#     lower = spectrum * qchisq(0.025, df = 2 * MY_K) / (2 * MY_K),
#     upper = spectrum * qchisq(0.975, df = 2 * MY_K) / (2 * MY_K),
#     sdish=sqrt((2 * spectrum^2)/MY_K)
#   )
# ##############plot for sanity
# ggplot() +
#   # 1. The Noisy Simulated Data (Black Points)
#   geom_point(data = filter(spectralEstDF,Date==unique(spectralEstDF$Date)[1], freq<F_MAX), aes(x = freq, y = spectrum),
#              color = "black", alpha = 0.5, size = 1)+
#   # Formatting
#   scale_x_log10() +
#   scale_y_log10() +
#   labs(
#     x = "Frequency (Hz)",
#     y = "Power Spectral Density"
#   ) +
#   theme_bw()


# 
# ############
# 
# 
# # Add Metadata & Clean
# spectralEstDF <- spectralEstDF %>%
#   mutate(
#     Ratio = str_extract(File, "(?<=spectralEstFor)[A-Za-z]+"),
#     Month = str_extract(File, "(?<=For[A-Za-z]{4})(\\d{2})"),
#     Day   = str_extract(File, "(?<=For[A-Za-z]{4}\\d{2})(\\d{2})"),
#     Date  = as.Date(paste0("2025-", Month, "-", Day)),
#     # Calculate approx SD for weighting
#     sdish = sqrt((1 * spectrum^2)/MY_K)
#   )

# # --- OPTION B: GENERATE SIMULATED DATA ---
# # filter(spectralEstDF,Date==unique(spectralEstDF$Date)[1], freq<F_MAX)
# # test=filter(spectralEstDF,Date==unique(spectralEstDF$Date)[1], freq<F_MAX)
# # table(diff(test$freq))
# # range(test$freq)
# # length(test$freq)
# # # posterior results from fit_AlYbApr15_Date_2025_02_27.rds
# # Parameter         Mean           SD
# # h0          h0 5.191944e-31 1.825889e-32
# # h_m1      h_m1 3.957758e-34 4.278109e-35
# # Kp          Kp 1.736019e-02 1.307403e-02
# # Ki          Ki 1.279804e-01 3.317706e-03
# # tau        tau 2.532619e+00 1.000327e-01
# # starting with these, might want something more different in the future
# #
dataName="ComplexModelSimDat"
f_sim <- seq(F_MIN,F_MAX,length.out=1000)#exp(seq(log(F_MIN), log(F_MAX), length.out=100))

#old # true_params <- list(h0=5.191944e-31, h_m1=1.617518e-33, Kp=1.736019e-02, Ki=1.279804e-01, tau=2.5)
true_params <- list(h0=5.580327e-31,h_m1= 1.013525e-33, Kp= 2.110357e-02, Ki= 1.212464e-01, tau= 3.112662e+00) ##TRY THESE TUESDAY
psd_true <- model_psd_r(f_sim, true_params$h0, true_params$h_m1,
                        true_params$Kp, true_params$Ki, true_params$tau, TP_VAL)
# Add multiplicative noise (Chi-squared like)
psd_noisy <- psd_true * rchisq(length(f_sim), df=2*MY_K) / (2*MY_K)

spectralEstDF <- data.frame(
  freq = f_sim, spectrum = psd_noisy,
  sdish = psd_noisy/sqrt(MY_K),
  Date = as.Date("2026-01-14"), File = "Simulated"
)

df_true <- data.frame(
  freq = f_sim,
  psd  = psd_true
)

ggplot() +
  # 1. The Noisy Simulated Data (Black Points)
  geom_point(data = spectralEstDF, aes(x = freq, y = spectrum),
             color = "black", alpha = 0.5, size = 1) +

  # 2. The True Physical Model (Red Line)
  geom_line(data = df_true, aes(x = freq, y = psd),
            color = "red", size = 1, linetype = "solid") +

  # Formatting
  scale_x_log10() +
  scale_y_log10() +
  labs(
    x = "Frequency (Hz)",
    y = "Power Spectral Density"
  ) +
  theme_bw()


# --- SAVE ---
saveRDS(spectralEstDF, file = DATA_PREP_FILE)
message(paste("Data prepared and saved to:", DATA_PREP_FILE))
