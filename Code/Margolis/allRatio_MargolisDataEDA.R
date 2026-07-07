library(Rmpfr)
library(rstan)
library(ggplot2)

# =====================================================================
# 1. GLOBAL SETTINGS & MAPPINGS
# =====================================================================
prec <- 128
scale_factor <- mpfr(1e15, precBits = prec)

# 2017 CIPM Recommended Standard Frequencies
nu_0 <- mpfr(c(
  "1267402452901050.0", # 1: 115In+ 
  "1233030706593509.0", # 2: H
  "1128575290808154.4", # 3: 199Hg
  "1121015393207851.0", # 4: 27Al+
  "1064721609899145.0", # 5: 199Hg+
  "688358979309308.3",  # 6: 171Yb+ (E2)
  "642121496772645.0",  # 7: 171Yb+ (E3)
  "518295836590863.6",  # 8: 171Yb
  "455986240494140.0",  # 9: 40Ca
  "444779044095486.5",  # 10: 88Sr+
  "429228066418007.0",  # 11: 88Sr
  "429228004229873.0",  # 12: 87Sr
  "411042129776399.8",  # 13: 40Ca+
  "6834682610.9043126"  # 14: 87Rb
), precBits = prec)

# Mapping indices to Clock Names
clock_names <- c(
  "nu1"  = "115In+",      "nu2"  = "H",           "nu3"  = "199Hg",
  "nu4"  = "27Al+",       "nu5"  = "199Hg+",      "nu6"  = "171Yb+ (E2)",
  "nu7"  = "171Yb+ (E3)", "nu8"  = "171Yb",       "nu9"  = "40Ca",
  "nu10" = "88Sr+",       "nu11" = "88Sr",        "nu12" = "87Sr",
  "nu13" = "40Ca+",       "nu14" = "87Rb"
)

# Mapping Clock Names to Noise Groups
group_map <- c(
  "115In+"      = "Ion Clocks",
  "H"           = "Legacy/Sparse Clocks",
  "199Hg"       = "Lattice Clocks",
  "27Al+"       = "Ion Clocks",
  "199Hg+"      = "Ion Clocks",
  "171Yb+ (E2)" = "Ion Clocks",
  "171Yb+ (E3)" = "Ion Clocks",
  "171Yb"       = "Lattice Clocks",
  "40Ca"        = "Legacy/Sparse Clocks",
  "88Sr+"       = "Ion Clocks",
  "88Sr"        = "Lattice Clocks",
  "87Sr"        = "Lattice Clocks",
  "40Ca+"       = "Ion Clocks",
  "87Rb"        = "Microwave Fountains"
)

# =====================================================================
# 2. DATA LOADING & PARSING
# =====================================================================
raw_data <- readLines("Data/Margolis2024DataMetrologia/ClockInputData2020_7.dat")

# Extract and clean measurement lines
q_lines <- raw_data[grep("^q\\d+", raw_data)]
parsed_list <- lapply(q_lines, function(line) strsplit(trimws(line), "\\s+")[[1]][1:4])
parsed_q <- do.call(rbind, parsed_list)

# Safely filter out q1 and q2 (instead of using hardcoded row indices)
parsed_q <- parsed_q[!parsed_q[,1] %in% c("q1", "q2"), ]

q_labels <- parsed_q[, 1]
q_types  <- parsed_q[, 2]
Q_mpfr   <- mpfr(parsed_q[, 3], precBits = prec)
U_mpfr   <- mpfr(parsed_q[, 4], precBits = prec)

N_meas <- length(Q_mpfr)
N_clocks <- length(nu_0)

# =====================================================================
# 3. MATHEMATICAL SCALING & IDENTIFICATION
# =====================================================================
ref_vals <- mpfr(rep(0, N_meas), precBits = prec)
is_ratio <- rep(0, N_meas)
clock_1  <- rep(0, N_meas)
clock_2  <- rep(0, N_meas)

for (i in 1:N_meas) {
  type_str <- q_types[i]
  if (grepl("_over_", type_str)) {
    is_ratio[i] <- 1
    parts <- strsplit(type_str, "_over_")[[1]]
    clock_1[i] <- as.numeric(gsub("nu", "", parts[1]))
    clock_2[i] <- as.numeric(gsub("nu", "", parts[2]))
    ref_vals[i] <- nu_0[clock_1[i]] / nu_0[clock_2[i]]
  } else {
    is_ratio[i] <- 0
    clock_1[i] <- as.numeric(gsub("nu", "", type_str))
    clock_2[i] <- 0
    ref_vals[i] <- nu_0[clock_1[i]]
  }
}

# Final scaled values for Stan
Y_stan <- as.numeric(((Q_mpfr - ref_vals) / ref_vals) * scale_factor)
u_scaled_stan <- as.numeric((U_mpfr / ref_vals) * scale_factor)
V_stan <- diag(u_scaled_stan^2)

# =====================================================================
# 4. DATAFRAME CONSTRUCTION (Master DataFrames)
# =====================================================================

# --- 4A. Master Measurement DataFrame ---
measurement_df <- data.frame(
  Q_ID = factor(q_labels, levels = q_labels[order(as.numeric(gsub("q", "", q_labels)))]),
  Y_Offset_Scaled = Y_stan,
  U_Unc_Scaled = u_scaled_stan,
  Nu_Type = q_types,
  Measurement_Type = ifelse(is_ratio == 1, "Frequency Ratios", "Absolute Frequencies"),
  Unc_Tier = factor(ifelse(u_scaled_stan > 5.0, "High Uncertainty", "Lower Uncertainty"), 
                    levels = c("Lower Uncertainty", "High Uncertainty"))
)

# Populate Transition Names & Noise Groups
measurement_df$Transition_Name <- sapply(1:N_meas, function(i) {
  if (is_ratio[i] == 1) {
    paste0(clock_names[paste0("nu", clock_1[i])], " / ", clock_names[paste0("nu", clock_2[i])])
  } else {
    clock_names[paste0("nu", clock_1[i])]
  }
})

measurement_df$Noise_Group <- sapply(1:N_meas, function(i) {
  if (is_ratio[i] == 1) {
    g1 <- group_map[ clock_names[paste0("nu", clock_1[i])] ]
    g2 <- group_map[ clock_names[paste0("nu", clock_2[i])] ]
    paste0(sort(c(g1, g2)), collapse = " / ")
  } else {
    group_map[ clock_names[paste0("nu", clock_1[i])] ]
  }
})

# --- 4B. Clock Participation DataFrame (Long Format) ---
clock_info_df <- data.frame(Clock_ID = character(), Meas_Type = character(), 
                            Uncertainty = numeric(), stringsAsFactors = FALSE)

for(i in 1:nrow(measurement_df)) {
  if(is_ratio[i] == 1) {
    clock_info_df <- rbind(clock_info_df, 
                           data.frame(Clock_ID = paste0("nu", clock_1[i]), Meas_Type = "Ratio", Uncertainty = u_scaled_stan[i]),
                           data.frame(Clock_ID = paste0("nu", clock_2[i]), Meas_Type = "Ratio", Uncertainty = u_scaled_stan[i]))
  } else {
    clock_info_df <- rbind(clock_info_df, 
                           data.frame(Clock_ID = paste0("nu", clock_1[i]), Meas_Type = "Absolute", Uncertainty = u_scaled_stan[i]))
  }
}

# Apply mappings to long dataframe
clock_info_df$Clock_ID <- factor(clock_names[clock_info_df$Clock_ID], levels = clock_names)
clock_info_df$Noise_Group <- factor(group_map[as.character(clock_info_df$Clock_ID)], 
                                    levels = c("Lattice Clocks", "Ion Clocks", "Microwave Fountains", "Legacy/Sparse Clocks"))


# =====================================================================
# 5. DATA VISUALIZATION
# =====================================================================
pd <- position_jitter(width = 0.15, height = 0, seed = 42)

# --- Plot 1: Input Data Offsets (Faceted by Noise Group Pairs) ---
ggplot(filter(measurement_df,Measurement_Type=="Absolute Frequencies"), aes(x = Transition_Name, y = Y_Offset_Scaled, color = Noise_Group)) +
  geom_point(size = 2, position = pd, alpha = 0.8) +
  geom_errorbar(aes(ymin = Y_Offset_Scaled - U_Unc_Scaled, ymax = Y_Offset_Scaled + U_Unc_Scaled), 
                width = 0.2, position = pd, alpha = 0.8) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 9, face = "bold"),
        legend.position = "right", strip.text = element_text(face = "bold", size = 11)) +
  labs(title = "Input Data Offsets Filtered by Noise Group",
       subtitle = "Ratio measurements are categorized by both clocks involved",
       x = "Transition / Clock Name", y = expression("Fractional Offset " ~ (Y %*% 10^{15})), color = "Measurement Type") 
  # facet_wrap(~ Noise_Group, scales = "free_x", ncol = 2)

# --- Plot 2: Data Volume by Clock ---
info_count_plot <- ggplot(clock_info_df, aes(x = Clock_ID, fill = Meas_Type)) +
  geom_bar(color = "black", alpha = 0.8) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold")) +
  labs(title = "Data Volume: Number of Measurements per Clock",
       subtitle = "Ratios count toward both participating clocks",
       x = "Clock Type", y = "Total Number of Measurements", fill = "Measurement Type") +
  scale_fill_manual(values = c("Absolute" = "#4e79a7", "Ratio" = "#f28e2b"))

# --- Plot 3: Network Density by Noise Group ---
group_vol_plot <- ggplot(clock_info_df, aes(x = Clock_ID, fill = Meas_Type)) +
  geom_bar(color = "black", alpha = 0.8) +
  facet_wrap(~ Noise_Group, scales = "free_x", nrow = 1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
        strip.text = element_text(face = "bold", size = 11), panel.spacing = unit(1, "lines")) +
  labs(title = "Network Density: Data Volume by Group",
       subtitle = "Visualizing the 'data-rich' vs. 'sparse' sectors of the clock network",
       x = "Clock Type", y = "Total Number of Measurements", fill = "Measurement Type") +
  scale_fill_manual(values = c("Absolute" = "#4e79a7", "Ratio" = "#f28e2b"))

# --- Plot 4: Noise Profiling by Group ---
group_unc_plot <- ggplot(clock_info_df, aes(x = Noise_Group, y = Uncertainty)) +
  geom_boxplot(outlier.shape = NA, color = "gray40", fill = "gray90", width = 0.5) +
  geom_jitter(aes(color = Meas_Type), position = pd, size = 2.5, alpha = 0.8) +
  scale_y_log10() +
  theme_bw() +
  theme(axis.text.x = element_text(size = 11, face = "bold"), legend.position = "right") +
  labs(title = "Noise Profiling: The Hybrid 4-Group Model",
       subtitle = "Comparing the spread of uncertainties across fundamentally different architectures",
       x = "Assigned Noise Group", y = "Scaled Uncertainty (Log10 Scale)", color = "Measurement Type") +
  scale_color_manual(values = c("Absolute" = "#4e79a7", "Ratio" = "#f28e2b"))

# Print the visualizations
print(info_count_plot)
print(group_vol_plot)
print(group_unc_plot)