# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(Rmpfr) 

# 1. Read the raw text file
file_path <- "Data/Margolis2024DataMetrologia/ClockInputData2020_7.dat"
raw_lines <- readLines(file_path)

# 2. Extract the main measurement data
# Filter for lines that start with 'q' followed by digits
q_lines <- raw_lines[grepl("^q\\d+", raw_lines)]

# Parse these lines into a Data Frame
clock_data <- read.table(text = q_lines, 
                         sep = "", 
                         header = FALSE, 
                         stringsAsFactors = FALSE,
                         colClasses = "character", # FORCE all columns to be read as text!
                         col.names = c("Measurement_ID", "Parameter", "Value", "Uncertainty", "Reference"))

# 3. Extract the correlation coefficients
r_lines <- raw_lines[grepl("^r\\(", raw_lines)]

correlations <- read.table(text = r_lines, 
                           sep = "", 
                           header = FALSE, 
                           stringsAsFactors = FALSE,
                           col.names = c("Measurement_Pair", "Correlation_Coeff"))

# 4. Clean up the Reference Column
clock_data_parsed <- clock_data %>%
  # Remove the opening '[' and closing ']'
  mutate(Reference_clean = gsub("^\\[|\\]$", "", Reference)) %>%
  
  # Extract Author, Year, and Notes. 
  extract(
    Reference_clean, 
    into = c("Author", "Year", "Notes"), 
    regex = "^(.*?)(19\\d{2}|20\\d{2})(.*)$", 
    remove = FALSE
  ) %>%
  
  # Clean up and catch the edge cases
  mutate(
    Author = ifelse(is.na(Author), Reference_clean, Author),
    Notes = gsub("^[,\\-]", "", Notes), 
    Notes = na_if(Notes, ""),           
    Year = as.numeric(Year)             
  )

# Update the 12thCCL and SYRTETAIData references
clock_data_parsed <- clock_data_parsed %>%
  mutate(
    Year = ifelse(Author == "12thCCL", 2005, Year),
    Author = ifelse(grepl("SYRTETAIData", Reference), "SYRTETAIData", Author),
    Year = ifelse(grepl("SYRTETAIData", Reference), 2017, Year),
    Notes = ifelse(grepl("SYRTETAIData", Reference), "MJD55954-57867(Jan12-Apr17)", Notes)
  )

# 5. Reconstruct the Published Uncertainties
clock_data_published <- clock_data_parsed %>%
  mutate(
    Published_Uncertainty = case_when(
      # Wrap Value and Uncertainty in as.numeric() because they are now text strings!
      Measurement_ID == "q31" ~ as.numeric(Value) * 5.34e-16, 
      Measurement_ID == "q73" ~ as.numeric(Value) * 1.5e-16,  
      Measurement_ID == "q98" ~ as.numeric(Value) * 1.3e-16,  
      
      # Reverse the basic multipliers
      grepl("Uncertx3", Notes) ~ as.numeric(Uncertainty) / 3,
      grepl("Uncertx1.5", Notes) ~ as.numeric(Uncertainty) / 1.5,
      grepl("Uncertx2", Notes) ~ as.numeric(Uncertainty) / 2,
      grepl("Uncertx6", Notes) ~ as.numeric(Uncertainty) / 6,
      grepl("Uncertx100", Notes) ~ as.numeric(Uncertainty) / 100, 
      
      # Default condition
      TRUE ~ as.numeric(Uncertainty)
    )
  )


# bacon2025_data <- data.frame(
#   Measurement_ID = c("bacon25_1", "bacon25_2", "bacon25_3"),
#   Parameter = c("nu4_over_nu12", "nu4_over_nu8", "nu8_over_nu12"),
#   Value = c("2.6117014317814627101", "2.1628871275166636674", "1.2075070393433377230"),
#   Uncertainty = c("58e-19", "70e-19", "37e-19"), 
#   Reference = "[Aeppli2025]",
#   Reference_clean = "Aeppli2025",
#   Author = "BACON Collab",
#   Year = 2025,
#   Notes = NA_character_,
#   Published_Uncertainty = c(58e-19, 70e-19, 37e-19),
#   stringsAsFactors = FALSE
# )
# 
# clock_data_published <- bind_rows(clock_data_published, bacon2025_data)

# 6. Calculate High-Precision Fractional Offsets
prec <- 128

# Reference Frequencies (nu1 to nu14)
nu_0 <- mpfr(c(
  "1267402452901050.0", # 1: 115In+ 
  "1233030706593509.0", # 2: H (1S-2S two-photon transition / 2)
  "1128575290808154.4", # 3: 199Hg
  "1121015393207851.0", # 4: 27Al+
  "1064721609899145.0", # 5: 199Hg+
  "688358979309308.3",  # 6: 171Yb+ (E2 Quadrupole)
  "642121496772645.0",  # 7: 171Yb+ (E3 Octupole)
  "518295836590863.6",  # 8: 171Yb
  "455986240494140.0",  # 9: 40Ca
  "444779044095486.5",  # 10: 88Sr+
  "429228066418007.0",  # 11: 88Sr
  "429228004229873.0",  # 12: 87Sr
  "411042129776399.8",  # 13: 40Ca+
  "6834682610.9043126"  # 14: 87Rb (Microwave transition)
), precBits = prec)

# Assign names to the mpfr array so we can look them up by parameter
names(nu_0) <- paste0("nu", 1:14)

# Helper function to grab the right reference value (handling both absolutes and ratios)
# Bulletproof helper function using numeric indexing for mpfr objects
get_ref_value <- function(param) {
  if (grepl("_over_", param)) {
    # Split "nu3_over_nu12" into "nu3" and "nu12"
    parts <- str_split(param, "_over_")[[1]]
    
    # Extract just the numbers to use as numeric indices (e.g., 3 and 12)
    idx1 <- as.numeric(str_extract(parts[1], "\\d+"))
    idx2 <- as.numeric(str_extract(parts[2], "\\d+"))
    
    return(nu_0[idx1] / nu_0[idx2])
    
  } else {
    # Extract the number for absolute frequencies (e.g., "nu1" -> 1)
    idx <- as.numeric(str_extract(param, "\\d+"))
    return(nu_0[idx])
  }
}
# Calculate exact offsets using rowwise operations
clock_data_comparable <- clock_data_published %>%
  rowwise() %>%
  mutate(
    # Get the high-precision reference value
    Ref_Value = list(get_ref_value(Parameter)),
    
    # Convert measurement to mpfr to preserve precision during subtraction
    Value_mpfr = list(mpfr(as.character(Value), precBits = prec)),
    
    # Calculate Offset: (Measurement - Reference) / Reference
    Fractional_Offset = as.numeric((Value_mpfr[[1]] - Ref_Value[[1]]) / Ref_Value[[1]]),
    
    # Calculate Fractional Uncertainty using the true original values
    Fractional_Uncertainty = as.numeric(
      mpfr(as.character(Published_Uncertainty), precBits = prec) / Value_mpfr[[1]]
    ),
    
    # Calculate original/inflated uncertainty bounds for future plotting
    Offset_Min = Fractional_Offset - Fractional_Uncertainty,
    Offset_Max = Fractional_Offset + Fractional_Uncertainty,
    
    # Also track the 'inflated' committee fractional bounds just in case
    Committee_Frac_Unc = as.numeric(
      mpfr(as.character(Uncertainty), precBits = prec) / Value_mpfr[[1]]
    ),
    Committee_Offset_Min = Fractional_Offset - Committee_Frac_Unc,
    Committee_Offset_Max = Fractional_Offset + Committee_Frac_Unc
  ) %>%
  ungroup()

# Preview the clean, comparable dataset
clock_data_comparable %>% 
  select(Measurement_ID, Parameter, Fractional_Offset, Fractional_Uncertainty) %>% 
  mutate(Fractional_Offset = sprintf("%e", Fractional_Offset)) 





# --------------------------------------------------------
# Now, you can proceed with Step 6 (Calculate High-Precision Fractional Offsets)
# and your plotting code exactly as they were!
### Plots!

library(ggplot2)
library(dplyr)
library(stringr)

# Define your exact list of Table 8 IDs
table_8_ids <- c("q1", "q31", "q51", "q52", "q73", "q74", "q78", "q88", "q98", "q105")

# 1. Categorize data, drop missing years, fix Factor levels, AND add the Table 8 Flag!
plot_data <- clock_data_comparable %>%
  filter(!is.na(Year)) %>%
  mutate(
    Fractional_Offset = as.numeric(Fractional_Offset),
    Offset_Min = as.numeric(Offset_Min),
    Offset_Max = as.numeric(Offset_Max),
    Committee_Offset_Min = as.numeric(Committee_Offset_Min),
    Committee_Offset_Max = as.numeric(Committee_Offset_Max),
    
    Measurement_Type = ifelse(grepl("_over_", Parameter), 
                              "Frequency Ratios", 
                              "Absolute Frequencies"),
    
    Parameter = factor(Parameter, levels = str_sort(unique(Parameter), numeric = TRUE)),
    
    # Check if the Measurement_ID is inside your list
    Table_8_Flag = ifelse(Measurement_ID %in% table_8_ids, 
                          "Adjusted", 
                          "Published")
  )

# Split the data into two separate dataframes
data_abs <- plot_data %>% filter(Measurement_Type == "Absolute Frequencies")
data_ratio <- plot_data %>% filter(Measurement_Type == "Frequency Ratios")

# custom color palette
high_contrast_colors <- c(
  "#e6194b", "#3cb44b", "#ffe119", "#4363d8", "#f58231", 
  "#911eb4", "#46f0f0", "#f032e6", "#bcf60c", "#fabebe", 
  "#008080", "#e6beff", "#9a6324", "#ede0a3", "#800000", 
  "#aaffc3", "#808000", "#ffd8b1", "#000075", "#808080"
)

plot_abs <- ggplot(data_abs, 
                   aes(x = factor(Year), 
                       y = Fractional_Offset * 1e15, 
                       fill = Parameter,
                       color = Parameter,
                       shape = Table_8_Flag, # Maps the shape!
                       group = Measurement_ID)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.8) +
  
  geom_errorbar(aes(ymin = Committee_Offset_Min * 1e15, ymax = Committee_Offset_Max * 1e15), 
                width = 0, alpha = 0.3, linewidth = 0.8, 
                position = position_dodge(width = 0.7)) +
  geom_errorbar(aes(ymin = Offset_Min * 1e15, ymax = Offset_Max * 1e15), 
                width = 0, alpha = 1, linewidth = 0.8, 
                position = position_dodge(width = 0.7)) +
  geom_point(color = "black", size = 3, stroke = 0.8, 
             position = position_dodge(width = 0.7)) +
  scale_fill_manual(values = high_contrast_colors) +
  scale_color_manual(values = high_contrast_colors) +
  scale_shape_manual(values = c("Published" = 21, 
                                "Adjusted" = 24)) +
  
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"),
    legend.position = "right",
    panel.grid.minor.x = element_blank(),
    
    # Make the legend symbols a bit larger so the colors are easier to see
    legend.key.size = unit(1.5, "lines") 
  ) +
  
  # Ensure the legend title applies to both the fill and the color
  labs(
    title = "Absolute Frequencies vs. 2017 References",
    subtitle = "Thick faded = Committee Uncertainty | Thin solid = Original Published",
    x = "Publication Year",
    y = expression("Fractional Offset (" %*% 10^{-15} ~ ")"),
    fill = "Absolute\nTransition",
    color = "Absolute\nTransition"
  )+
  # Add this to the bottom of BOTH of your ggplot codes!
  guides(
    # Force the Parameter legend to use a fillable shape (21) with a black outline
    fill = guide_legend(override.aes = list(shape = 21, color = "black")),
    color = guide_legend(override.aes = list(shape = 21, color = "black")),
    
    # Fill the 'Table 8' shape legend with grey so the symbols aren't empty/invisible
    shape = guide_legend(override.aes = list(fill = "grey50", size = 4))
  )+
  coord_cartesian(ylim = c(-15,15))

print(plot_abs)

# --- PLOT 2: FREQUENCY RATIOS ---
plot_ratio <- ggplot(data_ratio, 
                     aes(x = factor(Year), 
                         y = Fractional_Offset * 1e15, 
                         fill = Parameter,
                         color = Parameter,
                         shape = Table_8_Flag, # Maps the shape!
                         group = Measurement_ID)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.8) +
  
  geom_errorbar(aes(ymin = Committee_Offset_Min * 1e15, ymax = Committee_Offset_Max * 1e15), 
                width = 0, alpha = 0.3, linewidth = .8, 
                position = position_dodge(width = 0.7)) +
  geom_errorbar(aes(ymin = Offset_Min * 1e15, ymax = Offset_Max * 1e15), 
                width = 0, alpha = 1, linewidth = 0.8, 
                position = position_dodge(width = 0.7)) +
  geom_point(color = "black", size = 3, stroke = 0.8, 
             position = position_dodge(width = 0.7)) +
  scale_fill_manual(values = high_contrast_colors) +
  scale_color_manual(values = high_contrast_colors) +
  scale_shape_manual(values = c("Published" = 21, 
                                "Adjusted" = 24)) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"),
    legend.position = "right",
    panel.grid.minor.x = element_blank(),
    
    # Make the legend symbols a bit larger so the colors are easier to see
    legend.key.size = unit(1.5, "lines") 
  ) +
  
  # Ensure the legend title applies to both the fill and the color
  labs(
    title = "Frequency Ratios vs. 2017 References",
    subtitle = "Thick faded = Committee Uncertainty | Thin solid = Original Published",
    x = "Publication Year",
    y = expression("Fractional Offset (" %*% 10^{-15} ~ ")"),
    fill = "Frequency\nRatio",
    color = "Frequency\nRatio"
  )+
  # Add this to the bottom of BOTH of your ggplot codes!
  guides(
    # Force the Parameter legend to use a fillable shape (21) with a black outline
    fill = guide_legend(override.aes = list(shape = 21, color = "black")),
    color = guide_legend(override.aes = list(shape = 21, color = "black")),
    
    # Fill the 'Table 8' shape legend with grey so the symbols aren't empty/invisible
    shape = guide_legend(override.aes = list(fill = "grey50", size = 4))
  )
  # coord_cartesian(ylim = c(7.17,7.3))#zoom in on BACON


print(plot_ratio)


###organizing by nu

# --- PLOT 1: ABSOLUTE FREQUENCIES (Grouped by Nu Parameter) ---
plot_abs <- ggplot(data_abs, 
                   aes(x = Parameter, # <-- CHANGED TO PARAMETER
                       y = Fractional_Offset * 1e15, 
                       fill = Parameter,
                       color = Parameter,
                       shape = Table_8_Flag, 
                       group = Measurement_ID)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.8) +
  
  geom_errorbar(aes(ymin = Committee_Offset_Min * 1e15, ymax = Committee_Offset_Max * 1e15), 
                width = 0, alpha = 0.3, linewidth = 0.8, 
                position = position_dodge(width = 0.7)) +
  geom_errorbar(aes(ymin = Offset_Min * 1e15, ymax = Offset_Max * 1e15), 
                width = 0, alpha = 1, linewidth = 0.8, 
                position = position_dodge(width = 0.7)) +
  geom_point(color = "black", size = 3, stroke = 0.8, 
             position = position_dodge(width = 0.7)) +
  
  scale_fill_manual(values = high_contrast_colors) +
  scale_color_manual(values = high_contrast_colors) +
  scale_shape_manual(values = c("Published" = 21, 
                                "Adjusted" = 24)) +
  
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"),
    legend.position = "none", # <-- PRO TIP: We can hide the legend here since the x-axis IS the legend now!
    panel.grid.minor.x = element_blank()
  ) +
  
  labs(
    title = "Absolute Frequencies vs. 2017 References",
    subtitle = "Grouped by Transition | Triangles = Table 8 (Adjusted) | Circles = Published",
    x = "Absolute Transition (nu)",
    y = expression("Fractional Offset (" %*% 10^{-15} ~ ")")
  ) +
   coord_cartesian(ylim = c(-20,20))

print(plot_abs)


# --- PLOT 2: FREQUENCY RATIOS (Grouped by Ratio) ---
ggplot(data_ratio, 
                     aes(x = Parameter, # <-- CHANGED TO PARAMETER
                         y = Fractional_Offset * 1e15, 
                         fill = Parameter,
                         color = Parameter,
                         shape = Table_8_Flag, 
                         group = Measurement_ID)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.8) +
  
  geom_errorbar(aes(ymin = Committee_Offset_Min * 1e15, ymax = Committee_Offset_Max * 1e15), 
                width = 0, alpha = 0.3, linewidth = .8, 
                position = position_dodge(width = 0.7)) +
  geom_errorbar(aes(ymin = Offset_Min * 1e15, ymax = Offset_Max * 1e15), 
                width = 0, alpha = 1, linewidth = 0.8, 
                position = position_dodge(width = 0.7)) +
  geom_point(color = "black", size = 2, stroke = 0.8, 
             position = position_dodge(width = 0.7)) +
  
  scale_fill_manual(values = high_contrast_colors) +
  scale_color_manual(values = high_contrast_colors) +
  scale_shape_manual(values = c("Published" = 21, 
                                "Adjusted" = 24)) +
  
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"),
    legend.position = "right", 
    panel.grid.minor.x = element_blank(),
    legend.key.size = unit(1.5, "lines") 
  ) +
  
  labs(
    title = "Frequency Ratios vs. 2017 References",
    subtitle = "Grouped by Ratio | Triangles = Table 8 (Adjusted) | Circles = Published",
    x = "Frequency Ratio",
    y = expression("Fractional Offset (" %*% 10^{-15} ~ ")"),
    shape = "Data Source"
  ) +
  
  # Note: I removed the color/fill from the legend since they are on the x-axis now,
  # keeping only the shape legend!
  guides(
    fill = "none",
    color = "none",
    shape = guide_legend(override.aes = list(fill = "grey50", size = 4))
  )
  # coord_cartesian(ylim = c(-1,1))






