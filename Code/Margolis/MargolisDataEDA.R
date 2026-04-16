# install.packages("Rmpfr")
library(Rmpfr)
library(rstan)

# Set precision bits (128 bits gives ~38 decimal digits of precision, which is plenty)
prec <- 128

# 1. Define Nu_0 using STRINGS to prevent early truncation
nu_0 <- mpfr(c(
  "1128575290808154.4", # 1: Hg (nu3)
  "518295836590863.6",  # 2: Yb (nu8)
  "429228004229873.0"   # 3: Sr (nu12)
), precBits = prec)


# #2, matches figure
# q6	nu3		1128575290808155.1	6.4	[McFerran2012,corrected2015]
# q7	nu3		1128575290808154.62	0.41	[Tyumenev2016]
# 
# #10, matches figure
# q20	nu8		518295836590864.0	28.0	[Kohno2009]
# q21	nu8		518295836590863.1	2.0	[Yasuda2012]
# q22	nu8		518295836590865.2	0.7	[Lemke2009]
# q23	nu8		518295836590863.5	8.1	[Park2013]
# q24	nu8		518295836590863.59	0.31	[Pizzocaro2017]
# q25	nu8		518295836590863.38	0.57	[Kim2017]
# q70	nu8		518295836590863.30	0.378	[Luo2020]
# q75	nu8		518295836590863.71	0.11	[McGrew2019]
# q76	nu8		518295836590863.61	0.13	[Pizzocaro2020]
# q89	nu8		518295836590863.54	0.259	[Kobayashi2020]
# 
# #21, 20 in figure
# q36	nu12		429228004229874.0	1.1	[Boyd2007]
# q37	nu12		429228004229873.65	0.37	[Campbell2008]
# q38	nu12		429228004229873.6	1.1	[Baillard2008]
# q39	nu12		429228004229874.1	2.4	[Hong2009]
# q40	nu12		429228004229872.9	0.5	[Falke2011]
# q41	nu12		429228004229873.9	1.4	[Yamaguchi2012]
# q42	nu12		429228004229872.0	1.6	[Akamatsu2014b]
# q43	nu12		429228004229873.56	0.49	[Tanabe2015]
# q44	nu12		429228004229873.7	1.4	[Lin2015]
# q45	nu12		429228004229873.13	0.17	[Falke2014]
# q46	nu12		429228004229873.10	0.13	[LeTargat2013]
# q47	nu12		429228004229872.92	0.12	[Lodewyck2016]	
# q48	nu12		429228004229872.97	0.16	[Grebing2016(Oct14)]
# q49	nu12		429228004229873.04	0.11	[Grebing2016(Jun15)]
# q50	nu12		429228004229872.97	0.40	[Hachisu2017]
# q51	nu12		429228004229872.99	18.0	[Hachisu2017b,Uncertx100Exclude]
# q72	nu12		429228004229873.1	0.429	[Hobson2020]
# q73	nu12		429228004229873.00	0.0708	[Schwarz2020-1.65E-16unc]
# q90	nu12		429228004229873.082	0.0773	[Nemitz2020]
# q91	nu12		429228004229873.13	0.4	[Grotti2018]
# q96	nu12	429228004229873.19	1.50E-01	[Leopardi2020]
# 
# #2, matches figure
# q59	nu3_over_nu12	2.62931420989890960	2.2e-16	[Yamanake2015]
# q60	nu3_over_nu12	2.62931420989890915	4.6e-16	[Tyumenev2016]
# 
# #10, matches figure
# q64	nu8_over_nu12	1.2075070393433412	1.7e-15	[Akamatsu2014aErrata]
# q65	nu8_over_nu12	1.20750703934333776	2.9e-16	[Takamoto2015]
# q66	nu8_over_nu12	1.207507039343337749	5.5e-17	[Nemitz2016]
# q80	nu8_over_nu12	1.20750703934333790	7.00E-16	[Fujieda2018]
# q81	nu8_over_nu12	1.20750703934333841	3.40E-16	[Grotti2018]
# q82	nu8_over_nu12	1.20750703934333805	3.38E-16	[Pizzocaro2020b]
# q83	nu8_over_nu12	1.20750703934333738	3.02E-16	[Pizzocaro2020b]
# q93	nu8_over_nu12	1.20750703934333858	4.95E-16	[Hisai2020]
# q94	nu8_over_nu12	1.20750703934333782	7.49E-16	[Kobayashi2020]
# q102	nu8_over_nu12	1.2075070393433378482	8.20E-18	[Beloy2021]
# 
# #1, matches figure
# q79	nu3_over_nu8	2.17747319413456507	1.92E-16	[Ohmae2020]
# 

# 2. Define Q using STRINGS (All 46 Measurements)
Q <- mpfr(c(
  # --- Absolute Hg (nu3) [2 measurements] ---
  "1128575290808155.1",      # 1. q6   [McFerran2012,corrected2015]
  "1128575290808154.62",     # 2. q7   [Tyumenev2016]
  
  # --- Absolute Yb (nu8) [10 measurements] ---
  "518295836590864.0",       # 3. q20  [Kohno2009]
  "518295836590863.1",       # 4. q21  [Yasuda2012]
  "518295836590865.2",       # 5. q22  [Lemke2009]
  "518295836590863.5",       # 6. q23  [Park2013]
  "518295836590863.59",      # 7. q24  [Pizzocaro2017]
  "518295836590863.38",      # 8. q25  [Kim2017]
  "518295836590863.30",      # 9. q70  [Luo2020]
  "518295836590863.71",      # 10. q75 [McGrew2019]
  "518295836590863.61",      # 11. q76 [Pizzocaro2020]
  "518295836590863.54",      # 12. q89 [Kobayashi2020]
  
  # --- Absolute Sr (nu12) [21 measurements] ---
  "429228004229874.0",       # 13. q36 [Boyd2007]
  "429228004229873.65",      # 14. q37 [Campbell2008]
  "429228004229873.6",       # 15. q38 [Baillard2008]
  "429228004229874.1",       # 16. q39 [Hong2009]
  "429228004229872.9",       # 17. q40 [Falke2011]
  "429228004229873.9",       # 18. q41 [Yamaguchi2012]
  "429228004229872.0",       # 19. q42 [Akamatsu2014b]
  "429228004229873.56",      # 20. q43 [Tanabe2015]
  "429228004229873.7",       # 21. q44 [Lin2015]
  "429228004229873.13",      # 22. q45 [Falke2014]
  "429228004229873.10",      # 23. q46 [LeTargat2013]
  "429228004229872.92",      # 24. q47 [Lodewyck2016]
  "429228004229872.97",      # 25. q48 [Grebing2016(Oct14)]
  "429228004229873.04",      # 26. q49 [Grebing2016(Jun15)]
  "429228004229872.97",      # 27. q50 [Hachisu2017]
  "429228004229872.99",      # 28. q51 [Hachisu2017b]
  "429228004229873.1",       # 29. q72 [Hobson2020]
  "429228004229873.00",      # 30. q73 [Schwarz2020]
  "429228004229873.082",     # 31. q90 [Nemitz2020]
  "429228004229873.13",      # 32. q91 [Grotti2018]
  "429228004229873.19",      # 33. q96 [Leopardi2020]
  
  # --- Ratio Hg / Sr (nu3 / nu12) [2 measurements] ---
  "2.62931420989890960",     # 34. q59 [Yamanake2015]
  "2.62931420989890915",     # 35. q60 [Tyumenev2016]
  
  # --- Ratio Yb / Sr (nu8 / nu12) [10 measurements] ---
  "1.2075070393433412",      # 36. q64 [Akamatsu2014aErrata]
  "1.20750703934333776",     # 37. q65 [Takamoto2015]
  "1.207507039343337749",    # 38. q66 [Nemitz2016]
  "1.20750703934333790",     # 39. q80 [Fujieda2018]
  "1.20750703934333841",     # 40. q81 [Grotti2018]
  "1.20750703934333805",     # 41. q82 [Pizzocaro2020b]
  "1.20750703934333738",     # 42. q83 [Pizzocaro2020b]
  "1.20750703934333858",     # 43. q93 [Hisai2020]
  "1.20750703934333782",     # 44. q94 [Kobayashi2020]
  "1.2075070393433378482",   # 45. q102[Beloy2021]
  
  # --- Ratio Hg / Yb (nu3 / nu8) [1 measurement] ---
  "2.17747319413456507"      # 46. q79 [Ohmae2020]
), precBits = prec)


# 3. Define Uncertainties using STRINGS (All 46 Measurements)
U <- mpfr(c(
  # --- Absolute Hg (nu3) ---
  "6.4",         # 1. q6
  "0.41",        # 2. q7
  
  # --- Absolute Yb (nu8) ---
  "28.0",        # 3. q20
  "2.0",         # 4. q21
  "0.7",         # 5. q22
  "8.1",         # 6. q23
  "0.31",        # 7. q24
  "0.57",        # 8. q25
  "0.378",       # 9. q70
  "0.11",        # 10. q75
  "0.13",        # 11. q76
  "0.259",       # 12. q89
  
  # --- Absolute Sr (nu12) ---
  "1.1",         # 13. q36
  "0.37",        # 14. q37
  "1.1",         # 15. q38
  "2.4",         # 16. q39
  "0.5",         # 17. q40
  "1.4",         # 18. q41
  "1.6",         # 19. q42
  "0.49",        # 20. q43
  "1.4",         # 21. q44
  "0.17",        # 22. q45
  "0.13",        # 23. q46
  "0.12",        # 24. q47
  "0.16",        # 25. q48
  "0.11",        # 26. q49
  "0.40",        # 27. q50
  "18.0",        # 28. q51 
  "0.429",       # 29. q72
  "0.0708",      # 30. q73
  "0.0773",      # 31. q90
  "0.4",         # 32. q91
  "1.50e-01",    # 33. q96 
  
  # --- Ratio Hg / Sr (nu3 / nu12) ---
  "2.2e-16",     # 34. q59
  "4.6e-16",     # 35. q60
  
  # --- Ratio Yb / Sr (nu8 / nu12) ---
  "1.7e-15",     # 36. q64
  "2.9e-16",     # 37. q65
  "5.5e-17",     # 38. q66
  "7.00e-16",    # 39. q80
  "3.40e-16",    # 40. q81
  "3.38e-16",    # 41. q82
  "3.02e-16",    # 42. q83
  "4.95e-16",    # 43. q93
  "7.49e-16",    # 44. q94
  "8.20e-18",    # 45. q102
  
  # --- Ratio Hg / Yb (nu3 / nu8) ---
  "1.92e-16"     # 46. q79
), precBits = prec)

# Map the reference values to match the order of the 46 measurements
# Group 1: 2 Hg, Group 2: 10 Yb, Group 3: 21 Sr, 
# Group 4: 2 Hg/Sr, Group 5: 10 Yb/Sr, Group 6: 1 Hg/Yb
Ref_Vals <- c(
  rep(nu_0[1], 2),               # Abs Hg
  rep(nu_0[2], 10),              # Abs Yb
  rep(nu_0[3], 21),              # Abs Sr
  rep(nu_0[1] / nu_0[3], 2),     # Ratio Hg / Sr
  rep(nu_0[2] / nu_0[3], 10),    # Ratio Yb / Sr
  rep(nu_0[1] / nu_0[2], 1)      # Ratio Hg / Yb
)

scale_factor <- mpfr(1e15, precBits = prec)

# Now you can calculate all 46 offsets and uncertainties in one line
Y_mpfr <- ((Q - Ref_Vals) / Ref_Vals) * scale_factor
u_scaled_mpfr <- (U / Ref_Vals) * scale_factor


# Define the metadata vectors matching the exact order of your 46 measurements
q_labels <- c(
  "q6", "q7",                                                                          # Abs Hg
  "q20", "q21", "q22", "q23", "q24", "q25", "q70", "q75", "q76", "q89",                # Abs Yb
  "q36", "q37", "q38", "q39", "q40", "q41", "q42", "q43", "q44", "q45", "q46",         # Abs Sr (1-11)
  "q47", "q48", "q49", "q50", "q51", "q72", "q73", "q90", "q91", "q96",                # Abs Sr (12-21)
  "q59", "q60",                                                                        # Ratio Hg/Sr
  "q64", "q65", "q66", "q80", "q81", "q82", "q83", "q93", "q94", "q102",               # Ratio Yb/Sr
  "q79"                                                                                # Ratio Hg/Yb
)

nu_labels <- c(
  rep("nu3", 2),               # Abs Hg
  rep("nu8", 10),              # Abs Yb
  rep("nu12", 21),             # Abs Sr
  rep("nu3_over_nu12", 2),     # Ratio Hg / Sr
  rep("nu8_over_nu12", 10),    # Ratio Yb / Sr
  rep("nu3_over_nu8", 1)       # Ratio Hg / Yb
)

references <- c(
  # Abs Hg
  "[McFerran2012,corrected2015]", "[Tyumenev2016]",
  # Abs Yb
  "[Kohno2009]", "[Yasuda2012]", "[Lemke2009]", "[Park2013]", "[Pizzocaro2017]", 
  "[Kim2017]", "[Luo2020]", "[McGrew2019]", "[Pizzocaro2020]", "[Kobayashi2020]",
  # Abs Sr
  "[Boyd2007]", "[Campbell2008]", "[Baillard2008]", "[Hong2009]", "[Falke2011]", 
  "[Yamaguchi2012]", "[Akamatsu2014b]", "[Tanabe2015]", "[Lin2015]", "[Falke2014]", 
  "[LeTargat2013]", "[Lodewyck2016]", "[Grebing2016(Oct14)]", "[Grebing2016(Jun15)]", 
  "[Hachisu2017]", "[Hachisu2017b,Uncertx100Exclude]", "[Hobson2020]", 
  "[Schwarz2020-1.65E-16unc]", "[Nemitz2020]", "[Grotti2018]", "[Leopardi2020]",
  # Ratio Hg/Sr
  "[Yamanake2015]", "[Tyumenev2016]",
  # Ratio Yb/Sr
  "[Akamatsu2014aErrata]", "[Takamoto2015]", "[Nemitz2016]", "[Fujieda2018]", 
  "[Grotti2018]", "[Pizzocaro2020b]", "[Pizzocaro2020b]", "[Hisai2020]", 
  "[Kobayashi2020]", "[Beloy2021]",
  # Ratio Hg/Yb
  "[Ohmae2020]"
)

# Combine into a single data frame
measurement_df <- data.frame(
  Q_ID = q_labels,
  Nu_Type = nu_labels,
  Reference = references,
  # We convert to numeric here because standard R data frames do not natively 
  # support mpfr lists as columns without complex workarounds. 
  # Because the numbers are scaled to O(1), no precision is lost.
  Y_Offset_Scaled = as.numeric(Y_mpfr),
  U_Unc_Scaled = as.numeric(u_scaled_mpfr),
  stringsAsFactors = FALSE
)


library(ggplot2)

# 0. Re-create the Measurement_Type column
measurement_df$Measurement_Type <- ifelse(grepl("over", measurement_df$Nu_Type), 
                                          "Ratio", 
                                          "Absolute Frequency")
# Lock the order so Absolute appears on the left/top
measurement_df$Measurement_Type <- factor(measurement_df$Measurement_Type, 
                                          levels = c("Absolute Frequency", "Ratio"))

# 1. Set the order for your Nu_Types so they group logically
measurement_df$Nu_Type <- factor(measurement_df$Nu_Type, levels = c(
  "nu3", "nu8", "nu12",               # Absolutes first
  "nu3_over_nu8", "nu3_over_nu12", "nu8_over_nu12" # Ratios second
))

# 2. Sort the entire data frame by Nu_Type
measurement_df <- measurement_df[order(measurement_df$Nu_Type), ]

# 3. Lock the Q_ID factor levels to this new sorted order
measurement_df$Q_ID <- factor(measurement_df$Q_ID, levels = unique(measurement_df$Q_ID))

# 4. Plot it with the facet
ggplot(measurement_df, aes(x = Q_ID, y = Y_Offset_Scaled, color = Nu_Type)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = Y_Offset_Scaled - U_Unc_Scaled, 
                    ymax = Y_Offset_Scaled + U_Unc_Scaled), 
                width = 0.4) +
  
  # The vline was removed since the facet now handles the separation
  
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "right",
    strip.text = element_text(face = "bold", size = 11)
  ) +
  labs(
    x = "Measurement ID (Q_ID)",
    y = "Scaled Offset (Y x 1e15)",
    color = "Transition"
  ) +
  # scales = "free" unlocks BOTH the X and Y axes for each panel
  facet_wrap(~ Measurement_Type, scales = "free")




# 1. Create the Uncertainty Tier column
cutoff <- 5.0

measurement_df$Unc_Tier <- ifelse(
  measurement_df$U_Unc_Scaled > cutoff, 
  "High Uncertainty", 
  "High Precision"
)

# 2. Convert to a factor to control the order the panels appear in
measurement_df$Unc_Tier <- factor(
  measurement_df$Unc_Tier, 
  levels = c("High Precision", "High Uncertainty")
)

# 3. Plot with the multi-variable facet
ggplot(measurement_df, aes(x = Q_ID, y = Y_Offset_Scaled, color = Nu_Type)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = Y_Offset_Scaled - U_Unc_Scaled, 
                    ymax = Y_Offset_Scaled + U_Unc_Scaled), 
                width = 0.4) +
  
  geom_hline(yintercept = 0, color = "black", linetype = "dotted") +
  
  # Facet by BOTH the uncertainty tier and the measurement type
  # ncol = 2 organizes them nicely on the screen
  facet_wrap(~ Unc_Tier + Measurement_Type, scales = "free", ncol = 2) +
  
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "right",
    strip.text = element_text(face = "bold", size = 10)
  ) +
  labs(
    x = "Measurement ID (Q_ID)",
    y = "Scaled Offset (Y x 1e15)",
    color = "Transition"
  )
