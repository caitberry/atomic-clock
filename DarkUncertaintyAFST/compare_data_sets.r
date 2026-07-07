## code to create plots of simulated MB and REM data (for various c and xi values)
## and compare to BACON data

rm(list=ls())
library(readr)
library(data.table)
library(tidyverse)
library(gridExtra)

path = "/Users/smt3/Documents/GitHub/atomic-clock/"
simdatfolder = "DarkUncertaintyAFST/simulatedData/" #added to gitignore
figfolder = "DarkUncertaintyAFST/figures/" #added to gitignore

##BACON data
ratiolab = "YbSr"
ratiodf = read_csv(paste0(path, "Data/ClockComparison2025/BayesianAnalysisData/ErYb_",ratiolab,"_data.csv"))
ratiodf$offset_centered = ratiodf$offset - mean(ratiodf$offset)
bacon_measurements = ratiodf$offset
bacon_uncertainties = ratiodf$statistical_unc
summary(ratiodf)

set.seed(101)
mu = 0
k_cov_factor = 1.96 #1.96 for 95% intervals
N = length(bacon_measurements)
uncertainties = runif(N, min(bacon_uncertainties), max(bacon_uncertainties))

##---MB simulated data--------------------------
birg_constant = 2

measurements_MB = sapply(uncertainties, function(sd_term){rnorm(1, mean=mu, sd = sd_term*birg_constant)} )
sim_dat_MB = data.frame(day = 1:N, x = measurements_MB, u = uncertainties)
summary(sim_dat_MB)

##---REM simulated data-------------------------
xi_true = 3

lambda = rnorm(N, 0, xi_true)
epsilon = rnorm(N, 0, uncertainties)

measurements_REM = mu + lambda + epsilon
sim_dat_REM = data.frame(day=1:N, x=measurements_REM, u=uncertainties)
summary(sim_dat_REM)

##---Plots--------------------------------------

combined_df <- bind_rows(
  ratiodf %>% 
    transmute(x_val = as.character(date), y_val = offset_centered, uncertainty = statistical_unc, source = "BACON2"),
  sim_dat_MB %>% 
    transmute(x_val = as.character(ratiodf$date), y_val = x, uncertainty = u, source = "Simulated MB"),
  sim_dat_REM %>% 
    transmute(x_val = as.character(ratiodf$date), y_val = x, uncertainty = u, source = "Simulated REM")
)

p_combined_tidy <- ggplot(combined_df, aes(x = x_val, y = y_val, color = source)) +
  geom_point(size = 1) +
  geom_errorbar(aes(ymin = y_val - uncertainty, ymax = y_val + uncertainty), width = 0) +
  scale_color_manual(values = c(
    "BACON2"        = "black",
    "Simulated MB"  = "deepskyblue3",
    "Simulated REM" = "tomato"
  )) +
  theme_bw() +
  labs(
    title = "Combined Data Sets",
    x = "Day / Date",
    y = expression(x[i] %+-% u(x[i])),
    color = "Data Source"
  )

print(p_combined_tidy)
