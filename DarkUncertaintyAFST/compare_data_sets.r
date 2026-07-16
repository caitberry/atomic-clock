## code to create plots of simulated MB and REM data (for various c and xi values)

rm(list=ls())
library(readr)
library(data.table)
library(tidyverse)
library(gridExtra)

path = "/Users/smt3/Documents/GitHub/atomic-clock/"
simdatfolder = "DarkUncertaintyAFST/simulatedData/" #added to gitignore
figfolder = "DarkUncertaintyAFST/figures/" #added to gitignore

##---BACON data-------------------------------------
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

df1 <- ratiodf %>% 
    transmute(x_val = as.character(date), y_val = offset_centered, uncertainty = statistical_unc, source = "BACON2")
df2 <- sim_dat_MB %>% 
    transmute(x_val = as.character(ratiodf$date), y_val = x, uncertainty = u, source = "Simulated MB")
df3 <- sim_dat_REM %>% 
    transmute(x_val = as.character(ratiodf$date), y_val = x, uncertainty = u, source = "Simulated REM")

p1 <- ggplot(df1, aes(x = x_val, y = y_val)) +
  geom_point(size = 1) +
  geom_errorbar(aes(ymin = y_val - uncertainty, ymax = y_val + uncertainty), width = 0) +
  theme_bw() +
  ylim(-6,9) +
  labs(
    title = "BACON Data",
    x = "Day / Date",
    y = expression(x[i] %+-% u(x[i]))
  )

p2 <- ggplot(df2, aes(x = x_val, y = y_val)) +
  geom_point(size = 1) +
  geom_errorbar(aes(ymin = y_val - uncertainty, ymax = y_val + uncertainty), width = 0) +
  theme_bw() +
  ylim(-6,9) +
  labs(
    title = "MB Data",
    x = "Day / Date",
    y = expression(x[i] %+-% u(x[i]))
  )

p3 <- ggplot(df3, aes(x = x_val, y = y_val)) +
  geom_point(size = 1) +
  geom_errorbar(aes(ymin = y_val - uncertainty, ymax = y_val + uncertainty), width = 0) +
  theme_bw() +
  ylim(-6,9) +
  labs(
    title = "REM Data",
    x = "Day / Date",
    y = expression(x[i] %+-% u(x[i]))
  )


bacon_comparison <- grid.arrange(p1, p2, p3, nrow = 1)
ggsave(paste0("DarkUncertaintyAFST/figures/sim_data_comparison_N13.png"), bacon_comparison, device = "png")

##Note: currently not plotting N13 for c_val of 1.5 or 2.5 or xi of 1 or 10

##---Other N Simulated Data--------------------------
mu = 0
N_new = c(5, 13, 33, 100)
c_vals = c(1.5, 2, 3)
xi_vals = c(1, 2, 3)

MB_data = REM_data = list()
i_MB = i_REM = 1 

for(N in N_new){
  uncertainties <- runif(N, min(bacon_uncertainties), max(bacon_uncertainties))

  for(c_val in c_vals){
    measurements_MB <- sapply(uncertainties, function(sd_term){rnorm(1, mean=mu, sd = sd_term*c_val)} )
    sim_dat_MB <- data.frame(day = 1:N, x = measurements_MB, u = uncertainties)

    MB_data[[i_MB]] <- sim_dat_MB %>%
        transmute(x_val = as.character(sim_dat_MB$day), 
        y_val = x, 
        uncertainty = u, 
        source = "Simulated MB", 
        N = N,
        c = c_val)

    i_MB = i_MB + 1
  }

  for(xi in xi_vals){
    lambda <- rnorm(N, 0, xi)
    epsilon <- rnorm(N, 0, uncertainties)

    measurements_REM <- mu + lambda + epsilon
    sim_dat_REM <- data.frame(day=1:N, x=measurements_REM, u=uncertainties)

    REM_data[[i_REM]] <-  sim_dat_REM %>%
        transmute(x_val = as.character(sim_dat_REM$day), 
        y_val = x, 
        uncertainty = u, 
        source = "Simulated REM", 
        N = N,
        xi = xi)

    i_REM = i_REM + 1
  }
}
# length(REM_data)
# length(MB_data)

p_sim_MB = p_sim_REM = list()

for(i in 1:length(MB_data)){
  #set y-axis limits 
  if (MB_data[[i]]$c[1] < 2.5) {
    my_y_lim <- c(-5, 5)
  } 
  if ((MB_data[[i]]$c[1] > 2) & (REM_data[[i]]$xi[1] < 3.5)){
    my_y_lim <- c(-10, 10)
  } 
  if (MB_data[[i]]$c[1] > 3) {
    my_y_lim <- c(-40, 40)
  }
  p_sim_MB[[i]] <- ggplot(MB_data[[i]], aes(x = x_val, y = y_val)) +
    geom_point(size = 1) +
    geom_errorbar(aes(ymin = y_val - uncertainty, ymax = y_val + uncertainty), width = 0) +
    theme_bw() +
#    ylim(my_y_lim) + 
    ylim(-20,20) + 
    labs(
      title = paste0("MB Data; N = ", MB_data[[i]]$N[1], "; c = ", MB_data[[i]]$c[1]),
      x = "Day / Date",
      y = expression(x[i] %+-% u(x[i]))
    )
}

for(i in 1:length(REM_data)){
  #set y-axis limits 
  if (REM_data[[i]]$xi[1] < 3) {
    my_y_lim <- c(-5, 5)
  } 
  if ((REM_data[[i]]$xi[1] > 1) & (REM_data[[i]]$xi[1] < 10)){
    my_y_lim <- c(-10, 10)
  } 
  if (REM_data[[i]]$xi[1] > 3) {
    my_y_lim <- c(-40, 40)
  }

  p_sim_REM[[i]] <- ggplot(REM_data[[i]], aes(x = x_val, y = y_val)) +
    geom_point(size = 1) +
    geom_errorbar(aes(ymin = y_val - uncertainty, ymax = y_val + uncertainty), width = 0) +
    theme_bw() +
#    ylim(my_y_lim) +
    ylim(-20, 20) + 
    labs(
      title = paste0("REM Data; N = ", REM_data[[i]]$N[1], "; xi = ", REM_data[[i]]$xi[1]),
      x = "Day / Date",
      y = expression(x[i] %+-% u(x[i]))
    )
}


for (i in c(1,4,7,10)){
  combined_plot <- arrangeGrob(p_sim_MB[[i]], p_sim_MB[[i+1]], p_sim_MB[[i+2]],
                    p_sim_REM[[i]], p_sim_REM[[i+1]], p_sim_REM[[i+2]], nrow = 2)
  ggsave(paste0("DarkUncertaintyAFST/figures/sim_data_comparison_N", REM_data[[i]]$N[1], ".png"), combined_plot, device = "png")
}

i=10
grid.arrange(p_sim_MB[[i]], p_sim_MB[[i+1]], p_sim_MB[[i+2]],# p_sim_MB[[i+3]], 
            p_sim_REM[[i]], p_sim_REM[[i+1]], p_sim_REM[[i+2]], nrow = 2)
