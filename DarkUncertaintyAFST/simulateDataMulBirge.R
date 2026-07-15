############################################################################
### 
### Dark uncertainty for atomic clock data
### Simulate data assuming conditions necessary for multiplicative birge 
### 
### Suzanne Thornton
### March 2026
### Code adapted from Amanda Koepke and Angela Folz
### Code optimization assistance from Gemini
### 
### Assumption: each measurement x_i comes from a N(mu, c^2u_i^2) distribution
### using Angela's technique for generating u_i values
############################################################################
rm(list=ls())
library(data.table)
library(tidyverse)

path = "/Users/smt3/Documents/GitHub/atomic-clock/"
simdatfolder = "DarkUncertaintyAFST/simulatedData/" #added to gitignore
figfolder = "DarkUncertaintyAFST/figures/" #added to gitignore

##Supporter functions
source(paste0(path, "DarkUncertaintyAFST/", "MB_functions.R"))

##BACON data
ratiolab = "YbSr"
ratiodf = read_csv(paste0(path, "Data/ClockComparison2025/BayesianAnalysisData/ErYb_",ratiolab,"_data.csv"))
bacon_measurements = ratiodf$offset
bacon_uncertainties = ratiodf$statistical_unc

mu = 0
k_cov_factor = 1.96 #1.96 for 95% intervals





# --- Write CSV files containing MB simulated data 

N_new = c(5, 13, 33, 100)
true_c = c(1.5, 2, 2.5)
mu = 0
N_cov = 1E4

k = 1 #res_all index for bias/coverage plots
res_all = matrix(rep(NA, length(N_new)*length(true_c)*5), ncol = 5)
colnames(res_all) = c("N", "true_c", "bias", "unadj_coverage", "adj_coverage")

for(n in N_new){
  for(c_tmp in true_c){
    # 1. Vectorized generation of uncertainties and measurements
    uncertainties_cov = runif(n * N_cov, min(bacon_uncertainties), max(bacon_uncertainties))
    measurements_cov = rnorm(n * N_cov, mean = mu, sd = uncertainties_cov * c_tmp)
    
    # 2. Build a single long data table and calculate Birge results by group
    dt = data.table(
      run_id = rep(1:N_cov, each = n),
      day = rep(1:n, N_cov),
      u = uncertainties_cov,
      x = measurements_cov
    )
    
    # Calculate weighted mean and Birge errors for all runs at once
    # SD stands for subset of data, holds the actual data for the current group defined by "by"
    birge_res_dt = dt[, birge_wrapper(as.matrix(.SD)), by = run_id, .SDcols = c("x", "u")]
    
    # 3. Create bounds and merge back to format required for file output
    dt_merged = merge(dt, birge_res_dt, by = "run_id")
    # `:=` is data.table package command that saves memory and avoids creating multiple copies of large object
    dt_merged[, `:=`(
      lb = mean_birge - k_cov_factor * u_birge,
      ub = mean_birge + k_cov_factor * u_birge,
      lb_corrected = mean_birge - k_cov_factor * u_birge_corrected,
      ub_corrected = mean_birge + k_cov_factor * u_birge_corrected
    )]
    
    # .() is data.table package command that selects these specific columns
    dat_to_write_fin = dt_merged[, .(day, x, u, lb, ub, lb_corrected, ub_corrected)]
    fwrite(dat_to_write_fin, file = paste0(path, simdatfolder, "simDataMulBirge_N", n, "c", c_tmp, "_", N_cov, "iter_20260528.csv"))

    # 4. Vectorized calculation of bias and coverage
    cov_unadj = birge_res_dt[, mean((mean_birge - k_cov_factor * u_birge <= mu) & 
                                      (mean_birge + k_cov_factor * u_birge >= mu))]
    cov_adj = birge_res_dt[, mean((mean_birge - k_cov_factor * u_birge_corrected <= mu) & 
                                    (mean_birge + k_cov_factor * u_birge_corrected >= mu))]
    mean_bias = mean(birge_res_dt$mean_birge)
    
    res_all[k,] = c(n, c_tmp, mean_bias, cov_unadj, cov_adj)
    k = k + 1
  }
}


##---Single simulation----------------------------------
N = length(bacon_measurements) # 100
birg_constant = 2 
mu = 0

## Simulated data (Vectorized)
uncertainties = runif(N, min(bacon_uncertainties), max(bacon_uncertainties)) 
measurements = rnorm(N, mean = mu, sd = uncertainties * birg_constant)
sim_dat = data.frame(day = 1:N, x = measurements, u = uncertainties)
summary(sim_dat)

# FIX: Define the plotting variables that were missing in the original code
mu_hat_b = weighted_mean(sim_dat$x, sim_dat$u)
u_mu_hat_b = u_MB(sim_dat$x, sim_dat$u)
u_mu_hat_b_corrected = u_MB(sim_dat$x, sim_dat$u, correction = TRUE)


##---Plots--------------------------------------
# plot real BACON2 data
ggplot(ratiodf, aes(x=as.character(date), y=offset)) +
  geom_point(size=1) +
  geom_errorbar(aes(ymin=offset-statistical_unc, ymax=offset+statistical_unc), width=0) +
  ylim(c(-110,-90)) +
  theme_bw() +
  ggtitle("BACON2 data")

# plot simulated data set
# plot includes mu_hat with uncertainty bars showing k*u, where k=1 and true mu (red)
p1 <- ggplot(sim_dat, aes(x=day, y=x)) +
  geom_point(size=1) +
  geom_errorbar(aes(ymin=x-u, ymax=x+u), width=0) +
  geom_hline(aes(yintercept = mu, color = "true mu", linetype = "true mu")) + 
  geom_hline(aes(yintercept = mu_hat_b, color = "weighted mean", linetype = "weighted mean")) + 
  geom_hline(aes(yintercept = mu_hat_b + u_mu_hat_b, color = "unadjusted Birge bounds", linetype = "unadjusted Birge bounds")) + 
  geom_hline(aes(yintercept = mu_hat_b - u_mu_hat_b, color = "unadjusted Birge bounds", linetype = "unadjusted Birge bounds")) +
  geom_hline(aes(yintercept = mu_hat_b + u_mu_hat_b_corrected, color = "adjusted Birge bounds", linetype = "adjusted Birge bounds")) + 
  geom_hline(aes(yintercept = mu_hat_b - u_mu_hat_b_corrected, color = "adjusted Birge bounds", linetype = "adjusted Birge bounds")) + 
  # ylim(c(-6,6)) +
  # ylim(c(-110,-90)) + # if use mu=mean(BACON2 data)
  scale_color_manual(name = "Legend",
              breaks = c("true mu", 
               "weighted mean", 
               "unadjusted Birge bounds", 
               "adjusted Birge bounds"),
              values = c("true mu" = "red", 
                "weighted mean" = "purple", 
                "unadjusted Birge bounds" = "purple", 
                "adjusted Birge bounds" = "purple"),
              labels = c("true mu" = expression(mu), 
               "weighted mean" = expression(hat(mu)[WM]),
               "unadjusted Birge bounds" = "unadjusted Birge bounds",
               "adjusted Birge bounds" = "adjusted Birge bounds")) +
  scale_linetype_manual(name = "Legend",
              breaks = c("true mu", 
               "weighted mean", 
               "unadjusted Birge bounds", 
               "adjusted Birge bounds"),
              values = c("true mu" = 1, 
                "weighted mean" = 1, 
                "unadjusted Birge bounds" = 2, 
                "adjusted Birge bounds" = 3),
              labels = c("true mu" = expression(mu), 
               "weighted mean" = expression(hat(mu)[WM]),
               "unadjusted Birge bounds" = "unadjusted Birge bounds",
               "adjusted Birge bounds" = "adjusted Birge bounds")) +
  theme_bw() +
  ylab(expression(x[i] %+-% u(x[i]))) + xlab("Day") + 
  ggtitle(paste("Simulated data, N =", N, ", c =", birg_constant))

ggsave("DarkUncertaintyAFST/mul_birge_plot.png", plot = p1, width = 6, height = 4, units = "in", dpi = 300)


##---Investigate coverage------------------------------
N_cov = 1E4

uncertainties_cov = runif(N*N_cov, min(bacon_uncertainties), max(bacon_uncertainties))
measurements_cov = rnorm(N * N_cov, mean = mu, sd = uncertainties_cov * 2)

dt_single = data.table(
  run_id = rep(1:N_cov, each = N),
  day = rep(1:N, N_cov),
  u = uncertainties_cov,
  x = measurements_cov
)

birge_res_single = dt_single[, birge_wrapper(as.matrix(.SD)), by = run_id, .SDcols = c("x", "u")]

# Vectorized summary stats
single_stats = birge_res_single[, .(
  cov_unadj = mean((mean_birge - k_cov_factor * u_birge <= mu) & 
                     (mean_birge + k_cov_factor * u_birge >= mu)),
  cov_adj = mean((mean_birge - k_cov_factor * u_birge_corrected <= mu) & 
                   (mean_birge + k_cov_factor * u_birge_corrected >= mu)),
  mean_mu_hat = mean(mean_birge),
  mean_u_unadj = mean(u_birge),
  mean_u_adj = mean(u_birge_corrected)
)]

print(single_stats)