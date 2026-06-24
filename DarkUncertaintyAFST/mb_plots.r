# ---------------------------------------------------------
# --- Read CSV files and produce Bias/Coverage Plots
# --- for MB analysis of MB data and REM data
# ---------------------------------------------------------
rm(list=ls())
library(readr)
library(data.table)
library(tidyverse)
library(tidyr)

N_new = c(5, 13, 33, 100)
k_cov_factor = 1.96 #1.96 for 95% intervals
#mu = 0
#N_cov = 1E4 

path = "/Users/smt3/Documents/GitHub/atomic-clock/"
simdatfolder = "DarkUncertaintyAFST/simulatedData/" #added to gitignore
figfolder = "DarkUncertaintyAFST/figures/" #added to gitignore

##--- MB Analysis Functions -------------------

weighted_mean <- function(measurements, uncertainties){
    u_all = std_u(uncertainties)
    return(u_all^2 * sum(measurements/(uncertainties^2)))
}

#u(mu_hat_B) (from Toman et al. 2012) 
u_MB <- function(measurements, uncertainties, correction = FALSE){
  n = length(measurements)
  ave_WM = weighted_mean(measurements, uncertainties)
  chi_sq_obs = sum( ((ave_WM - measurements)^2) / ((uncertainties)^2))
  u = std_u(uncertainties)
  if(correction){  
    return(u*sqrt(chi_sq_obs/(n-3))) ##n-3 is correction suggested in Toman 2012
  } else { 
    return(u*sqrt(chi_sq_obs/(n-1)))
  }
}

birge_wrapper <- function(mtx){
  mean_birge = weighted_mean(mtx[,"x"], mtx[,"u"])
  u_birge = u_MB(mtx[,"x"], mtx[,"u"])
  u_birge_corrected = u_MB(mtx[,"x"], mtx[,"u"], correction = TRUE)
  return(data.frame(mean_birge = mean_birge, u_birge = u_birge, u_birge_corrected = u_birge_corrected))
}


##---REM Data----------------------------------
true_xi = c(1, 3, 10)
rem_data_names = "simDataRandomEffects_"

rem_data_mb_analysis = replicate(length(N_new)*length(true_xi),
                        list(N = NULL, xi = NULL, cov_un = NULL, cov = NULL, bias_both = NULL),
                        simplify = FALSE)
counter = 1

for (n in N_new){
  for (xi in true_xi){
    rem_data = read_csv(paste0(path, simdatfolder, rem_data_names, "N", n, "xi", xi, "_10000iter_20260326.csv"))
    rem_mu = rem_data$mu[1]

    ## Break data into a list of matricies with n rows
    chunk_size = rem_data$N[1]
    rem_data_list = split(rem_data, (seq_len(nrow(rem_data)) - 1) %/% chunk_size)
    rem_mtx_list = lapply(rem_data_list, as.matrix)
    rem_birge_res = lapply(rem_mtx_list, birge_wrapper) ##Note input must be a matrix

    coverage = 0
    coverage_corrected = 0
    mu_hat_un = NULL
    mu_hat_adj = NULL
    u_mu_hat_un = NULL 
    u_mu_hat_adj = NULL
    for (i in seq_along(rem_data_list)){
      tmp = rem_birge_res[[i]]
      if((tmp$mean_birge - k_cov_factor*tmp$u_birge <= rem_mu) & (tmp$mean_birge + k_cov_factor*tmp$u_birge >= rem_mu)){ coverage = coverage + 1}
      if((tmp$mean_birge - k_cov_factor*tmp$u_birge_corrected <= rem_mu) & (tmp$mean_birge + k_cov_factor*tmp$u_birge_corrected >= rem_mu)){ coverage_corrected = coverage_corrected + 1}

      mu_hat_un = c(mu_hat_un, tmp$mean_birge)
      mu_hat_adj = c(mu_hat_adj, tmp$mean_birge)

      u_mu_hat_un = c(u_mu_hat_un, tmp$u_birge)
      u_mu_hat_adj = c(u_mu_hat_adj, tmp$u_birge_corrected)
    }

    rem_data_mb_analysis[[counter]]$N = n
    rem_data_mb_analysis[[counter]]$xi = xi
    rem_data_mb_analysis[[counter]]$cov_un = coverage/length(rem_data_list)
    rem_data_mb_analysis[[counter]]$cov = coverage_corrected/length(rem_data_list)
    rem_data_mb_analysis[[counter]]$bias_both = mean(mu_hat_un)

    counter = counter + 1
  }
}

##---MB Data-----------------------------------
# method = c("unadjusted", "adjusted") 
true_c = c(1.5, 2, 2.5)
mb_data_names = "simDataMulBirge_"
mb_mu = 0

mb_data_mb_analysis = replicate(length(N_new)*length(true_c),
                        list(N = NULL, c = NULL, cov_un = NULL, cov = NULL, bias_both = NULL),
                        simplify = FALSE)
counter = 1

for (n in N_new){
  for (c in true_c){
    mb_data = read_csv(paste0(path, simdatfolder, mb_data_names, "N", n, "c", c, "_10000iter_20260528.csv"))
    
    coverage = 0
    coverage_corrected = 0
    mu_hat = NULL
    runs = length(mb_data) %% n

    for (i in seq_along(runs)){
      tmp = mb_data[i,]
      tmp_mu_hat = (tmp$ub/tmp$lb)/2
      ##double check: lb and ub already account for k? 
      if((tmp$lb <= mb_mu) & (tmp$ub >= mb_mu)){ coverage = coverage + 1}
      if((tmp$lb_corrected <= mb_mu) & (tmp$ub_corrected >= mb_mu)){ coverage_corrected = coverage_corrected + 1}

      mu_hat = c(mu_hat, tmp_mu_hat)
    }

    mb_data_mb_analysis[[counter]]$N = n
    mb_data_mb_analysis[[counter]]$c = c
    mb_data_mb_analysis[[counter]]$cov_un = coverage/runs
    mb_data_mb_analysis[[counter]]$cov = coverage_corrected/runs
    mb_data_mb_analysis[[counter]]$bias_both = mean(mu_hat)

    counter = counter + 1
  }
}

# --- Plots -------------------------------------
# coverage prob (y-axis)
# bias (y axis)
# N = c(5, 13, 33, 100) (x-axis)

# plot for mu bias 
#pdf(paste0(path, figfolder, "MBsim_MB_mubias_1E4iter.pdf"), width=8, height=6)

ggplot(data.frame(res_all), aes(x=N, y=bias, color=as.factor(true_c))) +
  geom_point(size=1) +
  geom_line() +
  theme_bw() +
  labs(color="true c") +
  ggtitle("MB est mu bias")

#dev.off()


# plot for mu coverage
#pdf(paste0(path, figfolder, "MBsim_MB_mucov_1E4iter.pdf"), width=8, height=6)

ggplot(res_all_long, aes(x=N, y=cov_prob, color=as.factor(true_c), linetype=as.factor(method))) +
  geom_point(size=1) +
  geom_line() +
  geom_hline(yintercept=0.95) +
  theme_bw() +
  labs(color="true c", linetype="method") +
  ggtitle("MB est mu coverage probability")

#dev.off()


