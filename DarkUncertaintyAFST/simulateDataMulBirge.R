############################################################################
### 
### Dark uncertainty for atomic clock data
### Simulate data assuming conditions necessary for multiplicative birge 
### 
### Suzanne Thornton
### March 2026
### Code adapted from Amanda Koepke and Angela Folz
### 
### Assumption: each measurement x_i comes from a N(mu, c^2u_i^2) distribution
### using Angela's technique for generating u_i values
############################################################################
rm(list=ls())
library(readr)
library(dplyr)
library(ggplot2)
library(metafor)

path = "/Users/smt3/Documents/GitHub/atomic-clock/"
simdatfolder = "DarkUncertaintyAFST/simulatedData/"

##BACON data 
ratiolab = "YbSr"
ratiodf = read_csv(paste0(path, "Data/ClockComparison2025/BayesianAnalysisData/ErYb_",ratiolab,"_data.csv"))
bacon_measurements = ratiodf$offset
bacon_uncertainties = ratiodf$statistical_unc
N = length(bacon_measurements) #100

##True parameters
birg_constant = 2 #0.5, 1.5
mu = 0

##---Supporting functions------------------------------
norm_gen <- function(sd_term){rnorm(1, mean=mu, sd = sd_term*birg_constant)}

#standard uncertainty of WM
std_u <- function(uncertainties){
    return(sum(1/(uncertainties^2))^(-1/2))
}

#mu_hat_B
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


##---Single simulation----------------------------------
##Simulated data
##these uncertainties are (likely) std deviations
uncertainties = runif(N, min(bacon_uncertainties), max(bacon_uncertainties)) # simulate "known" uncertainties for each day from a uniform, 
                                                         # with the lower and upper bound set to be what was observed 
                                                         # in the BACON2 data for a particular ratio
measurements = sapply(uncertainties, norm_gen)
sim_dat = data.frame(day = 1:N, x = measurements, u = uncertainties)
summary(sim_dat)


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
N_cov = 1E4  ## select number of runs

uncertainties_cov = runif(N*N_cov, min(bacon_uncertainties), max(bacon_uncertainties))
measurements_cov = sapply(uncertainties_cov, norm_gen)
day_cov = rep(1:N, N_cov)
sim_dat_long = cbind(day_cov, measurements_cov, uncertainties_cov) 

sim_dat_list = list()
run = 1
for(i in 1:N_cov){
  sim_dat_list[[i]] = sim_dat_long[(run):(N+run-1), ]
  colnames(sim_dat_list[[i]]) = c("day", "x", "u")
  run = run + N 
}

birge_res = lapply(sim_dat_list, birge_wrapper)

##Nominal coverage is for k=1 (68%)
k_cov_factor = 1.96 #1.96 for 95% intervals
coverage = 0
coverage_corrected = 0
for (i in 1:N_cov){
  tmp = birge_res[[i]]
  if((tmp$mean_birge - k_cov_factor*tmp$u_birge <= mu) & (tmp$mean_birge + k_cov_factor*tmp$u_birge >= mu)){ coverage = coverage + 1}
  if((tmp$mean_birge - k_cov_factor*tmp$u_birge_corrected <= mu) & (tmp$mean_birge + k_cov_factor*tmp$u_birge_corrected >= mu)){ coverage_corrected = coverage_corrected + 1}
}
coverage/N_cov
coverage_corrected/N_cov

##---Write CSV----------------------------------
## that contains: "Day","x","u","k", "lb","ub", "lb_corrected", "ub_corrected"
k_cov_factor = 1.96 #1.96 for 95% intervals
birge_CIs = list()
for(i in 1:N_cov){
  birge_CIs[[i]] = matrix(rep(
                    c(k_cov_factor,
                    birge_res[[i]]$mean_birge - k_cov_factor*birge_res[[i]]$u_birge,
                    birge_res[[i]]$mean_birge + k_cov_factor*birge_res[[i]]$u_birge,
                    birge_res[[i]]$mean_birge - k_cov_factor*birge_res[[i]]$u_birge_corrected,
                    birge_res[[i]]$mean_birge + k_cov_factor*birge_res[[i]]$u_birge_corrected),
                    N), ncol=5, byrow=TRUE)
  colnames(birge_CIs[[i]]) = c("k", "lb", "ub", "lb_corrected", "ub_corrected")
}

dat_to_write = list()
for(i in 1:N_cov){
  dat_to_write[[i]] = cbind(sim_dat_list[[i]], birge_CIs[[i]])
}

library(data.table)
dat_to_write_fin = rbindlist(lapply(dat_to_write, as.data.frame))
head(dat_to_write_fin,40)

# write.csv2(dat_to_write_fin, file = "DarkUncertaintyAFST/simulatedData/simDataMulBirge_N13mu0_10000iter_20260420.csv", row.names=FALSE)


##---Analyse REM Data----------------------------------
rem_data = read_csv(paste0(path, "DarkUncertaintyAFST/simulatedData/simDataRandomEffects_N13xi3_10000iter_20260326.csv"))
rem_mu = rem_data$mu[1]

## Break data into a list of matricies with 13 rows
chunk_size = rem_data$N[1]
rem_data_list = split(rem_data, (seq_len(nrow(rem_data)) - 1) %/% chunk_size)
rem_mtx_list = lapply(rem_data_list, as.matrix)
rem_birge_res = lapply(rem_mtx_list, birge_wrapper) ##Note input must be a matrix

k_cov_factor = 1.96 #1.96 for 95% intervals
coverage = 0
coverage_corrected = 0
for (i in seq_along(rem_data_list)){
  tmp = rem_birge_res[[i]]
  if((tmp$mean_birge - k_cov_factor*tmp$u_birge <= rem_mu) & (tmp$mean_birge + k_cov_factor*tmp$u_birge >= rem_mu)){ coverage = coverage + 1}
  if((tmp$mean_birge - k_cov_factor*tmp$u_birge_corrected <= rem_mu) & (tmp$mean_birge + k_cov_factor*tmp$u_birge_corrected >= rem_mu)){ coverage_corrected = coverage_corrected + 1}
}
coverage/length(rem_data_list)
coverage_corrected/length(rem_data_list)
