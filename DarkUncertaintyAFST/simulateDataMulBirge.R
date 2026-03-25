############################################################################
### 
### Dark uncertainty for atomic clock data
### Simulate data assuming conditions necessary for multiplicative birge 
### 
### Suzanne Thornton
### March 2026
### Code adapted from Amanda Koepke and Angela Folz
### 
############################################################################

library(readr)
library(dplyr)
library(ggplot2)
library(metafor)

path = "/Users/smt3/Documents/GitHub/atomic-clock/"
simdatfolder = "DarkUncertaintyAFST/simulatedData/"

## Assumption: each measurement x_i comes from a N(mu, c^2u_i^2) distribution
## use Angela's technique for generating u_i values

##BACON data 
ratiolab = "YbSr"
ratiodf = read_csv(paste0(path, "Data/ClockComparison2025/BayesianAnalysisData/ErYb_",ratiolab,"_data.csv"))
bacon_measurements = ratiodf$offset
bacon_uncertainties = ratiodf$statistical_unc
N = length(bacon_measurements)

##True parameters
birg_constant = 2
mu = 0

##Simulated data
##these uncertainties are (likely) std deviations
uncertainties = runif(N, min(bacon_uncertainties), max(bacon_uncertainties)) # simulate "known" uncertainties for each day from a uniform, 
                                                         # with the lower and upper bound set to be what was observed 
                                                         # in the BACON2 data for a particular ratio
norm_gen <- function(sd_term){rnorm(1, mean=mu, sd = sd_term*birg_constant)}
measurements = sapply(uncertainties, norm_gen)
sim_dat = data.frame(day = 1:N, x = measurements, u = uncertainties)
summary(sim_dat)

std_u <- function(uncertainties){
    return(sum(1/(uncertainties^2))^(-1/2))
}

#mu_hat_B
weighted_mean <- function(measurements, uncertainties){
    u_all = std_u(uncertainties)
    return(u_all^2 * sum(measurements/(uncertainties^2)))
}

#u(mu_hat_B) (from Toman et al. 2012) 
u_MB <- function(measurements, uncertainties){
    n = length(measurements)
    ave_WM = weighted_mean(measurements, uncertainties)
    chi_sq_obs = sum( ((ave_WM - measurements)^2) / ((uncertainties)^2))
    u = std_u(uncertainties)
    return(u*sqrt(chi_sq_obs/(n-1))) ##or use n-3 as in Toman 2012
}

## mu_hat_B and u(mu_hat_B) for simulated data
mu_hat_b = weighted_mean(sim_dat$x, sim_dat$u)
u_mu_hat_b = u_MB(sim_dat$x, sim_dat$u)
u_wm = std_u(sim_dat$u) ##unadjusted uncertainty


# plot real BACON2 data
ggplot(ratiodf, aes(x=as.character(date), y=offset)) +
  geom_point(size=1) +
  geom_errorbar(aes(ymin=offset-statistical_unc, ymax=offset+statistical_unc), width=0) +
  ylim(c(-110,-90)) +
  theme_bw() +
  ggtitle("BACON2 data")

# plot simulated data set
# TODO: add to this plot mu_hat with uncertainty bars showing k*u, where k=1 and plot true mu
ggplot(sim_dat, aes(x=day, y=x)) +
  geom_point(size=1) +
  geom_errorbar(aes(ymin=x-u, ymax=x+u), width=0) +
  geom_hline(yintercept = mu, color = "red") + 
  geom_hline(yintercept = mu_hat_b, color = "purple") + 
  geom_hline(yintercept = mu_hat_b + u_mu_hat_b, linetype = 2, color = "purple") + 
    geom_hline(yintercept = mu_hat_b - u_mu_hat_b, linetype = 2, color = "purple") + 
  ylim(c(-6,6)) +
  # ylim(c(-110,-90)) + # if use mu=mean(BACON2 data)
  theme_bw() +
  ggtitle(paste("Simulated data, N =", N, ", c =", birg_constant))


##TODO: generate many data sets to investigate coverage (i'm predicting over coverage) 