##---Supporting functions------------------------------
norm_gen <- function(sd_term){rnorm(1, mean=0, sd = sd_term*2)} #mu=0, birg_constant = 2
norm_gen_multiple <- function(sd_times_birge_c){rnorm(1, mean=0, sd = sd_times_birge_c)} #mu=0

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
