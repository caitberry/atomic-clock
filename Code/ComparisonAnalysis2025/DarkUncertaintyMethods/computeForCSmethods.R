adjust_uncertainty_birge <- function(measurements, uncertainties) {

  # Calculate initial Weighted Mean
  weights <- 1 / (uncertainties^2)
  w_mean <- sum(measurements * weights) / sum(weights)
  
  # Calculate initial Reduced Chi-Squared
  dof <- length(measurements) - 1
  chisq <- sum(((measurements - w_mean)^2) / (uncertainties^2))
  chisq_red <- chisq / dof
  
  # 4. Apply Logic
  if (chisq_red > 1) {
    birge_ratio <- sqrt(chisq_red)
    sigma_new <- uncertainties * birge_ratio
    
    # Return list with results
    return(list(
      method = "Birge Ratio (Multiplicative)",
      mu = w_mean,
      scaling_factor = birge_ratio,
      old_chisq_red = chisq_red,
      new_uncertainties = sigma_new
    ))
    
  } else {
    # No adjustment needed
    return(list(
      method = "Birge Ratio (Multiplicative)",
      mu = w_mean,
      scaling_factor = 1.0,
      old_chisq_red = chisq_red,
      new_uncertainties = uncertainties
    ))
  }
}



adjust_uncertainty_additive <- function(measurements, uncertainties) {
  
  n <- length(measurements)
  dof <- n - 1
  
  # Calculate initial stats to check if adjustment is needed
  weights_init <- 1 / (uncertainties^2)
  w_mean_init <- sum(measurements * weights_init) / sum(weights_init)
  chisq_init <- sum(((measurements - w_mean_init)^2) / (uncertainties^2))
  chisq_red_init <- chisq_init / dof
  
  # Define the target function for uniroot
  # We need (Reduced ChiSq - 1) to be 0
  target_func <- function(xi_a, x, sigma, df) {
    sigma_total_sq <- sigma^2 + xi_a^2
    w_new <- 1 / sigma_total_sq
    
    # Recalculate mean with new weights (critical step)
    mu_new <- sum(x * w_new) / sum(w_new)
    
    chisq_new <- sum(((x - mu_new)^2) / sigma_total_sq)
    return((chisq_new / df) - 1)
  }
  
  # 4. Apply Logic
  if (chisq_red_init > 1) {
    # Set search interval
    # Upper bound heuristic: 10x the max deviation from the mean
    max_search <- max(abs(measurements - w_mean_init)) * 10
    
    # Solve for xi_a
    solution <- uniroot(target_func, 
                        interval = c(0, max_search), 
                        x = measurements, 
                        sigma = uncertainties, 
                        df = dof)
    
    xi_a <- solution$root
    # Calculate final results
    sigma_new <- sqrt(uncertainties^2 + xi_a^2)
    weights_new <- 1 / (sigma_new^2)
    w_mean_new <- sum(measurements * weights_new) / sum(weights_new)
    
    return(list(
      method = "Dark Uncertainty (Additive)",
      mu = w_mean_new,
      additive_constant = xi_a,
      old_chisq_red = chisq_red_init,
      new_uncertainties = sigma_new
    ))
    
  } else {
    return(list(
      method = "Dark Uncertainty (Additive)",
      mu = w_mean_init,
      additive_constant = 0,
      old_chisq_red = chisq_red_init,
      new_uncertainties = uncertainties
    ))
  }
}

# --- Dummy Data ---
# A dataset where one point is very precise (0.01) but an outlier (12.1)
# # This highlights how Method B shifts the mean compared to Method A.
# vals <- c(10.5, 12.1, 9.8, 11.5, 10.0) 
# errs <- c(0.1,  0.01, 0.1, 0.1,  0.1) 

# # --- Run Method A ---
# multBirgeOut <- adjust_uncertainty_birge(vals, errs)
# addCSout <- adjust_uncertainty_additive(vals, errs)


analyzeData_TFmethods=function(data,trueTau,trueMu){

  multBirgeOut=adjust_uncertainty_birge(data$x,data$u)
  addCSout=adjust_uncertainty_additive(data$x,data$u)

  mu.est=c(multBirgeOut$mu, addCSout$mu)
  tau.est=c(NA,addCSout$additive_constant)
  muCI = data.frame(method=c(multBirgeOut$method,addCSout$method),
                    parameter="mu",
                    estimate=mu.est,
                    ci.lb=NA,
                    ci.ub=NA,
                    inInt=NA)
  tauCI = data.frame(method=c(multBirgeOut$method,addCSout$method),
                    parameter="tau",
                    estimate=tau.est,
                    ci.lb=NA,
                    ci.ub=NA,
                    inInt=NA)
  
  bothCIout=rbind(muCI,tauCI)
  placeholdersForBayesOut=data.frame(bad.rhat=0,
                                     min.neff=NA,
                                     num.Divergent=NA,
                                     tree.Depth.Ex=NA,
                                     min.bfmi=NA)
  return(cbind(bothCIout,placeholdersForBayesOut))
}
