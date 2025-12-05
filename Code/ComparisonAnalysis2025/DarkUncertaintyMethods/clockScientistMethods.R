# ---------------------------------------------------------
# 1. Setup / Dummy Data
# ---------------------------------------------------------
# Replace these with your actual data
# x: The measurements (e.g., frequency ratios)
# sigma: The statistical uncertainty for each day
# measurements <- df_AlSr$offset#c(1.000005, 1.000012, 0.999990, 1.000025, 0.999980)
# uncertainties <- df_AlSr$statistical_unc#c(0.000005, 0.000005, 0.000004, 0.000006, 0.000005)

# measurements <- df_AlYb$offset#c(1.000005, 1.000012, 0.999990, 1.000025, 0.999980)
# uncertainties <- df_AlYb$statistical_unc#c(0.000005, 0.000005, 0.000004, 0.000006, 0.000005)

measurements <- df_YbSr$offset#c(1.000005, 1.000012, 0.999990, 1.000025, 0.999980)
uncertainties <- df_YbSr$statistical_unc#c(0.000005, 0.000005, 0.000004, 0.000006, 0.000005)

# Calculate N (number of measurements) and Degrees of Freedom (dof)
N <- length(measurements)
dof <- N - 1

# ---------------------------------------------------------
# Part 1: Standard Weighted Mean & Initial Chi-Squared
# ---------------------------------------------------------

# Weights (w_i = 1 / sigma_i^2)
weights_initial <- 1 / (uncertainties^2)

# Weighted Mean
mean_val <- sum(measurements * weights_initial) / sum(weights_initial)

# Chi-Squared (Standard)
chisq <- sum(((measurements - mean_val)^2) / (uncertainties^2))


# Reduced Chi-Squared
chisq_red <- chisq / dof

print(paste("Initial Reduced Chi-Squared:", round(chisq_red, 4)))

# ---------------------------------------------------------
# Part 2: Method A - Multiplicative (Birge Ratio)
# ---------------------------------------------------------
# "Multiplying each day's statistical uncertainty by the Birge ratio"

if (chisq_red > 1) {
  # Calculate Birge Ratio
  birge_ratio <- sqrt(chisq_red)
  
  # Scale uncertainties
  sigma_birge <- uncertainties * birge_ratio
  
  # Note: The new Chi-sq red using sigma_birge will be exactly 1.0
  print(paste("Method 1 (Birge Ratio):", round(birge_ratio, 4)))
  print("New Uncertainties (Birge):")
  print(sigma_birge)
} else {
  print("Data is not overscattered (Chi_red <= 1). Birge scaling not applied.")
  sigma_birge <- uncertainties
}

# ---------------------------------------------------------
# Part 3: Method B - Additive Dark Uncertainty (Dorscher et al.)
# ---------------------------------------------------------
# "Adding a factor xi_a in quadrature... such that Chi_red = 1"

# We need to find xi_a such that the reduced chi-squared equals 1.
# Since the weighted mean ITSELF depends on xi_a (because weights change),
# this is a root-finding problem, not a simple algebraic solution.

target_chisq_function <- function(xi_a, x, sigma, dof) {
  # Add systematic in quadrature: sigma_total^2 = sigma^2 + xi_a^2
  sigma_total_sq <- sigma^2 + xi_a^2
  w_new <- 1 / sigma_total_sq
  
  # Recalculate mean with new weights
  mu_new <- sum(x * w_new) / sum(w_new)
  
  # Calculate new reduced chi-squared
  chisq_new <- sum(((x - mu_new)^2) / sigma_total_sq)
  chisq_red_new <- chisq_new / dof
  
  # We want (chisq_red_new - 1) to be 0
  return(chisq_red_new - 1)
}

if (chisq_red > 1) {
  # Use uniroot to find the xi_a that makes the function zero
  # Interval: search from 0 up to a large number (e.g., max deviation)
  max_search <- max(abs(measurements - mean_val)) * 10 
  
  solution <- uniroot(target_chisq_function, 
                      interval = c(0, max_search), 
                      x = measurements, 
                      sigma = uncertainties, 
                      dof = dof)
  
  xi_a <- solution$root
  
  # Calculate final total uncertainties for this method
  sigma_additive <- sqrt(uncertainties^2 + xi_a^2)
  
  print(paste("Method 2 (Additive Constant xi_a):", format(xi_a, scientific = TRUE)))
  print("New Uncertainties (Additive):")
  print(sigma_additive)
  
} else {
  print("Data is not overscattered. Additive systematic set to 0.")
  xi_a <- 0
  sigma_additive <- uncertainties
}

# ---------------------------------------------------------
# Summary Comparison
# ---------------------------------------------------------
results <- data.frame(
  Original_Sigma = uncertainties,
  Birge_Sigma    = sigma_birge,
  Additive_Sigma = sigma_additive
)

print("Comparison of Uncertainties:")
print(results)

xi_a
# mason: darkuncs=[.9e-18, 1e-18, 2.9e-18] for the AlYb, AlSr, and YbSr ratios respectively.
# me: darkuncs=[1.48 e-18, 1.72 e-18, 2.80 e-18] be the list of dark uncertainties for the AlYb, AlSr, and YbSr ratios respectively.

target_chisq_function(2.9,df_YbSr$offset,df_YbSr$statistical_unc,dim(df_YbSr)[1]-1)
target_chisq_function(2.8,df_YbSr$offset,df_YbSr$statistical_unc,dim(df_YbSr)[1]-1)


target_chisq_function(1,df_AlSr$offset,df_AlSr$statistical_unc,dim(df_AlSr)[1]-1)
target_chisq_function(1.72,df_AlSr$offset,df_AlSr$statistical_unc,dim(df_AlSr)[1]-1)


target_chisq_function(.9,df_AlYb$offset,df_AlYb$statistical_unc,dim(df_AlYb)[1]-1)
target_chisq_function(1.48,df_AlYb$offset,df_AlYb$statistical_unc,dim(df_AlYb)[1]-1)




# Weights (w_i = 1 / sigma_i^2)
weights_initial <- 1 / (uncertainties^2+)

# Weighted Mean
mean_val <- sum(measurements * weights_initial) / sum(weights_initial)

# Chi-Squared (Standard)
chisq <- sum(((measurements - mean_val)^2) / (uncertainties^2))


# Reduced Chi-Squared
chisq_red <- chisq / dof

