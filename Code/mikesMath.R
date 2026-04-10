# cleaned up final conclusions: 
# basic examples worked in more detail below

#################################################################
# Final Conclusion: Quadratic Form Distribution via Weighted Chi-Squares
#################################################################

library(tidyr)
library(ggplot2)
library(dplyr)
library(MASS) # Required for mvrnorm

# --- 1. Setup & Spectral Matrix Generation ---
N <- 50  
ar_coeffs <- c(0.5, -0.4, 0.3, -0.2)
mu <- rep(0, N)

# Create Covariance Matrix (Sigma) from AR(4)
s_acf <- stats::ARMAacf(ar = ar_coeffs, ma = numeric(0), lag.max = N-1)
Sigma <- stats::toeplitz(s_acf)

# Create Matrix A (Complex Demodulation / Multitaper Kernel)
# (Simplified example based on your Fourier frequency i=2)
freq_i <- 2/N 
t_vec  <- 1:N
# Complex projection matrix approach
U <- exp(-1i * 2 * pi * freq_i * t_vec) / sqrt(N)
A <- Re(U %*% Conj(t(U))) 

# --- 2. Eigen-Decomposition for Weights ---
# The weights (lambdas) of the quadratic form are the eigenvalues of A %*% Sigma
# We use the symmetric version: Sigma^(1/2) %*% A %*% Sigma^(1/2)
eig_Sig <- eigen(Sigma)
Sigma_half <- eig_Sig$vectors %*% diag(sqrt(eig_Sig$values)) %*% t(eig_Sig$vectors)

M <- Sigma_half %*% A %*% Sigma_half
lambdas <- eigen(M)$values
lambdas <- lambdas[lambdas > 1e-10] # Filter numerical noise

# --- 3. Calculate Approximation Parameters ---
s1 <- sum(lambdas)
s2 <- sum(lambdas^2)
s3 <- sum(lambdas^3)

# Satterthwaite (2-moment Gamma)
alpha_s <- s2 / s1
nu_s    <- s1^2 / s2

# Liu-Tang-Zhang (3-moment Shifted Gamma)
alpha_ltz <- s3 / s2
nu_ltz    <- s2^3 / s3^2
delta     <- s1 - (s2^2 / s3)

# --- 4. Simulation Comparison ---
n_sim <- 10000

# A. Ground Truth: Direct Quadratic Form
sim_direct <- replicate(n_sim, {
  x <- mvrnorm(1, mu, Sigma)
  as.numeric(t(x) %*% A %*% x)
})

# B. Theory: Weighted Sum of Independent Chi-Squares
# (Since mu=0, NCP is 0)
sim_weighted <- matrix(rchisq(n_sim * length(lambdas), 1), ncol=length(lambdas)) %*% lambdas

# C. Satterthwaite Approximation (Scaled Chi-Square)
sim_satt <- alpha_s * rchisq(n_sim, df = nu_s)

# D. Liu-Tang-Zhang (Shifted & Scaled Chi-Square)
sim_ltz  <- alpha_ltz * rchisq(n_sim, df = nu_ltz) + delta

# --- 5. Final Visualization ---
simres <- data.frame(
  values = c(as.vector(sim_direct), as.vector(sim_weighted), 
             as.vector(sim_satt), as.vector(sim_ltz)),
  type = rep(c("1. Direct (x'Ax)", "2. Weighted Sum", 
               "3. Satterthwaite", "4. Liu-Tang-Zhang"), each = n_sim)
)

ggplot(simres, aes(x = values, color = type)) +
  geom_density(linewidth = 1) +
  labs(title = "Quadratic Form Distribution Accuracy",
       subtitle = "Comparing Direct Simulation vs. Weighted Sum vs. Gamma Approximations",
       x = "Value", y = "Density", color = "Method") +
  theme_minimal() +
  scale_color_brewer(palette = "Set1")

# Numerical check of the means
cat("Theoretical Mean:", s1, "\n")
cat("Direct Mean:     ", mean(sim_direct), "\n")
cat("LTZ Mean:        ", mean(sim_ltz), "\n")







# # Load MASS for multivariate normal generation
# library(MASS)
# 
# 
# #### MIKE'S #1
# 
# # 1. Setup Interesting Data (Correlated)
# n <- 3
# mu <- c(10, 10, 10) # Constant mean vector
# 
# # Create a Covariance Matrix (Sigma)
# # High positive correlation (0.9) between variables
# Sigma <- matrix(c(1.0, 0.9, 0.9,
#                   0.9, 1.0, 0.9,
#                   0.9, 0.9, 1.0), nrow=3, ncol=3)
# 
# # Generate one random vector x from N(mu, Sigma)
# set.seed(123)
# x <- mvrnorm(n = 1, mu = mu, Sigma = Sigma)
# 
# # 2. Construct the Matrix A (Sample Variance Kernel)
# I <- diag(n)
# J <- matrix(1, nrow=n, ncol=n)
# C <- I - (1/n) * J
# A <- (1 / (n - 1)) * C
# 
# # 3. Compute the Quadratic Form
# # This calculates the scalar sample variance of our vector
# variance_quadratic <- t(x) %*% A %*% x
# variance_quadratic
# # 4. Compare with standard R function
# variance_standard <- var(x)
# variance_standard
# 
# # 5. The "Interesting" Part: The Trace Theorem
# # Theory says Expected Value E[S^2] = Trace(A * Sigma)
# # Let's calculate the theoretical expectation
# theoretical_expectation <- sum(diag(A %*% Sigma))
# 
# 
# #### MIKE'S #2
# 
# # 2. Perform Spectral Decomposition
# # The eigen() function returns values (lambda) and vectors (V)
# eigen_decomp <- eigen(Sigma)
# 
# # Extract components
# lambdas <- eigen_decomp$values   # Eigenvalues (vector)
# V       <- eigen_decomp$vectors  # Eigenvectors (matrix)
# 
# 
# 
# # 3. Construct the Matrix Square Root: Sigma^(1/2)
# # Formula: V * diag(sqrt(lambda)) * V'
# Lambda_sqrt <- diag(sqrt(lambdas))
# Sigma_sqrt  <- V %*% Lambda_sqrt %*% t(V)
# 
# # 4. Construct the Inverse Square Root: Sigma^(-1/2)
# # Formula: V * diag(1/sqrt(lambda)) * V'
# # This is the "Whitening Matrix" used to standardize data
# Lambda_inv_sqrt <- diag(1 / sqrt(lambdas))
# Sigma_inv_sqrt  <- V %*% Lambda_inv_sqrt %*% t(V)
# 
# # 5. Verification
# # Check if Sigma^(1/2) * Sigma^(1/2) equals the original Sigma
# Sigma_reconstructed <- Sigma_sqrt %*% Sigma_sqrt
# Sigma_reconstructed
# 
# 
# 
# #### MIKE'S #4
# 
# # 3. Calculate the Target Matrix M = Sigma^(1/2) * A * Sigma^(1/2)
# # This is the matrix that determines the distribution of the quadratic form
# M <- Sigma_sqrt %*% A %*% Sigma_sqrt
# 
# # 4. Output Results
# cat("Matrix M (Sigma^1/2 * A * Sigma^1/2):\n")
# print(round(M, 4))
# 
# # 5. Check Eigenvalues of M
# # The distribution of the quadratic form is a weighted sum of Chi-Squares:
# # sum(lambda_i * Chisq_1). If all non-zero eigenvalues are 1, it's a pure Chi-Square.
# eigen_M <- eigen(M)$values
# cat("\nEigenvalues of M (Weights of the Chi-Squares):\n")
# print(round(eigen_M, 5))
# 
# # Comparison: Check eigenvalues of A * Sigma
# # Theorem: The non-zero eigenvalues of (Sigma^1/2 A Sigma^1/2) are the same as (A Sigma)
# eigen_ASigma <- eigen(A %*% Sigma)$values
# cat("Eigenvalues of A * Sigma:\n")
# print(round(eigen_ASigma, 5))
# 
# 
# 
# 
# 
# 
# 
# 
# 



library(MASS) # For mvrnorm
library(ggplot2)

# for sample variance


# 1. Setup Parameters
n <- 3
# mu <- c(10, 10, 10) # Non-zero mean (crucial for Non-Centrality)
mu <- c(0, 0, 0) # Non-zero mean (crucial for Non-Centrality)

# Highly correlated covariance
Sigma <- matrix(c(1.0, 0.9, 0.9,
                  0.9, 1.0, 0.9,
                  0.9, 0.9, 1.0), nrow=3, ncol=3)

# Matrix A (Sample Variance Kernel)
I <- diag(n)
J <- matrix(1, n, n)
A <- (1 / (n - 1)) * (I - (1/n) * J)

# 2. Spectral Decomposition Helpers
# Get Sigma^(1/2) and Sigma^(-1/2)
eig_Sig <- eigen(Sigma)
V_Sig   <- eig_Sig$vectors
L_Sig   <- eig_Sig$values

Sigma_half     <- V_Sig %*% diag(sqrt(L_Sig)) %*% t(V_Sig)
Sigma_inv_half <- V_Sig %*% diag(1/sqrt(L_Sig)) %*% t(V_Sig)

# 3. Decompose the Core Matrix
# M = Sigma^(1/2) * A * Sigma^(1/2)
M <- Sigma_half %*% A %*% Sigma_half

eig_M   <- eigen(M)
lambdas <- eig_M$values   # These are the weights
P       <- eig_M$vectors  # Orthogonal rotation matrix

# 4. Calculate Non-Centrality Parameters (nu)
# Formula from image: nu = P^T * Sigma^(-1/2) * mu
nu_vec <- t(P) %*% Sigma_inv_half %*% mu

# The Non-Centrality Parameter (ncp) for Chi-Sq is nu^2
ncp_vec <- nu_vec^2

# 5. Simulation Verification
# We will generate 10,000 values using BOTH methods and compare means

n_sim <- 100000

# Method A: Direct Quadratic Form
# Generate X ~ N(mu, Sigma) and compute x'Ax
sim_direct <- replicate(n_sim, {
  x <- mvrnorm(1, mu, Sigma)
  as.numeric(t(x) %*% A %*% x)
})

# Method B: Weighted Sum of Non-Central Chi-Squares
# Sum of (lambda_i * Chisq(df=1, ncp=nu_i^2))
sim_weighted <- replicate(n_sim, {
  # Generate 3 independent chi-squares with different NCPs
  chisqs <- rchisq(n, df=1, ncp=ncp_vec)
  # Weighted sum
  sum(lambdas * chisqs)
})

# 6. Output Results
cat("--- 1. Decomposition Parameters ---\n")
cat("Weights (Eigenvalues of M):       ", round(lambdas, 4), "\n")
cat("Non-Centrality Params (nu^2):     ", round(as.numeric(ncp_vec), 4), "\n")

cat("\n--- 2. Simulation Comparison (Mean of Y) ---\n")
cat("Direct Quadratic Form Mean:       ", mean(sim_direct), "\n")
cat("Weighted Chi-Square Sum Mean:     ", mean(sim_weighted), "\n")

# Analytical Mean Check: Sum(lambda * (1 + nu^2))
analytical_mean <- sum(lambdas * (1 + ncp_vec))
cat("Analytical Expected Value:        ", analytical_mean, "\n")


simres=data.frame(values = c(sim_weighted, sim_direct),
                  type=rep(c("sim_weighted","sim_direct"),each=n_sim))

ggplot(simres,aes(values,col=type))+
  geom_density()
  

### compare with approx chi sqred
weights=lambdas
s1 <- sum(weights)
s2 <- sum(weights^2)
s3 <- sum(weights^3)

# 2. Satterthwaite Parameters (2-moment)

alpha_s <- s2 / s1
nu_s    <- s1^2 / s2

# 3. Liu-Tang-Zhang Parameters (3-moment)

# Formulas derived to match mean, variance, and skewness

alpha_ltz <- s3 / s2
nu_ltz    <- s2^3 / s3^2
delta     <- s1 - (s2^2 / s3)

# 4. Simulations

actual_Q <- matrix(rchisq(n_sim * length(weights), 1), ncol=length(weights)) %*% weights
sim_satt <- alpha_s * rchisq(n_sim, df = nu_s)
sim_ltz  <- alpha_ltz * rchisq(n_sim, df = nu_ltz) + delta


simres=data.frame(values = c(sim_weighted, sim_direct,sim_satt,sim_ltz),
                  type=rep(c("sim_weighted","sim_direct","sim_satt","sim_ltz"),each=n_sim))

ggplot(simres,aes(values,col=type))+
  geom_density()
  
ggplot(simres,aes(values,col=type))+
  geom_density()+
  coord_cartesian(ylim = c(7,8.5),xlim = c(0,.1))
  
###################################
###################################
### For spectrum
###################################
###################################


library(tidyr)
library(ggplot2)
library(dplyr)

mtse=modules::use("Functions.R")
source("Functions_SpectrumCovariance.R")

######################################
# 1. Setup Parameters
######################################

N <- 50  # length of the time series
# mu <- rep(0,n)
ar_coeffs <- c(0.5, -0.4, 0.3, -0.2)  # AR(4) coefficients

mu=rep(0,N)

t.vec <- 1:N

# # Simulate AR(4) process
x.t <- arima.sim(model = list(ar = ar_coeffs), n = N)
x.t=x.t-mean(x.t)

max.lag.acf=9
s_acf <- stats::ARMAacf(ar = ar_coeffs, ma = numeric(0), lag.max = max.lag.acf)#############REMOVE

# Create a Toeplitz matrix from the autocorrelation values
N.noNA=sum(!is.na(t.vec))

R_mat <- matrix(0, nrow = N.noNA, ncol = N.noNA)
R_mat <- stats::toeplitz(c(s_acf, rep(0, N.noNA - max.lag.acf - 1))) 

Sigma <- R_mat

# Matrix A (Sample Variance Kernel)
I <- diag(N)
J <- matrix(1, N, N)

myK=10
W1factor=12
V.mat3 <- mtse$get_tapers(t.vec, W = W1factor/N, K = myK)
taperMat=V.mat3$tapers
V.mat3$e.values

N.fourier=N/2+1

max.lag.acf=N-1

freq <- seq(0, 0.5, length.out = N.fourier)

#### do for one frequency
i=2
U <- taperMat*exp(-1i*2*pi*freq[i]*t.vec)/sqrt(myK)
U_star <- Conj(t(U))
A_f <- U %*% U_star

sum(Im(A_f))#hopefully all 0 so i can just take the real part below

A<-Re(A_f)
  
# 2. Spectral Decomposition Helpers
# Get Sigma^(1/2) and Sigma^(-1/2)
eig_Sig <- eigen(Sigma)
V_Sig   <- eig_Sig$vectors
L_Sig   <- eig_Sig$values

Sigma_half     <- V_Sig %*% diag(sqrt(L_Sig)) %*% t(V_Sig)
Sigma_inv_half <- V_Sig %*% diag(1/sqrt(L_Sig)) %*% t(V_Sig)

# 3. Decompose the Core Matrix
# M = Sigma^(1/2) * A * Sigma^(1/2)
M <- Sigma_half %*% A %*% Sigma_half

eig_M   <- eigen(M)
lambdas <- eig_M$values   # These are the weights
P       <- eig_M$vectors  # Orthogonal rotation matrix

# 4. Calculate Non-Centrality Parameters (nu)
# Formula from image: nu = P^T * Sigma^(-1/2) * mu
nu_vec <- t(P) %*% Sigma_inv_half %*% mu

# The Non-Centrality Parameter (ncp) for Chi-Sq is nu^2
ncp_vec <- nu_vec^2

# 5. Simulation Verification
# We will generate 10,000 values using BOTH methods and compare means

n_sim <- 5000

# Method A: Direct Quadratic Form
# Generate X ~ N(mu, Sigma) and compute x'Ax
sim_direct <- replicate(n_sim, {
  x <- mvrnorm(1, mu, Sigma)
  as.numeric(t(x) %*% A %*% x)
})

# Method B: Weighted Sum of Non-Central Chi-Squares
# Sum of (lambda_i * Chisq(df=1, ncp=nu_i^2))
sim_weighted <- replicate(n_sim, {
  # Generate independent chi-squares with different NCPs
  chisqs <- rchisq(N, df=1, ncp=ncp_vec)
  # Weighted sum
  sum(lambdas * chisqs)
})

# 6. Output Results
cat("--- 1. Decomposition Parameters ---\n")
cat("Weights (Eigenvalues of M):       ", round(lambdas, 4), "\n")
cat("Non-Centrality Params (nu^2):     ", round(as.numeric(ncp_vec), 4), "\n")

cat("\n--- 2. Simulation Comparison (Mean of Y) ---\n")
cat("Direct Quadratic Form Mean:       ", mean(sim_direct), "\n")
cat("Weighted Chi-Square Sum Mean:     ", mean(sim_weighted), "\n")

# Analytical Mean Check: Sum(lambda * (1 + nu^2))
analytical_mean <- sum(lambdas * (1 + ncp_vec))
cat("Analytical Expected Value:        ", analytical_mean, "\n")



### compare with approx chi sqred
weights=lambdas
s1 <- sum(weights)
s2 <- sum(weights^2)
s3 <- sum(weights^3)

# 2. Satterthwaite Parameters (2-moment)

alpha_s <- s2 / s1
nu_s    <- s1^2 / s2

# 3. Liu-Tang-Zhang Parameters (3-moment)

# Formulas derived to match mean, variance, and skewness

alpha_ltz <- s3 / s2
nu_ltz    <- s2^3 / s3^2
delta     <- s1 - (s2^2 / s3)

# 4. Simulations

actual_Q <- matrix(rchisq(n_sim * length(weights), 1), ncol=length(weights)) %*% weights
sim_satt <- alpha_s * rchisq(n_sim, df = nu_s)
sim_ltz  <- alpha_ltz * rchisq(n_sim, df = nu_ltz) + delta


simres=data.frame(values = c(sim_weighted, sim_direct,sim_satt,sim_ltz),
                  type=rep(c("sim_weighted","sim_direct","sim_satt","sim_ltz"),each=n_sim))

ggplot(simres,aes(values,col=type))+
  geom_density()










# 
# #################
# #################
# # with the spectrum
# #################
# #################
# 
# library(tidyr)
# library(ggplot2)
# library(dplyr)
# 
# mtse=modules::use("Functions.R")
# source("Functions_SpectrumCovariance.R")
# 
# N <- 2000  # length of the time series
# ar_coeffs <- c(0.5, -0.4, 0.3, -0.2)  # AR(4) coefficients
# 
# t.vec <- 1:N
# 
# # # Simulate AR(4) process
# x.t <- arima.sim(model = list(ar = ar_coeffs), n = N)
# x.t=x.t-mean(x.t)
# # 
# # omitted=c(100:300)
# # x.t[omitted]=NA
# # t.vec[which(is.na(x.t))] <- NA
# 
# myK=10
# W1factor=12
# V.mat3 <- mtse$get_tapers(t.vec, W = W1factor/N, K = myK) 
# taperMat=V.mat3$tapers
# V.mat3$e.values
# 
# N.fourier=N/2+1
# 
# max.lag.acf=N-1
# 
# freq <- seq(0, 0.5, length.out = N.fourier)
# 
# spec_est=numeric(N.fourier)
# for(i in 1:N.fourier){
#   if(i %% 100 == 0){print(paste(i," of ",N.fourier))}
#     # U (V) is the "taper matrix" from Bronez,
#     # exp(1i*2*pi*freq[i]*t.vec) is the Chave factor that has to be multiplied back in
#     
#     # U_star <- t(taperMat*exp(1i*2*pi*freq[i]*t.vec))/sqrt(myK)
#   
#     U <- taperMat*exp(-1i*2*pi*freq[i]*t.vec)/sqrt(myK)
#     U_star <- Conj(t(U))
#     # # 1. Construct the 'Middle Piece' matrix (N x N)
#     A_f <- U %*% U_star
#     # 
#     # # 2. Compute the Spectral Estimate (Scalar)
#     # # x is column vector of data (N x 1)
#     # # We use Re() because the result is theoretically real, 
#     # # but may have tiny imaginary artifacts due to precision.
#     x_vec <- as.matrix(x.t)
#     spec_est[i] <- Re(t(x_vec) %*% A_f %*% x_vec)
#     
# }
# 
# MTSE_full <- mtse$MT_spectralEstimate_fft(x.t, taperMat)
# 
# specEstsCait1=data.frame(freq = MTSE_full$freqs, 
#            spectrum = MTSE_full$spectrum,
#            type="CaitFFT")
# 
# # MTSE_full <- MT_spectralEstimate(x.t, taperMat)
# # 
# # specEstsCait2=data.frame(freq = MTSE_full$freqs, 
# #                          spectrum = MTSE_full$spectrum,
# #                          type="CaitOLD")
# 
# 
# specEstsQuad=data.frame(freq = freq, 
#                         spectrum = spec_est,
#                         type="Quad")
# 
# # Calculate the denominator of the AR spectral formula
# # We evaluate the polynomial at each frequency point
# complex_sum <- 0
# for (k in 1:length(ar_coeffs)) {
#   complex_sum <- complex_sum + ar_coeffs[k] * exp(-1i * 2 * pi * freq * k)
# }
# 
# # The True PSD
# true_psd <- 1 / (Mod(1 - complex_sum)^2)
# 
# specTruth=data.frame(freq = freq, 
#                      spectrum = true_psd,
#                      type="Truth")
# 
# allSpec=bind_rows(specEstsCait1,specEstsQuad,specTruth)
# 
# 
# ggplot(allSpec,aes(freq,spectrum,col=type))+
#   geom_point()
# 
# a=allSpec %>% group_by(type) %>%
#   filter(type == "CaitFFT")
# 
# b=allSpec %>% group_by(type) %>%
#   filter(type == "Quad")
# b$spectrum-a$spectrum
# # # Visualization
# # plot(freq, 10 * log10(true_psd), type = "l", col = "red", lwd = 2,
# #      xlab = "Frequency", ylab = "Power/Freq (dB)",
# #      main = "True Theoretical Spectrum of AR(4)")
# # # points(freq, 10 * log10(spec_est / N))
# # # Double the power for one-sided spectral density
# # # points(freq, 10 * log10(2 * spec_est / N))
# # points(freq, 10 * log10(spec_est))
# # points(freq, 10 * log10(MTSE_full$spectrum))
# # 
# # 
# # # Standard periodogram for comparison
# # pgram <- spec.pgram(x.t, plot = FALSE)
# # lines(pgram$freq, 10 * log10(pgram$spec), col = "gray")


# 
# 
# #### is the sum of weighted chi squared a scaled chisquared?
# 
# # 1. Define weights and calculate parameters
# weights <- c(0.5, 2, 5, 10)
# n_sim <- 100000
# 
# # Theoretical Moments
# mean_Q <- sum(weights)
# var_Q  <- 2 * sum(weights^2)
# 
# # Calculate Satterthwaite Parameters
# alpha <- var_Q / (2 * mean_Q)       # Multiplicative factor
# nu    <- (2 * mean_Q^2) / var_Q     # Effective degrees of freedom
# 
# # 2. Generate the ACTUAL weighted sum (Q)
# # Each column is a ChiSq(1), multiplied by its corresponding weight
# samples <- matrix(rchisq(n_sim * length(weights), df = 1), ncol = length(weights))
# Q_actual <- samples %*% weights
# 
# # 3. Generate the RESCALED Chi-Squared (alpha * ChiSq(nu))
# Q_approx <- alpha * rchisq(n_sim, df = nu)
# 
# # 4. Compare Distributions Visually
# par(mfrow=c(1,2)) # Side-by-side plots
# 
# # Plot 1: Density Comparison
# plot(density(Q_actual), col="blue", lwd=2, main="Density Comparison", 
#      xlab="Value", ylim=c(0, max(density(Q_actual)$y)))
# lines(density(Q_approx), col="red", lty=2, lwd=2)
# legend("topright", legend=c("Actual Weighted Sum", "Rescaled Chi-Sq"), 
#        col=c("blue", "red"), lty=c(1,2), lwd=2)
# 
# # Plot 2: Q-Q Plot (If they match, points will follow the y=x line)
# qqplot(Q_actual, Q_approx, main="Q-Q Plot", 
#        xlab="Actual Weighted Sum Quantiles", 
#        ylab="Rescaled Chi-Sq Quantiles", pch=20, col=rgb(0,0,0,0.1))
# abline(0, 1, col="red", lwd=2)
# 
# 
# 
# 
# 
# ###########
# # 1. Define weights and calculate cumulants (s_j = sum of weights^j)
# weights <- c(0.1, 0.5, 2, 10) # Highly unequal weights show the difference better
# n_sim <- 10^6
# 
# s1 <- sum(weights)
# s2 <- sum(weights^2)
# s3 <- sum(weights^3)
# 
# # 2. Satterthwaite Parameters (2-moment)
# alpha_s <- s2 / s1
# nu_s    <- s1^2 / s2
# 
# # 3. Liu-Tang-Zhang Parameters (3-moment)
# # Formulas derived to match mean, variance, and skewness
# alpha_ltz <- s3 / s2
# nu_ltz    <- s2^3 / s3^2
# delta     <- s1 - (s2^2 / s3)
# 
# # 4. Simulations
# actual_Q <- matrix(rchisq(n_sim * length(weights), 1), ncol=length(weights)) %*% weights
# sim_satt <- alpha_s * rchisq(n_sim, df = nu_s)
# sim_ltz  <- alpha_ltz * rchisq(n_sim, df = nu_ltz) + delta
# 
# # 5. Compare the "Tail" (95th to 99.9th percentiles)
# probs <- c(0.95, 0.99, 0.999)
# results <- rbind(
#   Actual = quantile(actual_Q, probs),
#   Satterthwaite = quantile(sim_satt, probs),
#   LTZ = quantile(sim_ltz, probs)
# )
# 
# print(round(results, 3))
# 
# # 6. Visualization: Q-Q Plot comparison
# par(mfrow=c(1,2))
# qqplot(actual_Q, sim_satt, main="Satterthwaite Q-Q", pch=".", col="blue")
# abline(0,1, col="red")
# 
# qqplot(actual_Q, sim_ltz, main="Liu-Tang-Zhang Q-Q", pch=".", col="darkgreen")
# abline(0,1, col="red")
# 
# 
# 
# ###########################
# 
# 
# # 1. Define weights and calculate cumulants (s_j = sum of weights^j)
# weights <- c(0.1, 0.5, 2, 10) # Highly unequal weights show the difference better
# n_sim <- 10^6
# 
# s1 <- sum(weights)
# s2 <- sum(weights^2)
# s3 <- sum(weights^3)
# 
# # 2. Satterthwaite Parameters (2-moment)
# alpha_s <- s2 / s1
# nu_s    <- s1^2 / s2
# 
# # 3. Liu-Tang-Zhang Parameters (3-moment)
# # Formulas derived to match mean, variance, and skewness
# alpha_ltz <- s3 / s2
# nu_ltz    <- s2^3 / s3^2
# delta     <- s1 - (s2^2 / s3)
# 
# # 4. Simulations
# actual_Q <- matrix(rchisq(n_sim * length(weights), 1), ncol=length(weights)) %*% weights
# sim_satt <- alpha_s * rchisq(n_sim, df = nu_s)
# sim_ltz  <- alpha_ltz * rchisq(n_sim, df = nu_ltz) + delta
# 
# 
# library(ggplot2)
# library(tidyr)
# 
# df <- data.frame(
#   Actual = as.vector(actual_Q),
#   Satterthwaite = as.vector(sim_satt),
#   LTZ = as.vector(sim_ltz)
# )
# 
# df_long <- pivot_longer(df, cols = everything(), 
#                         names_to = "Method", 
#                         values_to = "Value")
# 
# ggplot(df_long, aes(x = Value, color = Method, fill = Method)) +
#   geom_density(alpha = 0.1, linewidth = 1) +
#   labs(title = "Comparison of Quadratic Form Approximations",
#        subtitle = paste("Weights:", paste(weights, collapse = ", ")),
#        x = "Value of Q",
#        y = "Density") +
#   theme_minimal() +
#   scale_color_manual(values = c("Actual" = "black", 
#                                 "Satterthwaite" = "#0072B2", 
#                                 "LTZ" = "#D55E00")) +
#   scale_fill_manual(values = c("Actual" = "black", 
#                                "Satterthwaite" = "#0072B2", 
#                                "LTZ" = "#D55E00"))
# 
# # 4. Optional: Zoom into the tail to see the LTZ advantage
# # ggplot(df_long, aes(x = Value, color = Method)) +
# #   geom_density(linewidth = 1) +
# #   coord_cartesian(xlim = c(quantile(actual_Q, 0.9), quantile(actual_Q, 0.999))) +
# #   labs(title = "Tail Comparison (90th - 99.9th Percentile)") +
# #   theme_minimal()
# 
# # 5. Compare the "Tail" (95th to 99.9th percentiles)
# probs <- c(0.95, 0.99, 0.999)
# results <- rbind(
#   Actual = quantile(actual_Q, probs),
#   Satterthwaite = quantile(sim_satt, probs),
#   LTZ = quantile(sim_ltz, probs)
# )
# 
# print(round(results, 3))
# 
# # 6. Visualization: Q-Q Plot comparison
# par(mfrow=c(1,2))
# qqplot(actual_Q, sim_satt, main="Satterthwaite Q-Q", pch=".", col="blue")
# abline(0,1, col="red")
# 
# qqplot(actual_Q, sim_ltz, main="Liu-Tang-Zhang Q-Q", pch=".", col="darkgreen")
# abline(0,1, col="red")
# 
# 
# 
# # Plot 1: Density Comparison
# plot(density(actual_Q), col="blue", lwd=2, main="Density Comparison", 
#      xlab="Value", ylim=c(0, max(density(Q_actual)$y)))
# lines(density(Q_approx), col="red", lty=2, lwd=2)
# legend("topright", legend=c("Actual Weighted Sum", "Rescaled Chi-Sq"), 
#        col=c("blue", "red"), lty=c(1,2), lwd=2)
# 
# 
# library("fitdistrplus")
# 
# # Fit a Gamma distribution to your actual weighted data
# fit_gamma <- fitdist(as.vector(actual_Q), "gamma")
# 
# # Extract the fitted parameters
# fitted_shape <- fit_gamma$estimate["shape"]
# fitted_rate  <- fit_gamma$estimate["rate"]
# 
# cat("Fitted Shape:", fitted_shape, "\n")
# cat("Fitted Scale:", 1/fitted_rate, "\n")
# 
# # Compare to Satterthwaite nu/2 and 2*alpha
# cat("Satterthwaite Shape (nu/2):", nu_s / 2, "\n")
# cat("Satterthwaite Scale (2*alpha):", 2 * alpha_s, "\n")