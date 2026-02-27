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

N <- 100  # length of the time series
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
i=1
U <- taperMat*exp(-1i*2*pi*freq[i]*t.vec)/sqrt(myK)
U_star <- Conj(t(U))
A_f <- U %*% U_star

A<-A_f
  
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

library(ggplot2)
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
