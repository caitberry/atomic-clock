# =====================================================================
# 1. DATA PREPARATION & SCALING
# =====================================================================
# The 2017 reference anchor frequencies (nu_0)
# CHECK THESE!!!!!!!!!
nu_0 <- c(
  1128575290808154, # 1: Hg (nu3)
  518295836590864,  # 2: Yb (nu8)
  429228004229873   # 3: Sr (nu12)
)

# Our 6 Measurements (Q)
# q7	nu3		1128575290808154.62	0.41	[Tyumenev2016]
# q24	nu8		518295836590863.59	0.31	[Pizzocaro2017]
# q47	nu12		429228004229872.92	0.12	[Lodewyck2016]	
# q79	nu3_over_nu8	2.17747319413456507	1.92E-16	[Ohmae2020]
# q81	nu8_over_nu12	1.20750703934333841	3.40E-16	[Grotti2018]
# q59	nu3_over_nu12	2.62931420989890960	2.2e-16	[Yamanake2015]


Q <- c(
  1128575290808154.62, # 1. Absolute Hg (q7)
  518295836590863.59,  # 2. Absolute Yb (q24)
  429228004229872.92,  # 3. Absolute Sr (q47)
  2.17747319413456507,   # 4. Ratio Hg / Yb (q79)
  1.20750703934333841,   # 5. Ratio Yb / Sr (q81)
  2.62931420989890960    # 6. Ratio Hg / Sr (q59)
)

# Their absolute uncertainties (U)
# Corrected absolute uncertainties (U)
U <- c(
  0.41,       # 1. Hg absolute unc
  0.31,       # 2. Yb absolute unc
  0.12,       # 3. Sr absolute unc
  1.92e-16,   # 4. Hg/Yb unc 
  3.40e-16,   # 5. Yb/Sr unc 
  2.20e-16    # 6. Hg/Sr unc 
)

# Reference ratios (R_0) for measurements 4, 5, and 6
R_0 <- c(nu_0[1]/nu_0[2], nu_0[2]/nu_0[3], nu_0[1]/nu_0[3])

# Calculate scaled fractional offsets (Y matrix equivalent)
Y <- c(
  ((Q[1] - nu_0[1]) / nu_0[1]) * 1e15,
  ((Q[2] - nu_0[2]) / nu_0[2]) * 1e15,
  ((Q[3] - nu_0[3]) / nu_0[3]) * 1e15,
  ((Q[4] - R_0[1]) / R_0[1]) * 1e15,
  ((Q[5] - R_0[2]) / R_0[2]) * 1e15,
  ((Q[6] - R_0[3]) / R_0[3]) * 1e15
)

# Calculate scaled fractional uncertainties
u_scaled <- c(
  (U[1] / nu_0[1]) * 1e15,
  (U[2] / nu_0[2]) * 1e15,
  (U[3] / nu_0[3]) * 1e15,
  (U[4] / R_0[1]) * 1e15,
  (U[5] / R_0[2]) * 1e15,
  (U[6] / R_0[3]) * 1e15
)

# Construct Covariance Matrix (V)
V <- diag(u_scaled^2)

# Insert the 3 explicit off-diagonal correlations
V[1, 3] <- 0.438 * u_scaled[1] * u_scaled[3]; V[3, 1] <- V[1, 3]
V[2, 5] <- 0.088 * u_scaled[2] * u_scaled[5]; V[5, 2] <- V[2, 5]
V[4, 6] <- 0.826 * u_scaled[4] * u_scaled[6]; V[6, 4] <- V[4, 6]

# =====================================================================
# 2. LEAST-SQUARES ADJUSTMENT
# =====================================================================
# The Jacobian Design Matrix (A) [N x M]
# Columns represent the partial derivatives for x_Hg, x_Yb, x_Sr
A <- matrix(c(
  1,  0,  0,  # 1. d(Hg)/dx
  0,  1,  0,  # 2. d(Yb)/dx
  0,  0,  1,  # 3. d(Sr)/dx
  1, -1,  0,  # 4. d(Hg/Yb)/dx
  0,  1, -1,  # 5. d(Yb/Sr)/dx
  1,  0, -1   # 6. d(Hg/Sr)/dx
), nrow = 6, byrow = TRUE)

# Calculate Least-Squares Solution: X = (A^T V^-1 A)^-1 A^T V^-1 Y
V_inv <- solve(V)
At_Vinv <- t(A) %*% V_inv

# Covariance matrix of the adjusted variables
cov_X <- solve(At_Vinv %*% A) 

# Optimized scaled fractional offsets
X_hat <- cov_X %*% At_Vinv %*% Y 

# =====================================================================
# 3. SELF-CONSISTENCY CHECKS
# =====================================================================
# Calculate residuals and minimum S (Chi-squared)
residuals <- Y - (A %*% X_hat)
chi_squared <- as.numeric(t(residuals) %*% V_inv %*% residuals)

# Birge ratio: sqrt(chi_squared / degrees of freedom)
N <- 6
M <- 3
birge_ratio <- sqrt(chi_squared / (N - M))

# =====================================================================
# 4. HIGH-PRECISION RECONSTRUCTION
# =====================================================================
# Reconstruct absolute frequencies: nu = nu_0 * (1 + X_hat * 1e-15)
optimized_freqs <- nu_0 * (1 + as.vector(X_hat) * 1e-15)
optimized_unc <- nu_0 * (sqrt(diag(cov_X)) * 1e-15)

cat("--- LEAST SQUARES OPTIMIZED FREQUENCIES ---\n")
cat("Hg:", sprintf("%.3f", optimized_freqs[1]), "Hz  (Unc:", sprintf("%.2e", optimized_unc[1]), "Hz)\n")
cat("Yb:", sprintf("%.3f", optimized_freqs[2]), "Hz  (Unc:", sprintf("%.2e", optimized_unc[2]), "Hz)\n")
cat("Sr:", sprintf("%.3f", optimized_freqs[3]), "Hz  (Unc:", sprintf("%.2e", optimized_unc[3]), "Hz)\n\n")

cat("--- CONSISTENCY METRICS ---\n")
cat("Chi-squared:", round(chi_squared, 3), "\n")
cat("Birge Ratio:", round(birge_ratio, 3), "\n")

