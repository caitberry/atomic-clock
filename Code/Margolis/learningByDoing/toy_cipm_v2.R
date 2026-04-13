library(rstan)

# ---------------------------------------------------------
# 1. THE RAW DATA (Pulled from the CIPM 2021 dataset)
# ---------------------------------------------------------
# The 2017 reference anchor frequencies (nu_0) for Hg, Yb, Sr 
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

# Topography: Define what each measurement represents
# Index 4 represents Cesium
num_idx <- c(1, 2, 3, 1, 2, 1)
den_idx <- c(4, 4, 4, 2, 3, 3)

# ---------------------------------------------------------
# 2. TRANSFORMING THE DATA (Scaling for numerical stability)

# UNDERSTAND THIS PIECE NEXT!!!!!!!!!!!!
# ---------------------------------------------------------
Q_scaled <- numeric(6)
U_scaled <- numeric(6)

for (i in 1:6) {
  num_0 <- nu_0[num_idx[i]]
  
  if (den_idx[i] == 4) {
    # It's an absolute frequency measurement vs Cesium
    Q_scaled[i] <- ((Q[i] - num_0) / num_0) * 1e15
    U_scaled[i] <- (U[i] / num_0) * 1e15
  } else {
    # It's an optical-to-optical ratio
    den_0 <- nu_0[den_idx[i]]
    R_0 <- num_0 / den_0
    
    Q_scaled[i] <- ((Q[i] - R_0) / R_0) * 1e15
    U_scaled[i] <- (U[i] / R_0) * 1e15
  }
}


# Create the base diagonal covariance matrix
V_scaled <- diag(U_scaled^2)

# Add the off-diagonal covariance terms for our 6 measurements
# Index mapping: 1=q7, 2=q24, 3=q47, 4=q79, 5=q81, 6=q59

# r(q7, q47) = 0.438 -> indices 1 and 3
V_scaled[1, 3] <- 0.438 * U_scaled[1] * U_scaled[3]
V_scaled[3, 1] <- V_scaled[1, 3]

# r(q24, q81) = 0.088 -> indices 2 and 5
V_scaled[2, 5] <- 0.088 * U_scaled[2] * U_scaled[5]
V_scaled[5, 2] <- V_scaled[2, 5]

# r(q59, q79) = 0.826 -> indices 6 and 4
V_scaled[4, 6] <- 0.826 * U_scaled[4] * U_scaled[6]
V_scaled[6, 4] <- V_scaled[4, 6]

# ---------------------------------------------------------
# 3. RUN STAN AND RECONSTRUCT
# ---------------------------------------------------------
stan_data <- list(
  Q_scaled = Q_scaled,
  V_scaled = V_scaled,
  num_idx  = num_idx,
  den_idx  = den_idx
)

fit <- stan(
  file = "C:/Users/aak3/Documents/atomic-clock/Code/Margolis/toy_cipm.stan",
  data = stan_data,
  iter = 2000,
  warmup = 500,
  chains = 4
)

# Extract the fractional offsets Stan found
x_mean <- colMeans(rstan::extract(fit)$x)

# Un-scale them to get our optimized absolute frequencies
optimized_Hg <- nu_0[1] * (1 + x_mean[1] * 1e-15)
optimized_Yb <- nu_0[2] * (1 + x_mean[2] * 1e-15)
optimized_Sr <- nu_0[3] * (1 + x_mean[3] * 1e-15)

cat("\n--- TOY MODEL OPTIMIZED FREQUENCIES ---\n")
cat("Hg:", format(optimized_Hg, digits=16), "Hz\n")
cat("Yb:", format(optimized_Yb, digits=16), "Hz\n")
cat("Sr:", format(optimized_Sr, digits=16), "Hz\n")

