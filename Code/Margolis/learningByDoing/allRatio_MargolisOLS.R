library(Rmpfr)

# Set precision bits (128 bits gives ~38 decimal digits of precision)
prec <- 128

# =====================================================================
# 1. FILE PARSING & DATA PREPARATION
# =====================================================================
# Read the input file
file_path <- "Data/Margolis2024DataMetrologia/ClockInputData2020_7.dat"
lines <- readLines(file_path)

# Extract Number of Variables (M) and Measurements (N) dynamically or set manually
M <- 14
N <- 106

# Initialize a data frame to hold the measurement data
q_data <- data.frame(
  id = integer(), 
  type = character(), 
  val_str = character(), 
  unc_str = character(), 
  stringsAsFactors = FALSE
)

# Parse measurements (q1 to q106)
q_lines <- grep("^q\\d+\\s+", lines, value = TRUE)
for (line in q_lines) {
  parts <- unlist(strsplit(trimws(line), "\\s+"))
  id <- as.integer(sub("q", "", parts[1]))
  type <- parts[2]
  val_str <- parts[3]
  unc_str <- parts[4]
  
  q_data <- rbind(q_data, data.frame(
    id = id, 
    type = type, 
    val_str = val_str, 
    unc_str = unc_str, 
    stringsAsFactors = FALSE
  ))
}

# Ensure data is sorted by measurement ID
q_data <- q_data[order(q_data$id), ]

# Establish Reference Anchor Frequencies (nu_0)
# We will use the first absolute frequency measurement of each variable as its anchor. THIS ISN'T RIGHT
nu_0_str <- character(M)
for (i in 1:M) {
  target <- paste0("nu", i)
  match_row <- q_data[q_data$type == target, ][1, ]
  nu_0_str[i] <- match_row$val_str
}
nu_0 <- mpfr(nu_0_str, precBits = prec)

# Convert all measurements and uncertainties to mpfr
Q <- mpfr(q_data$val_str, precBits = prec)
U <- mpfr(q_data$unc_str, precBits = prec)

# =====================================================================
# 2. DESIGN MATRIX (A) & SCALING
# =====================================================================
A <- matrix(0, nrow = N, ncol = M)
R_0 <- mpfr(rep(NA, N), precBits = prec)

# Dynamically populate A and R_0 based on the 'type' column
for (i in 1:N) {
  type <- q_data$type[i]
  
  if (grepl("_over_", type)) {
    # It's a ratio measurement (e.g., "nu3_over_nu12")
    parts <- unlist(strsplit(type, "_over_"))
    num_idx <- as.integer(sub("nu", "", parts[1]))
    den_idx <- as.integer(sub("nu", "", parts[2]))
    
    R_0[i] <- nu_0[num_idx] / nu_0[den_idx]
    A[i, num_idx] <- 1
    A[i, den_idx] <- -1
  } else {
    # It's an absolute frequency measurement (e.g., "nu3")
    idx <- as.integer(sub("nu", "", type))
    
    R_0[i] <- nu_0[idx]
    A[i, idx] <- 1
  }
}

# Calculate scaled fractional offsets & uncertainties
scale_factor <- mpfr(1e15, precBits = prec)
Y_mpfr <- ((Q - R_0) / R_0) * scale_factor
u_scaled_mpfr <- (U / R_0) * scale_factor

# Cast back to standard double-precision for matrix algebra
Y <- as.numeric(Y_mpfr)
u_scaled <- as.numeric(u_scaled_mpfr)

# =====================================================================
# 3. COVARIANCE MATRIX (V) CONSTRUCTION
# =====================================================================
V <- diag(u_scaled^2)

# Extract and apply off-diagonal correlations
r_lines <- grep("^r\\(q\\d+,q\\d+\\)", lines, value = TRUE)
for (line in r_lines) {
  # RegEx to extract the two indices and the correlation coefficient
  matches <- regmatches(line, regexec("^r\\(q(\\d+),q(\\d+)\\)\\s+([-0-9.]+)", line))
  
  if (length(matches[[1]]) == 4) {
    idx1 <- as.integer(matches[[1]][2])
    idx2 <- as.integer(matches[[1]][3])
    r_val <- as.numeric(matches[[1]][4])
    
    # Apply to Covariance Matrix symmetrically
    cov_val <- r_val * u_scaled[idx1] * u_scaled[idx2]
    V[idx1, idx2] <- cov_val
    V[idx2, idx1] <- cov_val
  }
}

# =====================================================================
# 4. LEAST-SQUARES ADJUSTMENT
# =====================================================================
# Calculate Least-Squares Solution: X = (A^T V^-1 A)^-1 A^T V^-1 Y
V_inv <- solve(V)
cov_X <- solve(t(A) %*% V_inv %*% A) 
X_hat <- cov_X %*% t(A) %*% V_inv %*% Y

# =====================================================================
# 5. SELF-CONSISTENCY CHECKS
# =====================================================================
# Calculate residuals and Chi-squared
residuals <- Y - (A %*% X_hat)
chi_squared <- as.numeric(t(residuals) %*% V_inv %*% residuals)

# Birge ratio: sqrt(chi_squared / degrees of freedom)
birge_ratio <- sqrt(chi_squared / (N - M))

# =====================================================================
# 6. HIGH-PRECISION RECONSTRUCTION & OUTPUT
# =====================================================================
# Convert the offsets back to mpfr BEFORE any arithmetic involving the 1e-15 scale
X_hat_mpfr <- mpfr(as.vector(X_hat), precBits = prec)
unc_mpfr <- mpfr(sqrt(diag(cov_X)), precBits = prec)

# Reconstruct without adding to 1 (nu = nu_0 + nu_0 * X_hat * 1e-15)
optimized_freqs <- nu_0 + (nu_0 * X_hat_mpfr * 1e-15)
optimized_unc <- nu_0 * (unc_mpfr * 1e-15)

cat("--- LEAST SQUARES OPTIMIZED FREQUENCIES (nu1 to nu14) ---\n")
for (i in 1:M) {
  cat(sprintf("nu%d:\t", i), 
      format(optimized_freqs[i], nsmall=8, digits=24), "Hz\t",
      "(Unc:", sprintf("%.2e", as.numeric(optimized_unc[i])), "Hz)\n")
}

cat("\n--- CONSISTENCY METRICS ---\n")
cat("Degrees of Freedom (N - M):", N - M, "\n")
cat("Chi-squared:", round(chi_squared, 3), "\n")
cat("Birge Ratio:", round(birge_ratio, 3), "\n")
