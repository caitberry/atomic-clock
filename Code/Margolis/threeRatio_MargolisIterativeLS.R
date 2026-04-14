###### might need to try something like this
# library(Rmpfr)
# # Set precision to 128 bits (roughly 38 decimal digits)
# prec <- 128
# 
# # Define your frequency ratios as mpfr objects
# # Use character strings to prevent R from truncating them before conversion
# q_mpfr <- mpfr(c(
#   "1128575290808154.62", 
#   "2.17747319413456507"
# ), precBits = prec)
# Now, all arithmetic (+, -, *, /) and matrix operations 
# using q_mpfr will maintain this high precision.

options(digits = 22)
# 1. Define the measurements (q) and their covariance matrix (V)
# We assume the measurements are independent for simplicity, so V is diagonal.
q <- c(
  1128575290808154.62, # 1. Absolute Hg (q7)
  518295836590863.59,  # 2. Absolute Yb (q24)
  429228004229872.92,  # 3. Absolute Sr (q47)
  2.17747319413456507,   # 4. Ratio Hg / Yb (q79)
  1.20750703934333841,   # 5. Ratio Yb / Sr (q81)
  2.62931420989890960    # 6. Ratio Hg / Sr (q59)
)

q=q/c(9192631770, 9192631770, 9192631770, 1, 1, 1) # dividing by Cs

# Their absolute uncertainties (U)
# Corrected absolute uncertainties (U)
u_raw <- c(
  0.41,       # 1. Hg absolute unc
  0.31,       # 2. Yb absolute unc
  0.12,       # 3. Sr absolute unc
  1.92e-16,   # 4. Hg/Yb unc 
  3.40e-16,   # 5. Yb/Sr unc 
  2.20e-16    # 6. Hg/Sr unc 
)

u <- c(u_raw[1:3]/9192631770, #dividing by Cs
       u_raw[4]*q[4], u_raw[5]*q[5], u_raw[6]*q[6]) #convert relative uncertainties 

# 2. Covariance Matrix V (Assuming zero correlation for this example)
V <- diag(u^2)

# 3. Starting Values (Initial estimates for z1, z2, z3)
# z1 = Hg/Yb, z2 = Yb/Sr, z3 = Sr/Cs
s <- c(q[4]+.001, q[5]+.001, q[3]+.001) #starting a bit off to see it change in the iterations
# s <- c(q[4], q[5], q[3]) 

# 4. Iterative Least-Squares Adjustment
for (iteration in 1:5) {
  
  # Define functions f_i(z) for each measurement q_i
  # f1: Hg/Cs = z1 * z2 * z3
  # f2: Yb/Cs = z2 * z3
  # f3: Sr/Cs = z3
  # f4: Hg/Yb = z1
  # f5: Yb/Sr = z2
  # f6: Hg/Sr = z1 * z2
  
  f_s <- c(s[1]*s[2]*s[3], s[2]*s[3], s[3], s[1], s[2], s[1]*s[2])
  
  # Design Matrix A (Partial derivatives df_i / dz_j)
  A <- matrix(0, nrow=6, ncol=3)
  # row 1 (Hg/Cs)
  A[1,1] <- s[2]*s[3]; A[1,2] <- s[1]*s[3]; A[1,3] <- s[1]*s[2]
  # row 2 (Yb/Cs)
  A[2,2] <- s[3]; A[2,3] <- s[2]
  # row 3 (Sr/Cs)
  A[3,3] <- 1
  # row 4 (Hg/Yb)
  A[4,1] <- 1
  # row 5 (Yb/Sr)
  A[5,2] <- 1
  # row 6 (Hg/Sr)
  A[6,1] <- s[2]; A[6,2] <- s[1]
  
  # Vector Y = q - f(s)
  Y <- q - f_s
  
  # Solve: X = (A' V^-1 A)^-1 A' V^-1 Y
  V_inv <- solve(V)
  X <- solve(t(A) %*% V_inv %*% A) %*% t(A) %*% V_inv %*% Y
  
  # Update starting values
  s <- s + as.vector(X)
  
  cat(sprintf("Iteration %d: s1 = %.40f\n", iteration, s[1]))
}

# Final optimized frequencies
optimized_Sr <- s[3] * 9192631770
optimized_Yb <- s[2] * optimized_Sr
optimized_Hg <- s[1] * optimized_Yb

cat("\n--- Optimized Absolute Frequencies (Hz) ---\n")
cat(sprintf("Hg: %.4f\nYb: %.4f\nSr: %.4f\n", optimized_Hg, optimized_Yb, optimized_Sr))


# 
# 
# # After the loop has converged:
# 
# # 1. Final Covariance Matrix of the adjusted ratios (z)
# # Scale V_inv properly as you did in the solve step
# V_scaled_inv <- solve(V / c(9192631770^2, 9192631770^2, 9192631770^2, 1, 1, 1))
# cov_Z <- solve(t(A) %*% V_scaled_inv %*% A)
# 
# # 2. Extract standard uncertainties (square root of the diagonal)
# u_z <- sqrt(diag(cov_Z))
# 
# cat("\n--- Uncertainties of Adjusted Ratios ---\n")
# cat(sprintf("u(Hg/Yb): %.4e\nu(Yb/Sr): %.4e\nu(Sr/Cs): %.4e\n", u_z[1], u_z[2], u_z[3]))
# 
# # 3. Absolute frequency uncertainties (Propagation)
# # For Sr (simple, since q3 = z3)
# u_Sr <- u_z[3] * 9192631770
# 
# # For Yb (Yb = z2 * z3 * Cs)
# # u^2(Yb) = (dYb/dz2)^2 * u^2(z2) + (dYb/dz3)^2 * u^2(z3) + 2*(dYb/dz2)*(dYb/dz3)*cov(z2,z3)
# grad_Yb <- c(0, s[3]*9192631770, s[2]*9192631770)
# u_Yb <- sqrt(t(grad_Yb) %*% cov_Z %*% grad_Yb)
# 
# cat("\n--- Final Absolute Uncertainties (Hz) ---\n")
# cat(sprintf("Hg Uncertainty: (Needs similar propagation)\n"))
# cat(sprintf("Yb Uncertainty: %.4f Hz\n", u_Yb))
# cat(sprintf("Sr Uncertainty: %.4f Hz\n", u_Sr))