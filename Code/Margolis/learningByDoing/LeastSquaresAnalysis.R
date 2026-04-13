#' Least-Squares Analysis of Clock Frequency Comparison Data
#'
#' @param q Vector of N measured frequency ratios.
#' @param V Covariance matrix of the measured values (N x N).
#' @param s_initial Vector of M initial estimates for adjusted frequency ratios.
#' @param func_f A function that takes `s` (length M) and returns predicted `q` (length N).
#' @param func_A A function that takes `s` and returns the Jacobian matrix A (N x M).
#'               If NULL, you could implement a numerical derivative helper.
#' @param tol Convergence tolerance (fraction of uncertainty).
#' @param max_iter Maximum number of iterations.
#' @return A list containing the optimized ratios, covariance matrix, and consistency checks.

clock_least_squares <- function(q, V, s_initial, func_f, func_A, tol = 1e-12, max_iter = 100) {
  
  N <- length(q)
  M <- length(s_initial)
  
  # Ensure the system is over-determined or exactly determined
  if (N <= M) {
    stop("The number of measurements N must be strictly greater than M for Birge ratio calculation.")
  }
  
  # The paper notes that if V is well-conditioned, Cholesky decomposition can be used,
  # but orthogonal decomposition is preferred if poorly conditioned.
  # Here we use solve() which utilizes LU decomposition, but qr.solve() can be used for stability.
  V_inv <- solve(V) 
  
  s_current <- s_initial
  iteration <- 0
  converged <- FALSE
  
  while (iteration < max_iter && !converged) {
    iteration <- iteration + 1
    
    # 1. Calculate the predicted values F based on current estimates s
    F_val <- func_f(s_current)
    
    # 2. Calculate the difference Y = Q - F(s)
    Y <- q - F_val
    
    # 3. Calculate the Jacobian matrix A evaluated at s_current
    A <- func_A(s_current)
    
    # 4. Perform least-squares adjustment to find X that minimizes (Y-AX)^T V^-1 (Y-AX)
    # The normal equations solution is X = (A^T V^-1 A)^-1 A^T V^-1 Y
    At_Vinv <- t(A) %*% V_inv
    cov_X <- solve(At_Vinv %*% A)
    X_hat <- cov_X %*% At_Vinv %*% Y
    
    # 5. Update the adjusted frequency ratios
    s_next <- s_current + as.numeric(X_hat)
    
    # 6. Check for convergence (differ by a sufficiently small fraction)
    # Using the relative change in estimates as the convergence metric
    if (max(abs(X_hat / s_current)) < tol) {
      converged <- TRUE
    }
    
    s_current <- s_next
  }
  
  if (!converged) {
    warning("Least-squares adjustment did not converge within max_iter.")
  }
  
  # 7. Calculate final best estimates of measured quantities Q_hat
  # Q_hat = F + A * X_hat (evaluated at the final converged step)
  Q_hat <- func_f(s_current) 
  
  # 8. Perform self-consistency checks
  # Calculate minimum value of S (chi-squared)
  residuals <- q - Q_hat
  chi_squared <- as.numeric(t(residuals) %*% V_inv %*% residuals)
  
  # Calculate Birge ratio
  birge_ratio <- sqrt(chi_squared / (N - M))
  
  # Normalized residuals for individual measurements
  normalized_residuals <- residuals / sqrt(diag(V))
  
  # Return output list
  return(list(
    optimized_ratios_Z = s_current,
    covariance_matrix_Z = cov_X,
    calculated_measurements_Q = Q_hat,
    chi_squared = chi_squared,
    birge_ratio = birge_ratio,
    normalized_residuals = normalized_residuals,
    iterations = iteration
  ))
}


# 1. Define the measurements (q) and their covariance matrix (V)
# We assume the measurements are independent for simplicity, so V is diagonal.
q <- c(2.1, 4.95, 10.10)

# Uncertainties (standard deviations) for our 3 measurements
u <- c(0.1, 0.05, 0.10)
V <- diag(u^2) # Variances on the diagonal

# 2. Provide initial guesses for our adjusted ratios (z1 and z2)
# We can just use our direct measurements q1 and q2 as a starting point.
s_initial <- c(2.0, 5.0)

# 3. Define func_f: Express the measured quantities in terms of the adjusted ratios
# f1 = z1
# f2 = z2
# f3 = z1 * z2
func_f <- function(s) {
  z1 <- s[1]
  z2 <- s[2]
  return(c(
    z1, 
    z2, 
    z1 * z2
  ))
}

# 4. Define func_A: The Jacobian matrix (partial derivatives of func_f with respect to z1 and z2)
# Row 1 (df1/dz): [1, 0]
# Row 2 (df2/dz): [0, 1]
# Row 3 (df3/dz): [z2, z1]
func_A <- function(s) {
  z1 <- s[1]
  z2 <- s[2]
  A <- matrix(0, nrow = 3, ncol = 2)
  
  A[1, 1] <- 1;  A[1, 2] <- 0
  A[2, 1] <- 0;  A[2, 2] <- 1
  A[3, 1] <- z2; A[3, 2] <- z1
  
  return(A)
}

# 5. Run the least-squares adjustment (assuming the function from the previous response is loaded)
results <- clock_least_squares(
  q = q, 
  V = V, 
  s_initial = s_initial, 
  func_f = func_f, 
  func_A = func_A
)

# 6. View the optimized results
cat("Optimized Ratios (z1, z2):", results$optimized_ratios_Z, "\n")
cat("Birge Ratio:", results$birge_ratio, "\n")
cat("Iterations to converge:", results$iterations, "\n")