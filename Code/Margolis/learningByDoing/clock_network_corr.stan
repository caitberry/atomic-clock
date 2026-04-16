data {
  int<lower=1> N_clocks; 
  int<lower=1> N_meas;   
  
  vector[N_meas] Y;      
  matrix[N_meas, N_meas] V; // Pass the RAW covariance matrix V, not L_V
  
  array[N_meas] int<lower=0, upper=1> is_ratio;        
  array[N_meas] int<lower=1, upper=N_clocks> clock_1;  
  array[N_meas] int<lower=0, upper=N_clocks> clock_2;  
}

parameters {
  vector[N_clocks] x; 
  real<lower=0> sigma_dark; // The global dark uncertainty parameter (in parts per 10^15)
}

model {
  vector[N_meas] mu;
  matrix[N_meas, N_meas] V_total;
  
  // 1. Build the expected mean vector
  for (i in 1:N_meas) {
    if (is_ratio[i] == 0) {
      mu[i] = x[clock_1[i]];               
    } else {
      mu[i] = x[clock_1[i]] - x[clock_2[i]]; 
    }
  }

  // 2. Add dark uncertainty variance to the diagonal of the known covariance matrix
  V_total = V;
  for (i in 1:N_meas) {
    V_total[i, i] += square(sigma_dark);
  }

  // 3. Priors
  x ~ normal(0, 100); 
  
  // Half-normal prior for dark uncertainty (constrains it to be positive but prefers smaller values)
  // You can adjust the scale (e.g., 5) based on your domain knowledge of maximum plausible systematics
  sigma_dark ~ normal(0, 5); 

  // 4. Likelihood using the dynamically decomposed matrix
  Y ~ multi_normal_cholesky(mu, cholesky_decompose(V_total));
  // Y ~ multi_student_t_cholesky(3, mu, cholesky_decompose(V_total)); TRY THIS NEXT
  
}
