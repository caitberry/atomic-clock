data {
  int<lower=1> N_clocks; 
  int<lower=1> N_meas;   
  
  vector[N_meas] Y;      
  matrix[N_meas, N_meas] V; // Pass the RAW covariance matrix V
  
  array[N_meas] int<lower=0, upper=1> is_ratio;        
  array[N_meas] int<lower=1, upper=N_clocks> clock_1;  
  array[N_meas] int<lower=0, upper=N_clocks> clock_2;  
}

parameters {
  vector[N_clocks] x; 
  // Now a vector of dark uncertainties, one for each clock type!
  vector<lower=0>[N_clocks] sigma_dark; 
}

model {
  vector[N_meas] mu;
  matrix[N_meas, N_meas] V_total;
  
  // 1. Build the expected mean vector without a design matrix
  for (i in 1:N_meas) {
    if (is_ratio[i] == 0) {
      mu[i] = x[clock_1[i]];               
    } else {
      mu[i] = x[clock_1[i]] - x[clock_2[i]]; 
    }
  }

  // 2. Dynamically add clock-specific dark uncertainties to the diagonal
  V_total = V;
  for (i in 1:N_meas) {
    if (is_ratio[i] == 0) {
      // Absolute measurement: Just add the single clock's dark variance
      V_total[i, i] += square(sigma_dark[clock_1[i]]);
    } else {
      // Ratio measurement: Add BOTH clocks' dark variances in quadrature
      V_total[i, i] += square(sigma_dark[clock_1[i]]) + square(sigma_dark[clock_2[i]]);
    }
  }

  // 3. Priors
  x ~ normal(0, 100); 
  
  // Applies the prior to all elements in the sigma_dark vector
  sigma_dark ~ normal(0, 1); 

  // 4. Likelihood using the dynamically decomposed matrix
  // Y ~ multi_normal_cholesky(mu, cholesky_decompose(V_total));
  Y ~ multi_student_t_cholesky(3, mu, cholesky_decompose(V_total));

}