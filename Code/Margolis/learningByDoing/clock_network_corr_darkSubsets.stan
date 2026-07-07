data {
  int<lower=1> N_clocks; 
  int<lower=1> N_meas;   
  
  vector[N_meas] Y;      
  matrix[N_meas, N_meas] V; // Pass the RAW covariance matrix V
  
  array[N_meas] int<lower=0, upper=1> is_ratio;        
  array[N_meas] int<lower=1, upper=N_clocks> clock_1;
  array[N_meas] int<lower=0, upper=N_clocks> clock_2;  
  
  // NEW: Mapping for dark uncertainties
  int<lower=1> N_dark; // Number of unique dark uncertainty parameters
  array[N_clocks] int<lower=1, upper=N_dark> clock_to_dark_idx;
}

parameters {
  vector[N_clocks] x; 
  // Dark uncertainty vector is now dynamically sized based on the grouping
  vector<lower=0>[N_dark] sigma_dark; 
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

  // 2. Dynamically add mapped dark uncertainties to the diagonal
  V_total = V;
  for (i in 1:N_meas) {
    if (is_ratio[i] == 0) {
      // Absolute measurement: Add the mapped clock's dark variance
      V_total[i, i] += square(sigma_dark[clock_to_dark_idx[clock_1[i]]]);
    } else {
      // Ratio measurement: Add BOTH mapped clocks' dark variances in quadrature
      V_total[i, i] += square(sigma_dark[clock_to_dark_idx[clock_1[i]]]) + 
                       square(sigma_dark[clock_to_dark_idx[clock_2[i]]]);
    }
  }

  // 3. Priors
  x ~ normal(0, 100);
  sigma_dark ~ normal(0, 1);

  // 4. Likelihood using the dynamically decomposed matrix
  // Y ~ multi_normal_cholesky(mu, cholesky_decompose(V_total));
  Y ~ multi_student_t_cholesky(3, mu, cholesky_decompose(V_total)); 
}