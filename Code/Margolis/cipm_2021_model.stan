data {
  int<lower=1> N_measurements; 
  int<lower=1> N_transitions;  
  
  vector[N_measurements] Q_scaled;       // The scaled fractional offsets (~1.0 to 10.0)
  matrix[N_measurements, N_measurements] V_scaled; // The scaled covariance matrix
  
  array[N_measurements] int<lower=1, upper=N_transitions+1> num_idx; 
  array[N_measurements] int<lower=1, upper=N_transitions+1> den_idx;
}

parameters {
  // We are now estimating the scaled fractional offset from the 2017 CIPM values
  vector[N_transitions] x; 
}

transformed parameters {
  vector[N_measurements] F_pred;
  vector[N_transitions + 1] all_x;
  
  for (i in 1:N_transitions) {
    all_x[i] = x[i];
  }
  // Cesium's offset is exactly 0 by definition
  all_x[15] = 0.0; 
  
  // Predict the scaled fractional offset for each measurement
  for (i in 1:N_measurements) {
    if (den_idx[i] == 15) {
      // For absolute measurements, the offset is just the numerator's offset
      F_pred[i] = all_x[num_idx[i]];
    } else {
      // For ratios, the offset is the numerator offset minus the denominator offset
      F_pred[i] = all_x[num_idx[i]] - all_x[den_idx[i]];
    }
  }
}

model {
  // Weak prior: standard deviation of 100 in scaled space 
  // represents an enormous +/- 1e-13 fractional uncertainty 
  x ~ normal(0, 100);
  
  // The Likelihood is now numerically pristine
  Q_scaled ~ multi_normal(F_pred, V_scaled);
}