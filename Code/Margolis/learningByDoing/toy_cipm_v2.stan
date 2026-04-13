data {
  vector[6] Q_scaled;
  matrix[6, 6] V_scaled;    // Now contains the off-diagonal covariances
  array[6] int num_idx;
  array[6] int den_idx;
}

parameters {
  vector[3] x;
  
  // The additive excess standard deviation (noise)
  real<lower=0> sigma_excess; 
}

transformed parameters {
  vector[6] F_pred;
  vector[4] all_x;
  
  all_x[1:3] = x;
  all_x[4] = 0.0; 
  
  for (i in 1:6) {
    F_pred[i] = all_x[num_idx[i]] - all_x[den_idx[i]];
  }
}

model {
  // Priors
  x ~ normal(0, 100);
  
  // Weak half-normal prior for the excess noise.
  // Constrained by <lower=0> in the parameters block.
  sigma_excess ~ normal(0, 5); 
  
  // Create a modified covariance matrix
  matrix[6, 6] V_mod = V_scaled;
  
  // Add the excess variance to the diagonal (independent noise)
  for (i in 1:6) {
    V_mod[i, i] = V_scaled[i, i] + square(sigma_excess);
  }
  
  // Evaluate the likelihood using the dense, modified covariance matrix
  Q_scaled ~ multi_normal(F_pred, V_mod);
}
