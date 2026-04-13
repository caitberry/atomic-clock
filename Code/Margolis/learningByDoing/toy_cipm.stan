data {
  vector[6] Q_scaled;       // Our 6 scaled measurements
  matrix[6, 6] V_scaled;    // Our 6x6 diagonal covariance matrix
  
  array[6] int num_idx;     // Numerator indices (1=Hg, 2=Yb, 3=Sr)
  array[6] int den_idx;     // Denominator indices (4=Cesium)
}

parameters {
  vector[3] x; // The scaled fractional offsets for our 3 clocks
}

transformed parameters {
  vector[6] F_pred;
  vector[4] all_x;
  
  // Map our 3 clock parameters into the array
  all_x[1:3] = x;
  all_x[4] = 0.0; // Cesium offset is exactly 0
  
  // Calculate the predicted ratio/offset for each of the 6 measurements
  for (i in 1:6) {
    F_pred[i] = all_x[num_idx[i]] - all_x[den_idx[i]];
  }
}

model {
  // A wide, flat prior centered at 0
  x ~ normal(0, 100);
  
  // The Likelihood (evaluating our 6 predictions against the data)
  Q_scaled ~ multi_normal(F_pred, V_scaled);
}