functions {
  vector model_log_psd(int N, vector freq, real log_h0, real log_h_m1, 
                       real Kp, real Ki, real tau, real tp) {
    
    vector[N] log_y_model;
    
    real h0   = exp(log_h0);
    real h_m1 = exp(log_h_m1);
    
    complex J = to_complex(0, 1); // Imaginary unit 'i'
    real two_pi = 2.0 * pi();

    for (n in 1:N) {
      real f = freq[n];
      real w = two_pi * f;
      complex s = J * w;
      
      // --- THE PHYSICS (Simplified) ---
      
      // Calculate Open-Loop Gain G(s)
      // G(s) = [ (Kp + Ki/s) / (1 + s*tau) ] * e^(-s*tp/2)
      complex G = ( (Kp + Ki/s) / (1 + s*tau) ) * exp(-s * tp / 2.0);

      // Calculate Loop Shapes
      // abs(z) returns magnitude |z|
      real shape_h0  = abs(G / (1.0 + G))^2;  // | G / (1+G) |^2
      real shape_hm1 = abs(1.0 / (1.0 + G))^2; // | 1 / (1+G) |^2

      // Combine Noise Sources
      real y_val = shape_h0 * h0 + shape_hm1 * (h_m1 / f);
      
      // Store Log-Value
      log_y_model[n] = log(y_val);
    }
    return log_y_model;
  }
}

data {
  int<lower=1> N;
  vector<lower=0>[N] freq;        // Frequency (Hz)
  vector<lower=0>[N] y_obs;       // Observed PSD
  vector<lower=0>[N] sigma_log;   // Standard Deviation (on log scale)
  real<lower=0> tp;               // Probe time (s)
  
  real bias;                      // Bias correction term
}

parameters {
  // --- Noise Parameters (Log Scale) ---
  real log_h0; 
  
  // Constrained Log-Random Walk
  // log(1e-33) ~ -76.0, log(2e-33) ~ -75.3, bounds [-80, -70] are safe
  real<lower=-80, upper=-70> log_h_m1; 
  
  // --- Loop Parameters (Bounded) ---
  real<lower=0, upper=20> Kp;
  real<lower=0, upper=10> Ki;
  real<lower=1, upper=200> tau;
}

transformed parameters {
  // 1. Calculate the "True" Physical Model
  vector[N] log_y_hat = model_log_psd(N, freq, log_h0, log_h_m1, Kp, Ki, tau, tp);
  
  // 2. Add Bias for the Likelihood 
  log_y_hat = log_y_hat + bias;
}

model {
  // --- Priors ---
  log_h0   ~ normal(-72.0, 2.0); 
  log_h_m1 ~ normal(-75.63, 0.17); // Centered between 1e-33 and 2e-33
  tau ~ uniform(1, 200);
  Kp ~ lognormal(0, 4) T[, 20]; // keep an eye on this prior, the bounds Dave gave are too wide and influence the h0 est too much. Inform with science
  Ki  ~ uniform(0, 10);

  // --- Likelihood ---
  y_obs ~ lognormal(log_y_hat, sigma_log);
}

