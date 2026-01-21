//this is same as feedback_spectrum_model.stan, but using frequency instead of omega = 2 pi f

functions {
  // MODIFIED: Input is now 'vector freq' (Hz) instead of omega
  vector model_log_psd(int N, vector freq, real h0, real h_m1, 
                       real Kp, real Ki, real tau, real tp) {
    vector[N] log_y_model;
    
    real tp_half = tp / 2;
    real two_pi = 2.0 * pi();

    for (n in 1:N) {
      real f = freq[n];
      
      // CONVERSION: Physics (G) still needs omega (rad/s)
      real w = two_pi * f; 
      
      // Precomputed terms
      real angle = -w * tp_half;
      real cos_wtp = cos(angle);
      real sin_wtp = sin(angle);

      // --- Calculate G(w) ---
      // Note: We use 'w' here because transfer functions 
      // are defined in the s-domain (s = iw)
      
      // Numerator of G(w) = (Kp + Ki/iw) * e^(-iw*tp/2)
      // (Kp + Ki/iw) = (Kp - i*Ki/w)
      real term_re = Kp;
      real term_im = -Ki / w; 
      
      // Multiply by e^(-iw*tp/2) = cos + i*sin
      real re_num = term_re * cos_wtp - term_im * sin_wtp;
      real im_num = term_re * sin_wtp + term_im * cos_wtp;

      // Denominator of G(w) = 1 + iw*tau
      real re_den = 1.0;
      real im_den = w * tau;
      real denom_mag_sq = re_den^2 + im_den^2;

      // Complex division G = num / den
      real re_G = (re_num * re_den + im_num * im_den) / denom_mag_sq;
      real im_G = (im_num * re_den - re_num * im_den) / denom_mag_sq;

      // --- Calculate Loop Shapes ---
      
      // 1 + G(w)
      real re_1pG = 1.0 + re_G;
      real im_1pG = im_G;
      real one_plus_G_mag_sq = re_1pG^2 + im_1pG^2;

      // Shape 1: |G / (1+G)|^2 (Noise Transfer Function)
      real re_ratio = (re_G * re_1pG + im_G * im_1pG) / one_plus_G_mag_sq;
      real im_ratio = (im_G * re_1pG - re_G * im_1pG) / one_plus_G_mag_sq;
      real abs_ratio_sq = re_ratio^2 + im_ratio^2;

      // Shape 2: |1 / (1+G)|^2 (Sensitivity Function)
      real abs_inv_sq = 1.0 / one_plus_G_mag_sq;

      // --- Combine Noise Sources ---
      // SIMPLIFICATION: We use 'f' directly for the Random Walk term.
      // Previously: (2*pi*h_m1 / w). Now: (h_m1 / f).
      real y_val = abs_ratio_sq * h0 + abs_inv_sq * (h_m1 / f);

      log_y_model[n] = log(y_val);
    }
    return log_y_model;
  }
}

data {
  int<lower=1> N;
  vector<lower=0>[N] freq;        // MODIFIED: frequency in Hz
  vector<lower=0>[N] y_obs;       // observed PSD (per Hz)
  vector<lower=0>[N] rel_sd;      // relative sd
  real<lower=0> tp;               // probe time (s)
  
  real bias;
}

parameters {
  real<lower=0> h0;
  real<lower=1e-33, upper=2e-33> h_m1;
  real<lower=0,upper = 20> Kp;
  real<lower=0,upper=10> Ki;
  real<lower=1, upper =200> tau;
}

transformed parameters {
  // MODIFIED: Passing 'freq' instead of 'omega'
  vector[N] log_y_hat = model_log_psd(N, freq, h0, h_m1, Kp, Ki, tau, tp);
  
  // Shift the model down to match the biased data
  // The model predicts the "True" log-spectrum.
  // The likelihood compares it to the "Biased" data.
  // So we add the negative bias to the prediction here.
  log_y_hat = log_y_hat + bias;
}

model {
  // Priors 
  // h0 ~ lognormal(-72.78, 2.49);
  // h_m1 ~ lognormal(-75.79, 0.18);
  // Kp ~ normal(10, 5) T[0,]; 
  // Ki ~ lognormal(0, 1.17);
  // tau ~ lognormal(2.65, 1.35);


// --- Noise Parameters ---
  // h0   ~ lognormal(-72.78, 2.49); 
  h0   ~ lognormal(-72.78, 5); 
  // h_m1 ~ lognormal(-75.79, 1);
  h_m1 ~ uniform(1e-33,2e-33);
  // --- Loop Parameters (THE FIX) ---
  
  // 1. Kp: Was normal(10, 5). 
  // Change to Weakly Informative LogNormal centered near typical values (e.g., 0.1 to 10)
  // lognormal(-1, 2) covers 0.05 to 7.0 comfortably.
  // Kp ~ lognormal(-1, 2); 

  // 2. Ki: Was lognormal(0, 1.17) -> median 1.0. 
  // Your truth is 0.12. This is fine, but maybe widen it.
  // Ki ~ lognormal(-1, 2);

  // 3. Tau: Was lognormal(2.65, 1.35) -> median 14.
  // Your truth is 2.5. Let's center it closer to 5s or just make it very wide.
  // lognormal(1.0, 1.5) -> median 2.7, covers 0.1s to 50s.
  // tau ~ lognormal(1.0, 1.5);
  
  // Kp  ~ lognormal(0, 5);
  // Ki  ~ lognormal(0, 5);
  // tau ~ lognormal(0, 5);
  tau ~ uniform(1,200);
  Kp  ~ uniform(0, 20);
  Ki  ~ uniform(0, 10);
  
  // Likelihood
  y_obs ~ lognormal(log_y_hat, rel_sd);
}
