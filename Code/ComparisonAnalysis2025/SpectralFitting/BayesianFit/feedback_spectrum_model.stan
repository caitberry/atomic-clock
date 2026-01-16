functions {
  // MODIFICATION 1: Added 'real tp' to the function's argument list.
  vector model_log_psd(int N, vector omega, real h0, real h_m1,
                       real Kp, real Ki, real tau, real tp) {
    vector[N] log_y_model;
    
    // MODIFICATION 2: The hardcoded 'real tp = 1;' is now removed.
    // The function now uses the 'tp' value that was passed to it.
    real tp_half = tp / 2;

    for (n in 1:N) {
      real w = omega[n];
      real inv_w = 1.0 / w;
      real angle = -w * tp_half;

      // Precomputed trig
      real cos_wtp = cos(angle);
      real sin_wtp = sin(angle);

      // Numerator of G(w)
      real re_num = Kp * cos_wtp + Ki * inv_w * sin_wtp;
      real im_num = Kp * sin_wtp - Ki * inv_w * cos_wtp;

      // Denominator of G(w)
      real re_den = 1.0;
      real im_den = w * tau;
      real denom_mag_sq = re_den^2 + im_den^2;

      // Complex division: G = num / den
      real re_G = (re_num * re_den + im_num * im_den) / denom_mag_sq;
      real im_G = (im_num * re_den - re_num * im_den) / denom_mag_sq;

      // 1 + G(w)
      real re_1pG = 1.0 + re_G;
      real im_1pG = im_G;
      real one_plus_G_mag_sq = re_1pG^2 + im_1pG^2;

      // |G / (1 + G)|^2
      real re_ratio = (re_G * re_1pG + im_G * im_1pG) / one_plus_G_mag_sq;
      real im_ratio = (im_G * re_1pG - re_G * im_1pG) / one_plus_G_mag_sq;
      real abs_ratio_sq = re_ratio^2 + im_ratio^2;

      // |1 / (1 + G)|^2
      real abs_inv_sq = 1.0 / one_plus_G_mag_sq;

      // Model PSD (still linear)
      real y_val = abs_ratio_sq * h0 + abs_inv_sq * (2 * pi() * h_m1 * inv_w);

      // Take log
      log_y_model[n] = log(y_val);
    }
    return log_y_model;
  }
}

data {
  int<lower=1> N;
  vector<lower=0>[N] omega;       // angular frequencies (rad/s)
  vector<lower=0>[N] y_obs;       // observed PSD
  vector<lower=0>[N] rel_sd;      // relative sd on linear scale
  
  // MODIFICATION 3: Added 'tp' to the data block.
  real<lower=0> tp;               // probe time (s)
}

parameters {
  real<lower=0> h0;
  real<lower=0> h_m1;
  real<lower=0> Kp;
  real<lower=0> Ki;
  real<lower=0> tau;
}

transformed parameters {
  // MODIFICATION 4: Pass 'tp' from the data block to the function.
  vector[N] log_y_hat = model_log_psd(N, omega, h0, h_m1, Kp, Ki, tau, tp);
}

model {
  // Priors (these are unchanged)
  // h0 ~ lognormal(log(1e-31), 1.0);
  h0 ~ lognormal(-72.78, 2.49);
  // h_m1 ~ lognormal(log(1.5e-33), 0.5);
  h_m1 ~ lognormal(-75.79, 0.18);

  // Kp ~ normal(10, 5);
  Kp ~ normal(10, 5) T[0,]; 

  // Ki ~ normal(5, 2);
  Ki ~ lognormal(0, 1.17);
  // tau ~ normal(50, 30);
  tau ~ lognormal(2.65, 1.35);

  // Likelihood (this is unchanged)
  // The 'rel_sd' is calculated in R and passed in, so this remains correct.
  y_obs ~ lognormal(log_y_hat, rel_sd);
}

