data {
  int<lower=0> N_AlSr;
  int<lower=0> N_AlYb;
  int<lower=0> N_YbSr;
  array[N_AlSr] real x_j_data;
  array[N_AlYb] real y_j_data;
  array[N_YbSr] real z_j_data;
  array[N_AlSr] real sigma_x_j_data;
  array[N_AlYb] real sigma_y_j_data;
  array[N_YbSr] real sigma_z_j_data;
  array[N_AlSr] real a_x_j_data;
  array[N_AlYb] real a_y_j_data;
  array[N_YbSr] real a_z_j_data;
  array[N_AlSr] real b_x_j_data;
  array[N_AlYb] real b_y_j_data;
  array[N_YbSr] real b_z_j_data;
}

parameters {
  real mu_x;
  real mu_y;
  real mu_z;

  real<lower=0> xi_x;
  real<lower=0> xi_y;
  real<lower=0> xi_z;

  real alpha;
  real beta;
  real gamma;
  
  real eta_x;
  real eta_y;
  real eta_z;
}


model {
  // Priors
  mu_x ~ normal(0, 1e4);
  mu_y ~ normal(0, 1e4);
  mu_z ~ normal(0, 1e4);

  // xi_x ~ student_t(1, 0, 1000) T[0, ]; // df=1, location=0, scale=1000, truncated lower=0
  // xi_y ~ student_t(1, 0, 1000) T[0, ];
  // xi_z ~ student_t(1, 0, 1000) T[0, ];
  xi_x ~ normal(500, 1000) T[0, ];
  xi_y ~ normal(500, 1000) T[0, ];
  xi_z ~ normal(500, 1000) T[0, ];

  alpha ~ normal(0, 1);
  beta ~ normal(0, 1);
  gamma ~ normal(0, 1);

  // Priors for eta parameters
  eta_x ~ normal(mu_x, sqrt(0.2^2 + 0.4^2));
  eta_y ~ normal(mu_y, sqrt(0.2^2 + 0.3^2));
  eta_z ~ normal(mu_z, sqrt(0.2^2 + 0.4^2));

  // Likelihood
  for (j in 1:N_AlSr) {
    x_j_data[j] ~ normal(eta_x + alpha * a_x_j_data[j] - gamma * b_x_j_data[j],
                         sqrt(square(sigma_x_j_data[j]) + square(xi_x)));
  }
  for (j in 1:N_AlYb) {
    y_j_data[j] ~ normal(eta_y + alpha * a_y_j_data[j] - beta * b_y_j_data[j],
                         sqrt(square(sigma_y_j_data[j]) + square(xi_y)));
  }
  for (j in 1:N_YbSr) {
    z_j_data[j] ~ normal(eta_z + beta * a_z_j_data[j] - gamma * b_z_j_data[j],
                         sqrt(square(sigma_z_j_data[j]) + square(xi_z)));
  }
}
