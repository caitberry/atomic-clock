data {
  int<lower=1> N; // number of labs
  vector[N] x; // lab sample means 
  vector<lower=0>[N] u; // lab sample variances
  real<lower=0> tdf;
}

parameters {
  real mu; // overall mean of lab effects (lambda)
  real<lower=0> tau; // std deviation of lab effects (lambda)
  vector[N] lambda_std; // standardized lab effects
}

transformed parameters {
  vector[N] lambda; // lab effects

  lambda = mu + tau*lambda_std;
}

model {
  
  // priors
  mu ~ normal(0,10^5); // overall mean
  // tau ~ student_t(4,0,tdf); // std dev of true lab means
  tau ~ normal(500, 1000) T[0, ];
  
  lambda_std ~ normal(0.0, 1.0);  
  
  x ~ normal(lambda,u);   

}

