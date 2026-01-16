data {
  int<lower=1> N; // number of labs
  vector[N] x; // lab sample means 
  vector<lower=0>[N] u; // lab sample variances
  real<lower=0> tdf;
}

parameters {
  real mu; // overall mean of lab effects (lambda)
  real<lower=0> tau; // std deviation of lab effects (lambda)
  real lambda[N]; // lab effects
}


model {
  
  // priors
  mu ~ normal(0,10^5); // overall mean
  // tau ~ student_t(4,0,tdf); // std dev of true lab means
  tau ~ normal(500, 1000) T[0, ];

  lambda ~ normal(mu,tau); // true lab means
  
  x ~ normal(lambda,u);   

}

