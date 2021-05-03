//
// Mean trends by strata
// modeled on the State of Canada's Birds models
// 
//


data {
  int<lower=0> nstrata; // number of strata
  int<lower=0> N; // number of region by species combinations
  
  vector[N] t_hat; //vector of trends
  vector[N] sd_hat; //vector of sd estimates
  int<lower=0> strata[N]; //vector of strata
  
  
}


parameters {
  vector[nstrata] mu_strata_raw;
  real true_t[N];
  real<lower=0> epsilon;
  real<lower=0> sd_strata_raw;
  
  
}

transformed parameters {
  real mu_t[N];
  vector[nstrata] mu_strata;
  
  mu_strata = sd_strata_raw*mu_strata_raw; // non-centered prior
  
  for(i in 1:N){
  mu_t[i] = mu_strata[strata[i]];
  }
}


model {
// likelihood
// process model

  mu_strata_raw ~ normal(0,1);//regularizing prior
  sum(mu_strata_raw) ~ normal(0,0.001*nstrata); // zero-sum forcing
  sd_strata_raw ~ gamma(2,0.1); //zero-avoiding prior 
  epsilon ~ gamma(2,0.1); //zero-avoiding prior 
  
   true_t ~ normal(mu_t,epsilon);
   
// observation model   
for(i in 1:N){
   t_hat[i] ~ normal(true_t[i],sd_hat[i]);
    } 
  
}





