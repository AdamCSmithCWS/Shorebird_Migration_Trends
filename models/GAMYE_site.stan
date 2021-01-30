// This is a full hierarchical GAM time-series, with spatial gam parameters 
// as well as a gam-based Seasonal adjustment and sruvey-wide random year-effects


functions {
  real icar_normal_lpdf(vector bb, int nsites, int[] node1, int[] node2) {
    return -0.5 * dot_self(bb[node1] - bb[node2])
      + normal_lpdf(sum(bb) | 0, 0.001 * nsites); //soft sum to zero constraint on phi
 }
}


data {
  // int<lower=0> nstrata;
  int<lower=0> ncounts;
  int<lower=0> nsites;
  int<lower=0> nyears;
  
  // data for GAM s(date)
  int<lower=1> ndays;  // number of days in the basis function for season
  int<lower=1> nknots_season;  // number of knots in the basis function for season
  matrix[ndays, nknots_season] season_basispred; // basis function matrix
  
  // data for GAM s(year)
  int<lower=1> nknots_year;  // number of knots in the basis function for year
  matrix[nyears, nknots_year] year_basispred; // basis function matrix
  
  //data for spatial iCAR among sites (GAM parameters for s(year))
  int<lower=1> N_edges;
  int<lower=1, upper=nsites> node1[N_edges];  // node1[i] adjacent to node2[i]
  int<lower=1, upper=nsites> node2[N_edges];  // and node1[i] < node2[i]

  int<lower=0> count[ncounts];              // count observations
  //int<lower=1> strat[ncounts];              // strata indicators
  int<lower=1> site[ncounts];              // site indicators
  real year[ncounts];              // centered years
  int<lower=1> year_raw[ncounts]; // year index
  int<lower=1> date[ncounts];  // day indicator in the season
  
  
  vector[nsites] site_size; //log-scale size predictor for site level abundance
  // upper 98th quartile of all counts across species with similar flocking behaviour
  // conceptually, this should reflect the maximum size of flock that could be observed
  // at this site, acting like an offset for area surveyed (since the area is hard to define)
  // and generally not available for most sites
  
  
  //indices for re-scaling predicted counts within strata based on site-level intercepts
  // int<lower=1> max_sites; //dimension 1 of sites matrix
  // int<lower=0> sites[max_sites,nstrata]; //matrix of which sites are in each stratum
  // int<lower=1> nsites_strat[nstrata]; //number of unique sites in each stratum
}

parameters {
  vector[nsites] alpha_raw;             // intercepts
  vector[nyears] year_effect_raw;             // continental year-effects
  vector[ncounts] noise_raw;             // over-dispersion
  matrix[nsites,nknots_year] b_raw;         // spatial effect slopes (0-centered deviation from continental mean slope B)
  vector[nknots_year] B_raw;             // GAM coefficients year
  
  vector[nknots_season] B_season_raw;         // GAM coefficients
  
  real beta_size; //effect of site level predictor
 real<lower=0> sdnoise;    // sd of over-dispersion
 //real<lower=1> nu; 
  real<lower=0> sdalpha;    // sd of site effects
  real<lower=0> sdyear;    // sd of year effects
  real<lower=0> sdseason;    // sd of year effects
  real<lower=0> sdyear_gam;    // sd of GAM Hyperparameter coefficients
  real<lower=0> sdyear_gam_strat[nknots_year]; //sd of site level gams for each knot
  real ALPHA1; // overall intercept
}

transformed parameters { 
  vector[ncounts] E;           // log_scale additive likelihood
  vector[ndays] season_pred;
    matrix[nyears,nsites] year_pred;
  vector[nyears] Y_pred; 
  vector[nsites] alpha;
  vector[nknots_year] B;
   matrix[nsites,nknots_year] b;
  vector[ncounts] noise;             // over-dispersion
   vector[nyears] year_effect;             // continental year-effects
 
  season_pred = season_basispred*(sdseason*B_season_raw);
  alpha = sdalpha*alpha_raw + beta_size*site_size;
  B = sdyear_gam*B_raw;
  noise = sdnoise*noise_raw;
  year_effect = sdyear*year_effect_raw;
  
 
    for(k in 1:nknots_year){
    b[,k] = (sdyear_gam_strat[k] * b_raw[,k]) + B[k];
  }
  
  
  Y_pred = year_basispred * B; 
  
      for(s in 1:nsites){
     year_pred[,s] = Y_pred + (year_basispred * transpose(b[s,]));
}

  for(i in 1:ncounts){
    E[i] = ALPHA1 + year_pred[year_raw[i],site[i]] + alpha[site[i]] + year_effect[year_raw[i]] + season_pred[date[i]] + noise[i];
  }
  
  }
  
  
model { 
  sdnoise ~ std_normal(); //prior on scale of extra Poisson log-normal variance
  sdyear ~ std_normal(); //prior on scale of site level variation
  sdalpha ~ std_normal(); //prior on scale of site level variation
  sdyear_gam ~ normal(0,0.2); //prior on sd of gam hyperparameters
  sdyear_gam_strat ~ normal(0,0.05); // regularizing prior on variance of stratum level gam
 //nu ~ gamma(2,0.1); // prior on df for t-distribution of heavy tailed site-effects from https://github.com/stan-dev/stan/wiki/Prior-Choice-Recommendations#prior-for-degrees-of-freedom-in-students-t-distribution
  sdseason ~ std_normal();//variance of GAM parameters
  B_season_raw ~ std_normal();//GAM parameters
  ALPHA1 ~ std_normal();// overall species intercept 
  
  beta_size ~ normal(0,1);// effect of site-size predictor
  
  count ~ poisson_log(E); //vectorized count likelihood
  alpha_raw ~ std_normal(); // fixed site-effects
  noise_raw ~ std_normal(); //heavy tailed extra Poisson log-normal variance
  B_raw ~ std_normal();// prior on GAM hyperparameters
  year_effect_raw ~ std_normal(); //prior on ▲annual fluctuations
  sum(year_effect_raw) ~ normal(0,0.0001*nyears);//sum to zero constraint on year-effects
  sum(alpha_raw) ~ normal(0,0.001*nsites);//sum to zero constraint on site-effects
  sum(B_raw) ~ normal(0,0.001*nknots_year);//sum to zero constraint on GAM hyperparameters
  
  // site level spatial iCAR to define time-series GAM parameters (s(year))
    for(k in 1:nknots_year){
  b_raw[,k] ~ icar_normal_lpdf(nsites, node1, node2);
  }
  
 //  for(s in 1:nstrata){
 // b_raw[s,] ~ std_normal();
 //  sum(b_raw[s,]) ~ normal(0,0.001*nknots_year);//sum to zero constraint on GAM hyperparameters
 //  
 //  }
}

generated quantities {

  real<lower=0> N[nyears];
  real<lower=0> NSmooth[nyears];
  real<lower=0> NSmooth_n[nyears];
  real<lower=0> nsmooth[nsites,nyears];
    real seas_max = max(season_pred);
 
 // log_lik calculation for looic
      vector[ncounts] log_lik;
  for(i in 1:ncounts){
  log_lik[i] = poisson_log_lpmf(count[i] | E[i]);
  }
  
  
      for(s in 1:nsites){

  for(y in 1:nyears){

     nsmooth[s,y] = exp(ALPHA1 + year_pred[y,s] + seas_max + alpha[s] + 0.5*(sdnoise^2) );
    }
  }
  
    for(y in 1:nyears){

      N[y] = exp(ALPHA1 + Y_pred[y] + year_effect[y] + seas_max + 0.5*(sdalpha^2) + 0.5*(sdnoise^2) );
      NSmooth[y] = exp(ALPHA1 + Y_pred[y] + seas_max + 0.5*(sdalpha^2) + 0.5*(sdnoise^2) );
      
      NSmooth_n[y] = mean(nsmooth[,y]);
    }
    
}

