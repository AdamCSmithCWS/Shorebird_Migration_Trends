// This is a full hierarchical GAM time-series, with spatial gam parameters 
// as well as a gam-based Seasonal adjustment and sruvey-wide random year-effects
// same as GAMYE_strata but without the site-size predictor


functions {
  real icar_normal_lpdf(vector bb, int nsites, int[] node1, int[] node2) {
    return -0.5 * dot_self(bb[node1] - bb[node2])
      + normal_lpdf(sum(bb) | 0, 0.001 * nsites); //soft sum to zero constraint on phi
 }
}


data {
  int<lower=0> nstrata;
  int<lower=0> ncounts;
  int<lower=0> nsites;
  int<lower=0> nyears;
  
  // data for spline s(date)
  int<lower=1> ndays;  // number of days in the basis function for season
  int<lower=1> nknots_season;  // number of knots in the basis function for season
  matrix[ndays, nknots_season] season_basispred; // basis function matrix
  
  // data for spline s(year)
  int<lower=1> nknots_year;  // number of knots in the basis function for year
  matrix[nyears, nknots_year] year_basispred; // basis function matrix
  
  //data for spatial iCAR among strata
  int<lower=1> N_edges;
  int<lower=1, upper=nstrata> node1[N_edges];  // node1[i] adjacent to node2[i]
  int<lower=1, upper=nstrata> node2[N_edges];  // and node1[i] < node2[i]

  int<lower=0> count[ncounts];              // count observations
  int<lower=1> strat[ncounts];              // strata indicators
  int<lower=1> site[ncounts];              // site indicators
 // real year[ncounts];              // remove this
  int<lower=1> year_raw[ncounts]; // year index
  int<lower=1> date[ncounts];  // day indicator in the season
  //vector[nsites] site_size; //log-scale size predictor for site level abundance
  // upper 98th quartile of all counts across species with similar flocking behaviour
  // conceptually, this should reflect the maximum size of flock that could be observed
  // at this site, acting like an offset for area surveyed (since the area is hard to define)
  // and generally not available for most sites
  
  //indexes for re-scaling predicted counts within strata based on site-level intercepts
  int<lower=1> max_sites; //dimension 1 of sites matrix
  int<lower=0> sites[max_sites,nstrata]; //matrix of which sites are in each stratum
  int<lower=1> nsites_strat[nstrata]; //number of unique sites in each stratum
}

parameters {
  vector[nsites] alpha_raw;             // intercepts
  vector[nyears] year_effect_raw;             // continental year-effects
  vector[ncounts] noise_raw;             // over-dispersion
  matrix[nstrata,nknots_year] b_raw;         // spatial effect slopes (0-centered deviation from continental mean slope B)
  vector[nknots_year] B_raw;             // GAM coefficients year
  
  vector[nknots_season] B_season_raw;         // GAM coefficients
  //real beta_size; //effect of site level predictor
  
 real<lower=0> sdnoise;    // sd of over-dispersion
 //real<lower=2.1> nu; 
  real<lower=0> sdalpha;    // sd of site effects
  real<lower=0> sdyear;    // sd of year effects
  real<lower=0> sdseason;    // sd of year effects
  real<lower=0> sdyear_gam;    // sd of GAM coefficients
  real<lower=0> sdyear_gam_strat[nknots_year]; //sd of strata level gams for each knot
  real ALPHA1; // overall intercept
}

transformed parameters { 
  vector[ncounts] E;           // log_scale additive likelihood
  vector[ndays] season_pred = season_basispred*(sdseason*B_season_raw);
    matrix[nyears,nstrata] year_pred;
  vector[nyears] Y_pred; 
  vector[nsites] alpha;
  vector[nknots_year] B;
   matrix[nstrata,nknots_year] b;
  vector[ncounts] noise;             // over-dispersion
   vector[nyears] year_effect;             // continental year-effects
 
  alpha = sdalpha*alpha_raw;// + beta_size*site_size;
  B = sdyear_gam*B_raw;
  noise = sdnoise*noise_raw;
  year_effect = sdyear*year_effect_raw;
  
 
    for(k in 1:nknots_year){
    b[,k] = (sdyear_gam_strat[k] * b_raw[,k]) + B[k];
  }
  
  
  Y_pred = year_basispred * B; 
  
      for(s in 1:nstrata){
     year_pred[,s] = (year_basispred * transpose(b[s,]));
}

  for(i in 1:ncounts){
    E[i] = ALPHA1 + year_pred[year_raw[i],strat[i]] + alpha[site[i]] + year_effect[year_raw[i]] + season_pred[date[i]] + noise[i];
  }
  
  }
  
  
model { 
  sdnoise ~ std_normal(); //prior on scale of extra Poisson log-normal variance
  sdyear ~ normal(0,0.2); //prior on scale of annual fluctuations - 
  // above is informative so that 95% of the prior includes yearly fluctuations fall
  // between 33% decrease and a 50% increase
  sdalpha ~ std_normal(); //prior on scale of site level variation
  sdyear_gam ~ normal(0,1); //prior on sd of gam hyperparameters
  sdyear_gam_strat ~ gamma(2,4);//boundary avoiding prior (99% < 1.7) regularizing prior on variance of stratum level gam
 //nu ~ gamma(2,0.1); // prior on df for t-distribution of heavy tailed site-effects from https://github.com/stan-dev/stan/wiki/Prior-Choice-Recommendations#prior-for-degrees-of-freedom-in-students-t-distribution
  sdseason ~ std_normal();//variance of GAM parameters
  B_season_raw ~ std_normal();//GAM parameters
  ALPHA1 ~ std_normal();// overall species intercept 
 
   // beta_size ~ normal(0,1);// effect of site-size predictor

  count ~ poisson_log(E); //vectorized count likelihood
  alpha_raw ~ std_normal(); // fixed site-effects
  noise_raw ~ student_t(3,0,1);//std_normal(); // extra Poisson log-normal variance
  B_raw ~ std_normal();// prior on GAM hyperparameters
  year_effect_raw ~ std_normal(); //prior on ▲annual fluctuations
  sum(year_effect_raw) ~ normal(0,0.0001*nyears);//sum to zero constraint on year-effects
  sum(alpha_raw) ~ normal(0,0.0001*nsites);//sum to zero constraint on site-effects
  //sum(B_raw) ~ normal(0,0.001*nknots_year);//sum to zero constraint on GAM hyperparameters
  
    for(k in 1:nknots_year){
  b_raw[,k] ~ icar_normal(nstrata, node1, node2);
  }
  
 //  for(s in 1:nstrata){
 // b_raw[s,] ~ std_normal();
 //  sum(b_raw[s,]) ~ normal(0,0.001*nknots_year);//sum to zero constraint on GAM hyperparameters
 //  
 //  }
}

generated quantities {
  vector[nstrata] a;
  real<lower=0> N[nyears];
  real<lower=0> NSmooth[nyears];
  real<lower=0> n[nstrata,nyears];
  real<lower=0> nsmooth[nstrata,nyears];
    real seas_max = mean(season_pred);
 
 // log_lik calculation for looic
      vector[ncounts] log_lik;
  for(i in 1:ncounts){
  log_lik[i] = poisson_log_lpmf(count[i] | E[i]);
  }
  
  
      for(s in 1:nstrata){

        
  for(y in 1:nyears){
            real atmp[nsites_strat[s]];
            real atmp_smo[nsites_strat[s]];
            
        //a stratum-scaling component that tracks the alphas for sites in stratum
        for(j in 1:nsites_strat[s]){
          atmp[j] = exp(ALPHA1 + year_pred[y,s] + year_effect[y] + seas_max + 0.5*(sdnoise^2) + alpha[sites[j,s]]);
          atmp_smo[j] = exp(ALPHA1 + year_pred[y,s] + seas_max + 0.5*(sdnoise^2) + alpha[sites[j,s]]);
        }
        n[s,y] = mean(atmp);
        nsmooth[s,y] = mean(atmp_smo);
    }
  }
  
    for(y in 1:nyears){

      N[y] = exp(ALPHA1 + Y_pred[y] + year_effect[y] + seas_max + 0.5*(sdalpha^2) + 0.5*(sdnoise^2) );
      NSmooth[y] = exp(ALPHA1 + Y_pred[y] + seas_max + 0.5*(sdalpha^2) + 0.5*(sdnoise^2) );
      
    }
    
}

