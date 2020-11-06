// This is a full hierarchical GAM time-series with year-effects
// and a gam-based Seasonal adjustment


functions {
  real icar_normal_lpdf(vector bb, int nsites, int[] node1, int[] node2) {
    return -0.5 * dot_self(bb[node1] - bb[node2])
      + normal_lpdf(sum(bb) | 0, 0.001 * nsites); //soft sum to zero constraint on phi
 }
}


data {
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
  
 
  int<lower=1> N_edges;
  int<lower=1, upper=nsites> node1[N_edges];  // node1[i] adjacent to node2[i]
  int<lower=1, upper=nsites> node2[N_edges];  // and node1[i] < node2[i]

  int<lower=0> count[ncounts];              // count observations
  int<lower=1> site[ncounts];              // site indicators
  real year[ncounts];              // centered years
  int<lower=1> year_raw[ncounts]; // year index
  int<lower=1> date[ncounts];  // day indicator in the season
  

}

parameters {
  vector[nsites] alpha;             // intercepts
  vector[nyears] year_effect;             // continental year-effects
  vector[ncounts] noise;             // over-dispersion
  matrix[nsites,nknots_year] b_raw;         // spatial effect slopes (0-centered deviation from continental mean slope B)
  vector[nknots_year] B;             // GAM coefficients year
  
  vector[nknots_season] B_season;         // GAM coefficients
  
  real<lower=0> sigma[nknots_year];    // spatial standard deviation
  real<lower=0> sdnoise;    // sd of over-dispersion
 //real<lower=1> nu; 
  real<lower=0> sdalpha;    // sd of site effects
  real<lower=0> sdyear;    // sd of year effects
  real<lower=0> sdseason;    // sd of year effects
  real<lower=0> sdyear_gam;    // sd of GAM coefficients
  real ALPHA1; // overall intercept
}

transformed parameters { 
  matrix[nsites,nknots_year] b; //  
  vector[ncounts] E;           // log_scale additive likelihood
  vector[ndays] season_pred = season_basispred*B_season; //seasonal GAM
    matrix[nyears,nsites] year_pred;
  vector[nyears] Y_pred;  

  for(k in 1:nknots_year){
    b[,k] = (sigma[k] * b_raw[,k]) + B[k];
  }
  Y_pred = year_basispred * B; 
  
      for(s in 1:nsites){
     year_pred[,s] = year_basispred * transpose(b[s,]);
}

  for(i in 1:ncounts){
    E[i] = ALPHA1 + year_pred[year_raw[i],site[i]] + alpha[site[i]] + year_effect[year_raw[i]] + season_pred[date[i]] + noise[i];
  }
  
  }
model {
  sdnoise ~ normal(0,1); //prior on scale of extra Poisson log-normal variance
  sdyear ~ normal(0,1); //prior on scale of site level variation
  sdalpha ~ normal(0,1); //prior on scale of site level variation
  sdyear_gam ~ normal(0,0.2); //prior on sd of gam hyperparameters
  //nu ~ gamma(2,0.1); // prior on df for t-distribution of heavy tailed site-effects from https://github.com/stan-dev/stan/wiki/Prior-Choice-Recommendations#prior-for-degrees-of-freedom-in-students-t-distribution
  sdseason ~ normal(0,1);//variance of GAM parameters
  B_season ~ normal(0,sdseason);//GAM parameters
  ALPHA1 ~ normal(0,1);// overall species intercept
  
  count ~ poisson_log(E); //vectorized count likelihood
  alpha ~ normal(0,sdalpha); // fixed site-effects
  noise ~ normal(0,sdnoise); //heavy tailed extra Poisson log-normal variance
  B ~ normal(0,sdyear_gam);// prior on GAM hyperparameters
  sigma ~ normal(0,1); //prior on scale of spatial variation
  year_effect ~ normal(0,sdyear); //prior on scale of spatial variation
  sum(year_effect) ~ normal(0,0.0001*nyears);//sum to zero constraint on year-effects
  sum(alpha) ~ normal(0,0.001*nsites);//sum to zero constraint on site-effects
  sum(B) ~ normal(0,0.001*nknots_year);//sum to zero constraint on GAM hyperparameters
  sum(B_season) ~ normal(0,0.001*nknots_season);//sum to zero constraint on GAM seasonal parameters
  
  for(k in 1:nknots_year){
  b_raw[,k] ~ icar_normal_lpdf(nsites, node1, node2);
  }
}

generated quantities {
  real<lower=0> N[nyears];
  real<lower=0> NSmooth[nyears];
  real<lower=0> n[nsites,nyears];
  real<lower=0> nsmooth[nsites,nyears];
  
      for(s in 1:nsites){

  for(y in 1:nyears){

      n[s,y] = exp(ALPHA1 + year_pred[y,s] + year_effect[y] + season_pred[50] + alpha[s] + (0.5*(sdnoise^2)));
      nsmooth[s,y] = exp(ALPHA1 + year_pred[y,s] + season_pred[50] + alpha[s] + (0.5*(sdnoise^2)));
    }
  }
  
    for(y in 1:nyears){

      N[y] = exp(ALPHA1 + Y_pred[y] + year_effect[y] + season_pred[50] + (0.5*(sdnoise^2)) + (0.5*(sdalpha^2)));
      NSmooth[y] = exp(ALPHA1 + Y_pred[y] + season_pred[50] + (0.5*(sdnoise^2)) + (0.5*(sdalpha^2)));
      
    }
    
}

