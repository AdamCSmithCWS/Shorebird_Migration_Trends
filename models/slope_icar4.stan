// This is the same as slope_icar.stan but with fixed site effects
// integrate the annualSlopeSeasonalGAM_ZIP model with this iCar component on the slope parameters
// initially ignore season and the zip component
// just hte count model
// over dispersion
// site-effects

// next add

// year-effects (at what scale?) maybe continental since they're largely drive by annual recruitment and not clearly regional

// then
// season and potentially zip?


functions {
  real icar_normal_lpdf(vector b_raw, int nstrata, int[] node1, int[] node2) {
    return -0.5 * dot_self(b_raw[node1] - b_raw[node2])
      + normal_lpdf(sum(b_raw) | 0, 0.001 * nstrata); //soft sum to zero constraint on phi
 }
}
data {
  int<lower=0> nstrata;
  int<lower=0> ncounts;
  int<lower=0> nsites;
  int<lower=0> nyears;
 
  int<lower=0> N_edges;
  int<lower=1, upper=nstrata> node1[N_edges];  // node1[i] adjacent to node2[i]
  int<lower=1, upper=nstrata> node2[N_edges];  // and node1[i] < node2[i]

  int<lower=0> count[ncounts];              // count observations
  int<lower=0> strat[ncounts];              // strata indicators
  int<lower=0> site[ncounts];              // site indicators
  real year[ncounts];              // centered years
  int<lower=1> year_raw[ncounts]; // year index
}

parameters {
  vector[nsites] alpha;             // intercepts
  vector[nyears] year_effect;             // continental year-effects
  vector[ncounts] noise;             // over-dispersion
  vector[nstrata] b_raw;         // spatial effect slopes (0-centered deviation from continental mean slope B)
  real B;             // slope continental mean slope
  real<lower=0> sigma;    // spatial standard deviation
  real<lower=0> sdnoise;    // sd of over-dispersion
 //real<lower=1> nu; 
  real<lower=0> sdalpha;    // sd of site effects
  real<lower=0> sdyear;    // sd of year effects
}
transformed parameters { 
  vector[nstrata] b = (sigma * b_raw) + B; //  
  vector[ncounts] E;           // log_scale additive likelihood
  
  for(i in 1:ncounts){
    E[i] = b[strat[i]] * year[i] + alpha[site[i]] + year_effect[year_raw[i]] + noise[i];
  }
  
  }
model {
  sdnoise ~ normal(0,1); //prior on scale of extra Poisson log-normal variance
  sdyear ~ normal(0,1); //prior on scale of site level variation
  //boundary avoiding prior following Chung, Yeojin, Sophia Rabe-Hesketh, Vincent Dorie, Andrew Gelman, and Jingchen Liu. 2013. “A Nondegenerate Penalized Likelihood Estimator for Variance Parameters in Multilevel Models.” Psychometrika 78 (4): 685–709.
  sdalpha ~ normal(0,1); //prior on scale of site level variation
  //nu ~ gamma(2,0.1); // prior on df for t-distribution of heavy tailed site-effects from https://github.com/stan-dev/stan/wiki/Prior-Choice-Recommendations#prior-for-degrees-of-freedom-in-students-t-distribution
  count ~ poisson_log(E); //vectorized count likelihood
  alpha ~ normal(0,sdalpha); // fixed site-effects
  noise ~ normal(0,sdnoise); //heavy tailed extra Poisson log-normal variance
  B ~ normal(0,0.05);// prior on continental mean slope - informative because trends > 10%/year (abs(B) > 0.1) are very unlikely (> 2 SD from mean)
  sigma ~ normal(0,1); //prior on scale of spatial variation
  year_effect ~ normal(0,sdyear); //prior on scale of spatial variation
  sum(year_effect) ~ normal(0,0.0001*nyears);//sum to zero constraint on year-effects
  sum(alpha) ~ normal(0,0.001*nsites);//sum to zero constraint on site-effects
  
  b_raw ~ icar_normal_lpdf(nstrata, node1, node2);
}

