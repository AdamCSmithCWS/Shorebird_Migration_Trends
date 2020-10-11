// to do 
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
  real icar_normal_lpdf(vector phi, int nstrata, int[] node1, int[] node2) {
    return -0.5 * dot_self(phi[node1] - phi[node2])
      + normal_lpdf(sum(phi) | 0, 0.001 * nstrata); //soft sum to zero constraint on phi
 }
}
data {
  int<lower=0> nstrata;
  int<lower=0> ncounts;
  int<lower=0> nsites;
 
  int<lower=0> N_edges;
  int<lower=1, upper=nstrata> node1[N_edges];  // node1[i] adjacent to node2[i]
  int<lower=1, upper=nstrata> node2[N_edges];  // and node1[i] < node2[i]

  int<lower=0> count[ncounts];              // count observations
  int<lower=0> strat[ncounts];              // strata indicators
  int<lower=0> site[ncounts];              // site indicators
  real year[ncounts];              // centered years
}

parameters {
  vector[nsites] alpha;             // intercepts
  vector[ncounts] noise;             // over-dispersion
  real B;             // slope continental mean slope
  real<lower=0> sigma;    // spatial standard deviation
  real<lower=0> sdnoise;    // sd of over-dispersion
  vector[nstrata] b_raw;         // spatial effect slopes (0-centered deviation from continental mean slope B)
  real<lower=0> nu; 
  real<lower=0> sdalpha;    // sd of site effects
}
transformed parameters { 
  vector[nstrata] b = (sigma * b_raw) + B; //  
  vector[ncounts] E;           // log_scale additive likelihood
  
  for(i in 1:ncounts){
    E[i] = b[strat[i]] * year[i] + alpha[site[i]] + noise[i];
  }
  
  }
model {
  sdnoise ~ cauchy(0,5); //prior on scale of extra Poisson log-normal variance
  sdalpha ~ cauchy(0,5); //prior on scale of site level variation
  nu ~ normal(10,10); // prior on df for t-distribution of heavy tailed site-effects
  count ~ poisson_log(E); //vectorized count likelihood
  noise ~ normal(0,sdnoise); //extra Poisson log-normal variance
  alpha ~ student_t(nu,0,sdalpha); //heavy tailed site-effects
  
  target += -3*log(sigma) - 1/(sigma)^2;  // Stan equiv of BUGS model prior on tau
  b_raw ~ icar_normal_lpdf(nstrata, node1, node2);
}
generated quantities {
  real tau = sigma^-2; //variance on the spatial process
}
