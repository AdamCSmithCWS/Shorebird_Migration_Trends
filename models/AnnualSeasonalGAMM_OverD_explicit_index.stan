// generated with brms 2.10.0
functions {
}
data {
  //indicators
  int<lower=1> N;  // number of observations
  int<lower=1> nKnots_yr;  // number of observations
  int<lower=1> nKnots_date;  // number of observations
  int<lower=1> nStrata;  // number of geographic strata
   int<lower=1> nDates;  // number of days in the season
   int<lower=1> nYears;  // number of years
   int<lower=1> nSites[nStrata];  // number of sites in each stratum
  
  
  int count[N];  // response variable
  
  // predictors
    // basis function matrices for year (re-used for each stratum)
  matrix[nYears, nKnots_yr] basis_yr;
  
     // basis function matrices for date
  matrix[nDates, nKnots_date] basis_date;
 
 int strat[N]; //stratum indicators
  int site[N]; //site indicators - nested within stratum
 real yr_z[N]; //year predictor to match basis function matrix basis_yr
 real date_z[N]; //year predictor to match basis function matrix basis_date
 
  int prior_only;  // should the likelihood be ignored?
}
transformed data {
}
parameters {
  // temporary intercept for centered predictors
  real Intercept;
      // overdispersion variance
  real sd_noise;
  // spline coefficients year
  vector[nKnots_yr] bs_year;
    // standarized spline coefficients for year
  vector[nKnots_yr] zs_1_1;
  // standard deviations of the coefficients
  real<lower=0> sds_1_1;
  
      // standarized spline coefficients for year within strata
  matrix[nKnots_yr,nStrata] zs_2_1;
  // standard deviations of the spline coefficients for year within strata
  real<lower=0> sds_2_1;
  
  
  
    // spline coefficients date
  vector[nKnots_date] bs_date;
  
  // parameters for spline s(yr,k=13,id=1)

  // parameters for spline s(yr,by=strat,k=13,id=1)EastInland
  // standarized spline coefficients
  vector[knots_2[1]] zs_2_1;
  // standard deviations of the coefficients
  real<lower=0> sds_2_1;
  // parameters for spline s(yr,by=strat,k=13,id=1)Midcontinental
  // standarized spline coefficients
  vector[knots_3[1]] zs_3_1;
  // standard deviations of the coefficients
  real<lower=0> sds_3_1;
  // parameters for spline s(yr,by=strat,k=13,id=1)NortheastUSCoastal
  // standarized spline coefficients
  vector[knots_4[1]] zs_4_1;
  // standard deviations of the coefficients
  real<lower=0> sds_4_1;
  // parameters for spline s(yr,by=strat,k=13,id=1)Ontario
  // standarized spline coefficients
  vector[knots_5[1]] zs_5_1;
  // standard deviations of the coefficients
  real<lower=0> sds_5_1;
  // parameters for spline s(yr,by=strat,k=13,id=1)PacificandIntermountain
  // standarized spline coefficients
  vector[knots_6[1]] zs_6_1;
  // standard deviations of the coefficients
  real<lower=0> sds_6_1;
  // parameters for spline s(yr,by=strat,k=13,id=1)SoutheastCoastal
  // standarized spline coefficients
  vector[knots_7[1]] zs_7_1;
  // standard deviations of the coefficients
  real<lower=0> sds_7_1;
  // parameters for spline s(yr,by=strat,k=13,id=1)TexasCoastal
  // standarized spline coefficients
  vector[knots_8[1]] zs_8_1;
  // standard deviations of the coefficients
  real<lower=0> sds_8_1;
  // parameters for spline s(date,k=7)
  // standarized spline coefficients
  vector[knots_9[1]] zs_9_1;
  // standard deviations of the coefficients
  real<lower=0> sds_9_1;
    // means of Poisson including overdispersion
  vector[N] lambda;
}
transformed parameters {
  // actual spline coefficients
  vector[knots_1[1]] s_1_1 = sds_1_1 * zs_1_1;
  // actual spline coefficients
  vector[knots_2[1]] s_2_1 = sds_2_1 * zs_2_1;
  // actual spline coefficients
  vector[knots_3[1]] s_3_1 = sds_3_1 * zs_3_1;
  // actual spline coefficients
  vector[knots_4[1]] s_4_1 = sds_4_1 * zs_4_1;
  // actual spline coefficients
  vector[knots_5[1]] s_5_1 = sds_5_1 * zs_5_1;
  // actual spline coefficients
  vector[knots_6[1]] s_6_1 = sds_6_1 * zs_6_1;
  // actual spline coefficients
  vector[knots_7[1]] s_7_1 = sds_7_1 * zs_7_1;
  // actual spline coefficients
  vector[knots_8[1]] s_8_1 = sds_8_1 * zs_8_1;
  // actual spline coefficients
  vector[knots_9[1]] s_9_1 = sds_9_1 * zs_9_1;
  //add in the stratum by year predictions
}
model {
  // initialize linear predictor term
  vector[N] mu = Intercept + rep_vector(0, N) + Xs * bs + Zs_1_1 * s_1_1 + Zs_2_1 * s_2_1 + Zs_3_1 * s_3_1 + Zs_4_1 * s_4_1 + Zs_5_1 * s_5_1 + Zs_6_1 * s_6_1 + Zs_7_1 * s_7_1 + Zs_8_1 * s_8_1 + Zs_9_1 * s_9_1;
  // add the overdispersion variance
  for (n in 1:N)
  lambda[n] ~ normal(mu,sd_noise);
  // priors including all constants
  target += student_t_lpdf(Intercept | 3, -2, 10);
  target += normal_lpdf(zs_1_1 | 0, 1);
  target += student_t_lpdf(sds_1_1 | 3, 0, 10)
    - 1 * student_t_lccdf(0 | 3, 0, 10);
  target += normal_lpdf(zs_2_1 | 0, 1);
  target += student_t_lpdf(sds_2_1 | 3, 0, 10)
    - 1 * student_t_lccdf(0 | 3, 0, 10);
  target += normal_lpdf(zs_3_1 | 0, 1);
  target += student_t_lpdf(sds_3_1 | 3, 0, 10)
    - 1 * student_t_lccdf(0 | 3, 0, 10);
  target += normal_lpdf(zs_4_1 | 0, 1);
  target += student_t_lpdf(sds_4_1 | 3, 0, 10)
    - 1 * student_t_lccdf(0 | 3, 0, 10);
  target += normal_lpdf(zs_5_1 | 0, 1);
  target += student_t_lpdf(sds_5_1 | 3, 0, 10)
    - 1 * student_t_lccdf(0 | 3, 0, 10);
  target += normal_lpdf(zs_6_1 | 0, 1);
  target += student_t_lpdf(sds_6_1 | 3, 0, 10)
    - 1 * student_t_lccdf(0 | 3, 0, 10);
  target += normal_lpdf(zs_7_1 | 0, 1);
  target += student_t_lpdf(sds_7_1 | 3, 0, 10)
    - 1 * student_t_lccdf(0 | 3, 0, 10);
  target += normal_lpdf(zs_8_1 | 0, 1);
  target += student_t_lpdf(sds_8_1 | 3, 0, 10)
    - 1 * student_t_lccdf(0 | 3, 0, 10);
  target += normal_lpdf(zs_9_1 | 0, 1);
  target += student_t_lpdf(sds_9_1 | 3, 0, 10)
    - 1 * student_t_lccdf(0 | 3, 0, 10);
  // likelihood including all constants
  if (!prior_only) {
    target += poisson_log_lpmf(count | lambda);
  }
}
generated quantities {
  // actual population-level intercept
  real b_Intercept = Intercept;
}
