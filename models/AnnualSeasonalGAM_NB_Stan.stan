//
// This Stan program defines a single stratum temporal spline and seasonal spline for the shorebird migration monitoring data
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

// The input data is a vector 'y' of length 'N'.
data {
  int<lower=1> N;  // number of observations
  int Y[N];  // response variable
  // data for splines
  int Ks;  // number of linear effects
  matrix[N, Ks] Xs;  // design matrix for the linear effects
  
  
  // data for spline s(yr)
  int nb_1;  // number of bases
  int knots_1[nb_1];  // number of knots
  // basis function matrices
  matrix[N, knots_1[1]] Zs_1_1;
  
  
  // data for spline s(date)
  int nb_2;  // number of bases
  int knots_2[nb_2];  // number of knots
  // basis function matrices
  matrix[N, knots_2[1]] Zs_2_1;
  
  
  int prior_only;  // should the likelihood be ignored?
}

// The parameters accepted by the model. 
parameters {
 // temporary intercept for centered predictors
  real Intercept;
  // spline coefficients
  vector[Ks] bs;
  // parameters for spline s(yr)
  // standarized spline coefficients
  vector[knots_1[1]] zs_1_1;
  // standard deviations of the coefficients
  real<lower=0> sds_1_1;
  // parameters for spline s(date)
  // standarized spline coefficients
  vector[knots_2[1]] zs_2_1;
  // standard deviations of the coefficients
  real<lower=0> sds_2_1;

}

transformed parameters {
  // actual spline coefficients
  vector[knots_1[1]] s_1_1 = sds_1_1 * zs_1_1;
  // actual spline coefficients
  vector[knots_2[1]] s_2_1 = sds_2_1 * zs_2_1;
}

transformed data {
}

// The model to be estimated. 
model {
  // initialize linear predictor term
  vector[N] mu = Intercept + rep_vector(0, N) + Xs * bs + Zs_1_1 * s_1_1 + Zs_2_1 * s_2_1;
  // priors including all constants
  target += student_t_lpdf(Intercept | 3, 0, 10);
  target += normal_lpdf(zs_1_1 | 0, 1);
  target += student_t_lpdf(sds_1_1 | 3, 0, 10)
    - 1 * student_t_lccdf(0 | 3, 0, 10);
  target += normal_lpdf(zs_2_1 | 0, 1);
  target += student_t_lpdf(sds_2_1 | 3, 0, 10)
    - 1 * student_t_lccdf(0 | 3, 0, 10);
  // likelihood including all constants
  if (!prior_only) {
    target += poisson_log_lpmf(Y | mu);
  }
}

generated quantities {
  // actual population-level intercept
  real b_Intercept = Intercept;
}