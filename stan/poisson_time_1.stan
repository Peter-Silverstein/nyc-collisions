data {
  int<lower=0> N;
  array[N] int<lower=0> y; // count outcomes
  vector<lower=0>[N] E; // population
  int<lower=1> K; // num covariates
  matrix[N, K] xs; // design matrix
  
  // Census Tract Variables
  int <lower=1> I; // num tracts
  array[N] int<lower=1, upper=I> i_index; // tract index
  
  // Time Variables
  int <lower=1> J; // num time periods
  array[N] int<lower=1, upper=J> j_index; // time period index
}
transformed data {
  vector[N] log_E = log(E);
  
  // center continuous predictors 
  vector[K] means_xs;  // column means of xs before centering
  matrix[N, K] xs_centered;  // centered version of xs
  for (k in 1:K) {
    means_xs[k] = mean(xs[, k]);
    xs_centered[, k] = xs[, k] - means_xs[k];
  }
}
parameters {
  real beta0; // intercept
  vector[K] betas; // covariates
  vector[I] phi_1; // tract-level random effects
  vector[J] phi_2; // time-level random effects
  real<lower=0> sigma; // variance
  simplex[2] rho; // proportion of variance assigned to each source
}
transformed parameters {
  
  // combined temporal effects and random effects
  vector[N] gamma = (sqrt(rho[1]) * phi_1[i_index] + sqrt(rho[2]) * phi_2[j_index]) * sigma;
  
}
model {
  y ~ poisson_log(log_E + beta0 + xs_centered * betas + gamma);
  beta0 ~ normal(0, 0.5);
  betas ~ std_normal();
  phi_1 ~ std_normal();
  phi_2 ~ std_normal();
  sigma ~ normal(0, 5);
  rho ~ dirichlet(rep_vector(2, 2));
}
generated quantities {
  real beta_intercept = beta0 - dot_product(means_xs, betas);  // adjust intercept
  array[N] int y_rep;
  vector[N] log_lik;
  {
    vector[N] eta = log_E + beta0 + xs_centered * betas + gamma;   // centered data
    y_rep = max(eta) < 26 ? poisson_log_rng(eta) : rep_array(-1, N);
    for (n in 1:N) {
      log_lik[n] = poisson_log_lpmf(y[n] | eta[n]);
    }
  }
}
