data {
  int<lower=0> N;
  array[N] int<lower=0> y; // count outcomes
  vector<lower=0>[N] E; // exposure
  int<lower=1> K; // num covariates
  matrix[N, K] xs; // design matrix
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
  vector[N] theta; // heterogeneous random effects
  real<lower=0> sigma; // random effects variance
}
model {
  y ~ poisson_log(log_E + beta0 + xs_centered * betas + theta * sigma);
  beta0 ~ std_normal();
  betas ~ std_normal();
  theta ~ std_normal();
  sigma ~ normal(0, 5);
}
generated quantities {
  real beta_intercept = beta0 - dot_product(means_xs, betas);  // adjust intercept
  array[N] int y_rep;
  vector[N] log_lik;
  {
    vector[N] eta = log_E + beta0 + xs_centered * betas + theta * sigma;   // centered data
    y_rep = max(eta) < 26 ? poisson_log_rng(eta) : rep_array(-1, N);
    for (n in 1:N) {
      log_lik[n] = poisson_log_lpmf(y[n] | eta[n]);
    }
  }
}
