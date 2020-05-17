// regularied horseshoe for microbiome data
// prior structure adapted from Piironen and Vehtari, "Sparsity information and regularization in the horseshoe and other shrinkage priors" 2017
// model parameterization adapted from Betancourt https://betanalpha.github.io/assets/case_studies/bayes_sparse_regression.html

data {
  int < lower = 0 > n; // total observations
  int < lower = 0 > n_seqvar_id; // number of sequence variants
  int < lower = 0 > seqvar_id[n]; // sequence variant ids
  //real < lower = 0, upper = 1 > ET_Sa[n]; // proportional abundance SA in ET
  //real < lower = 0, upper = 1 > AN_asv[n]; // proportional abundance SA sequence variant in AN
  real ET_Sa[n]; // proportional abundance SA in ET
  real AN_asv[n]; // proportional abundance SA sequence variant in AN
}

transformed data {
  real m0 = 10;           // Expected number of large slopes
  real slab_scale = 2;    // Scale for large slopes (Betancourt = 3; brms = 2)
  real slab_scale2 = square(slab_scale);
  real slab_df = 4;      // Effective degrees of freedom for large slopes (Betancourt = 25; brms = 4)
  real half_slab_df = 0.5 * slab_df;
}

parameters {
  vector[n_seqvar_id] beta;
  vector<lower=0>[n_seqvar_id] lambda;
  real<lower=0> c2_tilde;
  real<lower=0> tau_tilde;
  real alpha;
  real<lower=0> sigma;
}

transformed parameters {
  vector[n_seqvar_id] b_seqvar;
  {
    real tau0 = (m0 / (n_seqvar_id - m0)) / sqrt(1.0 * n);
    real tau = tau0 * tau_tilde; // tau ~ cauchy(0, tau0)
    
    // c2 ~ inv_gamma(half_slab_df, half_slab_df * slab_scale2)
    // Implies that marginally beta ~ student_t(slab_df, 0, slab_scale)
    real c2 = slab_scale2 * c2_tilde;
    
    vector[n_seqvar_id] lambda_tilde =
      sqrt( c2 * square(lambda) ./ (c2 + square(tau) * square(lambda)) );
    
    // beta_seqvar ~ normal(0, tau * lambda_tilde)
    b_seqvar = tau * lambda_tilde .* beta;
  }
}

model {
  vector[n] mu;
  beta ~ normal(0, 1);
  lambda ~ cauchy(0, 1);
  tau_tilde ~ cauchy(0, 1);
  c2_tilde ~ inv_gamma(half_slab_df, half_slab_df);
  alpha ~ normal(0, 1);
  sigma ~ exponential( 1 );
  for ( i in 1:n ) {
    mu[i] = beta[seqvar_id[i]] * AN_asv[i] + alpha;
    }
  ET_Sa ~ normal( mu , sigma );
}
