// binomial likelihood with logit link
// regularied horseshoe for microbiome data
// model structure per Piironen and Vehtari, "Sparsity information and regularization in the horseshoe and other shrinkage priors" 2017
// model parameterization per Betancourt https://betanalpha.github.io/assets/case_studies/bayes_sparse_regression.html

data {
  int < lower = 0 > n; // total observations
  int < lower = 0 > n_seqvar_id; // number of sequence variants
  int < lower = 0 > seqvar_id[n]; // sequence variant ids
  real read_count[n]; // scaled read counts per sequence variant
  real specimen_read_total[n]; // scaled read counts per specimen
  int < lower = 0, upper = 1 > outcome[n]; // binomial outcome
  real asv_241b2b5f1ea0f1065ca7027534068198[n]; // scaled reads Staph aureus ASV of interest
  real asv_adcde8c78396022660fd6204471c65c8[n]; // scaled reads Staph aureus ASV of interest
  real asv_f011d78941f0c011aca75252ef604d17[n]; // scaled reads Staph aureus ASV of interest
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
  real b_total;
  real b_asv_241b;
  real b_asv_adcd;
  real b_asv_f011;
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
    
    b_seqvar = tau * lambda_tilde .* beta;
  }
}

model {
  vector[n] p;
  beta ~ normal(0, 1);
  lambda ~ cauchy(0, 1);
  tau_tilde ~ cauchy(0, 1);
  c2_tilde ~ inv_gamma(half_slab_df, half_slab_df);
  b_asv_241b ~ normal(0, 1);
  b_asv_adcd ~ normal(0, 1);
  b_asv_f011 ~ normal(0, 1);
  b_total ~ normal(0, 1);
  alpha ~ normal(0, 1);
  for ( i in 1:n ) {
        p[i] = beta[seqvar_id[i]] * read_count[i] + b_total * specimen_read_total[i] + b_asv_241b * asv_241b2b5f1ea0f1065ca7027534068198[i] + b_asv_adcd * asv_adcde8c78396022660fd6204471c65c8[i] + b_asv_f011 * asv_f011d78941f0c011aca75252ef604d17[i] + alpha;
    }
  outcome ~ binomial_logit(1, p);
}

