data{
    int<lower=1> N;
    int<lower=1> N_site_id;
    int<lower=1> N_spp_id;
    int obs_abd[N];
    int spp_id[N];
    int site_id[N];
    int<lower=2> D;              // number of latent dimensions 
}
transformed data {
  int<lower=1> M;
  // We calculate the number of below-zero loadings 
  // the second term is the number in the "square part" below diag
  // the first term is the "non square part" below that
  M  = D * (N_spp_id - D) + D * (D-1)/2;  // number of non-zero loadings
}
parameters{
    real inter;
    vector[N_site_id] site_intercept;
    real<lower=0> sitevar;
    vector[N_spp_id] spp_intercept;
    real<lower=0> sppvar;
    // adding this part: latent variables? is that you?
    matrix[N_site_id, D] latent_vars;
    // matrix[2, N_spp_id] spp_loadings;
    // blog parameters
    vector[M] spp_load_L_t;   // lower diagonal elements of L = spp_loadings
    vector<lower=0>[D] spp_load_L_d;   // diagonal elements of L = spp_loadings
   // vector<lower=0>[N_spp_id] psi;         // vector of variances
    // hyperparameters
    //real<lower=0>   mu_psi;
    //real<lower=0>  sigma_psi;
    // real   mu_lt;
    // real<lower=0>  sigma_lt;
}
transformed parameters{
  cholesky_factor_cov[N_spp_id,D] spp_loadings;  //lower triangular factor loadings Matrix 
  {
    int idx2;
    idx2 = 1;
    for(i in 1:N_spp_id){
      for(j in (i+1):D){
        spp_loadings[i,j] = 0; //constrain the upper triangular elements to zero 
      }
    }
    for (j in 1:D) {
      // add the diagonal elements
      spp_loadings[j,j] = spp_load_L_d[j];
      for (k in (j+1):N_spp_id) {
        // add the lower triangular elements
        spp_loadings[k,j] = spp_load_L_t[idx2];
        idx2 = idx2 + 1;
      } 
    }
  } 
}
model{
    vector[N] lamb;
    // the hyperpriors 
    // mu_psi ~ cauchy(0, 1);
    // sigma_psi ~ cauchy(0,1);
    // mu_lt ~ cauchy(0, 1);
    // sigma_lt ~ cauchy(0,1);
    // the priors 
    // note that spp_load_L_d is constrained with lower=0
    spp_load_L_d ~ cauchy(0,3);
    spp_load_L_t ~ cauchy(0, 4);
    // psi ~ cauchy(mu_psi,sigma_psi);
    sppvar ~ cauchy( 0 , 4 );
    spp_intercept ~ normal( 0 , sppvar );
    sitevar ~ cauchy( 0 , 4 );
    site_intercept ~ normal( 0 , sitevar );
    inter ~ normal( 0 , 5 );
    // latent variable priors??
    for (j in 1:N_site_id){
      for(k in 1:D ){
        latent_vars[j, k] ~ normal(0, 1);
      }
    }
    // the liklihood
    // add the row of 
    for ( i in 1:N ) {
        lamb[i] = inter + site_intercept[site_id[i]] + spp_intercept[spp_id[i]] + latent_vars[site_id[i],] * spp_loadings[spp_id[i],]';
    }
    obs_abd ~ poisson_log( lamb );
}
generated quantities{
    vector[N] lamb;
    matrix[N_spp_id,N_spp_id] Q;   //Covariance mat
    real dev;
    dev = 0;
    for ( i in 1:N ) {
        lamb[i] = inter + site_intercept[site_id[i]] + spp_intercept[spp_id[i]] + latent_vars[site_id[i],] * spp_loadings[spp_id[i],]';
    }
    dev = dev + (-2)*poisson_log_lpmf( obs_abd | lamb );
    //add the correlation matrix in here! Q = L*L'
    Q = spp_loadings * spp_loadings'; //+ diag_matrix(psi);
    // OK this gives me an error when I declare that Q needs to be positive definite, by using cov_matrix[N_spp_id]
}
