data{
    int<lower=1> N;
    int<lower=1> N_site_id;
    int<lower=1> N_spp_id;
    int obs_abd[N];
    int spp_id[N];
    int site_id[N];
}
parameters{
    real inter;
    vector[N_site_id] site_intercept;
    real<lower=0> sitevar;
    vector[N_spp_id] spp_intercept;
    real<lower=0> sppvar;
    // adding this part: latent variables? is that you?
    matrix[N_site_id, 2] latent_vars;
    matrix[2, N_spp_id] spp_loadings;
}
model{
    vector[N] lamb;
    sppvar ~ cauchy( 0 , 4 );
    spp_intercept ~ normal( 0 , sppvar );
    sitevar ~ cauchy( 0 , 4 );
    site_intercept ~ normal( 0 , sitevar );
    inter ~ normal( 0 , 5 );
    // latent variable priors??
    for (j in 1:N_site_id){
      for(k in 1:2 ){
        latent_vars[j, k] ~ normal(0, 1);
      }
    }
    // and the species loadings?
    for (j in 1:N_spp_id){
      for(k in 1:2 ){
        spp_loadings[k, j] ~ normal(0, 3);
      }
    }
    // add the row of 
    for ( i in 1:N ) {
        lamb[i] = inter + site_intercept[site_id[i]] + spp_intercept[spp_id[i]] + latent_vars[site_id[i],] * spp_loadings[,spp_id[i]];
    }
    obs_abd ~ poisson_log( lamb );
}
generated quantities{
    vector[N] lamb;
    real dev;
    dev = 0;
    for ( i in 1:N ) {
        lamb[i] = inter + site_intercept[site_id[i]] + spp_intercept[spp_id[i]];
    }
    dev = dev + (-2)*poisson_log_lpmf( obs_abd | lamb );
}
