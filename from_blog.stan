data {
  int<lower=1> N;                // number of 
  int<lower=1> P;                // number of 
  matrix[N,P] Y;                 // data matrix of order [N,P]
  int<lower=1> D;              // number of latent dimensions 
}
transformed data {
  int<lower=1> M;
  vector[P] mu;
  // We calculate the number of below-zero loadings 
  // the second term is the number in the "square part" below diag
  // the first term is the "non square part" below that
  M  = D * (P - D) + D * (D-1)/2;  // number of non-zero loadings
  mu = rep_vector(0.0,P);
}
parameters {    
  vector[M] L_t;   // lower diagonal elements of L
  vector<lower=0>[D] L_d;   // diagonal elements of L
  vector<lower=0>[P] psi;         // vector of variances
  real<lower=0>   mu_psi;
  real<lower=0>  sigma_psi;
  real   mu_lt;
  real<lower=0>  sigma_lt;
}
transformed parameters{
  cholesky_factor_cov[P,D] L;  //lower triangular factor loadings Matrix 
  cov_matrix[P] Q;   //Covariance mat
  {
    int idx2;
    idx2 = 1;
    for(i in 1:P){
      for(j in (i+1):D){
        L[i,j] = 0; //constrain the upper triangular elements to zero 
      }
    }
    for (j in 1:D) {
      L[j,j] = L_d[j];
      for (k in (j+1):P) {
        L[k,j] = L_t[idx2];
        idx2 = idx2 + 1;
      } 
    }
  } 
  Q = L*L' + diag_matrix(psi); 
}
model {
// the hyperpriors 
mu_psi ~ cauchy(0, 1);
sigma_psi ~ cauchy(0,1);
mu_lt ~ cauchy(0, 1);
sigma_lt ~ cauchy(0,1);
// the priors 
L_d ~ cauchy(0,3);
L_t ~ cauchy(mu_lt,sigma_lt);
psi ~ cauchy(mu_psi,sigma_psi);
//The likelihood
for( j in 1:N)
Y[j] ~ multi_normal(mu,Q); 
}
