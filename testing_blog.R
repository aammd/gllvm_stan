library("MASS")
set.seed(42)
D <-3
P <- 10 
N <-7

mu_theta <-rep(0,D) # the mean of eta
mu_epsilon<-rep(0,P) # the mean of epsilon
Phi<-diag(rep(1,D))
Psi <- diag(c(0.2079, 0.19, 0.1525, 0.20, 0.36, 0.1875, 0.1875, 1.00, 0.27, 0.27))
l1 <- c(0.99, 0.00, 0.25, 0.00, 0.80, 0.00, 0.50, 0.00, 0.00, 0.00)
l2 <- c(0.00, 0.90, 0.25, 0.40, 0.00, 0.50, 0.00, 0.00, -0.30, -0.30)
l3<-  c(0.00, 0.00, 0.85, 0.80, 0.00, 0.75, 0.75, 0.00, 0.80, 0.80)
L <-cbind(l1,l2,l3) # the loading matrix

Theta <-mvrnorm(N, mu_theta, Phi) # sample factor scores
Epsilon <-mvrnorm(N, mu_epsilon, Psi) # sample error vector
Y<-Theta%*%t(L)+Epsilon# generate observable data
pairs(Y)


### fit model


# a function to generate intial values that are slightly jittered for each chain.

init_fun <- function(){
  init_values <- list(L_t = rep(0,24),
                      L_d = rep(.5,D),
                      psi = rep(.2,P),
                      sigma_psi=0.15,
                      mu_psi=0.2,
                      sigma_lt=0.5,
                      mu_lt=0.0)
  return(init_values)
}



library(rstan)

blog_model <- stan_model("from_blog.stan")

fa.data <-list(P=P,N=N,Y=Y,D=D)

blog_results <- sampling(blog_model, chains = 1, init = init_fun)


# two ways of getting size of lower corner:

(5 * (5 - 1)) / 2

D*(P-D) + D*(D-1)/2