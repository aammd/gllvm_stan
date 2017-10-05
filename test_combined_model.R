library(rstan)
library(rethinking)
library(tidyverse)
library(vegan)
library(rstanarm)
library(ggplot2)


number_of_species <- 3

set.seed(4812)
# species correlation
sppcors <- rethinking::rlkjcorr(1, K = number_of_species, eta = 1.5)

# variance of each species
# sppvars <- runif(5, min = 1, max = 1.5)
# actually not using this; this is a complication that can be added later if I want. 

# lognormal abundances
sppmeans <- rlnorm(number_of_species, log(17))

## sample the species -- using a constant variace in each species of 0.7
spp_lambdas <- MASS::mvrnorm(n = 100, mu = sppmeans,
                             Sigma = 
                               diag(0.7, number_of_species) %*% 
                               sppcors %*% 
                               diag(0.7, number_of_species))

# reshaping the data: the dataset needs a bit of massaging to be useful later
species_sites <- spp_lambdas %>% 
  as.data.frame %>%
  set_names(paste0("sp", LETTERS[1:number_of_species])) %>% 
  rownames_to_column("site") %>% 
  gather(spp, true_mean, -site) %>% 
  arrange(site) %>% 
  # I are using a multivariate normal distribution to sample "average abundances". Actually sampling abundances from a poisson like normal.  I suppose these will be overdispersed, in fact, because the average that they are based on is itself variable. 
  mutate(obs_abd = map_dbl(true_mean, rpois, n = 1)) %>% 
  # add index for species
  mutate(spp_id  = rethinking::coerce_index(spp),
         site_id = rethinking::coerce_index(site))

knitr::kable(head(species_sites))




data_for_stan <- list(obs_abd  = species_sites$obs_abd,
                      site_id  = species_sites$site_id,
                      spp_id   = species_sites$spp_id,
                      N = nrow(species_sites),
                      N_site_id = length(unique(species_sites$site_id)),
                      N_spp_id = length(unique(species_sites$spp_id)),
                      D = 2)

model_lv_stan <- stan_model("combined_blog_rethinking_model.stan")


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


results_lv_stan <- sampling(model_lv_stan,data = data_for_stan, chains = 1, init = 0)


# , data = data_for_stan, chains = 2

saveRDS(results_lv_stan, file = "results_lv_stan.rds")

results_lv_stan <- readRDS("results_lv_stan.rds")

list_of_draws <- rstan::extract(results_lv_stan)

str(list_of_draws$spp_loadings)

median_loadings <- apply(list_of_draws$spp_loadings, c(2,3), median)

## i think this is supposed to be an estimate of the correlation matrix, 
##.. but the diagonal is not even close to 1

median_loadings %*% t(median_loadings)
## compare with the real data..
sppcors
# well... the signs are in the right direction anyways??

median_latent <- apply(list_of_draws$latent_vars, c(2,3), median)

species_obs_matrix <- species_sites %>% 
  dplyr::select(site, spp, obs_abd) %>% 
  spread(spp, obs_abd) %>% 
  select(-site)

species_obs_matrix %>% 
  prcomp %>%
  biplot

# via vegan
library(vegan)
species_obs_matrix %>% 
  decorana() %>% 
  plot(display = "sites")

decor_sites <- species_obs_matrix %>% 
  decorana()

scores(decor_sites)

post_crustes <- procrustes(decor_sites, median_latent, scale = TRUE)

plot(post_crustes)

par(mfrow = c(1, 2))
plot(decor_sites, display = "sites")
title("DCA of simulated data")
plot(post_crustes$Yrot[,1:2])
title("a best guess??")


