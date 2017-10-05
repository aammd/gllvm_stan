# simple case where species are correlated
library(rethinking)


# normal distribution -- its biomass or something!

sppcors <- rethinking::rlkjcorr(1, K = 5, eta = 1.5)

sppvars <- runif(5, min = 1, max = 1.5)

sppmeans <- rlnorm(5, log(15))

spp_lambdas <- MASS::mvrnorm(n = 10, mu = sppmeans, Sigma = diag(sppvars) %*% sppcors %*% diag(sppvars))

library(tidyverse)

species_sites <- spp_lambdas %>% 
  as.data.frame %>%
  set_names(paste0("sp", LETTERS[1:5])) %>% 
  rownames_to_column("site") %>% 
  gather(spp, true_mean, -site) %>% 
  arrange(site) %>% 
  mutate(obs_abd = map_dbl(true_mean, rpois, n = 1))


library(vegan)

sol_spp <- species_sites %>% 
  select(site, spp, obs_abd) %>% 
  spread(spp, obs_abd) %>% 
  remove_rownames %>% 
  column_to_rownames("site") %>% 
  as.matrix %>% 
  decostand("hellinger") %>% 
  metaMDS()

plot(sol_spp, type = "t")

comm_data <- species_sites %>% 
  select(site, spp, obs_abd) %>% 
  as.data.frame() %>% 
  mutate(site_id = as.numeric(site),
         spp_id = rethinking::coerce_index(spp))


# nah i don’t get this part -----------------------------------------------

# could also add a species intercept??
corr_mod <- alist(
  obs_abd   ~  dpois(lamb),
  log(lamb) <- inter + spp_eff[spp_id] + site_eff[site_id] + a[spp_id] * b[site_id] + a2[spp_id] * b2[site_id] ,# + n1 * h1 + n2 * h2,
  spp_eff[spp_id] ~ dnorm(0, spvar),
  site_eff[site_id] ~ dnorm(0, sitevar),
  a[spp_id]~dnorm(0,5),
  b[site_id]~dnorm(0,5),
  a2[spp_id]~dnorm(0,5),
  b2[site_id]~dnorm(0,5),
  spvar     ~ dcauchy(0, 2),
  sitevar   ~ dcauchy(0, 2),
  inter     ~ dnorm(0, 5)
  # c(n1, n2, h1, h2) ~ dnorm(0, 1)
  )

library(rstan)
test_model <- rethinking::map2stan(corr_mod, comm_data)

rethinking::precis(test_model, depth = 2)

test_samples <- rethinking::extract.samples(test_model)

test_samples %>% str

species_names <- comm_data %>% 
  select(spp, spp_id) %>% 
  distinct %>% 
  {.[["spp"]][order(.[["spp_id"]])]}
site_names <- comm_data %>% 
  select(site, site_id) %>% 
  distinct %>% glimpse %>% 
  {.[["site"]][order(.[["site_id"]])]}

test_medians <- test_samples[c("a", "b", "a2", "b2")] %>% 
  map(apply, 2, median)

site_data <- test_medians[c("b", "b2")] %>% 
  as_data_frame() %>%
  cbind(site_names)

site_data %>% 
  ggplot(aes(x = b, y = b2)) + geom_label(aes(label = site_names))

spp_data <- test_medians[c("a", "a2")] %>% 
  as_data_frame() %>%
  cbind(species_names)

site_data %>% 
  ggplot(aes(x = b, y = b2)) + 
  geom_label(aes(label = site_names)) + 
  geom_label(aes(x = a, y = a2, label = species_names), data = spp_data)

# does not look like a coplot -- because they are on different scales!

head(site_data)
head(spp_data)

site_data2 <- site_data[c("b", "b2")] * 
  matrix(sqrt(colSums(spp_data[c("a", "a2")] ^ 2)) / 
           sqrt(colSums(site_data[c("b", "b2")] ^ 2)), 
         nrow(site_data), 2, byrow = TRUE)

spp_data2 <- spp_data[c("a", "a2")] * 
  matrix(sqrt(colSums(site_data2[c("b", "b2")] ^ 2)) / 
           sqrt(colSums(spp_data[c("a", "a2")] ^ 2)), 
         nrow(spp_data), 2, byrow = TRUE)


site_data2 %>% 
  add_column(site_names = site_data$site_names) %>% 
  ggplot(aes(x = b, y = b2)) + 
  geom_label(aes(label = site_names)) + 
  geom_label(aes(x = a, y = a2, label = species_names),
             data = spp_data2 %>% 
               add_column(species_names = spp_data$species_names))

lv_coefs <- do.call(cbind, test_medians[c("a", "a2")])


# why are these numbers so much smaller than the original??
corrplot::corrplot(cov2cor(lv_coefs %*% t(lv_coefs) + diag(5)))



# OK how about the VERY simpliest form -- multivariate random numbers:
