# in overdispersed count data, every observation has its own count number:

library(tidyverse)

mydat <- data_frame(
  samplenumber = letters[1:20],
  onemean      = 16,
  pois_onemean = rpois(20, onemean),
  manymeans    = rnorm(20, mean = 16, sd = 4),
  pois_lotmean = rpois(20, manymeans)
)
mydat_long <- mydat %>% 
  select(samplenumber, pois_onemean, pois_lotmean) %>% 
  gather(meantype, val, -samplenumber)

mydat_long %>% 
  ggplot(aes(x  = val, fill = meantype)) + geom_histogram(binwidth = 3) + facet_wrap(~meantype, ncol = 1)

mydat_long
library(lme4)

vanilla_poisson <- glm(val ~ 1, data = filter(mydat_long, meantype == "pois_onemean"), family = "poisson")

summary(vanilla_poisson)
exp(2.7408)


vanilla_poisson_many <- glm(val ~ 1,
                            data = filter(mydat_long,
                                          meantype == "pois_lotmean"),
                            family = "poisson")

summary(vanilla_poisson_many)


vanilla_poisson_many_mix <- glmer(val ~ 1 + (1|samplenumber),
                            data = filter(mydat_long,
                                          meantype == "pois_lotmean"),
                            family = "poisson")
summary(vanilla_poisson_many_mix)

AIC(vanilla_poisson_many_mix, vanilla_poisson_many)

exp(2.65 + 0.12)

negbin_many_ <- glmer(val ~ 1 + (1|samplenumber),
                      data = filter(mydat_long,
                                    meantype == "pois_lotmean"),
                      family = "poisson")



# offsets -----------------------------------------------------------------


#' we are interested in the occurance of a strange kind of toadstool. They are 
#' found in farmer's fields. They are found at a constant density per hectare --
#' unfortunately, we don't have data at the hectare level, just at the level of
#' _farms_

farms <- data_frame(farmsize = rep(c(5, 15, 20), each = 5))

farms_mushrooms <- farms %>% 
  mutate(mushrooms = map_dbl(farmsize, ~sum(rpois(.x, lambda = 3))))

avg_shrooms <- glm(mushrooms ~ 1, family = poisson, data = farms_mushrooms)

summary(avg_shrooms)

exp(1.317)
#' so this approach gives us a mean of about 3 mushrooms per farm -- but it is badly overdispersed
#' 
mean(farms_mushrooms$mushrooms)
#' but not wrong! This gives us the average of that column. But we don't want
#' the average _per farm_, we want the average _per hectare_.
#' 

farms_mushrooms <- farms_mushrooms %>% 
  mutate(farmsize_log = log(farmsize),
         farmsize_scale = farmsize_log - mean(farmsize_log))

farms_mushrooms %>% 
  ggplot(aes(x = farmsize, y = mushrooms)) + geom_point(position = position_jitter(width = 1))

# this works fine! gets close to the answer
avg_shrooms_offset <- glm(mushrooms ~ 1 + offset(farmsize_log), family = poisson, data = farms_mushrooms)
summary(avg_shrooms_offset)

# Not what you want! gives the average response for the average group.
avg_shrooms_offset_mod <- glm(mushrooms ~ 1 + offset(farmsize_scale), family = poisson, data = farms_mushrooms)
summary(avg_shrooms_offset_mod)

#' See this [perfectly good SO answer](https://stats.stackexchange.com/questions/237963/how-to-formulate-the-offset-of-a-glm)

farms <- data_frame(farmsize = rep(c(5, 15, 20), each = 5))

farms_mushrooms <- farms %>% 
  mutate(mushrooms = map_dbl(farmsize, ~sum(rpois(.x, lambda = 3))),
         farmsize_log = log(farmsize))

farms_mushrooms

# avg_shrooms <- glm(mushrooms ~ 1, family = poisson, data = farms_mushrooms)
# 
# summary(avg_shrooms)
# this works fine! gets close to the answer
avg_shrooms_offset <- glm(mushrooms ~ 1+ offset(log(farmsize)), family = poisson, data = farms_mushrooms)
summary(avg_shrooms_offset)
exp(coef(avg_shrooms_offset))

# if you screw up and place the offset in the model, the coefficient it gets is close to 1 -- it is set to 1 automatically.
