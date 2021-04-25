### estimating the mean relative trend among spatial units (strata and regions)


# PACKAGES ----------------------------------------------------------------


library(tidyverse)
library(tidybayes)
library(rstan)
rstan_options(auto_write = TRUE, javascript = FALSE)
library(shinystan)
library(sf)
library(spdep)
library(ggforce)
library(GGally)

source("Functions/Utility_functions.R")

# DATA LOAD --------------------------------------------------

# extract the centered trends
# run the comparison model on the centered trends
# map the mean trend in space
# plot the mean trends by region

trends_all <- read.csv("trends/All_strata_level_trends.csv")

trends_all$centered_log_trend_sd = ((trends_all$centered_log_trend_uci-trends_all$centered_log_trend_lci)/(1.96*2))

t_types = unique(trends_all$trend_type)

# Simple strata-only model ------------------------------------------------

for(time in t_types){
## in this model each species is treated as an independent measure within each
## stratum - 
## this is reasonable given that each species estimates are centered on the species
## mean trend

  trends <- filter(trends_all,trend_type == time)
trends$strataF <- factor(trends$region)
trends$strat <- as.integer(trends$strataF)
nstrata = max(trends$strat)
strata = trends$strat
t_hat = trends$centered_log_trend
sd_hat = trends$centered_log_trend_sd
N = nrow(trends)

# generate stan data ------------------------------------------------------

stan_data <- list(nstrata = nstrata,
                 N = N,
                 t_hat = t_hat,
                 sd_hat = sd_hat,
                 strata = strata)

mod.file = "models/trend_comparison.stan"



parms = c("mu_t",
          "mu_strata",
          "epsilon",
          "sd_strata_raw",
          "epsilon")


## compile model
model = stan_model(file=mod.file)

## run sampler on model, data
stanfit <- sampling(model,
                    data=stan_data,
                    verbose=TRUE, refresh=200,
                    chains=4, iter=2500,
                    warmup=1000,
                    cores = 4,
                    pars = parms,
                    control = list(adapt_delta = 0.8,
                                   max_treedepth = 15))






save(list = c("stanfit","stan_data","trends"),
     file = paste0("output/Strata_comparison",time,"_all_species.RData"))

#launch_shinystan(stanfit) 

}


# loading the base hexagon map
load("data/hexagon_grid.RData")
tp_out <- vector(mode = "list",length = length(t_types))
names(tp_out) <- t_types
mp_out <- tp_out


for(time in t_types){
  
  nyrs <- 15
  if(time == "Long-term"){
    nyrs <- (2019-1980)
  }
load(paste0("output/Strata_comparison",time,"_all_species.RData"))
mu_strata <- gather_draws(stanfit,mu_strata[strat]) %>% 
  mutate(pos = ifelse(.value > 0,TRUE,FALSE)) %>% 
  group_by(strat) %>% 
  summarise(mean = mean(.value),
            lci = quantile(.value,0.05),
            uci = quantile(.value,0.95),
            lqrt = quantile(.value,0.25),
            uqrt = quantile(.value,0.75),
            p_pos = sum(pos)/n(),
            p_neg = 1-(sum(pos)/n())) %>% 
  mutate(sig = ifelse(lci > 0 | uci < 0,TRUE,FALSE),
         sigq = ifelse(lqrt > 0 | uqrt < 0,TRUE,FALSE),
         trend = texp(exp(mean),ny = nyrs),
         trend_lci = texp(exp(lci),ny = nyrs),
         trend_uci = texp(exp(uci),ny = nyrs),
         trend_lqrt = texp(exp(lqrt),ny = nyrs),
         trend_uqrt = texp(exp(uqrt),ny = nyrs),
         p_not_zero = ifelse(p_pos > 0.5,p_pos,p_neg))

strat_sums <- trends %>% 
  select(strat,region,species) %>% 
  distinct() %>% 
  group_by(strat,region) %>% 
  summarise(nspecies = n()) %>% 
  left_join(.,mu_strata,by = "strat")


tp <- ggplot(data = strat_sums)+
  geom_point(aes(x = strat,y = trend,size = p_not_zero,colour = sig))+
  geom_errorbar(aes(x = strat,y = trend,ymin = trend_lci,ymax = trend_uci),
                alpha = 0.2,width = 0)+
geom_errorbar(aes(x = strat,y = trend,ymin = trend_lqrt,ymax = trend_uqrt),
              alpha = 0.2,width = 0,size = 1.3)+
geom_abline(slope = 0,intercept = 0)+
  labs(title = time)

tp_out[[time]] <- tp


# Mapping mean trends -----------------------------------------------------

mp <- trend_map_composite(trends = strat_sums,
                            map.file = "BBS_ProvState_strata",
                            hex_map = poly_grid,
                            #size_value = "Species Count",
                          size_value = "Probability > or < zero")


mp_out[[time]] <- mp


}#end time loop

pdf(file = "Figures/Comparison_graphs_all_species.pdf",
    width = 11,
    height = 8.5)
for(jj in 1:length(t_types)){
  print(tp_out[[jj]])
}
dev.off()


pdf(file = "Figures/Comparison_maps_all_species.pdf",
    width = 11,
    height = 8.5)
for(jj in 1:length(t_types)){
  print(mp_out[[jj]])
}
dev.off()





