# SPECIES MCMC data-prep -------------------------------------------------------
library(tidyverse)
library(rstan)
rstan_options(auto_write = TRUE, javascript = FALSE)
library(shinystan)
library(sf)
library(spdep)
library(ggforce)
library(tidybayes)


# library(doParallel)
# library(foreach)

#load data
load("data/allShorebirdPrismFallCounts.RData")
grid_spacing <- 300000  # size of squares, in units of the CRS (i.e. meters for lae)

 
FYYYY = 1980

w_cosewic = sps[c(2:4,7,10,12:20,22,11,25)]

# for(sp in w_cosewic){
#   
#   if(file.exists(paste0("output/",sp,"_GAMYE_strat_simple",grid_spacing/1000,".RData"))){w_cosewic <- w_cosewic[-which(w_cosewic == sp)]}
# 
# }  

sps_remain = sps[-which(sps %in% w_cosewic)]



 for(sp in sps_remain[1:6]){
  
   if(file.exists(paste0("output/",sp,"_GAMYE_strat_simple",grid_spacing/1000,".RData"))){next}
   
   
    load(paste0("data/data",sp,"_GAMYE_strat_simple",grid_spacing/1000,".RData"))



## compile model
slope_icar_model = stan_model(file=mod.file)

print(sp)
## run sampler on model, data
slope_icar_stanfit <- sampling(slope_icar_model,
                               data=stan_data,
                               verbose=TRUE, refresh=100,
                               chains=4, iter=1800,
                               warmup=1200,
                               cores = 4,
                               pars = c(parms,"nu"),
                               control = list(adapt_delta = 0.9,
                                              max_treedepth = 14))





save(list = c("slope_icar_stanfit",
              "stan_data",
              "dts",
              "real_grid",
              "strats_dts",
              "strat_regions",
              "mod.file",
              "parms"),
     file = paste0("output/",sp,"_GAMYE_strat_simple",grid_spacing/1000,".RData"))


}#end modeling loop



















launch_shinystan(slope_icar_stanfit) 



 

