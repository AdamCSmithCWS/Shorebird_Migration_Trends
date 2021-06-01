# SPECIES MCMC data-prep -------------------------------------------------------
library(tidyverse)
library(cmdstanr)
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



 for(sp in sps_remain[4:6]){
  
   if(file.exists(paste0("output/",sp,"_GAMYE_fit_new.RData"))){next}
   
   
    load(paste0("data/data",sp,"_GAMYE_strat_simple_new.RData"))

    mod.file = "models/GAMYE_strata_two_season_gammaprior_beta_normal.stan"
    prior = "gamma"
    noise_dist = "normal"
    
    mod.file = "models/GAMYE_strata_two_season_gammaprior_beta.stan"
    prior = "gamma"
    noise_dist = "t"
    
    # mod.file = "models/GAMYE_strata_two_season_tprior_beta_normal.stan"
    # prior = "t"
    # 
    
    # mod.file = "models/GAMYE_strata_tprior_beta_normal.stan"
    # prior = "t"
    # 
    # mod.file = "models/GAMYE_strata_tprior_beta.stan"
    # prior = "t"
    
    mod.file = "models/GAMYE_strata_gammaprior_beta_normal.stan"
    prior = "gamma"
    noise_dist = "normal"
    
    mod.file = "models/GAMYE_strata_gammaprior_beta.stan"
    prior = "gamma"
    noise_dist = "t"
    
    
    ## compile model
    modl = cmdstan_model(stan_file=mod.file)

print(sp)
## run sampler on model, data
# slope_icar_stanfit <- sampling(slope_icar_model,
#                                data=stan_data,
#                                verbose=TRUE, refresh=100,
#                                chains=4, iter=1800,
#                                warmup=1200,
#                                cores = 4,
#                                pars = c(parms),
#                                control = list(adapt_delta = 0.9,
#                                               max_treedepth = 14))

cmdstanfit<- modl$sample(data=stan_data,
               refresh=100,
               chains=4, iter_sampling =800,
               iter_warmup=1000,
               parallel_chains = 4,
               max_treedepth = 15,
               adapt_delta = 0.95)

spf = gsub(sp,pattern = " ",replacement = "_")

cmdstanfit$save_output_files(dir = "output",
                             basename = paste0(spf,"-",prior),
                             timestamp = FALSE,
                             random = FALSE)
csvfl = paste0(getwd(),"/output/",spf,"-",prior,"-",1:4,".csv")


if(shiny_explore){
  sl_rstan <- rstan::read_stan_csv(csvfl)
  launch_shinystan(as.shinystan(sl_rstan))
}


cmdstanfit$save_object(file = paste0("output/",spf,"_",prior,".RDS"))


save(list = c("stan_data",
              "dts",
              "real_grid",
              "strats_dts",
              "strat_regions",
              "mod.file",
              "parms",
              "sp",
              "spf"),
     file = paste0("output/",spf,"_",prior,"_GAMYE_fit_new.RData"))


}#end modeling loop


# save csv files ----------------------------------------------------------






# Alternate model for SESA ------------------------------------------------

# 
# 
# sp <- w_cosewic[14]
#    
#    #if(file.exists(paste0("output/",sp,"_GAMYE_strat_simple",grid_spacing/1000,".RData"))){next}
#    
#    
#    load(paste0("data/data",sp,"_GAMYE_strat_simple",grid_spacing/1000,".RData"))
#    
#    mod.file <- "models/GAMYE_strata_simple_boundary_avoid.stan"
#    
#    ## compile model
#    slope_icar_model = stan_model(file=mod.file)
#    
#    print(sp)
#    ## run sampler on model, data
#    slope_icar_stanfit <- sampling(slope_icar_model,
#                                   data=stan_data,
#                                   verbose=TRUE, refresh=100,
#                                   chains=4, iter=1800,
#                                   warmup=1200,
#                                   cores = 4,
#                                   pars = c(parms,"seas_max","seas_max2"),
#                                   control = list(adapt_delta = 0.9,
#                                                  max_treedepth = 15))
#    
#    
#    
#    
#    
#    save(list = c("slope_icar_stanfit",
#                  "stan_data",
#                  "dts",
#                  "real_grid",
#                  "strats_dts",
#                  "strat_regions",
#                  "mod.file",
#                  "parms"),
#         file = paste0("output/",sp,"_GAMYE_strat_simple_boundary_avoid",grid_spacing/1000,".RData"))
#    
# 
# 
# 
# 
# 
# 
# 
# 
# 
# loo_out = loo(slope_icar_stanfit)
# 
# 
# 
# launch_shinystan(slope_icar_stanfit) 
# 
# 

 

