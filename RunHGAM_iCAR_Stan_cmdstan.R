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
output_dir <- "g:/Shorebird_Migration_Trends/output"
#output_dir2 <- "g:/Shorebird_Migration_Trends/output"

 for(sp in sps){
   load(paste0("data/data",sp,"_cmdstanr_data.RData"))
   spf = gsub(sp,pattern = " ",replacement = "_")
   spf = gsub(pattern = "\'",replacement = "",
              x = spf)
   sp_file_name <- paste0(spf,"-",prior,"-",noise_dist2)
   
   if(file.exists(paste0(output_dir,"/",sp_file_name,"-",1,".csv"))){next}
   
   
 

    
    ## compile model
    modl = cmdstan_model(stan_file=mod.file2)

print(sp)


cmdstanfit<- modl$sample(data=stan_data,
               refresh=100,
               chains=4, iter_sampling =800,
               iter_warmup=1000,
               parallel_chains = 4,
               max_treedepth = 15,
               adapt_delta = 0.8,
               init = init_def)


cmdstanfit$save_output_files(dir = output_dir,
                             basename = sp_file_name,
                             timestamp = FALSE,
                             random = FALSE)

csvfl = paste0(output_dir,"/",sp_file_name,"-",1:4,".csv")
#cmdstanfit$save_object(file = paste0(output_dir,"/",sp_file_name,".RDS"))

shiny_explore <- FALSE
if(shiny_explore){
   #load(paste0(output_dir,"/",sp_file_name,"_fit_add.RData"))
  sl_rstan <- rstan::read_stan_csv(csvfl)
  launch_shinystan(as.shinystan(sl_rstan))
}




save(list = c("cmdstanfit",
              "stan_data",
              "dts",
              "real_grid",
              "strats_dts",
              "strat_regions",
              "parms",
              "sp",
              "spf",
              "sp_file_name",
              "output_dir",
              "prior",
              "mod.file1",
              "noise_dist1",
              "csvfl"),
     file = paste0(output_dir,"/",sp_file_name,"_fit_add.RData"))


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

 

