# SPECIES MCMC data-prep -------------------------------------------------------
library(tidyverse)
library(cmdstanr)
library(shinystan)


sps <- readRDS("data/species_vector.rds")
grid_spacing <- 300000  # size of squares, in units of the CRS (i.e. meters for lae)


FYYYY = 1980

output_dir <- "output"
#output_dir2 <- "g:/Shorebird_Migration_Trends/output"

 for(sp in sps){
   load(paste0("data_local/data",sp,"_cmdstanr_data.RData"))
   spf = gsub(sp,pattern = " ",replacement = "_")
   spf = gsub(pattern = "\'",replacement = "",
              x = spf)
   sp_file_name <- paste0(spf,"-gamma-t")
   
   #if(file.exists(paste0(output_dir,"/",sp_file_name,"-",1,".csv"))){next}
   
   

    ## compile model
    modl = cmdstan_model(stan_file=mod.file2, stanc_options = list("Oexperimental"))

print(sp_file_name)


cmdstanfit<- modl$sample(data=stan_data,
               refresh=200,
               chains=4, iter_sampling =1000,
               iter_warmup=1000,
               parallel_chains = 4,
               max_treedepth = 13,
               adapt_delta = 0.8,
               init = init_def)



#csvfl = paste0(output_dir,"/",sp_file_name,"-",1:4,".csv")
#cmdstanfit$save_object(file = paste0(output_dir,"/",sp_file_name,".RDS"))

shiny_explore <- FALSE
if(shiny_explore){
   #load(paste0(output_dir,"/",sp_file_name,"_fit_add.RData"))
  sl_rstan <- rstan::read_stan_csv(csvfl)
  launch_shinystan(as.shinystan(sl_rstan))
}


sumt = cmdstanfit$summary()

save(list = c("sumt",
              "cmdstanfit",
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
              "mod.file2"),
     file = paste0(output_dir,"/",sp_file_name,"_fit_add.RData"))


}#end modeling loop


# save csv files ----------------------------------------------------------








 

