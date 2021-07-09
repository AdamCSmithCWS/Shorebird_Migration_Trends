library(tidyverse)
library(sf)
library(spdep)
library(ggforce)
library(tidybayes)
library(GGally)
library(posterior)
library(cmdstanr)
source("functions/utility_functions.R")
source("functions/posterior_summary_functions.R")

library(loo)
load("data/allShorebirdPrismFallCounts.RData")
source("Functions/palettes.R")


#lists for stored figures

blank_list <- vector(mode = "list",length = length(sps))
names(blank_list) <- sps


indices_out <- NULL


# Loading gen times from Bird et al 2020 --------
gens = read.csv("data/cobi13486-sup-0004-tables4.csv")
fullgensnames = read.csv("data/cobi13486-sup-0001-tables1.csv")

fullgensnames <- fullgensnames %>% select(Scientific_name,Common_name)
gens <- gens %>% select(Scientific_name,
                        GenLength) %>% 
  left_join(.,fullgensnames)

gens[which(gens$Common_name == "American Golden Plover"),"Common_name"] <- "American Golden-Plover"
gens[which(gens$Common_name == "Grey Plover"),"Common_name"] <- "Black-bellied Plover"

sps[-which(sps %in% gens$Common_name)]
 
  gens <- gens %>% filter(Common_name %in% sps)



grid_spacing <- 300000  # size of squares, in units of the CRS (i.e. meters for lae)


FYYYY = 1980

t1 = Sys.time()


w_cosewic = sps[c(2:4,7,10,12:20,22,11,25)]

# Species loop ------------------------------------------------------------
output_dir <- "g:/Shorebird_Migration_Trends/output"

for(sp in sps[c(21,22)]){
  #if(sp == "Semipalmated Sandpiper"){next}
  spf = gsub(sp,pattern = " ",replacement = "_")
  spf = gsub(pattern = "\'",replacement = "",
             x = spf)
  
  load(paste0("data/data",sp,"_cmdstanr_data.RData"))
  
  noise_dist_sel <- noise_dist2
  
  sp_file_name <- paste0(spf,"-",prior,"-",noise_dist_sel)
  # paste0(output_dir,"/",sp_file_name,".RDS")
  # 
  #paste0(output_dir,"/",sp_file_name,"_fit_add.RData")

  if(file.exists(paste0(output_dir,"/",sp_file_name,"_fit_add.RData"))){
    load(paste0(output_dir,"/",sp_file_name,"_fit_add.RData"))
    if(length(csvfl) > 4){csvfl <- csvfl[1:4]}
    cmdstanfit <- rstan::read_stan_csv(csvfl) 
    #<- readRDS(paste0(output_dir,"/",sp_file_name,".RDS"))
    #load(paste0(output_dir,"/",spf,"_fit_add.RData"))
    three_gen <- max(10,ceiling(gens[which(gens$Common_name == sp),"GenLength"]*3))
    #Three generation assessment time in COSEWIC report
    y3g <- 2019-three_gen
      
      # CONSIDERATIONS ----------------------------------------------------------
    
    ## map the looic values to see if some areas are being poorly predicted
    
    ## add some variation to the plotting of the observed means (box and whisker info)
    
    #################################
    ## plot the changing mean alpha through time - are smaller sites being surveyed
    ## more in recent years?
    #################################
    
    
    ## generate some fake data to ensure there isn't a negative bias
    
    
     ## export a table with an n-surveys by site by year matrix for the data holders to confirm
    ## send the model summary to Paul - seasonal-split, maps of regions, etc.
    ## rationalize the 1980 start date (see above table), plot number of surveys and sites by time
    
    # Calculate Annual indices using samples ----------------------------------
    
    syear = min(dts$YearCollected)
    
    # gather_draws2 <- funtion(model,
    #                          vars = "",
    #                          dims = 1){
    #   
    #   cmt = paste0("(",paste(vars,sep = "|"),")")
    #   dmt = paste(rep("[:digit:]",times = dims),sep = ",")
    #   nmpt = paste0("(",cmt,"\\[",dmt,"\\]",")")
    #   plong <- draws %>% pivot_longer(
    #     cols = matches(cmt),
    #     names_pattern = regex(nmpt),
    #     names_to = c("variable","group"),
    #     values_to = ".value"
    #   )
    # }
    # 
    


    Nsamples <- posterior_samples(fit = cmdstanfit,
                 parm = "N",
                 dims = c("y"))
      
      
   
    Nsamples$year <- Nsamples$y + (syear-1)
    
    NSmoothsamples <- posterior_samples(fit = cmdstanfit,
                                  parm = "NSmooth",
                                  dims = c("y"))
    NSmoothsamples$year <- NSmoothsamples$y + (syear-1)
    
   
    
   
    
    indicesN <- index_summary(samples = Nsamples,
                              parm = "N",
                              dims = "year")
    
    
    indicesNSmooth <- index_summary(samples = NSmoothsamples,
                                    parm = "NSmooth",
                                    dims = "year")
    
    
    
    
    
    
    indices = bind_rows(indicesN,indicesNSmooth)
    indices$species <- sp
    
    indices_out <- bind_rows(indices_out,indices)
    
  }
  print(sp)
}#end of species loop


write.csv(indices_out,paste0("trends/Indices_continental_gamma_",noise_dist_sel,".csv"))



