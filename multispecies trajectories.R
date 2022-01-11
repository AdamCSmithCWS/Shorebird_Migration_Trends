## multi species trajectory plots

library(tidyverse)
library(posterior)
library(cmdstanr)
source("functions/utility_functions.R")
source("functions/posterior_summary_functions.R")
source("Functions/palettes.R")

load("data/allShorebirdPrismFallCounts.RData")
source("Functions/palettes.R")

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



# Species loop ------------------------------------------------------------
output_dir <- "g:/Shorebird_Migration_Trends/output"



for(sp in sps[-1]){
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
    
    
    
    cmdstanfit <- as_cmdstan_fit(csvfl) 
    three_gen <- max(10,ceiling(gens[which(gens$Common_name == sp),"GenLength"]*3))
    
    
    syear = min(dts$YearCollected)
    
    stg <- c(rep("recent",three_gen+1),
                 rep("previous",three_gen))
    nadd = length(syear:2019)-length(stg) 
    if(nadd > 0){
      stg <- c(stg,rep("earlier",nadd))
    }
    if(nadd < 0){
      stg <- stg[1:length(syear:2019)]
    }
    
    
    stag <- data.frame(year = syear:2019,
                      stage = rev(stg))
    
    
    
    
    Nsamples <- posterior_samples(fit = cmdstanfit,
                                  parm = "N",
                                  dims = c("y")) %>% 
      mutate(year = y + (syear-1),
             parm = "Full")
    
    
    
   
    NSmoothsamples <- posterior_samples(fit = cmdstanfit,
                                        parm = "NSmooth",
                                        dims = c("y"))%>% 
      mutate(year = y + (syear-1))
    


    indicesN <- posterior_sums(Nsamples,
                               quantiles = NULL,
                               dims = "y") %>% 
      mutate(parm = "Full")
      
    indicesNSmooth <- posterior_sums(NSmoothsamples,
                               quantiles = NULL,
                               dims = "y") %>% 
      mutate(parm = "Smooth")
    
      
    indices = bind_rows(indicesN,indicesNSmooth) %>% 
      mutate(species = sp) %>% 
      left_join(.,stag,by = "year")
    
    
    indices_out <- bind_rows(indices_out,indices)
    
  }
  
}


save(list = "indices_out",
     file = "Trends/all_survey_wide_indices.RData")




# Plotting ----------------------------------------------------------------


load("Trends/all_survey_wide_indices.RData")

inds_sm <- indices_out %>% 
  filter(parm = "Smooth")





