

### checking convergence on models


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

# for(sp in w_cosewic){
# 
#   if(file.exists(paste0("output/",sp,"_GAMYE_strat_simple",grid_spacing/1000,".RData"))){w_cosewic <- w_cosewic[-which(w_cosewic == sp)]}
# 
# }
output_dir <- "g:/Shorebird_Migration_Trends/output"
#output_dir2 <- "g:/Shorebird_Migration_Trends/output"

sum_output <- NULL

for(sp in sps){
  load(paste0("data/data",sp,"_cmdstanr_data.RData"))
  spf = gsub(sp,pattern = " ",replacement = "_")
  spf = gsub(pattern = "\'",replacement = "",
             x = spf)
  sp_file_name <- paste0(spf,"-",prior,"-",noise_dist2)
  

  load(paste0(output_dir,"/",sp_file_name,"_fit_add.RData"))

print(sp)


tmpsum <- (cmdstanfit$summary())
tmpsum <- tmpsum %>% 
  mutate(species = sp)
sum_output <- bind_rows(sum_output,tmpsum)
}


write.csv(sum_output,paste0("convergence/convergence_summary",prior,"-",noise_dist2,".csv"),
          row.names = FALSE)



betas <- sum_output %>% 
  filter(grepl("b[",variable,fixed = TRUE))

betaplot = ggplot(data = betas,aes(x = rhat,group = species))+
  geom_histogram()+
  labs(title = "beta[]")+
  facet_wrap(~species,nrow = 5,ncol = 6,scales = "free")
print(betaplot)





Bs <- sum_output %>% 
  filter(grepl("B[",variable,fixed = TRUE))

Bplot = ggplot(data = Bs,aes(x = rhat,group = species))+
  geom_histogram()+
  labs(title = "B[]")+
  facet_wrap(~species,nrow = 5,ncol = 6,scales = "free")
print(Bplot)




nsmooths <- sum_output %>% 
  filter(grepl("nsmooth[",variable,fixed = TRUE))

nsmoothplot = ggplot(data = nsmooths,aes(x = rhat,group = species))+
  geom_histogram()+
  labs(title = "nsmooth[]")+
  facet_wrap(~species,nrow = 5,ncol = 6,scales = "free")
print(nsmoothplot)




ns <- sum_output %>% 
  filter(grepl("n[",variable,fixed = TRUE))

nsplot = ggplot(data = ns,aes(x = rhat,group = species))+
  geom_histogram()+
  labs(title = "n[]")+
  facet_wrap(~species,nrow = 5,ncol = 6,scales = "free")
print(nsplot)


pdf(file = paste0("figures/Rhat_distributions",prior,"-",noise_dist2,".pdf"),
    width = 11,
    height = 8.5)
print(nsplot)
print(nsmoothplot)
print(Bplot)
print(betaplot)

dev.off()



