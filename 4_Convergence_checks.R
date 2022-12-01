

### checking convergence on models


# SPECIES MCMC data-prep -------------------------------------------------------
library(tidyverse)
library(cmdstanr)
library(shinystan)



# library(doParallel)
# library(foreach)

#load data
load("data/full_observation_dataset.Rdata")
#load the hexagon map
load( "data/hexagon_grid.RData")
sps <- readRDS("data/species_vector.rds")



FYYYY = 1980


output_dir <- "output"
#output_dir2 <- "g:/Shorebird_Migration_Trends/output"

sum_output <- NULL

for(sp in sps){
  load(paste0("data/data",sp,"_cmdstanr_data.RData"))
  spf = gsub(sp,pattern = " ",replacement = "_")
  spf = gsub(pattern = "\'",replacement = "",
             x = spf)
  sp_file_name <- paste0(spf,"gamma-t")
  

  load(paste0(output_dir,"/",sp_file_name,"_fit_add.RData"))

print(sp)


tmpsum <- (cmdstanfit$summary())
tmpsum <- tmpsum %>% 
  mutate(species = sp)
sum_output <- bind_rows(sum_output,tmpsum)
}


write.csv(sum_output,paste0("convergence/convergence_summarygamma-t.csv"),
          row.names = FALSE)




# sum_output = read.csv(paste0("convergence/convergence_summarygamma-t.csv"))

noise <- sum_output  %>% 
  filter(grepl("noise[",variable,fixed = TRUE))

noise_comp <- NULL
for(sp in sps){
  load(paste0("data/data",sp,"_cmdstanr_data.RData"))
  tmp <- noise %>% filter(species == sp)
  tmp <- bind_cols(tmp,dts)
  noise_comp <- bind_rows(noise_comp,tmp)
}


noise_by <- ggplot(data = noise_comp,aes(x = count+1,y = rhat))+
  geom_point(aes(colour = year))+
  scale_x_continuous(trans = "log10")+
  facet_wrap(~species,scales = "free")

print(noise_by)



sdnoise <- sum_output  %>% 
  filter(grepl("sdnoise",variable,fixed = TRUE)) 

sdnoise_by = ggplot(data = sdnoise,aes(x = mean,y = rhat))+
  geom_point(aes(size = ess_bulk))+
  ggrepel::geom_text_repel(data = sdnoise,aes(label = species))
print(sdnoise_by)



sdalpha <- sum_output  %>% 
  filter(grepl("sdalpha",variable,fixed = TRUE)) 

sdalpha_by = ggplot(data = sdalpha,aes(x = mean,y = rhat))+
  geom_point(aes(size = ess_bulk))+
  ggrepel::geom_text_repel(data = sdalpha,aes(label = species))
print(sdalpha_by)



ALPHA <- sum_output  %>% 
  filter(grepl("ALPHA",variable,fixed = TRUE)) 

ALPHA_by = ggplot(data = ALPHA,aes(x = mean,y = rhat))+
  geom_point(aes(size = ess_bulk))+
  ggrepel::geom_text_repel(data = ALPHA,aes(label = species))
print(ALPHA_by)



alpha <- sum_output  %>% 
  filter(grepl("alpha[",variable,fixed = TRUE),
         !grepl("sdalpha",variable,fixed = TRUE)) 

alpha_by = ggplot(data = alpha,aes(x = mean,y = rhat))+
  geom_point(aes(alpha = ess_bulk))+
  facet_wrap(~species,scales = "free")
print(alpha_by)




infff = sum_output %>% filter(!is.finite(rhat))

allplot = ggplot(data = sum_output,aes(x = rhat,group = species))+
  geom_histogram()+
  labs(title = "all_parameters")+
  facet_wrap(~species,nrow = 5,ncol = 6,scales = "free")
print(allplot)

all_notnoise <- sum_output %>% 
  filter(!grepl("noise",variable,fixed = TRUE))

all_notnoiseplot = ggplot(data = all_notnoise,aes(x = rhat,group = species))+
  geom_histogram()+
  labs(title = "all_parameters_except_noise")+
  facet_wrap(~species,nrow = 5,ncol = 6,scales = "free")
print(all_notnoiseplot)


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





sdyear_gam_strats <- sum_output %>% 
  filter(grepl("sdyear_gam_strat",variable,fixed = TRUE))

sdyear_gam_stratplot = ggplot(data = sdyear_gam_strats,aes(x = rhat,group = species))+
  geom_histogram()+
  labs(title = "sdyear_gam_strat")+
  facet_wrap(~species,nrow = 5,ncol = 6,scales = "free")
print(sdyear_gam_stratplot)


sdyear_gam <- sum_output %>% 
  filter(grepl("sdyear_gam",variable,fixed = TRUE),
         !grepl("sdyear_gam_",variable,fixed = TRUE))

sdyear_gamplot = ggplot(data = sdyear_gam,aes(x = species,y = mean))+
  geom_errorbar(aes(ymin = q5,ymax = q95),alpha = 0.2,width = 0)+
  geom_point()+
  coord_flip()+
  labs(title = "sdyear_gam")
print(sdyear_gamplot)



sdalpha <- sum_output %>% 
  filter(grepl("sdalpha",variable,fixed = TRUE))



ALPHA1 <- sum_output %>% 
  filter(grepl("ALPHA1",variable,fixed = TRUE))


alpha <- sum_output %>% 
  filter(grepl("alpha",variable,fixed = TRUE),
         !grepl("alpha_raw",variable,fixed = TRUE))

alphaplot = ggplot(data = alpha,aes(x = rhat,group = species))+
  geom_histogram()+
  labs(title = "alpha")+
  facet_wrap(~species,nrow = 5,ncol = 6,scales = "free")
print(alphaplot)

alpha_fail <- alpha %>% filter(rhat > 1.1)



pdf(file = paste0("figures/Rhat_distributionsgamma-t.pdf"),
    width = 11,
    height = 8.5)
print(allplot)
print(all_notnoiseplot)

print(nsplot)
print(nsmoothplot)
print(Bplot)
print(betaplot)
print(sdyear_gam_stratplot)
print(sdyear_gamplot)
print(noise_by)
print(sdnoise_by)
print(sdalpha_by)
print(ALPHA_by)
print(alpha_by)


dev.off()



