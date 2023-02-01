## multi species trajectory plots

library(tidyverse)
library(posterior)
library(cmdstanr)
source("functions/utility_functions.R")
source("functions/posterior_summary_functions.R")
source("Functions/palettes.R")

source("Functions/palettes.R")

#load observation data
load("data/full_observation_dataset.Rdata")
#load the hexagon map
load( "data/hexagon_grid.RData")
sps <- readRDS("data/species_vector.rds")
sp_groups <- read.csv("data/seasons_by_species.csv")

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
output_dir <- "output"



for(sp in sps){
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
      mutate(year = y + (syear-1))
    
    
    
   
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
indices_out2 <- indices_out

save(list = "indices_out2",
     file = "Trends/all_survey_wide_indices.RData")






load("Trends/all_survey_wide_indices.RData")



# exporting indices for distribution --------------------------------------

inds_distribution <- indices_out2 %>% 
  filter(parm == "Full") %>% 
  select(species,
         year,
         median,
         lci,
         uci) %>% 
  rename(predicted_mean_abundance = median,
         lower_95percent_CL = lci,
         upper_95percent_CL = uci)

write.csv(inds_distribution,
          "trends/Survey_wide_annual_indices_shorebird.csv",
          row.names = FALSE)


inds_distribution_all <- indices_out2 %>% 
  select(species,
         parm,
         year,
         median,
         lci,
         uci) %>% 
  rename(smooth_or_full_indices = parm,
         predicted_mean_abundance = median,
         lower_95percent_CL = lci,
         upper_95percent_CL = uci)

write.csv(inds_distribution_all,
          "trends/Survey_wide_smooth_and_full_annual_indices_shorebird.csv",
          row.names = FALSE)


# Export for State of Canada's Birds publication --------------------------
source("functions/loess_func.r")
inds_socb <- inds_distribution_all %>% 
  filter(smooth_or_full_indices == "Full") %>% 
  select(-smooth_or_full_indices) %>% 
  rename(Index = predicted_mean_abundance,
         Year = year,
         lower_ci = lower_95percent_CL,
         upper_ci = upper_95percent_CL) %>% 
  group_by(species) %>% 
mutate(LOESS_index = loess_func(Index,Year))

add_sm <- inds_distribution_all %>% 
  filter(smooth_or_full_indices == "Smooth") %>% 
  select(species,year,predicted_mean_abundance) %>% 
  rename(Alternate_LOESS_index_oringinal_GAM = predicted_mean_abundance,
         Year = year)

inds_socb <- inner_join(inds_socb,add_sm,
                        by = c("species","Year"))


write.csv(inds_socb, "Trends/Shorebird_migration_2019_annual_indices_for_socb.csv")




# Plotting ----------------------------------------------------------------

swp <- function(x){
  nx <- gsub(pattern = "change",
             replacement = "recent",
             x)
  nx2 <- gsub(pattern = "full",
             replacement = "remainder",
             nx)
  return(nx2)
}





inds_sm1 <- indices_out2 %>% 
  filter(parm == "Smooth") %>% 
  mutate(stage = swp(stage))

inds_sm <- inds_sm1 %>% 
  filter(stage != "remainder")

inds_smr <- inds_sm1 %>% 
  filter(stage == "remainder")


inds_f <- indices_out2 %>% 
  filter(parm == "Full")


cuts <- indices_out2 %>% 
  mutate(stage = swp(stage)) %>% 
  group_by(species,stage) %>% 
  summarise(fyr = min(year)) %>% 
  filter(stage != "full")




np <- ggplot(data = inds_f,
             aes(x = year,y = median))+
  geom_ribbon(aes(ymin = lci,ymax = uci),
              alpha = 0.1)+
  geom_line()+
  geom_line(data = inds_sm, aes(x = year,y = median,
                                colour = stage))+
  geom_line(data = inds_smr, aes(x = year,y = median),
            alpha = 0.6)+
  geom_vline(data = cuts,aes(xintercept = fyr),
             colour = grey(0.3))+
  scale_y_continuous(labels = scales::comma,
                     trans = "log10")+
  xlab("")+
  ylab("Modeled annual abundance (mean count/survey)")+
  my_col2_traj+
  theme_bw()+
  labs(caption = str_wrap("Figure S5. Smoothed and unsmoothed population trajectories (displayed as esimated mean count per survey on a log scale) for 28 species of North American shorebirds 1980-2019. The most recent three-generation trends, the previous three-generation trends, and any prior years are separated by vertical bars, and displayed in different colours.",110))+
  theme(text = element_text(family = "serif"),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        plot.caption = element_text(size = 10,hjust = 0),
        plot.margin = margin(0.75,0.75,0.75,0.75,"in"))+
  facet_wrap(~species,nrow = 7,ncol = 4,
             scales = "free_y")

pdf(file = "Figures/Figure_S5_temp.pdf",
    width = 8.5,
    height = 11)
print(np)
dev.off()



