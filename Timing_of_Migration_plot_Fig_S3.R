### script to prepare Figure S3 

library(tidyverse)
library(posterior)
library(cmdstanr)
library(lubridate)
source("functions/utility_functions.R")
source("functions/posterior_summary_functions.R")
source("Functions/palettes.R")



#load observation data
load("data/full_observation_dataset.Rdata")
#load the hexagon map
load( "data/hexagon_grid.RData")
sps <- readRDS("data/species_vector.rds")
sp_groups <- read.csv("data/seasons_by_species.csv")

#lists for stored figures

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

season_out<- NULL
seasonEffect_out <- NULL


# Species loop ------------------------------------------------------------
output_dir <- "output"


for(sp in sps){
  spf = gsub(sp,pattern = " ",replacement = "_")
  spf = gsub(pattern = "\'",replacement = "",
             x = spf)
  
  load(paste0("data/data",sp,"_cmdstanr_data.RData"))
  
  sp_file_name <- paste0(spf,"-","gamma-t")

  if(file.exists(paste0(output_dir,"/",sp_file_name,"_fit_add.RData"))){
    load(paste0(output_dir,"/",sp_file_name,"_fit_add.RData"))
    if(length(csvfl) > 4){csvfl <- csvfl[1:4]}
    cmdstanfit <- as_cmdstan_fit(csvfl) 
    #<- readRDS(paste0(output_dir,"/",sp_file_name,".RDS"))
    #load(paste0(output_dir,"/",spf,"_fit_add.RData"))
    three_gen <- max(10,ceiling(gens[which(gens$Common_name == sp),"GenLength"]*3))
    #Three generation assessment time in COSEWIC report
    y3g <- 2019-three_gen
    
   
    syear = min(dts$YearCollected)
    
 
    alpha_samples <- posterior_samples(fit = cmdstanfit,
                                       parm = "alpha",
                                       dims = c("s"))

    sdalpha_samples <- posterior_samples(fit = cmdstanfit,
                                         parm = "sdalpha")
    sdalpha = sdalpha_samples %>%
      summarise(mean = mean((.value)),
                lci = quantile((.value),0.025),
                uci = quantile((.value),0.975))


    sites_strat = (stan_data$sites)
    nstrata = stan_data$nstrata
    nsites_strat = stan_data$nsites_strat
    
    sitesbystrat = NULL
    for(st in 1:stan_data$nstrata){
      tmp = data.frame(strat = st,
                       site = sites_strat[1:nsites_strat[st],st])
      sitesbystrat <- bind_rows(sitesbystrat,tmp)
    }
    sitesbystrat <- arrange(sitesbystrat,site)
    #
    alphas = alpha_samples %>% group_by(s) %>%
      summarise(mean = mean(exp(.value)),
                lci = quantile(exp(.value),0.025),
                uci = quantile(exp(.value),0.975)) %>%
      mutate(site = s) %>% left_join(.,sitesbystrat)
    
    ALm <- posterior_samples(fit = cmdstanfit,
                             parm = "ALPHA1",
                             dims = NULL) %>%
      summarise(mean = mean(.value)) %>%
      as.numeric()


    sdnm <- posterior_samples(fit = cmdstanfit,
                              parm = "sdnoise",
                              dims = NULL) %>%
      summarise(mean = mean(.value)) %>%
      as.numeric()

    scale_adj <- ((sdnm^2)*0.5 + ALm)
    # 
    
    n_countsby_site <- dts %>% group_by(site) %>%
      summarise(n_c = n())
    # 
    if(grepl(x = mod.file1,pattern = "two_season")){
      
      strat_season_strat <- dts %>% distinct(seas_strat,strat,hex_name)
      
      alphas <- left_join(alphas,strat_season_strat,by = "strat")
      
      
      strat_offs <- alphas %>% left_join(.,n_countsby_site,by = "site") %>% 
        group_by(strat,seas_strat) %>% 
        summarise(adjs = mean(mean))
      
      
      
      
      season_samples <- posterior_samples(fit = cmdstanfit,
                                          parm = "season_pred",
                                          dims = c("d","s"))
      
      seasonEffect2 = season_samples %>% group_by(d,s) %>% 
        summarise(mean = mean(exp(.value)),
                  lci = quantile(exp(.value),0.025),
                  uci = quantile(exp(.value),0.975)) %>% 
        mutate(day = d,
               seas_strat = s,
               species = sp) 
      
      seasonEffect_out <- bind_rows(seasonEffect_out,seasonEffect2)
      
      seasonEffectT = season_samples %>% group_by(d,s) %>% 
        summarise(mean = mean(exp(.value+scale_adj)),
                  lci = quantile(exp(.value+scale_adj),0.025),
                  uci = quantile(exp(.value+scale_adj),0.975)) %>% 
        mutate(day = d,
               seas_strat = s) 
      
      obs_season <- dts %>% group_by(date,seas_strat) %>% 
        summarise(mean = mean(count),
                  median = median(count),
                  lqrt = quantile(count,0.05),
                  uqrt = quantile(count,0.95)) %>% 
        mutate(day = date)
      
      sse <- seasonEffectT %>% group_by(seas_strat) %>% 
        group_split()
      
      #seasonEffect <- expand_grid(seasonEffect,strat = c(1:nstrata))
      tout <- NULL
      for(j in 1:nstrata){
        wg = as.integer(strat_offs[which(strat_offs$strat == j),"seas_strat"])
        tmp <- sse[[wg]]
        tmp$strat <- j
        tout <- bind_rows(tout,tmp)
      }
      
      seasonEffect <- left_join(tout,strat_offs,by = "strat") %>% 
        mutate(mean = mean*adjs,
               lci = lci*adjs,
               uci = uci*adjs)
      
      
      seas_strat_offs <- alphas %>% left_join(.,n_countsby_site,by = "site") %>% 
        group_by(seas_strat) %>% 
        summarise(adjs = mean(mean),
                  adjs2 = median(mean))
      
      seasonEffect_plot <- left_join(seasonEffectT,seas_strat_offs,by = "seas_strat") %>% 
        mutate(mean = mean*adjs + 1,
               lci = lci*adjs + 1,
               uci = uci*adjs + 1,
               species = sp)
      
      season_out <- bind_rows(season_out,seasonEffect_plot)
      
      
      
    }else{
      strat_offs <- alphas %>% group_by(strat) %>% 
        summarise(adjs = mean(mean))
      
      
      season_samples <- posterior_samples(fit = cmdstanfit,
                                          parm = "season_pred",
                                          dims = c("d"))
      
      
      #season_samples <- slope_icar_stanfit %>% gather_draws(season_pred[d])
      seasonEffectT = season_samples %>% group_by(d) %>% 
        summarise(mean = mean(exp(.value+scale_adj)),
                  lci = quantile(exp(.value+scale_adj),0.025),
                  uci = quantile(exp(.value+scale_adj),0.975)) %>% 
        mutate(day = d) 
      
      
      seasonEffect2 = season_samples %>% group_by(d) %>% 
        summarise(mean = mean(exp(.value)),
                  lci = quantile(exp(.value),0.025),
                  uci = quantile(exp(.value),0.975)) %>% 
        mutate(day = d,
               species = sp) 
      
      seasonEffect_out <- bind_rows(seasonEffect_out,seasonEffect2)
      
      
      obs_season <- dts %>% group_by(date) %>% 
        summarise(mean = mean(count),
                  median = median(count),
                  lqrt = quantile(count,0.05),
                  uqrt = quantile(count,0.95)) %>% 
        mutate(day = date)
      
      seasonEffect <- expand_grid(seasonEffectT,strat = c(1:nstrata))
      seasonEffect <- left_join(seasonEffect,strat_offs) %>% 
        mutate(mean = mean*adjs,
               lci = lci*adjs,
               uci = uci*adjs)
      
      
      
      seasonEffect_plot <- seasonEffectT %>% 
        mutate(mean = mean*mean(alphas$mean)+1,
               lci = lci*mean(alphas$mean)+1,
               uci = uci*mean(alphas$mean)+1,
               species = sp)
      
      season_out <- bind_rows(season_out,seasonEffect_plot)
      
    }
    
save(list = c("seasonEffect_out",
              "season_out"),
     file = "Data_local/Season_effects.RData")
    
print(sp)
  }
}





# Plotting ----------------------------------------------------------------



load("Data_local/Season_effects.RData")

d1 <- yday(as.Date.character("2019-6-30","%Y-%m-%d"))
od <- function(x){
  y <- as.Date(paste0("2019-",x + d1),"%Y-%j")
}
season_out <- season_out %>% 
  mutate(n_season = factor(ifelse(is.na(seas_strat),"Single","Split"),
                           levels = c("Single","Split"),
                           ordered = TRUE),
         Season_Region = factor(ifelse(is.na(seas_strat)|seas_strat == 1,
                                       "Most or all of survey area",
                                       "Southern portion")),
         date = od(day)) %>% 
  arrange(n_season)
spl <- unique(season_out$species)
season_out <- season_out %>% 
  mutate(species = factor(species,levels = spl,
                          ordered = TRUE))

season_spag2 <- ggplot(data = season_out,aes(x = date,y = mean,
                                             group = Season_Region,
                                             colour = Season_Region))+
  geom_ribbon(aes(x = date,y = mean,
                  ymin = lci,ymax = uci,fill = Season_Region),
              inherit.aes = FALSE,alpha = 0.2)+
  geom_line(alpha = 0.5)+
  my_col2_traj+
  guides(fill = guide_legend(title = "Region of survey area"),
         colour = guide_legend(title = "Region of survey area"))+
  scale_y_continuous(trans = "log10")+
  scale_x_date(date_breaks = "1 month",
               #date_minor_breaks = "1 week",
               date_labels =  "%b %d")+
  theme_bw()+
  ylab("Mean count across all sites")+
  #xlab("day")+
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 7.5),
        axis.text.x=element_text(angle=45, hjust=1),
        legend.position = "top")+
  facet_wrap(vars(species),
             nrow = 5,
             ncol = 6,
             scales = "free_y")

print(season_spag2)


pdf("Figures/Seasonal_means.pdf",
    width = 10,
    height = 7)
print(season_spag2)
dev.off()












