### estimating the mean relative trend among spatial units (strata and regions)

### also produces supplementary figures of maps of the strata neighbourhood relationships for each species
### as well as the supplementary figures showing site locations

# PACKAGES ----------------------------------------------------------------


library(tidyverse)
library(tidybayes)
library(rstan)
rstan_options(auto_write = TRUE, javascript = FALSE)
library(shinystan)
library(sf)
library(spdep)
library(ggforce)
library(GGally)
library(patchwork)


source("Functions/Utility_functions.R")
# loading the base hexagon map
load("data/hexagon_grid.RData")

# DATA LOAD --------------------------------------------------

# extract the centered trends
# run the comparison model on the centered trends
# map the mean trend in space
# plot the mean trends by region

trends_all <- read.csv("trends/All_strata_gamma_t_level_trends.csv")

trends_all$centered_log_trend_sd = ((trends_all$centered_log_trend_uci-trends_all$centered_log_trend_lci)/(1.96*2))

t_types = unique(trends_all$trend_type)

# Simple strata-only model ------------------------------------------------

for(time in t_types){
## in this model each species is treated as an independent measure within each
## stratum - 
## this is reasonable given that each species estimates are centered on the species
## mean trend

  trends <- filter(trends_all,trend_type == time)
trends$strataF <- factor(trends$region)
trends$strat <- as.integer(trends$strataF)
nstrata = max(trends$strat)
strata = trends$strat
t_hat = trends$centered_log_trend
sd_hat = trends$centered_log_trend_sd
N = nrow(trends)

n_species <- length(unique(trends$species))
# generate stan data ------------------------------------------------------

stan_data <- list(nstrata = nstrata,
                 N = N,
                 t_hat = t_hat,
                 sd_hat = sd_hat,
                 strata = strata)

mod.file = "models/trend_comparison.stan"



parms = c("mu_t",
          "mu_strata",
          "epsilon",
          "sd_strata_raw",
          "epsilon")


## compile model
model = stan_model(file=mod.file)

## run sampler on model, data
stanfit <- sampling(model,
                    data=stan_data,
                    verbose=TRUE, refresh=200,
                    chains=4, iter=2500,
                    warmup=1000,
                    cores = 4,
                    pars = parms,
                    control = list(adapt_delta = 0.8,
                                   max_treedepth = 15))






save(list = c("stanfit",
              "stan_data",
              "trends",
              "n_species"),
     file = paste0("output/Strata_comparison",time,"_all_species.RData"))

#launch_shinystan(stanfit) 

}



tp_out <- vector(mode = "list",length = length(t_types))
names(tp_out) <- t_types
mp_out <- tp_out


for(time in t_types){
  

load(paste0("output/Strata_comparison",time,"_all_species.RData"))
mu_strata <- gather_draws(stanfit,mu_strata[strat]) %>% 
  mutate(pos = ifelse(.value > 0,TRUE,FALSE)) %>% 
  group_by(strat) %>% 
  summarise(mean = mean(.value),
            lci = quantile(.value,0.05),
            uci = quantile(.value,0.95),
            lqrt = quantile(.value,0.25),
            uqrt = quantile(.value,0.75),
            p_pos = sum(pos)/n(),
            p_neg = 1-(sum(pos)/n())) %>% 
  mutate(sig = ifelse(lci > 0 | uci < 0,TRUE,FALSE),
         sigq = ifelse(lqrt > 0 | uqrt < 0,TRUE,FALSE),
         trend = (exp(mean)-1)*100,
         trend_lci = (exp(lci)-1)*100,
         trend_uci = (exp(uci)-1)*100,
         trend_lqrt = (exp(lqrt)-1)*100,
         trend_uqrt = (exp(uqrt)-1)*100,
         p_not_zero = ifelse(p_pos > 0.5,p_pos,p_neg))

strat_sums <- trends %>% 
  select(strat,region,species) %>% 
  distinct() %>% 
  group_by(strat,region) %>% 
  summarise(nspecies = n()) %>% 
  left_join(.,mu_strata,by = "strat")


tp <- ggplot(data = strat_sums)+
  geom_point(aes(x = strat,y = trend,size = p_not_zero,colour = sig))+
  geom_errorbar(aes(x = strat,y = trend,ymin = trend_lci,ymax = trend_uci),
                alpha = 0.2,width = 0)+
geom_errorbar(aes(x = strat,y = trend,ymin = trend_lqrt,ymax = trend_uqrt),
              alpha = 0.2,width = 0,size = 1.3)+
geom_abline(slope = 0,intercept = 0)+
  labs(title = time)

tp_out[[time]] <- tp


# Mapping mean trends -----------------------------------------------------

mp <- trend_map_composite(trends = strat_sums,
                            map.file = "BBS_ProvState_strata",
                            hex_map = poly_grid,
                            #size_value = "Species Count",
                          size_value = "Probability > or < zero",
                          tlab = paste(time,n_species,"species"))


mp_out[[time]] <- mp


}#end time loop

pdf(file = "Figures/Comparison_graphs_all_species.pdf",
    width = 11,
    height = 8.5)
for(jj in 1:length(t_types)){
  print(tp_out[[jj]])
}
dev.off()


pdf(file = "Figures/Comparison_maps_all_species.pdf",
    width = 11,
    height = 8.5)
for(jj in 1:length(t_types)){
  print(mp_out[[jj]])
}
dev.off()




# publication maps --------------------------------------------------------


# library(gridExtra)


  time = "Recent-three-generation"

load(paste0("output/Strata_comparison",time,"_all_species.RData"))

mu_strata <- gather_draws(stanfit,mu_strata[strat]) %>% 
  mutate(pos = ifelse(.value > 0,TRUE,FALSE)) %>% 
  group_by(strat) %>% 
  summarise(mean = mean(.value),
            lci = quantile(.value,0.05),
            uci = quantile(.value,0.95),
            lqrt = quantile(.value,0.25),
            uqrt = quantile(.value,0.75),
            p_pos = sum(pos)/n(),
            p_neg = 1-(sum(pos)/n())) %>% 
  mutate(sig = ifelse(lci > 0 | uci < 0,TRUE,FALSE),
         sigq = ifelse(lqrt > 0 | uqrt < 0,TRUE,FALSE),
         trend = (exp(mean)-1)*100,
         trend_lci = (exp(lci)-1)*100,
         trend_uci = (exp(uci)-1)*100,
         trend_lqrt = (exp(lqrt)-1)*100,
         trend_uqrt = (exp(uqrt)-1)*100,
         p_not_zero = ifelse(p_pos > 0.5,p_pos,p_neg))

strat_sums <- trends %>% 
  select(strat,region,species) %>% 
  distinct() %>% 
  group_by(strat,region) %>% 
  summarise(nspecies = n()) %>% 
  left_join(.,mu_strata,by = "strat")


# Mapping mean trends -----------------------------------------------------
  


  mp1 <- trend_map_composite_simple(trends = strat_sums,
                            map.file = "BBS_ProvState_strata",
                            hex_map = poly_grid,
                            #size_value = "Species Count",
                            size_value = "Probability > or < zero",
                            tlab = "A",
                            add_legend = TRUE)
  
  
  # print(mp1)
  # leg <- get_legend(mp1)
  # 
  # mp1 <- mp1 +
  #   theme(legend.position = "none")
  # 
  
  
  time = "Previous-three-generation"
  
  load(paste0("output/Strata_comparison",time,"_all_species.RData"))
  
  mu_strata <- gather_draws(stanfit,mu_strata[strat]) %>% 
    mutate(pos = ifelse(.value > 0,TRUE,FALSE)) %>% 
    group_by(strat) %>% 
    summarise(mean = mean(.value),
              lci = quantile(.value,0.05),
              uci = quantile(.value,0.95),
              lqrt = quantile(.value,0.25),
              uqrt = quantile(.value,0.75),
              p_pos = sum(pos)/n(),
              p_neg = 1-(sum(pos)/n())) %>% 
    mutate(sig = ifelse(lci > 0 | uci < 0,TRUE,FALSE),
           sigq = ifelse(lqrt > 0 | uqrt < 0,TRUE,FALSE),
           trend = (exp(mean)-1)*100,
           trend_lci = (exp(lci)-1)*100,
           trend_uci = (exp(uci)-1)*100,
           trend_lqrt = (exp(lqrt)-1)*100,
           trend_uqrt = (exp(uqrt)-1)*100,
           p_not_zero = ifelse(p_pos > 0.5,p_pos,p_neg))
  
  strat_sums <- trends %>% 
    select(strat,region,species) %>% 
    distinct() %>% 
    group_by(strat,region) %>% 
    summarise(nspecies = n()) %>% 
    left_join(.,mu_strata,by = "strat")
  
  
  # Mapping mean trends -----------------------------------------------------
  
  

  mp2 <- trend_map_composite_simple(trends = strat_sums,
                                    map.file = "BBS_ProvState_strata",
                                    hex_map = poly_grid,
                                    #size_value = "Species Count",
                                    size_value = "Probability > or < zero",
                                    tlab = "B",
                                    add_legend = FALSE)
  
  
  
  # layt = "
  # AAAC
  # AAAC
  # BBBC
  # BBBC
  # "
 mpboth <- mp1 / mp2 + 
   plot_layout(guides = "collect")
  
  print(mpboth)
  
  
  pdf(file = "Figures/trend_maps_combined.pdf",
      width = 6,
      height = 6)
  print(mpboth)
  dev.off()
  

  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
# demo hexagon map with sites, hexagons, and political jurisdictions --------

species <- "Least Sandpiper"

map.file = "BBS_ProvState_strata"
hex_map = poly_grid
tlab = time

laea = st_crs("+proj=laea +lat_0=40 +lon_0=-95") # Lambert equal area coord reference system

locat = system.file("maps",
                    package = "bbsBayes")

strata_map = read_sf(dsn = locat,
                     layer = map.file)
strata_map = st_transform(strata_map,crs = laea)

load(paste0("data_local/arch/data",species,"_cmdstanr_data0.5_10_2.RData"))


centres = suppressWarnings(st_centroid(real_grid_regs))

coords = st_coordinates(centres)

nb_l <- spdep::nb2listw(nb_db)
nt = length(attributes(nb_l$neighbours)$region.id)
DA = data.frame(
  from = rep(1:nt,sapply(nb_l$neighbours,length)),
  to = unlist(nb_l$neighbours)
)
DA = cbind(DA,coords[DA$from,c("X","Y")],coords[DA$to,c("X","Y")])
colnames(DA)[3:6] = c("long","lat","long_to","lat_to")


box <- st_as_sfc(st_bbox(vintj))
xb = range(st_coordinates(box)[,"X"])
yb = range(st_coordinates(box)[,"Y"])



ggp <- ggplot(data = centres) + 
  geom_sf(data = strata_map,#alpha = 0,
          fill = grey(0.95),
          colour = "white",
          inherit.aes = FALSE)+ 
#  geom_sf(data = vintj,alpha = 0,colour = grey(0.9),inherit.aes = FALSE)+ 
  geom_sf(data = real_grid_regs,alpha = 0,colour = grey(0.85),inherit.aes = FALSE)+
  geom_segment(data=DA,aes(x = long, y = lat,xend=long_to,yend=lat_to),inherit.aes = FALSE,
               colour = "black",linewidth=0.3,alpha=0.1) +
  xlab("")+
  ylab("")+
  geom_sf(size = 0.9)+
  theme_bw() +
   theme(rect = element_blank(),
         panel.grid.major = element_line(color = "white"),
         axis.text = element_text(size = rel(0.8)))+
  coord_sf(xlim = xb,ylim = yb)+
  theme(legend.position = "none")

print(ggp)

pdf("Figures/Demo_strata_figure.pdf",
    width = 3.5,
    height = 2.5)
print(ggp)
dev.off()





# pdf("Figures/Demo_strata_figure_all_species.pdf",
#     width = 10,
#     height = 8)

maps_out <- vector("list",length = length(sps))
names(maps_out) <- sps
for(species in sps){
  map.file = "BBS_ProvState_strata"
  hex_map = poly_grid
  tlab = time
  
  laea = st_crs("+proj=laea +lat_0=40 +lon_0=-95") # Lambert equal area coord reference system
  
  locat = system.file("maps",
                      package = "bbsBayes")
  
  strata_map = read_sf(dsn = locat,
                       layer = map.file)
  strata_map = st_transform(strata_map,crs = laea)
 #
  load(paste0("data_local/arch/data",species,"_cmdstanr_data0.5_10_2.RData"))
  
  
 
  centres = suppressWarnings(st_centroid(real_grid_regs))
  
  coords = st_coordinates(centres)
  
  nb_l <- spdep::nb2listw(nb_db)
  nt = length(attributes(nb_l$neighbours)$region.id)
  DA = data.frame(
    from = rep(1:nt,sapply(nb_l$neighbours,length)),
    to = unlist(nb_l$neighbours)
  )
  DA = cbind(DA,coords[DA$from,c("X","Y")],coords[DA$to,c("X","Y")])
  colnames(DA)[3:6] = c("long","lat","long_to","lat_to")
  
  
  box <- st_as_sfc(st_bbox(vintj))
  xb = range(st_coordinates(box)[,"X"])
  yb = range(st_coordinates(box)[,"Y"])
  
  
  
  ggp <- ggplot(data = centres) + 
    geom_sf(data = strata_map,#alpha = 0,
            fill = grey(0.95),
            colour = "white",
            inherit.aes = FALSE)+ 
    #geom_sf(data = vintj,alpha = 0,colour = grey(0.9),inherit.aes = FALSE)+ 
    geom_sf(data = real_grid_regs,alpha = 0,colour = grey(0.85),inherit.aes = FALSE)+
    geom_segment(data=DA,aes(x = long, y = lat,xend=long_to,yend=lat_to),inherit.aes = FALSE,
                 colour = "black",linewidth=0.3,alpha=0.1) +
    xlab("")+
    ylab("")+
    geom_sf(size = 0.9)+
    theme_bw() +
    labs(title = paste(species))+
    theme(rect = element_blank(),
          panel.grid.major = element_line(color = "white"),
          axis.text = element_text(size = rel(0.8)))+
    coord_sf(xlim = xb,ylim = yb)+
    theme(legend.position = "none")
  
  maps_out[[species]] <- ggp
  
  rm(list = c("real_grid_regs","nb_db",
              "DA"))
}


pdf("Figures/Figure_S4_temp.pdf",
    width = 8.5,
    height = 11)
 
  plott <- (maps_out[[1]] +
    maps_out[[2]] +
    maps_out[[3]] +
    maps_out[[4]] +
    maps_out[[5]] +
    maps_out[[6]] +
    maps_out[[7]] +
    maps_out[[8]] +
    maps_out[[9]] +
    maps_out[[10]] +
    maps_out[[11]] +
    maps_out[[12]]+
      maps_out[[15]]+
      maps_out[[14]] ) + 
    plot_layout(ncol = 3, nrow = 5) +
    plot_annotation(caption = str_wrap("Figure S4. Neighbourhood relationships for the 28 shorebird species included in analyses. The hexagons represent the spatial strata, within which individual survey sites were grouped. The points are centroids of the hexagons, and the lines connecting points represent the strata that are considered neighbours in the intrinsic Conditional Autoregressive component of the model",75),
                    theme = theme(plot.caption = element_text(size = 12,hjust = 0),
                                  plot.margin = margin(0.75,0.75,0.75,0.75,"in")))
  
  print(plott)
  
  
  plott <- (maps_out[[15]] +
              maps_out[[16]] +
              maps_out[[17]] +
              maps_out[[18]] +
              maps_out[[19]] +
              maps_out[[20]] +
              maps_out[[21]] +
              maps_out[[22]] +
              maps_out[[23]] +
              maps_out[[24]] +
              maps_out[[25]] +
              maps_out[[26]]+
              maps_out[[27]]+
              maps_out[[28]] ) + 
    plot_layout(ncol = 3, nrow = 5) +
    plot_annotation(caption = "Figure S4-Cont. Neighbourhood relationships for the 28 shorebird species included in analyses",
                    theme = theme(plot.caption = element_text(size = 12,hjust = 0),
                                  plot.margin = margin(0.75,0.75,0.75,0.75,"in")))
  

  print(plott)
dev.off()



# map of site locations by span -------------------------------------------



map.file = "BBS_ProvState_strata"
hex_map = poly_grid
tlab = time

laea = st_crs("+proj=laea +lat_0=40 +lon_0=-95") # Lambert equal area coord reference system

locat = system.file("maps",
                    package = "bbsBayes")

strata_map = read_sf(dsn = locat,
                     layer = map.file)
strata_map = st_transform(strata_map,crs = laea)

## read in the map of all sites

load("Data/site_map.RData")
load("Data/full_observation_dataset.Rdata")

source("functions/palettes.R")
sample_sum <- ssData %>% 
  filter(YearCollected >= 1980) %>% 
  group_by(SamplingEventIdentifier,
           YearCollected,
           SiteName,
           hex_name,
           SurveyAreaIdentifier) %>% 
  summarise(n = n())

nsurveys_y <- sample_sum %>% 
  group_by(SurveyAreaIdentifier,
           YearCollected,
           hex_name) %>% 
  summarise(nsurveys = n())

nyrs_span <- nsurveys_y %>% 
  group_by(SurveyAreaIdentifier,
           hex_name) %>% 
  summarise(nyears = n(),
            fyr = min(YearCollected),
            lyr = max(YearCollected),
            span = 1+lyr - fyr)

iss_sites_lcc_map <- iss_sites_lcc %>% 
  inner_join(.,nyrs_span,by = "SurveyAreaIdentifier")

iss_sites_lcc_map3 <- iss_sites_lcc_map %>% 
  filter(span <10)
box <- st_as_sfc(st_bbox(iss_sites_lcc_map))
xb = range(st_coordinates(box)[,"X"])
yb = range(st_coordinates(box)[,"Y"])

 
  
sites_map <- ggplot(data = iss_sites_lcc_map3)+
  geom_sf(data = strata_map,#
          alpha = 0,
          #fill = grey(0.95),
          colour = grey(0.85),#= "white",
          inherit.aes = FALSE)+ 
  geom_sf(aes(size = span),
          alpha = 0.3,
          colour = my_simpl_col)+
  #labs(title = "B")+
  scale_radius(
    #max_size = 3,
                  range = c(0,0.5),
                  breaks = c(2,5,8,10),
                  trans = "identity",
                  guide = guide_legend(title = "Years"))+
  # scale_colour_viridis_c(aesthetics = c("colour"),
  #                        begin = 0.25,
  #                        end = 0.9,
  #                        guide = guide_legend(title = "Years"),
  #                        #guide = "none",
  #                        direction = 1)+
  theme_bw() +
  theme(axis.text = element_text(size = rel(0.8)))+
   coord_sf(xlim = xb,ylim = yb)#+
  # theme(rect = element_blank(),
  #       panel.grid.major = element_line(color = "white"),
  #       axis.text = element_text(size = rel(0.8)))



# same but with only sites > 10 year span ---------------------------------



nyrs_span2 <- nsurveys_y %>% 
  group_by(SurveyAreaIdentifier,
           hex_name) %>% 
  summarise(nyears = n(),
            fyr = min(YearCollected),
            lyr = max(YearCollected),
            span = 1+lyr - fyr) %>% 
  filter(span >= 10)

iss_sites_lcc_map2 <- iss_sites_lcc %>% 
  inner_join(.,nyrs_span2,by = "SurveyAreaIdentifier")





sites_map2 <- ggplot(data = iss_sites_lcc_map2)+
  geom_sf(data = strata_map,#
          alpha = 0,
          #fill = grey(0.95),
          colour = grey(0.85),#= "white",
          inherit.aes = FALSE)+ 
  geom_sf(aes(size = span),
          alpha = 0.2,
          colour = my_simpl_col)+
  # geom_sf(data = strata_map,#
  #         alpha = 0,
  #         #fill = grey(0.95),
  #         colour = grey(0.85),#= "white",
  #         inherit.aes = FALSE)+ 
  #labs(title = "A")+
  scale_radius(
    #max_size = 3,
    range = c(0.5,3),
    trans = "identity",
    guide = guide_legend(title = "Years"))+
  # scale_colour_viridis_c(aesthetics = c("colour"),
  #                        begin = 0.25,
  #                        end = 0.9,
  #                        guide = guide_legend(title = "Years"),
  #                        #guide = "none",
  #                        direction = 1)+
  theme_bw()+
  theme(axis.text = element_text(size = rel(0.8)))+
  coord_sf(xlim = xb,ylim = yb)#+
# theme(rect = element_blank(),
#       panel.grid.major = element_line(color = "white"),
#       axis.text = element_text(size = rel(0.8)))


pdf("Figures/Sites_figure_supplement.pdf",
    width = 3.5,
    height = 3)
print(sites_map)
dev.off()


pdf("Figures/Figure1.pdf",
    width = 3.5*1.5,
    height = 3*1.5)
print(sites_map2)
dev.off()

