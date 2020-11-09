# SPECIES MCMC data-prep -------------------------------------------------------
library(tidyverse)
library(rstan)
library(shinystan)
library(sf)
library(spdep)
library(ggforce)
library(tidybayes)


# library(doParallel)
# library(foreach)

#load data
load("data/allShorebirdPrismFallCounts.RData")

#load spatial function to transform GeoBugs format to Stan
source("functions/mungeCARdata4stan.R")
 
# removing Alaska, NWT, and Hawaii ----------------------------------------

 
 ssData <- filter(ssData,StateProvince != "US-HI")
 ssData <- filter(ssData,StateProvince != "US-AK")
 ssData <- filter(ssData,StateProvince != "NT")
 
 # generate equal-area grid stratification and neighbourhoods --------------------------------------
 
 ### site-level model would be useful, but it implies a tremendous number of "trend" parameters that may cause some significant
 ### computational limits
 ### instead, etsablish an equal-area grid as a geographic stratification
 
 ssData$DecimalLatitude = round(ssData$DecimalLatitude,4)
 ssData$DecimalLongitude = round(ssData$DecimalLongitude,4)
 
 all_sites = ssData %>% distinct(SurveyAreaIdentifier,DecimalLatitude,DecimalLongitude)
 # dupl_sites = all_sites[duplicated(all_sites$SurveyAreaIdentifier),"SurveyAreaIdentifier"]
 # dupl_coords = all_sites[which(all_sites$SurveyAreaIdentifier %in% dupl_sites$SurveyAreaIdentifier),]
 # dupl_coords <- dupl_coords[order(dupl_coords$SurveyAreaIdentifier),]
 #write.csv(dupl_coords,"dupl_coords.csv")
 
 laea = st_crs("+proj=laea +lat_0=40 +lon_0=-95") # Lambert equal area
 
 
 
 locat = system.file("maps",
                     package="bbsBayes")
 map.file = "BBS_ProvState_strata"
 
 prov_state = sf::read_sf(dsn = locat,
                              layer = map.file)
 
 # bbs_strata_proj = st_transform(bbs_strata_map,crs = 102003) #USA Contiguous albers equal area projection stanadard from ESRI - https://mgimond.github.io/Spatial/coordinate-systems-in-r.html
 
 prov_state = st_transform(prov_state,crs = laea)
 
 
 iss_sites = st_as_sf(all_sites,coords = c("DecimalLongitude","DecimalLatitude"), crs = 4326)
 iss_sites_lcc <- st_transform(iss_sites, laea)
 bb = st_bbox(iss_sites_lcc) %>% 
   st_as_sfc()

 grid_spacing <- 50000  # size of squares, in units of the CRS (i.e. meters for lae)
 
 poly_grid <- st_make_grid(bb, square = F, cellsize = c(grid_spacing, grid_spacing)) %>% # the grid, covering bounding box
   st_sf() # 
 crd <- st_coordinates(st_centroid(poly_grid))
 poly_grid$hex_name <- paste(round(crd[,"X"]),round(crd[,"Y"]),sep = "_")
 
 iss_sites_lcc <- st_join(iss_sites_lcc, poly_grid, join = st_nearest_feature)
 
 strats <- iss_sites_lcc
 st_geometry(strats) <- NULL
 ssData <- left_join(ssData,strats,by = "SurveyAreaIdentifier")
 ### hex_name is now the new stratification
 
 
 source("functions/GAM_basis_function.R")

# n_cores <- 4
# cluster <- makeCluster(n_cores, type = "PSOCK")
# registerDoParallel(cluster)
# 
# 
# 
# fullrun <- foreach(sp = sps[c(11,25,3,8)],
#                    .packages = c("jagsUI","tidyverse","ggmcmc"),
#                    .inorder = FALSE,
#                    .errorhandling = "pass") %dopar%
#   {
 # sp = sps[25]  
    #load("data/allShorebirdPrismFallCounts.RData")
    
  #  source("functions/GAM_basis_function.R")
    

#for(sp in sps){
sp = sps[11]
FYYYY = 1974
dts <- filter(ssData,CommonName == sp,
              YearCollected >= FYYYY)
nyrs_study <- 2019-FYYYY #length of the time-series being modeled


dts$present <- FALSE

dts[which(dts$ObservationCount > 0),"present"] <- TRUE

# #number of non-zero observations by region
# nobs_reg <- dts %>% 
#   filter(present == TRUE) %>% 
#   group_by(Region) %>% 
#   summarise(nobs = n())

#number of non-zero observations by site
nobs_site <- dts %>% 
  filter(present == TRUE) %>% 
  group_by(Region,SurveyAreaIdentifier) %>% 
  summarise(nobs = n())

yspn = function(x){
  diff(range(x))
}
#number of years with non-zero observations, by site
nyrs_site <- dts %>% 
  filter(present == TRUE) %>%  
  group_by(Region,SurveyAreaIdentifier,YearCollected) %>% 
  summarise(nobs = n()) %>% 
  group_by(Region,SurveyAreaIdentifier) %>% 
  summarise(nyears = n(),
            span_years = yspn(YearCollected),
            fyear = min(YearCollected),
            lyear = max(YearCollected))

#number of sites with > 5 year span by region
nsites_w5 <- nyrs_site %>% 
  filter(span_years > 5) %>%  
  group_by(Region) %>% 
  summarise(nobs = n()) 


#sites with > 5 years of non-zero, observations
sites_keep <- nyrs_site[which(nyrs_site$span_years >= 10),"SurveyAreaIdentifier"]

# drop sites with <5 year span of data ------------------------------------


dts <- filter(dts,SurveyAreaIdentifier %in% sites_keep$SurveyAreaIdentifier) 


# #number of years with non-zero observations, by region
# nyrs_region <- dts %>% 
#   filter(present == TRUE) %>%  
#   group_by(hex_name,YearCollected) %>% 
#   summarise(nobs = n()) %>% 
#   group_by(hex_name) %>% 
#   summarise(nyears = n(),
#             span_years = yspn(YearCollected),
#             fyear = min(YearCollected),
#             lyear = max(YearCollected))
# 
# #strats with 7 or more years of non-zero, observations - species has to be observed in a region in at least 2/3 of the years in the time-series
# regions_keep <- nyrs_region[which(nyrs_region$span_years >= nyrs_study*0.5),"hex_name"]
# 

# drop strata with < 7 years of non-zero observations ---------------------


# dts <- filter(dts,hex_name %in% regions_keep$hex_name) 

real_grid <- iss_sites_lcc %>% filter(SurveyAreaIdentifier %in% sites_keep$SurveyAreaIdentifier) 

strats_dts <- data.frame(SurveyAreaIdentifier = real_grid$SurveyAreaIdentifier,
                         stratn = 1:length(real_grid$SurveyAreaIdentifier))

dts <- left_join(dts,strats_dts,by = "SurveyAreaIdentifier")

# generate neighbourhoods -------------------------------------------------

coords = st_coordinates(st_centroid(real_grid))

nb_db <- spdep::tri2nb(coords,row.names = real_grid$SurveyAreaIdentifier)

nb_info = spdep::nb2WB(nb_db)

### re-arrange GEOBUGS formated nb_info into appropriate format for Stan model
car_stan_dat <- mungeCARdata4stan(adjBUGS = nb_info$adj,
                                  numBUGS = nb_info$num)



# DOY season definition ---------------------------------------------------


fday = min(dts$doy)-1

syear = min(dts$YearCollected)

dts <- dts %>% mutate(count = as.integer(ObservationCount),
                      year = as.integer(YearCollected),
                      yr = as.integer(year-(syear-1)),
                      strat = stratn,
                      date = doy-fday,
                      site = as.integer(factor(SurveyAreaIdentifier))) 

nstrata = max(dts$strat)
nsites = max(dts$site)





ncounts = nrow(dts)

# dmin_pred <- tapply(dts$date,dts$strat,min)
# dmax_pred <- tapply(dts$date,dts$strat,max)
# 
# dmax = max(dts$date)
nyears = max(dts$yr)
midyear = floor(nyears/2)

# GAM seasonal basis function ---------------------------------------------

nKnots_season = 10
basis_season <- gam.basis.func(orig.preds = as.integer(unlist(dts[,"date"])),
                               nknots = nKnots_season,
                               standardize = "z",
                               random = F,
                               npredpoints = max(dts$date),
                               even_gaps = FALSE,
                               sm_name = "season")

ndays <- basis_season$npredpoints_season
# 

# GAM year basis function ---------------------------------------------

nKnots_year = 13
basis_year <- gam.basis.func(orig.preds = as.integer(unlist(dts[,"yr"])),
                               nknots = nKnots_year,
                               standardize = "z",
                               random = F,
                               npredpoints = max(dts$yr),
                               even_gaps = TRUE,
                               sm_name = "year")

# 




# Explore site-level trajectories of observed means ------------------------

# site_means <- dts %>% group_by(year,SurveyAreaIdentifier) %>%
#   summarise(means = log(mean(count,na.rm = T)+1))
# nrow(site_means)/nyears
# smp = ggplot(data = site_means,aes(x = year,y = means,colour = SurveyAreaIdentifier))+
#   geom_point(alpha = 0.05)+
#   geom_smooth(method = "lm",se = FALSE)+
#   #scale_y_continuous(trans = "log10")+
#   theme(legend.position = "none")
# print(smp)#   

# prepare stan data -------------------------------------------------------

stan_data <- list(count = as.integer(unlist(dts$count)),
                  year = as.integer(unlist(dts$yr-midyear)),
                  year_raw = as.integer(unlist(dts$yr)),
                  site = as.integer(unlist(dts$site)),
                  #strat = as.integer(unlist(dts$strat)),
                  date = as.integer(unlist(dts$date)),
                  
                  nyears = nyears,
                  #nstrata = nstrata,
                  nsites = nsites,
                  ncounts = ncounts,
                  ndays = ndays,
                  
                  #nsites_strat = nsites_strat,
                  #max_sites = max_sites,
                  #sites = sites,
                  
                  # season_basis = basis_season$season_basis,
                   season_basispred = basis_season$season_basispred,
                   nknots_season = basis_season$nknots_season,
                  
                  year_basispred = basis_year$year_basispred,
                  nknots_year = basis_year$nknots_year,
                  
                  #midyear = midyear,
                  
                  N_edges = car_stan_dat$N_edges,
                  node1 = car_stan_dat$node1,
                  node2 = car_stan_dat$node2)



parms = c("sdnoise",
          #"nu", #
          "sdalpha",
          "b",
          "B",
          "alpha",
          "ALPHA1",
          "sigma",
          "sdyear",
          "year_effect",
          "sdseason",
          "B_season",
          "season_pred",
          "n",
          "nsmooth",
          "N",
          "NSmooth")

mod.file = "models/slope_icar_GAM_site.stan"

## compile model
slope_icar_model = stan_model(file=mod.file)

## run sampler on model, data
stime = system.time(slope_icar_stanfit <-
                      sampling(slope_icar_model,
                               data=stan_data,
                               verbose=TRUE, refresh=50,
                               chains=5, iter=3000,
                               warmup=2000,
                               cores = 5,
                               pars = parms,
                               control = list(adapt_delta = 0.8,
                                              max_treedepth = 15)))


 launch_shinystan(slope_icar_stanfit) 



save(list = c("slope_icar_stanfit",
              "stan_data"),
     file = paste0(sp,"_Stan_icar_GAM_site_level.RData"))

source("functions/utility_functions.R")


# Extract Annual Indices ------------------------------------------------



indicesN <- index_summary(parm = "N",
                          dims = "year")


indicesNSmooth <- index_summary(parm = "NSmooth",
                                dims = "year")


  
  indices = bind_rows(indicesN,indicesNSmooth)
indices$year = indices$year + (syear-1)


N_gg = ggplot(data = indices,aes(x = year, y = PI50,fill = parm))+
  geom_ribbon(aes(ymin = PI2_5,ymax = PI97_5),alpha = 0.2)+
  geom_line(aes(colour = parm))+
  geom_point(aes(y = obsmean),colour = grey(0.5),alpha = 0.3)

print(N_gg)



indicesnsmooth <- index_summary(parm = "nsmooth",
                                dims = c("site","year"),
                                probs = c(0.05,0.5,0.95))

indicesn <- index_summary(parm = "n",
                          dims = c("site","year"),
                          probs = c(0.05,0.5,0.95))

indices_strat = bind_rows(indicesn,indicesnsmooth)
indices_strat$year = indices_strat$year + (syear-1)
sites_dts <- distinct(dts,site,SurveyAreaIdentifier)
indices_strat <- left_join(indices_strat,sites_dts, by = "site")

pdf(file = paste0("figures/", sp," Site trajectories.pdf"),
    width = 8.5,
    height = 11)
for(jj in 1:ceiling(nsites/12)){
n_gg = ggplot(data = indices_strat,aes(x = year, y = PI50,fill = parm))+
  geom_ribbon(aes(ymin = PI5,ymax = PI95),alpha = 0.2)+
  geom_line(aes(colour = parm))+
  geom_point(aes(y = obsmean),colour = grey(0.5),alpha = 0.1)+
  facet_wrap_paginate(facets = ~SurveyAreaIdentifier,page = jj,nrow = 4, ncol = 3,scales = "free")
print(n_gg)
}
  dev.off()




# Calculate Annual indices using samples ----------------------------------

Nsamples <- slope_icar_stanfit %>% gather_draws(N[y])
Nsamples$year <- Nsamples$y + (syear-1)

NSmoothsamples <- slope_icar_stanfit %>% gather_draws(NSmooth[y])
NSmoothsamples$year <- NSmoothsamples$y + (syear-1)



# calculate trends continent --------------------------------------------------------

t_NSmooth <- ItoT(inds = NSmoothsamples,
            start = 1980,
            end = 2019,
            regions = FALSE,
            qs = 95,
            trend_type = "endpoint",
            index_type = "smooth",
            retransformation_type = "standard")

t_NSmooth_short <- ItoT(inds = NSmoothsamples,
                  start = 2009,
                  end = 2019,
                  regions = FALSE,
                  qs = 95,
                  index_type = "smooth")




# # map the trend estimates -------------------------------------------------
# 
# slopes = as.data.frame(summary(slope_icar_stanfit,
#                  pars = c("b"),
#                  probs = c(0.025,0.5,0.975))$summary)
# 
# myrename = function(fit){
#   rename_with(fit,~ paste0("PI",gsub(".","_",gsub("%", "", .x, fixed = TRUE), fixed = TRUE)),ends_with("%"))
#   
# }
# slopes = myrename(slopes)
# slopes$stratn = 1:nrow(slopes)
# slopes$trend = (exp(slopes$mean)-1)*100
# 
# breaks <- c(-7, -4, -2, -1, -0.5, 0.5, 1, 2, 4, 7)
# labls = c(paste0("< ",breaks[1]),paste0(breaks[-c(length(breaks))],":", breaks[-c(1)]),paste0("> ",breaks[length(breaks)]))
# labls = paste0(labls, " %")
# 
# slopes$Tplot <- cut(slopes$trend,breaks = c(-Inf, breaks, Inf),labels = labls)
# 
# map_palette <- c("#a50026", "#d73027", "#f46d43", "#fdae61", "#fee090", "#ffffbf",
#                  "#e0f3f8", "#abd9e9", "#74add1", "#4575b4", "#313695")
# names(map_palette) <- labls
# 
# hex_map <- left_join(real_grid,strats_dts)
# 
# hex_map <- left_join(hex_map,slopes,by = "stratn")
# 
# trend_map <- ggplot()+
#   geom_sf(data = prov_state,colour = grey(0.8))+
#   geom_sf(data = hex_map,aes(fill = Tplot,colour = Tplot))+
#   coord_sf(ylim = c(-1509123,1890567),expand = TRUE)+
#   scale_colour_manual(values = map_palette, aesthetics = c("fill","colour"),
#                       guide = ggplot2::guide_legend(reverse=TRUE),
#                       name = "Trends")#paste0("Trend\n",fyr,"-",lyr))
# print(trend_map)
# 
# 
# 
# 
# 
# save(list = c("stan_data",
#               #"basis_season",
#               "dts",
#               "slope_icar_model",
#               "hex_map",
#               "strats_dts"),
#      file = paste0("output/",sp,"slope_iCAR_results.RData"))





}#end species loops


stopCluster(cl = cluster)





# Plotting ----------------------------------------------------------------

library(tidybayes)

load("data/allShorebirdPrismFallCounts.RData")
source("functions/Utility_functions.R")

for(sp in sps){
  
  if(file.exists(paste0("output/",sp,"slope_results.RData"))){
  
  load(paste0("output/",sp,"slope_ZIP_results.RData"))
  
    fyear = (min(dts$YearCollected))
    
  strats = unique(dts[,c("strat","Region")])
  strats = rename(strats,s = strat)
  
  
  
  sums = data.frame(out2$summary)
  names(sums) <- c("mean","sd","lci","lqrt","median","uqrt","uci","Rhat","n.eff","overlap0","f")
  sums$Parameter = row.names(sums)
  
  # compiling indices -------------------------------------
  
  n_inds <- extr_inds(param = "n_s")
  N_inds <- extr_inds(param = "N",regions = FALSE)
 
 
  
  sdnoise_st = extr_sum(param = "sdnoise",
                    index = c("s"),
                    log_retrans = F) 
  sdnoise_st <- left_join(sdnoise_st,strats,by = "s")
  
  nu_st = extr_sum(param = "nu",
                        index = c("s"),
                        log_retrans = F) 
  nu_st <- left_join(nu_st,strats,by = "s")
  
  
  
  psi_st = extr_sum(param = "psi",
                       index = c("j","s"),
                       log_retrans = F) 
 
  psiSamples <- out2$samples %>% gather_draws(psi)
  sdnoiseSamples <- out2$samples %>% gather_draws(sdnoise[s])
  
# extracting the seasonal smooth ------------------------------------------

  season_sm = extr_sum(param = "vis.sm_season",
                       index = c("day","s"),
                       log_retrans = TRUE) 
  
  season_sm <- left_join(season_sm,strats)
  pp <- ggplot()+
    geom_line(data = season_sm,aes(x = day,y = mean,colour = Region,group = Region))+
    geom_ribbon(data = season_sm,aes(x = day,ymax = uci,ymin = lci),alpha = 0.2)+
    ylab("")+
    xlab("Days since July 1")+
    facet_wrap(facets = ~Region,ncol = 3,scales = "free")
  
  
  
  pdf(file = paste0("Figures/",sp,"_Season_slope.pdf"),
      width = 8.5,
      height = 11)
  print(pp)
  dev.off()
  
  # calculating trends  -----------------------------------------------------
  
  NSamples <- out2$samples %>% gather_draws(N[y])
  NSamples$year <- NSamples$y + fyear-1
  

  
  N_compSamples <- out2$samples %>% gather_draws(N_comp[y])
  N_compSamples$year <- N_compSamples$y + fyear-1
  
  n_sSamples <- out2$samples %>% gather_draws(n_s[s,y])
  n_sSamples$year <- n_sSamples$y + fyear-1
  n_sSamples <- left_join(n_sSamples,strats,by = "s")
  
 

  t_n_s <- ItoT(inds = n_sSamples,regions = TRUE)

  
  
  
  t_n_sS <- ItoT_slope(inds = n_sSamples,regions = TRUE)

  
  # t_n_s_a1 <- ItoT(inds = n_s_a1Samples,regions = TRUE,retransformation_type = "lognormal_only")
  # 
  # 
  # t_n_s_a2 <- ItoT(inds = n_s_a2Samples,regions = TRUE,retransformation_type = "none")
  # 
  # 
  # 
  t_N <- ItoT(inds = NSamples,regions = FALSE)

  
  t_n_s_15 <- ItoT(inds = n_sSamples,regions = TRUE,start= 2004)

  
  t_n_sS_15 <- ItoT_slope(inds = n_sSamples,regions = TRUE,start= 2004)

  
  t_N_15 <- ItoT(inds = NSamples,regions = FALSE,start= 2004)

  t_NS <- ItoT_slope(inds = NSamples,regions = FALSE)
  t_NS_15 <- ItoT_slope(inds = NSamples,regions = FALSE,start= 2004)
  
  
  
  t_N_comp <- ItoT(inds = N_compSamples,regions = FALSE,start = 1978)
  
  t_N_comp_15 <- ItoT(inds = N_compSamples,regions = FALSE,start= 2004)
  
  t_N_compS <- ItoT_slope(inds = N_compSamples,regions = FALSE)
  t_N_compS_15 <- ItoT_slope(inds = N_compSamples,regions = FALSE,start= 2004)
  
  

  trend_out <- bind_rows(t_N,
                         t_N_15,
                         t_NS,
                         t_NS_15,
                         t_N_comp,
                         t_N_comp_15,
                         t_N_compS,
                         t_N_compS_15,
                         t_n_s,
                         t_n_sS,
                         t_n_s_15,
                         t_n_sS_15)
  
  write.csv(trend_out,file = paste0("Trends/trends_slope_",sp,".csv"),row.names = F)
  
  # plotting indices --------------------------------------------------------
  
  
  # plot_Hyper <- plot_ind(inds = N_inds,
  #                        #smooth_inds = ,
  #                        raw = dts,
  #                        add_observed = TRUE,
  #                        add_samplesize = TRUE,
  #                        species = sp,
  #                        regions = FALSE,
  #                        title_size = 20,
  #                        axis_title_size = 18,
  #                        axis_text_size = 16)  
  # 
  # pdf(file = paste0("Figures/",sp,"_Hyperparameter_slope.pdf"),
  #     width = 8.5,
  #     height = 11)
  # print(plot_Hyper)
  # dev.off()
  
  
  
  
#   plot_by_st <- plot_ind(inds = n_inds_a2,
#                          #smooth_inds = ,
#                          raw = dts,
#                          add_observed = TRUE,
#                          add_samplesize = TRUE,
#                          species = sp,
#                          regions = TRUE,
#                          title_size = 20,
#                          axis_title_size = 18,
#                          axis_text_size = 16)  
#   
#   pdf(file = paste0("Figures/",sp,"_A2_slope.pdf"),
#       width = 8.5,
#       height = 11)
# print(plot_by_st)
# dev.off()


plot_by_st <- plot_ind(inds = n_inds,
                       smooth_inds = NULL,
                       raw = dts,
                       add_observed = TRUE,
                       add_samplesize = TRUE,
                       species = sp,
                       regions = TRUE,
                       title_size = 20,
                       axis_title_size = 18,
                       axis_text_size = 16)  

pdf(file = paste0("Figures/",sp,"_A1_slope.pdf"),
    width = 8.5,
    height = 11)
print(plot_by_st)
dev.off()




}# end if output exists

}

