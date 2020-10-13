# SPECIES MCMC data-prep -------------------------------------------------------
library(tidyverse)
library(ggmcmc)
library(rstan)
library(shinystan)
library(sf)
library(spdep)


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
 all_sites = ssData %>% distinct(SurveyAreaIdentifier,DecimalLatitude,DecimalLongitude)
 
 
 
 
 iss_sites = st_as_sf(all_sites,coords = c("DecimalLongitude","DecimalLatitude"), crs = 4326)
 laea = st_crs("+proj=laea +lat_0=40 +lon_0=-95") # Lambert equal area
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
 
 
# source("functions/GAM_basis_function.R")

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
sp = sps[25]

dts <- filter(ssData,CommonName == sp,
              YearCollected > 1977)
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
            span_years = yspn(YearCollected))

#number of sites with > 5 year span by region
nsites_w5 <- nyrs_site %>% 
  filter(span_years > 5) %>%  
  group_by(Region) %>% 
  summarise(nobs = n()) 


#sites with > 5 years of non-zero, observations
sites_keep <- nyrs_site[which(nyrs_site$span_years >= 5),"SurveyAreaIdentifier"]

# drop sites with <5 year span of data ------------------------------------


dts <- filter(dts,SurveyAreaIdentifier %in% sites_keep$SurveyAreaIdentifier) 


#number of years with non-zero observations, by region
nyrs_region <- dts %>% 
  filter(present == TRUE) %>%  
  group_by(hex_name,YearCollected) %>% 
  summarise(nobs = n()) %>% 
  group_by(hex_name) %>% 
  summarise(nyears = n(),
            span_years = yspn(YearCollected))

#strats with 7 or more years of non-zero, observations - species has to be observed in a region in at least 2/3 of the years in the time-series
regions_keep <- nyrs_region[which(nyrs_region$nyears >= 7),"hex_name"]


# drop strata with < 7 years of non-zero observations ---------------------


dts <- filter(dts,hex_name %in% regions_keep$hex_name) 

real_grid <- poly_grid %>% filter(hex_name %in% regions_keep$hex_name) 

strats_dts <- data.frame(hex_name = real_grid$hex_name,
                         stratn = 1:length(real_grid$hex_name))

dts <- left_join(dts,strats_dts,by = "hex_name")

# generate neighbourhoods -------------------------------------------------

coords = st_coordinates(st_centroid(real_grid))

nb_db <- spdep::tri2nb(coords,row.names = real_grid$hex_name)

nb_info = spdep::nb2WB(nb_db)

### re-arrange GEOBUGS formated nb_info into appropriate format for Stan model
car_stan_dat <- mungeCARdata4stan(adjBUGS = nb_info$adj,
                                  numBUGS = nb_info$num)



# DOY season definition ---------------------------------------------------


fday = min(dts$doy)-1

syear = min(dts$YearCollected)

dts <- dts %>% mutate(count = as.integer(ObservationCount),
                      year = as.integer(YearCollected),
                      yr = as.integer(year-syear),
                      strat = stratn,
                      date = doy-fday,
                      site = as.integer(factor(SurveyAreaIdentifier))) 

nstrata = max(dts$strat)
nsites = max(dts$site)

### commented out site-identification by  regions - shouldn't be necessary with such small strata
# sByReg = unique(dts[,c("SurveyAreaIdentifier","strat")])
# sByReg <- arrange(sByReg,strat,SurveyAreaIdentifier)
# 
# nsites <- table(sByReg$strat)
# 
# site <- NULL
# for(s in 1:nstrata){
#   site <- c(site,1:nsites[s]) 
# }
# sByReg$site <- site
# sByReg <- select(sByReg,-strat)
# dts <- left_join(dts,sByReg,by = "SurveyAreaIdentifier")



ncounts = nrow(dts)

# dmin_pred <- tapply(dts$date,dts$strat,min)
# dmax_pred <- tapply(dts$date,dts$strat,max)
# 
# dmax = max(dts$date)
nyears = max(dts$yr)
midyear = floor(nyears/2)

# # GAM seasonal basis function ---------------------------------------------
# 
# nKnots_season = 5
# basis_season <- gam.basis.func(orig.preds = as.integer(unlist(dts[,"date"])),
#                                nknots = nKnots_season,
#                                standardize = "z",
#                                random = F,
#                                npredpoints = max(dts$date),
#                                even_gaps = FALSE,
#                                sm_name = "season")
# 
# ndays <- basis_season$npredpoints_season
# 



# prepare stan data -------------------------------------------------------




stan_data <- list(count = as.integer(unlist(dts$count)),
                  year = as.integer(unlist(dts$yr-midyear)),
                  site = as.integer(unlist(dts$site)),
                  strat = as.integer(unlist(dts$strat)),
                  
                  #nyears = nyears,
                  nstrata = nstrata,
                  nsites = nsites,
                  ncounts = ncounts,
                  #ndays = ndays,
                  
                  # season_basis = basis_season$season_basis,
                  # season_basispred = basis_season$season_basispred,
                  # nknots_season = basis_season$nknots_season,
                  
                  #midyear = midyear,
                  
                  N_edges = car_stan_dat$N_edges,
                  node1 = car_stan_dat$node1,
                  node2 = car_stan_dat$node2)



parms = c("sdnoise",
          #"nu", #
          "b",
          "sdalpha")
mod.file = "models/slope_iCAR.stan"

## compile model
slope_icar_model = stan_model(file=mod.file)

## run sampler on model, data
stime = system.time(slope_icar_stanfit <-
                      sampling(slope_icar_model,seed=12345,
                               data=stan_data,
                               verbose=TRUE, refresh=100,
                               chains=1, iter=4000,
                               cores = 1,
                               pars = parms,
                               control = list(adapt_delta = 0.99,
                                              max_treedepth = 18)))

launch_shinystan(stime) 







save(list = c("jags_data",
              "basis_season",
              "dts",
              "t2",
              "t1",
              "out2"),
     file = paste0("output/",sp,"slope_ZIP_results.RData"))





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

