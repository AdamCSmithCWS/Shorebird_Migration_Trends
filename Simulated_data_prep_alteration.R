# Simulated MCMC data-prep -------------------------------------------------------
library(tidyverse)
library(rstan)
rstan_options(auto_write = TRUE, javascript = FALSE)
library(sf)
library(spdep)
library(sgt) #includes the generalized t-distribution

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
 laea = st_crs("+proj=laea +lat_0=40 +lon_0=-95") # Lambert equal area
 
 
 
 locat = system.file("maps",
                     package="bbsBayes")
 map.file = "BBS_ProvState_strata"
 
 prov_state = sf::read_sf(dsn = locat,
                              layer = map.file)
 
 orig_ss_regions = sf::read_sf(dsn = "data",
                               layer = "region_polygons")
 orig_ss_regions <- select(orig_ss_regions,Region)
 orig_ss_regions = st_transform(orig_ss_regions,crs = laea)
 
 # bbs_strata_proj = st_transform(bbs_strata_map,crs = 102003) #USA Contiguous albers equal area projection stanadard from ESRI - https://mgimond.github.io/Spatial/coordinate-systems-in-r.html
 
 prov_state = st_transform(prov_state,crs = laea)
 
 
 iss_sites = st_as_sf(all_sites,coords = c("DecimalLongitude","DecimalLatitude"), crs = 4326)
 iss_sites_lcc <- st_transform(iss_sites, laea)
 bb = st_bbox(iss_sites_lcc) %>% 
   st_as_sfc()

 grid_spacing <- 300000  # size of squares, in units of the CRS (i.e. meters for lae)
 
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


 #### 
 

# Simulation Characteristics ----------------------------------------------
# 
#  sp = "Red Knot" #example species with moderately large spatial distribution
#  # and reasonably difficult overdispersion
#  
#  FYYYY = 1980
#  
# load(paste0("data/data",sp,"_GAMYE_strat_simple",grid_spacing/1000,".RData"))
# load(paste0("output/",sp,"_GAMYE_strat_simple",grid_spacing/1000,".RData"))
# 
# full_sum <- summary(slope_icar_stanfit)$summary
# 
# 
# 
#  # stable population with no year-effects
#  # steadily increasing population with specific year effects
#  # stable population with change point in 2000 to increasing
#  # Simulate data using estimated values from a species analysis and plug into model (REKN?)
#  
# 
#  # convincing to simulate a population with the same abundance distribution in space, but with the opposite temporal pattern to what we've seen in many of the results
#  ALPHA1 <- full_sum["ALPHA1","mean"]
#  sdalpha <- full_sum["sdalpha","mean"]
#  alphas <- full_sum[paste0("alpha[",1:stan_data$nsites,"]"),"mean"]
#  
#  set_seed(2019)
#  year_effect <- rep(c(-1,1),length = stan_data$nyears)*abs(rnorm(stan_data$nyears,0,full_sum[paste0("sdyear"),"mean"]))
#  year_effect[c(1,stan_data$nyears)] <- 0
#  ## year-effects are 0 in first and last years, and alternate positive-negative for all other years
#  season_pred1 <- full_sum[paste0("season_pred[",1:stan_data$ndays,",",1,"]"),"mean"]
#  season_pred2 <- full_sum[paste0("season_pred[",1:stan_data$ndays,",",2,"]"),"mean"]
#  
#  season_pred <- matrix(c(season_pred1,season_pred2),
#                        ncol = 2)
#  sdnoise <- full_sum["sdnoise","mean"]
#  #nu <- 4#full_sum["nu","mean"]
#   

 
 load(paste0("data/data_simulated_stable_uptick_GAMYE_strat_simple",grid_spacing/1000,".RData"))
 sdnoise <- 1.970721
 noise = rnorm(stan_data$ncounts,0,sdnoise)
 year_pred <- matrix(NA,nrow = stan_data$nyears,ncol = stan_data$nstrata)
 
 B1 = 0 #initial stable linear trend
 BC = 0.04 #change in trend in year 20
 sdBC = 0.005 #sd on change in trend
 BCs = rnorm(stan_data$nstrata,BC,sdBC)
 bx = 0.005 #implies a 2%/year range in trends by longitude
 by = 0.01 #implies a 2%/year range in trends by latitude
 
 strats_xy <- strats_dts %>% 
   mutate(x = as.numeric(str_split(hex_name,pattern = "_",simplify = TRUE)[,1]),
          y = as.numeric(str_split(hex_name,pattern = "_",simplify = TRUE)[,2]),
          x_scale = (x-mean(x))/(0.5*diff(range(x))),#scaled x coordinate to create a longitudinal gradient in trends
          y_scale = (y-mean(y))/(0.5*diff(range(y))),#scaled y coordinate to create a latitudinal gradient in trends
          b1 = x_scale*(bx+B1) + y_scale*(by+B1), #initial log-linear slopes for each stratum with spatial gradients
          b2 = b1+BCs) #final log-linear slope with spatial gradients
 
midyear = 20
b1 = strats_xy$b1 #initial slopes for each stratum
b2 = strats_xy$b2 #second slopes for each stratum

#generating a linear breakpoint trajectory for each stratum
for(s in 1:stan_data$nstrata){
year_pred[1:midyear,s] <- b1[s]*((1:midyear)-midyear)
year_pred[(midyear+1):stan_data$nyears,s] <- b2[s]*((midyear+1):stan_data$nyears-midyear)
}

#noise = rsgt(stan_data$ncounts,mu = 0,sigma = sdnoise,q = 4,p = 4) #p*q = df

### for the uptick scenario, must reduce the alphas, otherwise the increase
### leads to incredibly unrealistic counts
#alphas <- (alphas-3)

stan_data_sim = stan_data
for(i in 1:stan_data_sim$ncounts){
stan_data_sim$count[i] = rpois(n = 1,exp(ALPHA1 + year_pred[stan_data_sim$year_raw[i],stan_data_sim$strat[i]] + 
  alphas[stan_data_sim$site[i]] + year_effect[stan_data_sim$year_raw[i]] + 
  season_pred[stan_data_sim$date[i],stan_data_sim$seas_strat[i]] + noise[i]));

}
 
plot(x = strats_xy$y_scale,strats_xy$b2) 
plot(x = strats_xy$x_scale,strats_xy$b2) 

t2_map <- ggplot(data = strats_xy,aes(x = x_scale,y = y_scale))+
   geom_point(aes(colour = b1+b2))+
   scale_colour_viridis_c()
print(t2_map)

stan_data <- stan_data_sim


save(list = c("stan_data",
                   "dts",
                   "real_grid",
                   "real_grid_regs",
                   "strats_dts",
                   "strat_regions",
                   "mod.file",
                   "parms",
              "strats_xy",
              "season_pred",
              "alphas",
              "ALPHA1",
              "year_pred",
              "year_effect",
              "noise"),
          file = paste0("data/data_simulated_stable_upturn_GAMYE_strat_simple.RData"))

 
sites_by_yr <- data.frame(year = stan_data$year_raw,
                          site = stan_data$site,
                          strat = stan_data$strat)

sites_by_yr$alpha <- alphas[sites_by_yr$site]

mn_sites_yr = sites_by_yr %>% 
   group_by(year,strat) %>% 
   summarise(mn = mean(alpha))

 
sites_yr_plot = ggplot(data = mn_sites_yr,aes(x = year,y = mn))+
   geom_point()+
   geom_smooth(method = "lm")+
   facet_wrap(~strat,nrow = 5,ncol = 6)
   
 print(sites_yr_plot)
 
 slp <- function(x,y){
    df = data.frame(x = x,y = y)
    df = na.omit(df)
    m = lm(y~x,data = df)
    return(as.numeric(m$coefficients["x"]))
 }
 
slope_sites <- mn_sites_yr %>% group_by(strat) %>% 
   summarise(slp = slp(y = mn,x = year)) %>% 
   left_join(.,strats_xy,by = c("strat" = "stratn"))

plot(slope_sites$slp,slope_sites$b2)
plot(slope_sites$slp,slope_sites$b1+slope_sites$b2)




# Balanced simulated dataset ----------------------------------------------

stan_data_bal <- stan_data_sim

sites_strats <- unique(data.frame(strat = stan_data_sim$strat,
                                  site = stan_data_sim$site,
                                  seas_strat = stan_data_sim$seas_strat))

sample_days = floor(seq(6,stan_data_sim$ndays-6,length = 12))

bal_days <- NULL
for(j in 1:nrow(sites_strats)){
   tmp = sites_strats[j,]
   tmp2 = sample_days + round(runif(length(sample_days),-5,5))
   tmp = expand_grid(tmp,date = as.integer(tmp2))
   
   bal_days <- bind_rows(bal_days,tmp)
}

bal_dat <- expand_grid(bal_days,year_raw = c(1:stan_data_sim$nyears))


stan_data_bal$ncounts <- nrow(bal_dat)
stan_data_bal$year_raw <- bal_dat$year_raw
stan_data_bal$site <- bal_dat$site
stan_data_bal$strat <- bal_dat$strat
stan_data_bal$date <- bal_dat$date
stan_data_bal$seas_strat <- bal_dat$seas_strat
noise = rnorm(stan_data_bal$ncounts,0,sdnoise)

stan_data_bal$count <- vector(mode = "integer",length = stan_data_bal$ncounts)
for(i in 1:stan_data_bal$ncounts){
   t1 = rpois(n = 1,exp(ALPHA1 + year_pred[stan_data_bal$year_raw[i],stan_data_bal$strat[i]] + 
                                               alphas[stan_data_bal$site[i]] + year_effect[stan_data_bal$year_raw[i]] + 
                                               season_pred[stan_data_bal$date[i],stan_data_bal$seas_strat[i]] + noise[i]))
   if(is.na(t1)){break}
   stan_data_bal$count[i] = t1
}

stan_data <- stan_data_bal


save(list = c("stan_data",
              "dts",
              "real_grid",
              "real_grid_regs",
              "strats_dts",
              "strat_regions",
              "mod.file",
              "parms",
              "strats_xy",
              "season_pred",
              "alphas",
              "ALPHA1",
              "year_pred",
              "year_effect",
              "noise"),
     file = paste0("data/data_simulated_stable_Balanced_upturn_GAMYE_strat_simple.RData"))


