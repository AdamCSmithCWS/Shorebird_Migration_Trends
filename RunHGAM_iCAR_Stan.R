# SPECIES MCMC data-prep -------------------------------------------------------
library(tidyverse)
library(rstan)
rstan_options(auto_write = TRUE, javascript = FALSE)
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
 
 
 # generate site-size predictors for species groups ------------------------
 
 sp_groups = read.csv("Flock_Sizes.csv") 
 
 # sp_groups$two_seasons <- FALSE
 # sp_groups$two_seasons[which(sp_groups$Species %in% c("Willet",
 #                                                      "Short-billed Dowitcher",
 #                                                      "Sanderling",
 #                                                      "Black-bellied Plover",
 #                                                      "Red Knot",
 #                                                      "Piping Plover",
 #                                                      "Rudy Turnstone",
 #                                                      "Hudsonian Godwit"))] <- TRUE
 # 
 # sp_groups$regions_w_alt_season <- NA
 # sp_groups$regions_w_alt_season[which(sp_groups$Species %in% c("Willet",
 #                                                      "Short-billed Dowitcher",
 #                                                      "Sanderling",
 #                                                      "Red Knot",
 #                                                      "Piping Plover",
 #                                                      "Rudy Turnstone"))] <- "Southeast Coastal - Texas Coastal"
 # 
 # sp_groups$regions_w_alt_season[which(sp_groups$Species %in% c( "Black-bellied Plover"))] <- "Southeast Coastal - Texas Coastal - Pacific and Intermountain"
 # sp_groups$regions_w_alt_season[which(sp_groups$Species %in% c(  "Hudsonian Godwit"))] <- "East Inland - Ontario"
 # 
 #
 # write.csv(sp_groups,"Flock_Sizes.csv",row.names = FALSE)
 
 
 
 sp_groups$flock <- ifelse(is.na(sp_groups$Large_flocks),"small","large")
 
 ssData <- inner_join(ssData,sp_groups[,c("Species","flock")],
                      by = c("CommonName" = "Species"))
 
 q_max_flock <- ssData %>% group_by(SurveyAreaIdentifier,flock) %>% 
   summarise(max_count = max(ObservationCount),
             Q98_count = quantile(ObservationCount,0.98),
             dif = max_count-Q98_count) %>% 
   mutate(size = ifelse(Q98_count == 0,log(max_count/2),log(Q98_count))) 
 # 98th percentile of all counts for species in size-group OR half of the maximum count if 98th percentile == 0
 # half of max limit only relevant for a handful of sites
 
 q_m_large = filter(q_max_flock,flock == "large",
                    max_count > 0)
 
 
 q_m_small = filter(q_max_flock,flock == "small",
                    max_count > 0)
 
 species_large <- sp_groups[which(sp_groups$flock == "large"),"Species"]
 species_small <- sp_groups[which(sp_groups$flock == "small"),"Species"]
 
 
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
 mean_counbts_doy_out <- vector(mode = "list",length = length(sps))
 names(mean_counbts_doy_out) <- sps
 ggp_out <- vector(mode = "list",length = length(sps))
 names(ggp_out) <- sps
 mean_counbts_year_out <- vector(mode = "list",length = length(sps))
 names(mean_counbts_year_out) <- sps
 
 
for(sp in sps){
#sp = sps[11]
FYYYY = 1980


if(sp %in% species_large){site_sizes <- q_m_large}
if(sp %in% species_small){site_sizes <- q_m_small}

if(sp_groups[which(sp_groups$Species == sp),"two_seasons"]){
  two_seasons <- TRUE
  regs_alt_season <- as.character(str_split(sp_groups[which(sp_groups$Species == sp),"regions_w_alt_season"],pattern = " - ",simplify = TRUE))
}



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

# drop sites with <10 year span of data ------------------------------------


dts <- filter(dts,SurveyAreaIdentifier %in% sites_keep$SurveyAreaIdentifier) 


#number of years with non-zero observations, by region
nyrs_region <- dts %>% 
  filter(present == TRUE) %>%  
  group_by(hex_name,YearCollected) %>% 
  summarise(nobs = n()) %>% 
  group_by(hex_name) %>% 
  summarise(nyears = n(),
            span_years = yspn(YearCollected),
            fyear = min(YearCollected),
            lyear = max(YearCollected))

#strats with 7 or more years of non-zero, observations - species has to be observed in a region in at least 2/3 of the years in the time-series
p_time_series = 0.5 #only required to have dat that span 20% of the time-series in a region.
regions_keep <- nyrs_region[which(nyrs_region$span_years >= nyrs_study*p_time_series),"hex_name"]


# drop strata with < 7 years of non-zero observations ---------------------


dts <- filter(dts,hex_name %in% regions_keep$hex_name) 

real_grid <- poly_grid %>% filter(hex_name %in% regions_keep$hex_name) 

dom_value <- function(x){
  tt = table(x)
  tt = sort(tt)
  return(names(tt)[1])
}
hex_by_reg <- dts %>% group_by(hex_name) %>% 
  summarise(Region = dom_value(Region))

real_grid_regs <- left_join(real_grid, hex_by_reg,by = "hex_name")

dts <- dts %>% select(-Region) %>% 
  left_join(.,hex_by_reg)

strats_dts <- data.frame(hex_name = real_grid$hex_name,
                         stratn = 1:length(real_grid$hex_name))

dts <- left_join(dts,strats_dts,by = "hex_name")
real_grid <- inner_join(real_grid,strats_dts)
real_grid_regs <- inner_join(real_grid_regs,strats_dts)

# Add site size offsets ---------------------------------------------------


dts <- left_join(dts,site_sizes[,c("SurveyAreaIdentifier","size")])



# DOY season definition ---------------------------------------------------


fday = min(dts$doy)-1

syear = min(dts$YearCollected)

dts <- dts %>% mutate(count = as.integer(ObservationCount),
                      year = as.integer(YearCollected),
                      yr = as.integer(year-(syear-1)),
                      strat = stratn,
                      date = doy-fday,
                      site = as.integer(factor(SurveyAreaIdentifier)),
                      size_cent = as.numeric(scale(size))) 

nstrata = max(dts$strat)
nsites = max(dts$site)

sizes_by_site <- unique(dts[,c("site","size_cent")])
sizes_by_site <- sizes_by_site[order(sizes_by_site$site),]


## indexing of sites by strata for annual index calculations
sByReg = unique(dts[,c("site","strat")])
sByReg <- arrange(sByReg,strat,site)
# 
nsites_strat <- table(sByReg$strat)
max_sites <- max(nsites_strat)
sites <- matrix(data = 0,nrow = max_sites,ncol = nstrata)
for(j in 1:nstrata){
  sites[1:nsites_strat[j],j] <- as.integer(unlist(sByReg[which(sByReg$strat == j),"site"]))
}



# generate neighbourhoods -------------------------------------------------

# Voronoi polygons from strata-centroids -----------------------------------
# voronoi polygons ensures all strata have neighbours
reg_bounds <- st_union(orig_ss_regions)
reg_bounds_buf = st_buffer(reg_bounds,dist = grid_spacing)

centres = st_centroid(real_grid)
#centres_buf <- st_buffer(centres, dist=100)
coords = st_coordinates(centres)


box <- st_as_sfc(st_bbox(centres))

cun = st_union(centres)

v <- try(st_cast(st_voronoi(cun, envelope = box)),silent = TRUE)

# Adding random noise to hex centres to avoid geometry errors -------------


while(class(v)[1] == "try-error"){
 
  #tmp = st_geometry(centres)
  centres = st_centroid(real_grid)
  for(i in 1:nrow(centres)){
    centres$geometry[[i]] <- centres$geometry[[i]] + st_point(rnorm(2,0,1))
  }
  #centres_buf <- st_buffer(centres, dist=100)
  coords = st_coordinates(centres)
  
  
  box <- st_as_sfc(st_bbox(centres))
  
  cun = st_union(centres)
  
  v <- try(st_cast(st_voronoi(cun, envelope = box)),silent = TRUE)
  
  
  
}

vint = try(st_sf(st_cast(st_intersection(v,reg_bounds_buf),"POLYGON")),silent = TRUE)

# Adding random noise to hex centres to avoid geometry errors -------------

j = 0
while(class(vint)[1] == "try-error" & j < 10){

  #tmp = st_geometry(centres)
  centres = st_centroid(real_grid)
  for(i in 1:nrow(centres)){
    centres$geometry[[i]] <- centres$geometry[[i]] + st_point(rnorm(2,0,1))
  }
  #centres_buf <- st_buffer(centres, dist=100)
  coords = st_coordinates(centres)
  
  
  box <- st_as_sfc(st_bbox(centres))
  
  cun = st_union(centres)
  
  v <- try(st_cast(st_voronoi(cun, envelope = box)),silent = TRUE)
  
  
  reg_bounds_buf = st_buffer(reg_bounds_buf,dist = 100)
  vint = try(st_sf(st_cast(st_intersection(v,reg_bounds_buf),"POLYGON")),silent = TRUE)
  j = j+1
  
  }

if(j == 10){stop(paste("Geometry ERRORS check vintj for",sp))}

# tmp = ggplot()+
#   geom_sf(data = cun)+
#   geom_sf(data = v,alpha = 0.1)
#   #geom_sf(data = centres)
# print(tmp)

vintj = st_join(vint,centres,join = st_contains)
vintj = arrange(vintj,stratn)

nb_db = poly2nb(vintj,row.names = vintj$stratn,queen = FALSE)


### currently using 2 nearest neighbours to define the spatial relationships
## many regions may  have > 2 neighbours because of the symmetry condition
# nb_db <- spdep::knn2nb(spdep::knearneigh(coords,k = 4),row.names = route_map$route,sym = TRUE)
cc = st_coordinates(st_centroid(vintj))
#

ggp = ggplot(data = real_grid_regs)+
  geom_sf(data = vintj,alpha = 0.3)+ 
  geom_sf(aes(col = Region))+
  geom_sf_text(aes(label = stratn),size = 3,alpha = 0.3)+
  labs(title = sp)
pdf(file = paste0("Figures/",sp,"strata_connections ",grid_spacing/1000,".pdf"))
plot(nb_db,cc,col = "pink")
text(labels = rownames(cc),cc ,pos = 2)
print(ggp)
dev.off()

ggp_out[[sp]] <- ggp

# wca = which(grepl(route_map$strat,pattern = "-CA-",fixed = T))
# wak = which(grepl(route_map$strat,pattern = "-AK-",fixed = T))
# 
# nb2[[wak]]



mean_counbts_doy = ggplot(data = dts,aes(x = doy,y = count+1,colour = Region))+
  scale_y_log10()+
  geom_point(position = position_jitter(width = 0.1,height = 0),alpha = 0.3)+
  geom_smooth()+
  labs(title = sp)+
  facet_wrap(facets = ~strat,nrow = 8,ncol = 5,scales = "free_y")

pdf(paste0("Figures/",sp,"seasonal_counts ",grid_spacing/1000,".pdf"),
    width = 11,height = 8.5)
print(mean_counbts_doy)
dev.off()
mean_counbts_doy_out[[sp]] <- mean_counbts_doy



nb_info = spdep::nb2WB(nb_db)



### re-arrange GEOBUGS formated nb_info into appropriate format for Stan model
car_stan_dat <- mungeCARdata4stan(adjBUGS = nb_info$adj,
                                  numBUGS = nb_info$num)







mean_counbts_year = ggplot(data = dts,aes(x = year,y = count+1,colour = Region))+
  scale_y_log10()+
  geom_point(position = position_jitter(width = 0.1,height = 0),alpha = 0.3)+
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs",k = 5))+
  labs(title = sp)+
  facet_wrap(facets = ~strat,nrow = 8,ncol = 5,scales = "free_y")

pdf(paste0("Figures/",sp,"annual_counts ",grid_spacing/1000,".pdf"),
    width = 11,height = 8.5)
print(mean_counbts_year)
dev.off()
mean_counbts_year_out[[sp]] <- mean_counbts_year



# Identifying the strata and region combinations --------------------------
strat_regions <- unique(dts[,c("Region","strat")])
strat_regions <- strat_regions[order(strat_regions$strat),]




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

nKnots_year = ceiling(nyears/3)
basis_year <- gam.basis.func(orig.preds = as.integer(unlist(dts[,"yr"])),
                               nknots = nKnots_year,
                               standardize = "z",
                               random = F,
                               npredpoints = max(dts$yr),
                               even_gaps = TRUE,
                               sm_name = "year")

# 


if(two_seasons){
  
  strat_regions$seas_strat <- 1
  strat_regions$seas_strat[which(strat_regions$Region %in% regs_alt_season)] <- 2
  dts <- left_join(dts,strat_regions,by = c("Region","strat"))

  seasons <- as.matrix(strat_regions[,c("strat","seas_strat")])
  
  stan_data <- list(count = as.integer(unlist(dts$count)),
                    year = as.integer(unlist(dts$yr-midyear)),
                    year_raw = as.integer(unlist(dts$yr)),
                    site = as.integer(unlist(dts$site)),
                    strat = as.integer(unlist(dts$strat)),
                    date = as.integer(unlist(dts$date)),
                    seas_strat = as.integer(unlist(dts$seas_strat)),
                    
                    nyears = nyears,
                    nstrata = nstrata,
                    nsites = nsites,
                    ncounts = ncounts,
                    ndays = ndays,
                    
                    nsites_strat = nsites_strat,
                    max_sites = max_sites,
                    sites = sites,
                    seasons = seasons,
                    
                    site_size = sizes_by_site$size_cent,
                    
                    # season_basis = basis_season$season_basis,
                    season_basispred = basis_season$season_basispred,
                    nknots_season = basis_season$nknots_season,
                    
                    year_basispred = basis_year$year_basispred,
                    nknots_year = basis_year$nknots_year,
                    
                    #midyear = midyear,
                    
                    N_edges = car_stan_dat$N_edges,
                    node1 = car_stan_dat$node1,
                    node2 = car_stan_dat$node2)
  mod.file = "models/GAMYE_strata_two_season.stan"
  
  parms = c("sdnoise",
            #"nu", #
            "sdalpha",
            "b",
            "B",
            "alpha",
            "ALPHA1",
            "sdyear",
            "sdyear_gam_strat",
            "sdyear_gam",
            "year_effect",
            "sdseason",
            "B_season_raw1",
            "B_season_raw2",
            "season_pred",
            "n",
            "nsmooth",
            "N",
            "NSmooth",
            "log_lik",
            "beta_size")
  
}else{
  stan_data <- list(count = as.integer(unlist(dts$count)),
                    year = as.integer(unlist(dts$yr-midyear)),
                    year_raw = as.integer(unlist(dts$yr)),
                    site = as.integer(unlist(dts$site)),
                    strat = as.integer(unlist(dts$strat)),
                    date = as.integer(unlist(dts$date)),
                    
                    nyears = nyears,
                    nstrata = nstrata,
                    nsites = nsites,
                    ncounts = ncounts,
                    ndays = ndays,
                    
                    nsites_strat = nsites_strat,
                    max_sites = max_sites,
                    sites = sites,
                    
                    site_size = sizes_by_site$size_cent,
                    
                    # season_basis = basis_season$season_basis,
                    season_basispred = basis_season$season_basispred,
                    nknots_season = basis_season$nknots_season,
                    
                    year_basispred = basis_year$year_basispred,
                    nknots_year = basis_year$nknots_year,
                    
                    #midyear = midyear,
                    
                    N_edges = car_stan_dat$N_edges,
                    node1 = car_stan_dat$node1,
                    node2 = car_stan_dat$node2)
  
  mod.file = "models/GAMYE_strata.stan"
  
  parms = c("sdnoise",
            #"nu", #
            "sdalpha",
            "b",
            "B",
            "alpha",
            "ALPHA1",
            "sdyear",
            "sdyear_gam_strat",
            "sdyear_gam",
            "year_effect",
            "sdseason",
            "B_season_raw",
            "season_pred",
            "n",
            "nsmooth",
            "N",
            "NSmooth",
            "log_lik",
            "beta_size")
}

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

save(list = c("stan_data",
              "dts",
              "real_grid",
              "strats_dts",
              "strat_regions",
              "mod.file",
              "parms"),
     file = paste0("data/data",sp,"_GAMYE_strat_offset_",grid_spacing/1000,".RData"))
}


 pdf(paste0("Figures/","All_seasonal_counts ",grid_spacing/1000,".pdf"),
     width = 11,height = 8.5)
 for(sp in sps){
   print(mean_counbts_doy_out[[sp]]+
           labs(title = sp))
   
 }
 dev.off()
 
 pdf(file = paste0("Figures/","ALL_Strata_",grid_spacing/1000,".pdf"))
 
 for(sp in sps){
   print(ggp_out[[sp]]+
           labs(title = sp))
   
 }
 
 dev.off()
 
 pdf(paste0("Figures/","All_annual_counts ",grid_spacing/1000,".pdf"),
     width = 11,height = 8.5)
 for(sp in sps){
   print(mean_counbts_year_out[[sp]]+
           labs(title = sp))
   
 }
 dev.off()
 
 
 for(sp in sps){
  
   
    load(paste0("data/data",sp,"_GAMYE_strat_offset_",grid_spacing/1000,".RData"))



## compile model
slope_icar_model = stan_model(file=mod.file)

print(sp)
## run sampler on model, data
slope_icar_stanfit <- sampling(slope_icar_model,
                               data=stan_data,
                               verbose=TRUE, refresh=100,
                               chains=4, iter=1500,
                               warmup=1000,
                               cores = 4,
                               pars = parms,
                               control = list(adapt_delta = 0.85,
                                              max_treedepth = 14))





save(list = c("slope_icar_stanfit",
              "stan_data",
              "dts",
              "real_grid",
              "strats_dts",
              "strat_regions",
              "mod.file",
              "parms"),
     file = paste0("output/",sp,"_GAMYE_strat_offset",grid_spacing/1000,".RData"))


}#end modeling loop

launch_shinystan(slope_icar_stanfit) 

load(paste0(sp,"_GAMYE_strat_offset",grid_spacing/1000,".RData"))


# CONSIDERATIONS ----------------------------------------------------------

## map the looic values to see if some areas are being poorly predicted

## add some variation to the plotting of the observed means (box and whisker info)

## assess the fit of the site adjuestment model and the non-adjusted model



# Calculate Annual indices using samples ----------------------------------


source("functions/utility_functions.R")


Nsamples <- slope_icar_stanfit %>% gather_draws(N[y])
Nsamples$year <- Nsamples$y + (syear-1)

NSmoothsamples <- slope_icar_stanfit %>% gather_draws(NSmooth[y])
NSmoothsamples$year <- NSmoothsamples$y + (syear-1)



# calculate trends continent --------------------------------------------------------


t_NSmooth_90 <- ItoT(inds = NSmoothsamples,
                  start = 1990,
                  end = 2019,
                  regions = NULL,
                  qs = 95,
                  trend_type = "endpoint",
                  index_type = "smooth",
                  retransformation_type = "standard")
t_NSmooth_80 <- ItoT(inds = NSmoothsamples,
                  start = 1980,
                  end = 2019,
                  regions = NULL,
                  qs = 95,
                  trend_type = "endpoint",
                  index_type = "smooth",
                  retransformation_type = "standard")

t_NSmooth_07 <- ItoT(inds = NSmoothsamples,
                        start = 2007,
                        end = 2019,
                        regions = NULL,
                        qs = 95,
                        index_type = "smooth")

anot_funct <- function(x){
  ant = paste(signif(x$percent_change,3),
        "% since",
        x$start_year,
        "[",
        signif(x$p_ch_lci,3),
        " : ",
        signif(x$p_ch_uci,3),
        "]")
}

anot_90 = anot_funct(t_NSmooth_90)
anot_80 = anot_funct(t_NSmooth_80)
anot_07 = anot_funct(t_NSmooth_07)


# Extract Annual Indices ------------------------------------------------


indicesN <- index_summary(parm = "N",
                          dims = "year",
                          site_scale = FALSE,
                          season_scale = TRUE)


indicesNSmooth <- index_summary(parm = "NSmooth",
                                dims = "year",
                                site_scale = FALSE,
                                season_scale = TRUE)





  
  indices = bind_rows(indicesN,indicesNSmooth)
indices$year = indices$year + (syear-1)
yup = max(max(indices$PI97_5)*1.5,quantile(indices$obsmean,0.95))

N_gg = ggplot(data = indices,aes(x = year, y = PI50,fill = parm))+
  geom_ribbon(aes(ymin = PI2_5,ymax = PI97_5),alpha = 0.2)+
  geom_line(aes(colour = parm))+
  labs(title = paste(sp,"GAMYE continental traj (full and smooth) with obs means"))+
  theme(legend.position = "none")+
  annotate("text", x = 1997, y = yup*0.9, label = anot_80)+
  annotate("text", x = 1997, y = yup*0.8, label = anot_90)+
  annotate("text", x = 1997, y = yup*0.7, label = anot_07)+
  coord_cartesian(ylim = c(0,yup))+
  geom_point(aes(y = obsmean),colour = grey(0.5),alpha = 0.3)


pdf(paste0("figures/",sp,FYYYY,"_GAMYE_survey_wide_trajectory",grid_spacing/1000,".pdf"))
print(N_gg)
dev.off()

indicesnsmooth <- index_summary(parm = "nsmooth",
                                dims = c("stratn","year"),
                                season_scale = FALSE,
                                site_scale = FALSE)

indicesn <- index_summary(parm = "n",
                          dims = c("stratn","year"),
                          season_scale = FALSE,
                          site_scale = FALSE)

indices_strat = bind_rows(indicesn,indicesnsmooth)
indices_strat$year = indices_strat$year + (syear-1)
indices_strat <- left_join(indices_strat,strats_dts, by = "stratn")

pdf(file = paste0("figures/", sp,FYYYY,"_GAMYE_Strata_trajectories",grid_spacing/1000,".pdf"),
    width = 8.5,
    height = 11)
print(N_gg)
for(jj in 1:ceiling(nstrata/12)){
n_gg = ggplot(data = indices_strat,aes(x = year, y = PI50,fill = parm))+
  geom_ribbon(aes(ymin = PI2_5,ymax = PI97_5),alpha = 0.2)+
  geom_line(aes(colour = parm))+
  geom_point(aes(y = obsmean),colour = grey(0.5),alpha = 0.1)+
  facet_wrap_paginate(facets = ~hex_name,page = jj,nrow = 4, ncol = 3,scales = "free")
print(n_gg)
}
  dev.off()





# combine the trajectories within the original regional strata ------------




# map the trend estimates -------------------------------------------------

  # real_grid
  

  reg_strats <- data.frame(hex_name = real_grid_regs$hex_name,
                           Region = real_grid_regs$Region)
  
 strats_dts <- left_join(strats_dts,reg_strats,by = "hex_name")
  
 
  nsmoothsamples <- slope_icar_stanfit %>% gather_draws(nsmooth[s,y])
  nsmoothsamples$year <- nsmoothsamples$y + (syear-1)
  nsmoothsamples <- left_join(nsmoothsamples,strats_dts,by = c("s" = "stratn"))
 # st_drop <- "2305040_1068494"
  #nsmoothsamples <- nsmoothsamples %>% filter(hex_name != st_drop)
  
  
  t_nsmooth_strat_07 <- ItoT(inds = nsmoothsamples,
                       start = 2007,
                       end = 2019,
                       regions = "hex_name",
                       qs = 95,
                       index_type = "smooth")
  
  t_nsmooth_reg_07 <- ItoT(inds = nsmoothsamples,
                       start = 2007,
                       end = 2019,
                       regions = "Region",
                       qs = 95,
                       index_type = "smooth")
  
  t_nsmooth_reg_80 <- ItoT(inds = nsmoothsamples,
                           start = 1980,
                           end = 2019,
                           regions = "Region",
                           qs = 95,
                           index_type = "smooth")
  
  
  t_nsmooth_reg_00 <- ItoT(inds = nsmoothsamples,
                           start = 2000,
                           end = 2019,
                           regions = "Region",
                           qs = 95,
                           index_type = "smooth")
  
  
  
  
  
  nsamples <- slope_icar_stanfit %>% gather_draws(n[s,y])
  nsamples$year <- nsamples$y + (syear-1)
  nsamples <- left_join(nsamples,strats_dts,by = c("s" = "stratn"))
  
  st_drop <- "2305040_1068494"
  nsamples <- nsamples %>% filter(hex_name != st_drop)
  
  indsn_region = ItoI(inds = nsmoothsamples,
               regions = "Region")
  
  ind_fc = ggplot(data = indsn_region,aes(x = year,y = median,group = region))+
    geom_ribbon(aes(ymin = lci,ymax = uci),alpha = 0.2)+
    geom_line()+
    coord_cartesian(ylim = c(0,NA))+
    facet_wrap(~region,scales = "free")
  
print(ind_fc)  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
slopes = as.data.frame(summary(slope_icar_stanfit,
                 pars = c("b"),
                 probs = c(0.025,0.5,0.975))$summary)

myrename = function(fit){
  rename_with(fit,~ paste0("PI",gsub(".","_",gsub("%", "", .x, fixed = TRUE), fixed = TRUE)),ends_with("%"))

}
slopes = myrename(slopes)
slopes$stratn = 1:nrow(slopes)
slopes$trend = (exp(slopes$mean)-1)*100

breaks <- c(-7, -4, -2, -1, -0.5, 0.5, 1, 2, 4, 7)
labls = c(paste0("< ",breaks[1]),paste0(breaks[-c(length(breaks))],":", breaks[-c(1)]),paste0("> ",breaks[length(breaks)]))
labls = paste0(labls, " %")

slopes$Tplot <- cut(slopes$trend,breaks = c(-Inf, breaks, Inf),labels = labls)

map_palette <- c("#a50026", "#d73027", "#f46d43", "#fdae61", "#fee090", "#ffffbf",
                 "#e0f3f8", "#abd9e9", "#74add1", "#4575b4", "#313695")
names(map_palette) <- labls

hex_map <- left_join(real_grid,strats_dts)

hex_map <- left_join(hex_map,slopes,by = "stratn")

trend_map <- ggplot()+
  geom_sf(data = prov_state,colour = grey(0.8))+
  geom_sf(data = hex_map,aes(fill = Tplot,colour = Tplot))+
  coord_sf(ylim = c(-1509123,1890567),expand = TRUE)+
  scale_colour_manual(values = map_palette, aesthetics = c("fill","colour"),
                      guide = ggplot2::guide_legend(reverse=TRUE),
                      name = "Trends")#paste0("Trend\n",fyr,"-",lyr))
print(trend_map)





save(list = c("stan_data",
              #"basis_season",
              "dts",
              "slope_icar_model",
              "hex_map",
              "strats_dts"),
     file = paste0("output/",sp,"slope_iCAR_results.RData"))





}#end species loops


stopCluster(cl = cluster)


# 
# 
# 
# # Plotting ----------------------------------------------------------------
# 
# library(tidybayes)
# 
# load("data/allShorebirdPrismFallCounts.RData")
# source("functions/Utility_functions.R")
# 
# for(sp in sps){
#   
#   if(file.exists(paste0("output/",sp,"slope_results.RData"))){
#   
#   load(paste0("output/",sp,"slope_ZIP_results.RData"))
#   
#     fyear = (min(dts$YearCollected))
#     
#   strats = unique(dts[,c("strat","Region")])
#   strats = rename(strats,s = strat)
#   
#   
#   
#   sums = data.frame(out2$summary)
#   names(sums) <- c("mean","sd","lci","lqrt","median","uqrt","uci","Rhat","n.eff","overlap0","f")
#   sums$Parameter = row.names(sums)
#   
#   # compiling indices -------------------------------------
#   
#   n_inds <- extr_inds(param = "n_s")
#   N_inds <- extr_inds(param = "N",regions = FALSE)
#  
#  
#   
#   sdnoise_st = extr_sum(param = "sdnoise",
#                     index = c("s"),
#                     log_retrans = F) 
#   sdnoise_st <- left_join(sdnoise_st,strats,by = "s")
#   
#   nu_st = extr_sum(param = "nu",
#                         index = c("s"),
#                         log_retrans = F) 
#   nu_st <- left_join(nu_st,strats,by = "s")
#   
#   
#   
#   psi_st = extr_sum(param = "psi",
#                        index = c("j","s"),
#                        log_retrans = F) 
#  
#   psiSamples <- out2$samples %>% gather_draws(psi)
#   sdnoiseSamples <- out2$samples %>% gather_draws(sdnoise[s])
#   
# # extracting the seasonal smooth ------------------------------------------
# 
#   season_sm = extr_sum(param = "vis.sm_season",
#                        index = c("day","s"),
#                        log_retrans = TRUE) 
#   
#   season_sm <- left_join(season_sm,strats)
#   pp <- ggplot()+
#     geom_line(data = season_sm,aes(x = day,y = mean,colour = Region,group = Region))+
#     geom_ribbon(data = season_sm,aes(x = day,ymax = uci,ymin = lci),alpha = 0.2)+
#     ylab("")+
#     xlab("Days since July 1")+
#     facet_wrap(facets = ~Region,ncol = 3,scales = "free")
#   
#   
#   
#   pdf(file = paste0("Figures/",sp,"_Season_slope.pdf"),
#       width = 8.5,
#       height = 11)
#   print(pp)
#   dev.off()
#   
#   # calculating trends  -----------------------------------------------------
#   
#   NSamples <- out2$samples %>% gather_draws(N[y])
#   NSamples$year <- NSamples$y + fyear-1
#   
# 
#   
#   N_compSamples <- out2$samples %>% gather_draws(N_comp[y])
#   N_compSamples$year <- N_compSamples$y + fyear-1
#   
#   n_sSamples <- out2$samples %>% gather_draws(n_s[s,y])
#   n_sSamples$year <- n_sSamples$y + fyear-1
#   n_sSamples <- left_join(n_sSamples,strats,by = "s")
#   
#  
# 
#   t_n_s <- ItoT(inds = n_sSamples,regions = TRUE)
# 
#   
#   
#   
#   t_n_sS <- ItoT_slope(inds = n_sSamples,regions = TRUE)
# 
#   
#   # t_n_s_a1 <- ItoT(inds = n_s_a1Samples,regions = TRUE,retransformation_type = "lognormal_only")
#   # 
#   # 
#   # t_n_s_a2 <- ItoT(inds = n_s_a2Samples,regions = TRUE,retransformation_type = "none")
#   # 
#   # 
#   # 
#   t_N <- ItoT(inds = NSamples,regions = FALSE)
# 
#   
#   t_n_s_15 <- ItoT(inds = n_sSamples,regions = TRUE,start= 2004)
# 
#   
#   t_n_sS_15 <- ItoT_slope(inds = n_sSamples,regions = TRUE,start= 2004)
# 
#   
#   t_N_15 <- ItoT(inds = NSamples,regions = FALSE,start= 2004)
# 
#   t_NS <- ItoT_slope(inds = NSamples,regions = FALSE)
#   t_NS_15 <- ItoT_slope(inds = NSamples,regions = FALSE,start= 2004)
#   
#   
#   
#   t_N_comp <- ItoT(inds = N_compSamples,regions = FALSE,start = 1978)
#   
#   t_N_comp_15 <- ItoT(inds = N_compSamples,regions = FALSE,start= 2004)
#   
#   t_N_compS <- ItoT_slope(inds = N_compSamples,regions = FALSE)
#   t_N_compS_15 <- ItoT_slope(inds = N_compSamples,regions = FALSE,start= 2004)
#   
#   
# 
#   trend_out <- bind_rows(t_N,
#                          t_N_15,
#                          t_NS,
#                          t_NS_15,
#                          t_N_comp,
#                          t_N_comp_15,
#                          t_N_compS,
#                          t_N_compS_15,
#                          t_n_s,
#                          t_n_sS,
#                          t_n_s_15,
#                          t_n_sS_15)
#   
#   write.csv(trend_out,file = paste0("Trends/trends_slope_",sp,".csv"),row.names = F)
#   
#   # plotting indices --------------------------------------------------------
#   
#   
#   # plot_Hyper <- plot_ind(inds = N_inds,
#   #                        #smooth_inds = ,
#   #                        raw = dts,
#   #                        add_observed = TRUE,
#   #                        add_samplesize = TRUE,
#   #                        species = sp,
#   #                        regions = FALSE,
#   #                        title_size = 20,
#   #                        axis_title_size = 18,
#   #                        axis_text_size = 16)  
#   # 
#   # pdf(file = paste0("Figures/",sp,"_Hyperparameter_slope.pdf"),
#   #     width = 8.5,
#   #     height = 11)
#   # print(plot_Hyper)
#   # dev.off()
#   
#   
#   
#   
# #   plot_by_st <- plot_ind(inds = n_inds_a2,
# #                          #smooth_inds = ,
# #                          raw = dts,
# #                          add_observed = TRUE,
# #                          add_samplesize = TRUE,
# #                          species = sp,
# #                          regions = TRUE,
# #                          title_size = 20,
# #                          axis_title_size = 18,
# #                          axis_text_size = 16)  
# #   
# #   pdf(file = paste0("Figures/",sp,"_A2_slope.pdf"),
# #       width = 8.5,
# #       height = 11)
# # print(plot_by_st)
# # dev.off()
# 
# 
# plot_by_st <- plot_ind(inds = n_inds,
#                        smooth_inds = NULL,
#                        raw = dts,
#                        add_observed = TRUE,
#                        add_samplesize = TRUE,
#                        species = sp,
#                        regions = TRUE,
#                        title_size = 20,
#                        axis_title_size = 18,
#                        axis_text_size = 16)  
# 
# pdf(file = paste0("Figures/",sp,"_A1_slope.pdf"),
#     width = 8.5,
#     height = 11)
# print(plot_by_st)
# dev.off()
# 
# 
# 
# 
# }# end if output exists
# 
# }

