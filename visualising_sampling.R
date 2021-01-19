library(tidyverse)
library(sf)
library(spdep)
library(ggforce)
library(gganimate)
library(magick)
library(transformr)
library(gifski)
library(patchwork)


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



# animated gif ------------------------------------------------------------


# visualise the sampling over time
# 

nevents_by_site_yr <- ssData %>% distinct(SurveyAreaIdentifier,YearCollected,doy) %>% 
  group_by(SurveyAreaIdentifier,YearCollected) %>% 
  summarise(nobs = n(),
            log_nobs = log(n(),10))

map_events <- inner_join(iss_sites_lcc,nevents_by_site_yr,by = c("SurveyAreaIdentifier"))
map_events$YearCollected <- as.integer(map_events$YearCollected)

mp.plot = ggplot2::ggplot()+
  ggplot2::geom_sf(data = map_events,ggplot2::aes(fill = log_nobs,colour = log_nobs,group = YearCollected),size = 2)+
  ggplot2::theme_minimal()+
  ggplot2::ylab("")+
  ggplot2::xlab("")+
  #ggplot2::labs(title = ptit)+
  # ggplot2::theme(legend.position = "right", line = ggplot2::element_line(size = 0.4),
  #                rect = ggplot2::element_rect(size = 0.1),
  #                axis.text = ggplot2::element_blank(),
  #                axis.line = ggplot2::element_blank())+
  ggplot2::scale_colour_viridis_c(aesthetics = c("fill","colour"),
                                  direction = -1,
                                guide = ggplot2::guide_legend(reverse=TRUE),
                                name = paste0("Number Surveys year (log_10)"))+
  gganimate::transition_time(YearCollected)+
  ggplot2::labs(title = paste("log_n"," : {frame_time}"))



panim = animate(mp.plot, nframes = 200, fps = 4, end_pause = 15, rewind = FALSE,height = 800,width = 800)#,renderer = magick_renderer())

anim_save(filename = paste0("animated_sampling_map_noNA.gif"),animation = panim,path = "Figures")




sites_yr = expand(ssData,SurveyAreaIdentifier,YearCollected)
iss_sites_lcc_yr = right_join(iss_sites_lcc,sites_yr,by = "SurveyAreaIdentifier")



nevents_by_site_yr <- ssData %>% distinct(SurveyAreaIdentifier,YearCollected,doy) %>% 
  group_by(SurveyAreaIdentifier,YearCollected) %>% 
  summarise(nobs = n(),
            log_nobs = log(n(),10))

map_events <- left_join(iss_sites_lcc_yr,nevents_by_site_yr,by = c("SurveyAreaIdentifier","YearCollected"))
map_events$YearCollected <- as.integer(map_events$YearCollected)


mp.plot = ggplot2::ggplot()+
  ggplot2::geom_sf(data = map_events,ggplot2::aes(fill = log_nobs,colour = log_nobs,group = YearCollected),size = 2)+
  ggplot2::theme_minimal()+
  ggplot2::ylab("")+
  ggplot2::xlab("")+
  #ggplot2::labs(title = ptit)+
  # ggplot2::theme(legend.position = "right", line = ggplot2::element_line(size = 0.4),
  #                rect = ggplot2::element_rect(size = 0.1),
  #                axis.text = ggplot2::element_blank(),
  #                axis.line = ggplot2::element_blank())+
  ggplot2::scale_colour_viridis_c(aesthetics = c("fill","colour"))+#,
  # guide = ggplot2::guide_legend(reverse=TRUE),
  # name = paste0("Number of sites by year"))+
  gganimate::transition_time(YearCollected)+
  ggplot2::labs(title = paste("log_n"," : {frame_time}"))



panim = animate(mp.plot, nframes = 100, fps = 7, end_pause = 15, rewind = FALSE,height = 800,width = 800)#,renderer = magick_renderer())

anim_save(filename = paste0("animated_sampling_map.gif"),animation = panim,path = "Figures")



# static plot -------------------------------------------------------------





nyears_by_site <- ssData %>%  
  distinct(SurveyAreaIdentifier,YearCollected) %>% 
  mutate(decade = 10*floor(YearCollected/10)) %>% 
  distinct(SurveyAreaIdentifier,decade) %>%
  group_by(SurveyAreaIdentifier) %>% 
  summarise(nobs = n(),
            nobs_alpha = n()/5) 

map_events <- inner_join(iss_sites_lcc,nyears_by_site,by = c("SurveyAreaIdentifier"))


mp.plot = ggplot2::ggplot()+
  geom_sf(data = prov_state)+
  ggplot2::geom_sf(data = map_events,ggplot2::aes(fill = nobs,colour = nobs,size = nobs,alpha = nobs_alpha))+
  ggplot2::theme_minimal()+
  ggplot2::ylab("")+
  ggplot2::xlab("")+
  ggplot2::labs(title = "Number of decades with at least one observation")+
  ggplot2::theme(legend.position = "right", line = ggplot2::element_line(size = 0.4),
                 rect = ggplot2::element_rect(size = 0.1),
                 axis.text = ggplot2::element_blank(),
                 axis.line = ggplot2::element_blank())+
  ggplot2::scale_colour_viridis_c(aesthetics = c("fill","colour"),
                                  direction = -1,
   guide = ggplot2::guide_legend(reverse=TRUE),
   name = paste0("Number of decades"))+
  guides(alpha = "none",size = "none")


# first 2 decades ---------------------------------------------------------



nyears_by_site <- ssData %>% 
  filter(YearCollected < 1986) %>% 
  distinct(SurveyAreaIdentifier,YearCollected) %>% 
  mutate(decade = 1970) %>% 
  distinct(SurveyAreaIdentifier,decade) %>%
  group_by(SurveyAreaIdentifier) %>% 
  summarise(nobs = n(),
            nobs_alpha = n()/2) 

map_events <- inner_join(iss_sites_lcc,nyears_by_site,by = c("SurveyAreaIdentifier"))


mp.plot2 = ggplot2::ggplot()+
  geom_sf(data = prov_state)+
  ggplot2::geom_sf(data = map_events,ggplot2::aes(fill = nobs,colour = nobs,alpha = nobs_alpha))+
  ggplot2::theme_minimal()+
  ggplot2::ylab("")+
  ggplot2::xlab("")+
  ggplot2::labs(title = "Sites with observations 1974 - 1985")+
  ggplot2::theme(legend.position = "none", line = ggplot2::element_line(size = 0.4),
                 rect = ggplot2::element_rect(size = 0.1),
                 axis.text = ggplot2::element_blank(),
                 axis.line = ggplot2::element_blank())+
  ggplot2::scale_colour_viridis_c(aesthetics = c("fill","colour"),
                                  direction = -1,
                                  guide = ggplot2::guide_legend(reverse=TRUE),
                                  name = paste0("Number of decades"))+
  guides(alpha = "none")


 print(mp.plot2)


# last decade -------------------------------------------------------------


 
 nyears_by_site <- ssData %>% 
   filter(YearCollected > 2010) %>% 
   distinct(SurveyAreaIdentifier,YearCollected) %>% 
   mutate(decade = 10*floor(YearCollected/10)) %>% 
   distinct(SurveyAreaIdentifier,decade) %>%
   group_by(SurveyAreaIdentifier) %>% 
   summarise(nobs = n(),
             nobs_alpha = n()/2) 
 
 map_events <- inner_join(iss_sites_lcc,nyears_by_site,by = c("SurveyAreaIdentifier"))
 
 
 mp.plot3 = ggplot2::ggplot()+
   geom_sf(data = prov_state)+
   ggplot2::geom_sf(data = map_events,ggplot2::aes(fill = nobs,colour = nobs,alpha = nobs_alpha))+
   ggplot2::theme_minimal()+
   ggplot2::ylab("")+
   ggplot2::xlab("")+
   ggplot2::labs(title = "Sites with observations post 2010")+
   ggplot2::theme(legend.position = "none", line = ggplot2::element_line(size = 0.4),
                  rect = ggplot2::element_rect(size = 0.1),
                  axis.text = ggplot2::element_blank(),
                  axis.line = ggplot2::element_blank())+
   ggplot2::scale_colour_viridis_c(aesthetics = c("fill","colour"),
                                   direction = -1,
                                   guide = ggplot2::guide_legend(reverse=TRUE),
                                   name = paste0("Number of decades"))+
   guides(alpha = "none")
 
 
 print(mp.plot3)
 
 
 
 
 # first and last decades -------------------------------------------------------------
 
 
 
 nyears_by_site <- ssData %>% 
   filter(YearCollected > 2010 | YearCollected < 1986) %>% 
   distinct(SurveyAreaIdentifier,YearCollected) %>% 
   mutate(decade = ifelse(YearCollected > 2000,1970,2019)) %>% 
   distinct(SurveyAreaIdentifier,decade) %>%
   group_by(SurveyAreaIdentifier) %>% 
   summarise(nobs = n()) %>% 
   filter(nobs > 1)
 
 map_events <- inner_join(iss_sites_lcc,nyears_by_site,by = c("SurveyAreaIdentifier"))
 
 
 mp.plot4 = ggplot2::ggplot()+
   geom_sf(data = prov_state)+
   ggplot2::geom_sf(data = map_events,ggplot2::aes(fill = nobs,colour = nobs))+
   ggplot2::theme_minimal()+
   ggplot2::ylab("")+
   ggplot2::xlab("")+
   ggplot2::labs(title = "Sites with observations in 1974-1985 and 2010s")+
   ggplot2::theme(legend.position = "none", line = ggplot2::element_line(size = 0.4),
                  rect = ggplot2::element_rect(size = 0.1),
                  axis.text = ggplot2::element_blank(),
                  axis.line = ggplot2::element_blank())+
   ggplot2::scale_colour_viridis_c(aesthetics = c("fill","colour"),
                                   direction = -1,
                                   guide = ggplot2::guide_legend(reverse=TRUE),
                                   name = paste0("Number of decades"))+
   guides(alpha = "none",size = "none")
 
 
 print(mp.plot4)
 
 
 
 nyears_by_site <- ssData %>% 
   filter(YearCollected > 2010 | YearCollected < 1990) %>% 
   distinct(SurveyAreaIdentifier,YearCollected) %>% 
   mutate(decade = ifelse(YearCollected > 2000,1970,2019)) %>% 
   distinct(SurveyAreaIdentifier,decade) %>%
   group_by(SurveyAreaIdentifier) %>% 
   summarise(nobs = n()) %>% 
   filter(nobs > 1)
 
 map_events <- inner_join(iss_sites_lcc,nyears_by_site,by = c("SurveyAreaIdentifier"))
 
 
 mp.plot5 = ggplot2::ggplot()+
   geom_sf(data = prov_state)+
   ggplot2::geom_sf(data = map_events,ggplot2::aes(fill = nobs,colour = nobs))+
   ggplot2::theme_minimal()+
   ggplot2::ylab("")+
   ggplot2::xlab("")+
   ggplot2::labs(title = "Sites with observations in 1974-1989 and 2010s")+
   ggplot2::theme(legend.position = "none", line = ggplot2::element_line(size = 0.4),
                  rect = ggplot2::element_rect(size = 0.1),
                  axis.text = ggplot2::element_blank(),
                  axis.line = ggplot2::element_blank())+
   ggplot2::scale_colour_viridis_c(aesthetics = c("fill","colour"),
                                   direction = -1,
                                   guide = ggplot2::guide_legend(reverse=TRUE),
                                   name = paste0("Number of decades"))+
   guides(alpha = "none",size = "none")
 
 

pdf("Figures/Number of decades at sites.pdf",
    width = 11,
    height = 8.5)
print(mp.plot)
print(mp.plot2)
print(mp.plot3)
print(mp.plot4)
print(mp.plot5)
dev.off()


