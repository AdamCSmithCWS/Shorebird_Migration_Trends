### extracting the ISS data from eBird

library(auk)
library(tidyverse)
library(lubridate)

# untar("E:/eBird_all/ebd_sampling_relJun-2020.tar",exdir = "E:/eBird_all")
# untar("E:/eBird_all/ebd_relJun-2020.tar",exdir = "E:/eBird_all")
### setting path to eBird full dataset
#auk::auk_set_ebd_path("E:/eBird_all",overwrite = T)

sps1 = read.csv("data/all shorebird migration Shorebird continental annual indices 1970-2016 gam.csv",
                stringsAsFactors = FALSE)
sps = unique(sps1$species)


ebd <- auk_ebd("E:/eBird_all/ebd_relJun-2020.txt", 
               file_sampling = "E:/eBird_all/ebd_sampling_relJun-2020.txt")

# 
# 
# ebd_filters <- ebd %>% 
#    auk_species(species = sps) %>% 
#   # # southeastern coastal plain bcr
#   auk_country(country = c("US","CA")) %>% 
#   # fall surveys only
#   auk_date(date = c("*-07-01", "*-12-31")) %>% 
#   # restrict to the ISS protocol
#   auk_protocol(protocol = c("International Shorebird Survey (ISS)")) 
# 


# post filtering ----------------------------------------------------------


data_dir <- "data"
if (!dir.exists(data_dir)) {
  dir.create(data_dir)
}
f_ebd <- file.path(data_dir, "ebd_ISS.txt")
f_sampling <- file.path(data_dir, "ebd_ISS_sampling.txt")

# only run if the files don't already exist
if (!file.exists(f_ebd)) {
  auk_filter(ebd_filters, file = f_ebd, file_sampling = f_sampling)
}

iss = read_ebd(f_ebd)

iss_s = read_sampling(f_sampling)


# Zero-fill manual --------------------------------------------------------
iss_s2 <- expand_grid(iss_s,sps) #making complete list of checklists and species
iss_s2 <- mutate(iss_s2,common_name = sps)

iss_m <- full_join(iss[,c("checklist_id","common_name","observation_count")],iss_s2) #merging species observations

iss_m[which(is.na(iss_m$observation_count)),"observation_count"] <- "0"

# function to convert time observation to hours since midnight
time_to_decimal <- function(x) {
  x <- hms(x, quiet = TRUE)
  hour(x) + minute(x) / 60 + second(x) / 3600
}

# clean up variables
iss_m1 <- iss_m %>% 
  mutate(
    # convert X to NA
    observation_count = if_else(observation_count == "X", 
                                NA_character_, observation_count),
    observation_count = as.integer(observation_count),
    # effort_distance_km to 0 for non-travelling counts
    effort_distance_km = if_else(protocol_type != "Traveling", 
                                 0, effort_distance_km),
    # convert time to decimal hours since midnight
    time_observations_started = time_to_decimal(time_observations_started),
    # split date into year and day of year
    year = year(observation_date),
    day_of_year = yday(observation_date)
  )




