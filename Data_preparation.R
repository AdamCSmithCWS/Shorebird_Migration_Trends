### extracting the ISS data from eBird

library(auk)
library(tidyverse)
library(lubridate)
library(sf)


# untar("E:/eBird_all/ebd_sampling_relJun-2020.tar",exdir = "E:/eBird_all")
# untar("E:/eBird_all/ebd_relJun-2020.tar",exdir = "E:/eBird_all")
### setting path to eBird full dataset
#auk::auk_set_ebd_path("E:/eBird_all",overwrite = T)

sps1 = read.csv("data/all shorebird migration Shorebird continental annual indices 1970-2016 gam.csv",
                stringsAsFactors = FALSE)
sps = unique(sps1$species)


#ebd <- auk_ebd("E:/eBird_all/ebd_relJun-2020.txt", 
 #              file_sampling = "E:/eBird_all/ebd_sampling_relJun-2020.txt")

# 
# 
# ebd_filters <- ebd %>% 
#    auk_species(species = sps) %>% 
#   # # southeastern coastal plain bcr
#   auk_country(country = c("US","CA")) %>% 
#   # fall surveys only
#   auk_date(date = c("*-07-01", "*-11-30")) %>% 
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

iss_samp = read_sampling(f_sampling)






# fitting to strata - geographic overlay ----------------------------------
map <- read_sf(dsn = "data",
               layer ="region_polygons")

st_crs(map)

iss_sites = unique(iss_samp[,c("locality_id",
                        "latitude",
                        "longitude")])




iss_sites = st_as_sf(iss_sites,coords = c("longitude","latitude"), crs = 4326)


iss_sites_regs <- st_join(iss_sites, map, join = st_nearest_feature)

iss_samp <- left_join(iss_samp,iss_sites_regs[,c("locality_id","Region","Region_FR")])


# Zero-fill manual --------------------------------------------------------
iss_s2 <- expand_grid(iss_samp,sps) #making complete list of checklists and species
iss_s2 <- mutate(iss_s2,common_name = sps)

iss_m <- full_join(iss[,c("checklist_id","common_name","observation_count")],iss_s2) #merging species observations

iss_m[which(is.na(iss_m$observation_count)),"observation_count"] <- "0"


# function to convert time observation to hours since midnight
time_to_decimal <- function(x) {
  x <- hms(x, quiet = TRUE)
  hour(x) + minute(x) / 60 + second(x) / 3600
}

# clean up variables and further filtering to just US data and just years 1976 - 2019
iss_m1 <- iss_m %>% 
  mutate(
    # convert X to NA
    observation_count = if_else(observation_count == "X", 
                                NA_character_, observation_count),
    observation_count = as.integer(observation_count),
    # convert time to decimal hours since midnight
    time_observations_started = time_to_decimal(time_observations_started),
    # split date into year and day of year
    year = year(observation_date),
    day_of_year = yday(observation_date),
  ) %>% 
  filter(.,year > 1973 & year < 2020,
         country == "United States")



### effort info in ISS is only common in the last ~10 years.
wEMin = which(!is.na(iss_m1$duration_minutes))
wEArea = which(!is.na(iss_m1$effort_area_ha))
wEdist = which(!is.na(iss_m1$effort_distance_km))

wEffort = unique(c(wEMin,wEArea,wEdist))
eff_y = table(iss_m1[wEffort,"year"])
all_y = table(iss_m1$year)
plot(eff_y/all_y)
eff_y/all_y
# 1976        1977        1978        1979        1980        1981        1982        1983        1984        1985        1986        1987        1988        1989        1990        1991 
# 0.011201629 0.014334862 0.008196721 0.015767132 0.019313305 0.022562450 0.015394089 0.029808774 0.019551050 0.026139410 0.031914894 0.047206924 0.033912324 0.032198712 0.031645570 0.025196850 
# 1992        1993        1994        1995        1996        1997        1998        1999        2000        2001        2002        2003        2004        2005        2006        2007 
# 0.033626902 0.018443804 0.018308081 0.020177563 0.022222222 0.020176545 0.021428571 0.015455066 0.018724400 0.015025042 0.038336933 0.132538105 0.108346293 0.095898004 0.608996540 0.696969697 
# 2008        2009        2010        2011        2012        2013        2014        2015        2016        2017        2018        2019 
# 0.799697657 0.805788982 0.802350427 0.964864865 0.979458450 0.945194599 0.977719528 0.990298507 0.993333333 0.982922201 1.000000000 1.000000000
all_y  # what happened in 2006?
# 1976  1977  1978  1979  1980  1981  1982  1983  1984  1985  1986  1987  1988  1989  1990  1991  1992  1993  1994  1995  1996  1997  1998  1999  2000  2001  2002  2003  2004  2005  2006  2007 
# 27496 48832 47824 46172 39144 34748 45472 49784 38668 41776 42112 35588 33852 30436 35392 35560 34972 48580 44352 34692 35280 44408 43120 48916 47852 50316 51856 42252 54012 50512  8092 29568 
# 2008  2009  2010  2011  2012  2013  2014  2015  2016  2017  2018  2019 
# 37044 29988 26208 20720 29988 35252 42728 37520 29400 29512 21840 29092

#iss_m1$locality_id #this is the unique site identifier

iss_m1 <- iss_m1[which(iss_m1$day_of_year < 334),] ### dropping the observations outside of the July-November window

# ISS data to bind --------------------------------------------------------

iss_full <- iss_m1[,c("checklist_id",
                   "common_name",
                   "observation_count",
                   "country",
                   "state_code",
                   "locality_id",
                   "latitude",
                   "longitude",
                   "Region",
                   "year",
                   "day_of_year")]

#renaming to match AKN headers
iss_full <- rename(iss_full,
                   SamplingEventIdentifier = checklist_id,
                   CommonName = common_name,
                   ObservationCount = observation_count,
                   Country = country,
                   StateProvince = state_code,
                   SurveyAreaIdentifier = locality_id,
                   DecimalLatitude = latitude,
                   DecimalLongitude = longitude,
                   YearCollected = year,
                   doy = day_of_year)

# Nature Counts data ------------------------------------------------------




# OSS data ---------------------------------------------------------------

# oss <- read.delim("data/OSS_8July2020.txt",stringsAsFactors = F,nrows = 48006)
## not sure why, but the original file won't load properly - 
# Warning message:
#   In scan(file = file, what = what, sep = sep, quote = quote, dec = dec,  :
#             EOF within quoted string
# after warning, it only loads ~6000 rows and even those rows are not properly loaded
# Solution was to open it in excel (tab delimited), save as csv, then import csv below        
oss <- read.csv("data/OSS_8July2020.csv",stringsAsFactors = F)



# ACSS data ---------------------------------------------------------------

acss <- read.delim("data/ACSS_8July2020.txt",stringsAsFactors = F)





# Filtering columns -------------------------------------------------------


#names(oss)
#c("GlobalUniqueIdentifier","DateLastModified","BasisOfRecord","InstitutionCode","CollectionCode","CatalogNumber","ScientificName","HigherTaxon","Kingdom","Phylum","Class","Order","Family","Genus","SpecificEpithet","InfraspecificRank","InfraspecificEpithet","ScientificNameAuthor","IdentificationQualifier","HigherGeography","Continent","WaterBody","IslandGroup","Island","Country","StateProvince","County","Locality","MinimumElevationInMeters","MaximumElevationInMeters","MinimumDepthInMeters","MaximumDepthInMeters","DecimalLatitude","DecimalLongitude","GeodeticDatum","CoordinateUncertaintyInMeters","YearCollected","MonthCollected","DayCollected","TimeCollected","JulianDay","Collector","Sex","LifeStage","ImageURL","RelatedInformation","CollectorNumber","FieldNumber","FieldNotes","OriginalCoordinatesSystem","LatLongComments","GeoreferenceMethod","GeoreferenceReferences","GeoreferenceVerificationStatus","Remarks","FootprintWKT","FootprintSRS","ProjectCode","ProtocolType","ProtocolCode","ProtocolSpeciesTargeted","ProtocolReference","ProtocolURL","SurveyAreaIdentifier","SurveyAreaSize","SurveyAreaPercentageCovered","SurveyAreaShape","SurveyAreaLongAxisLength","SurveyAreaShortAxisLength","SurveyAreaLongAxisOrientation","CoordinatesScope","SamplingEventIdentifier","SamplingEventStructure","RouteIdentifier","TimeObservationsStarted","TimeObservationsEnded","DurationInHours","TimeIntervalStarted","TimeIntervalEnded","TimeIntervalsAdditive","NumberOfObservers","EffortMeasurement1","EffortUnits1","EffortMeasurement2","EffortUnits2","EffortMeasurement3","EffortUnits3","EffortMeasurement4","EffortUnits4","EffortMeasurement5","EffortUnits5","EffortMeasurement6","EffortUnits6","EffortMeasurement7","EffortUnits7","EffortMeasurement8","EffortUnits8","EffortMeasurement9","EffortUnits9","EffortMeasurement10","EffortUnits10","EffortMeasurement11","EffortUnits11","EffortMeasurement12","EffortUnits12","EffortMeasurement13","EffortUnits13","EffortMeasurement14","EffortUnits14","EffortMeasurement15","EffortUnits15","EffortMeasurement16","EffortUnits16","EffortMeasurement17","EffortUnits17","EffortMeasurement18","EffortUnits18","NoObservations","DistanceFromObserver","DistanceFromObserverMin","DistanceFromObserverMax","DistanceFromStart","BearingInDegrees","SpecimenDecimalLatitude","SpecimenDecimalLongitude","SpecimenGeodeticDatum","SpecimenUTMZone","SpecimenUTMNorthing","SpecimenUTMEasting","ObservationCount","ObservationDescriptor","ObservationCount2","ObservationDescriptor2","ObservationCount3","ObservationDescriptor3","ObservationCount4","ObservationDescriptor4","ObservationCount5","ObservationDescriptor5","ObservationCount6","ObservationDescriptor6","ObsCountAtLeast","ObsCountAtMost","ObservationDate","DateUncertaintyInDays","AllIndividualsReported","AllSpeciesReported","UTMZone","UTMNorthing","UTMEasting","CoordinatesUncertaintyInDecimalDegrees","CommonName","RecordPermissions","MultiScientificName1","MultiScientificName2","MultiScientificName3","MultiScientificName4","MultiScientificName5","MultiScientificName6","TaxonomicAuthorityAuthors","TaxonomicAuthorityVersion","TaxonomicAuthorityYear","SpeciesCode","TaxonConceptID","BreedingBirdAtlasCode","HabitatDescription","Remarks2","LastModifiedAction","RecordReviewStatus")
nc_cols <- c("CatalogNumber","ScientificName","Genus","SpecificEpithet",
  "Country","StateProvince","County","Locality","DecimalLatitude","DecimalLongitude","GeodeticDatum",
  "YearCollected","MonthCollected","DayCollected","JulianDay",
  "CollectorNumber",
  "ProjectCode","ProtocolType","SurveyAreaIdentifier",
  "SurveyAreaShape",
  "SamplingEventIdentifier",
  "DurationInHours","NumberOfObservers",
  "EffortMeasurement1","EffortUnits1",
  "EffortMeasurement2","EffortUnits2",
  "NoObservations",
  "ObservationCount","ObservationDescriptor",
  "ObservationDate",
  "AllIndividualsReported","AllSpeciesReported",
  "UTMZone","UTMNorthing","UTMEasting",
  "CommonName","TaxonomicAuthorityAuthors","TaxonomicAuthorityVersion","TaxonomicAuthorityYear","SpeciesCode")


oss <- oss[,nc_cols]
acss <- acss[,nc_cols]

oss[,"SurveyAreaIdentifier"] <- as.character(oss[,"SurveyAreaIdentifier"])
acss[,"SurveyAreaIdentifier"] <- as.character(acss[,"SurveyAreaIdentifier"])

css <- bind_rows(acss,oss)
# css$SurveyAreaIdentifier <- paste(css$ProjectCode,css$SurveyAreaIdentifier,sep = "_")
#above is unecessary column is unique across both datasets

## drops or fixes the observations with no valid dates and outside of the fall migration period
css <- css[which(!is.na(css$MonthCollected) & css$MonthCollected > 6 & css$MonthCollected < 12),] 
css <- css[which(css$YearCollected > 1973),] 

css[which(is.na(css$DayCollected)),"DayCollected"] <- 15 #this places the observation in teh middle of the month

css$ObservationDate <- lubridate::ymd(paste(css$YearCollected,css$MonthCollected,css$DayCollected,sep = "/"))
css$doy <- lubridate::yday(css$ObservationDate)


## drops all sites with no coordinates (there are ~2000 observations and this is about 1% of the total ACSS and OSS combined) most have no information on province or county as well
css <- css[which(!is.na(css$DecimalLatitude)),]


# identifying unique sites and sampling events ------------------------------------------------

css_samp <- unique(css[,c("SamplingEventIdentifier",
                          "Country",
                          "StateProvince",
                          "SurveyAreaIdentifier",
                          "DecimalLatitude",
                          "DecimalLongitude",
                          "YearCollected",
                          "doy")])


# overlaying with regions -------------------------------------------------

css_sites = unique(css_samp[,c("SurveyAreaIdentifier",
                               "DecimalLatitude",
                               "DecimalLongitude")])

css_sites = st_as_sf(css_sites,coords = c("DecimalLongitude","DecimalLatitude"), crs = 4326)


css_sites_regs <- st_join(css_sites, map, join = st_nearest_feature)

css_samp <- left_join(css_samp,data.frame(css_sites_regs[,c("SurveyAreaIdentifier","Region")]))




#full sampling event by species matrix
css_full <- expand_grid(css_samp,sps)
css_full <- rename(css_full,
                   CommonName = sps)


#confirming common species names are the same
css_sps <- unique(css$CommonName)
#sps[-which(sps %in% css_sps)]
#character(0)


# dropping species not included, extra columns, and pres/abs --------------------------------
css <- css[which(css$CommonName %in% sps &
                   css$ObservationDescriptor != "Presence/Absence"),c("SamplingEventIdentifier",
                                                                 "CommonName",
                                            "ObservationCount")]


# zero fill css data  -----------------------------------------------------


css_full <- left_join(css_full,css,by = c("SamplingEventIdentifier","CommonName"))
css_full[which(is.na(css_full$ObservationCount)),"ObservationCount"] <- 0

css_full <- select(css_full,-geometry)
css_full$Country <- "Canada" #some country values are missing from the ACSS file





# Combining ISS and CSS data ----------------------------------------------
# includes observations, zero-filled, and all sampling events




ssData = bind_rows(css_full,iss_full)



# adding national divisions to regions ------------------------------------


ssData[which(ssData$Region == "East Inland" &
               ssData$StateProvince == "ON"),"Region"] <- "Ontario"
ssData[which(ssData$Region == "Northeast Coastal" &
               ssData$Country == "Canada"),"Region"] <- "Atlantic Canada"
ssData[which(ssData$Region == "Northeast Coastal" &
               ssData$Country == "United States"),"Region"] <- "Northeast US Coastal"




save(list = c("ssData",
              "sps",
              "map"),file = "data/allShorebirdPrismFallCounts.RData")











