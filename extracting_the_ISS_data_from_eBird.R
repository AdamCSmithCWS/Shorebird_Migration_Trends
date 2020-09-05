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
    # convert time to decimal hours since midnight
    time_observations_started = time_to_decimal(time_observations_started),
    # split date into year and day of year
    year = year(observation_date),
    day_of_year = yday(observation_date)
  ) %>% 
  filter(.,year > 1975)


hist(iss_m1$observation_count)
hist(iss_m1$effort_area_ha)
hist(iss_m1$effort_distance_km)
hist(iss_m1$duration_minutes)








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
  "TimeObservationsStarted","TimeObservationsEnded","DurationInHours","NumberOfObservers",
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










