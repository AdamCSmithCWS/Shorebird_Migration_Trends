library(tidyverse)

source("functions/utility_functions.R")
source("Functions/palettes.R")

load("data/allShorebirdPrismFallCounts.RData")
# removing Alaska, NWT, and Hawaii ----------------------------------------


ssData <- filter(ssData,StateProvince != "US-HI")
ssData <- filter(ssData,StateProvince != "US-AK")
ssData <- filter(ssData,StateProvince != "NT")


#lists for stored figures


grid_spacing <- 300000  # size of squares, in units of the CRS (i.e. meters for lae)


FYYYY = 1980


length(unique(ssData$SurveyAreaIdentifier))

length(unique(ssData$CommonName))

length(unique(ssData$SamplingEventIdentifier))

sum(ssData$ObservationCount)


# Species loop ------------------------------------------------------------
output_dir <- "g:/Shorebird_Migration_Trends/output"

full_long_data <- NULL

for(sp in sps){
  #if(sp == "Semipalmated Sandpiper"){next}
  spf = gsub(sp,pattern = " ",replacement = "_")
  spf = gsub(pattern = "\'",replacement = "",
             x = spf)
  
  load(paste0("data/data",sp,"_cmdstanr_data.RData"))
  
 full_long_data <- bind_rows(full_long_data,dts) 
  
  
}


length(unique(full_long_data$SurveyAreaIdentifier))

length(unique(full_long_data$CommonName))
length(unique(full_long_data$SiteName))

sp_sums <- full_long_data %>% 
  group_by(CommonName) %>% 
  summarise(sum_count = sum(count),
            n_sites = max(site),
            .groups = "keep") %>% 
  arrange(-sum_count)

write.csv(sp_sums,
          file = "output/Species_raw_data_summary.csv")

sum(sp_sums$sum_count)
