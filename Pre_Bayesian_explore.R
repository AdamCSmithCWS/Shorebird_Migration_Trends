# SPECIES MCMC data-prep -------------------------------------------------------

library(jagsUI)
library(ggmcmc)
library(gamm4)

fyear <- 1978 # first year in which there are sufficient data in all regions


 load("data/allShorebirdPrismFallCounts.RData")
# source("functions/GAM_basis_function.R")

# n_cores <- 6
# cluster <- makeCluster(n_cores, type = "PSOCK")
# registerDoParallel(cluster)
# 
# 
# 
# fullrun <- foreach(sp = sps[c(3,4,6,8,21,22)],
#                    .packages = c("jagsUI","tidyverse","ggmcmc"),
#                    .inorder = FALSE,
#                    .errorhandling = "pass") %dopar%
#   {
    
    #load("data/allShorebirdPrismFallCounts.RData")
    
    sp = sps[25]

dts <- filter(ssData,CommonName == sp,
              YearCollected > fyear-1)
dts$present <- FALSE
dts[which(dts$ObservationCount > 0),"present"] <- TRUE

#number of non-zero observations by region
nobs_reg <- dts %>% 
  filter(present == TRUE) %>% 
  group_by(Region) %>% 
  summarise(nobs = n())

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
sites_keep <- nyrs_site[which(nyrs_site$span_years > 5),"SurveyAreaIdentifier"]

dts <- filter(dts,SurveyAreaIdentifier %in% sites_keep$SurveyAreaIdentifier) 


fday = min(dts$doy)-1


dts <- dts %>% mutate(count = as.integer(ObservationCount),
                      year = as.integer(YearCollected),
                      yr = as.integer(year-(fyear-1)),
                      strat = as.integer(factor(Region)),
                      date = doy-fday) 

nstrata = max(dts$strat)

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

dtst <- filter(dts,Region == "Atlantic Canada")

dts$Region <- factor(dts$Region)
dts$site <- factor(dts$SurveyAreaIdentifier)


M <- gam(count ~ yr +s(yr,Region,bs = "re") + Region + s(date, k = 7,m = 2) + s(date,Region, k = 7,m = 2, bs = "fs") + s(site,Region, bs = "re"),
           data = dts,
         family = "poisson")
           #family = negbin(theta = 1))

plot.gam(M)





#  }#end species loops


#stopCluster(cl = cluster)






