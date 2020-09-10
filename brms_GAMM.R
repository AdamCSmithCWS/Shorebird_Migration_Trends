#### Stan GAMM

library(brms)
library(tidyverse)
library(rstan)
options(mc.cores = 3)


load("data/allShorebirdPrismFallCounts.RData")

#source("functions/GAM_basis_function.R")

#for(sp in sps){
sp = sps[25]

dts <- filter(ssData,CommonName == sp)
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
                      yr = as.integer(year-1973),
                      strat = as.integer(factor(Region)),
                      date = doy-fday) 

nstrata = max(dts$strat)

sByReg = unique(dts[,c("SurveyAreaIdentifier","strat")])
sByReg <- arrange(sByReg,strat,SurveyAreaIdentifier)

nsites <- table(sByReg$strat)

site <- NULL
for(s in 1:nstrata){
  site <- c(site,1:nsites[s]) 
}
sByReg$site <- site
sByReg <- select(sByReg,-strat)
dts <- left_join(dts,sByReg,by = "SurveyAreaIdentifier")

dtst = filter(dts,strat == 1)

nyears = max(dts$yr)
# brms GAMM definition year----------------------------------------------------

dts$strat <- factor(dts$Region)

# mod = make_stancode(formula = count ~ s(yr, k = 13,id = 1) + s(yr,by = strat, k = 13,id = 1) + s(date,k = 7),
#                     random = ~(1 | strat/site),
#                     data = dts,
#                     family = poisson)
# cat(mod,file = "models/AnnualSeasonalGAMM.stan")

dat = make_standata(formula = count ~ s(yr, k = 13,id = 1) + s(yr,by = strat, k = 13,id = 1) + s(date,k = 7),
                    random = ~(1 | strat/site),
                    data = dts,
                    family = poisson)

M = stan(file = "models/AnnualSeasonalGAMM_OverD.stan",
         data = dat,
         chains = 1,
         pars = c("sd_noise",
                  "s_1_1"))









# brms GAMM definition year----------------------------------------------------

yr_dat <- data.frame(yr = 1:nyears,
                     count = rpois(nyears,3))
mod_yr = make_stancode(formula = count ~ s(yr, k = 13),
                    data = yr_dat,
                    family = poisson)

yr_b_dat <- make_standata(formula = count ~ s(yr, k = 13),
                          data = yr_dat,
                          family = poisson)

basis_yr <- yr_b_dat$Zs_1_1
yr_pred = data.frame(yr = 1:nyears,
                     yr_z = yr_b_dat$Xs[1:nyears,1]) # to be joined to original data frame


ndoy = max(dts$date)

date_dat = data.frame(date = 1:ndoy,
                      count = rpois(ndoy,3))

date_b_dat <- make_standata(formula = count ~ s(date, k = 7),
                          data = date_dat,
                          family = poisson)

basis_date <- date_b_dat$Zs_1_1

date_pred = data.frame(date = 1:ndoy,
                     date_z = date_b_dat$Xs[1:ndoy,1]) # to be joined to original data frame











