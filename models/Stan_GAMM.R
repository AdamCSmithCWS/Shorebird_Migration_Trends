#### Stan GAMM

#install.packages("rstanarm")

library(rstan)
library(rstanarm)
options(mc.cores = 3)


load("data/allShorebirdPrismFallCounts.RData")

source("functions/GAM_basis_function.R")

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




test2 <- stan_gamm4(formula = count ~ s(yr),
                   random = ~(1 | site),
                   data = dts,
                   family = neg_binomial_2,
                   #QR = TRUE, # only if multiple predictors
                   cores = 3,
                   chains = 3,
                   adapt_delta = 0.99,
                   seed = 12345)


test <- stan_gamm4(formula = count ~ s(yr),
                   random = ~(1 | site),
                   data = dtst,
                   family = neg_binomial_2,
                   #QR = TRUE, # only if multiple predictors
                   cores = 3,
                   chains = 3,
                   adapt_delta = 0.99,
                   seed = 12345)

launch_shinystan(test)
rstan::get_stanmodel(test$stanfit)




















# GAM seasonal basis function ---------------------------------------------

nKnots_season = 5
basis_season <- gam.basis.func(orig.preds = as.integer(unlist(dts[,"date"])),
                               nknots = nKnots_season,
                               standardize = "z",
                               random = F,
                               npredpoints = max(dts$date),
                               even_gaps = FALSE,
                               sm_name = "season")




# GAM annual basis function ---------------------------------------------

nKnots_year = 13
basis_year <- gam.basis.func(orig.preds = as.integer(unlist(dts[,"yr"])),
                             nknots = nKnots_year,
                             standardize = "z",
                             random = F,
                             npredpoints = max(dts$yr),
                             even_gaps = FALSE,
                             sm_name = "year")





jags_data <- list(count = as.integer(unlist(dts$count)),
                  yr = as.integer(unlist(dts$yr)),
                  site = as.integer(unlist(dts$site)),
                  strat = as.integer(unlist(dts$strat)),
                  
                  nyears = nyears,
                  nstrata = nstrata,
                  nsites = nsites,
                  ncounts = ncounts,
                  
                  season_basis = basis_season$season_basis,
                  season_basispred = basis_season$season_basispred,
                  nknots_season = basis_season$nknots_season,
                  
                  year_basispred = basis_year$year_basispred,
                  nknots_year = basis_year$nknots_year)








