# SPECIES MCMC data-prep -------------------------------------------------------

library(jagsUI)
library(ggmcmc)

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



ncounts = nrow(dts)

dmin_pred <- tapply(dts$date,dts$strat,min)
dmax_pred <- tapply(dts$date,dts$strat,max)

dmax = max(dts$date)
nyears = max(dts$yr)

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



mod.file = "models/AnnualSeasonalGAM.R"




parms = c("sdnoise",
          # "nu", #if optional heavy-tailed noise
          "sdgam_season",
          "sdgam_year",
          "sdgam_year_b",
          "beta_season",
          "b_year",
          "B_year",
          "beta_year",
          "sdsite",
          "N",
          "n_s",
          "alpha")


#adaptSteps = 200              # Number of steps to "tune" the samplers.
burnInSteps = 5000            # Number of steps to "burn-in" the samplers.
nChains = 1                   # Number of chains to run.
numSavedSteps=1000          # Total number of steps in each chain to save.
thinSteps=10                   # Number of steps to "thin" (1=keep every step).
nIter = ceiling( ( (numSavedSteps * thinSteps )+burnInSteps)) # Steps per chain.

t1 = Sys.time()




# MCMC sampling -----------------------------------------------------------



out2 = jagsUI(data = jags_data,
              parameters.to.save = parms,
              n.chains = 3,
              n.burnin = burnInSteps,
              n.thin = thinSteps,
              n.iter = nIter,
              parallel = T,
              modules = NULL,
              model.file = mod.file)


t2 = Sys.time()

out2$n.eff
out2$rhat

gg = ggs(out2$samples)

ggy = ggs(out2$samples,family = "B_year")
bby2 = ggs(out2$samples,family = "beta_year")
gga = ggs(out2$samples,family = "alpha")
ggall = rbind(ggy,gga)
ggmcmc(ggall,file = paste0("output/mcmc_",sp,".pdf"))
ggmcmc(bby2,file = paste0("output/mcmc_bby2_",sp,".pdf"))

ggsd = ggs(out2$samples,family = "sd")
ggmcmc(ggsd,file = paste0("output/mcmc_ggsd_",sp,".pdf"))

