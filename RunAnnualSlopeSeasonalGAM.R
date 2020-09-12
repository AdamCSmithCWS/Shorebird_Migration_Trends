# SPECIES MCMC data-prep -------------------------------------------------------

library(jagsUI)
library(ggmcmc)


library(doParallel)
library(foreach)


 load("data/allShorebirdPrismFallCounts.RData")
# source("functions/GAM_basis_function.R")

n_cores <- 3
cluster <- makeCluster(n_cores, type = "PSOCK")
registerDoParallel(cluster)



fullrun <- foreach(sp = sps[c(25,11,13)],
                   .packages = c("jagsUI","tidyverse","ggmcmc"),
                   .inorder = FALSE,
                   .errorhandling = "pass") %dopar%
  {
 # sp = sps[25]  
    load("data/allShorebirdPrismFallCounts.RData")
    
    source("functions/GAM_basis_function.R")
    

#for(sp in sps){
#sp = sps[25]

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
midyear = floor(nyears/2)

# GAM seasonal basis function ---------------------------------------------

nKnots_season = 5
basis_season <- gam.basis.func(orig.preds = as.integer(unlist(dts[,"date"])),
                               nknots = nKnots_season,
                               standardize = "z",
                               random = F,
                               npredpoints = max(dts$date),
                               even_gaps = FALSE,
                               sm_name = "season")







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
                  
                  midyear = midyear)



mod.file = "models/AnnualSlopeSeasonalGAM.R"




parms = c("sdnoise",
          # "nu", #if optional heavy-tailed noise
          "sdgam_season",
          "B",
          "b",
          "sdsite",
          "N",
          "n_s",
          "alpha",
          "vis.sm_season")


#adaptSteps = 200              # Number of steps to "tune" the samplers.
burnInSteps = 10000            # Number of steps to "burn-in" the samplers.
nChains = 1                   # Number of chains to run.
numSavedSteps=2000          # Total number of steps in each chain to save.
thinSteps=20                   # Number of steps to "thin" (1=keep every step).
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
              #modules = NULL,
              model.file = mod.file)


t2 = Sys.time()

out2$n.eff
out2$Rhat
save(list = c("jags_data",
              "basis_season",
              "dts",
              "t2",
              "t1",
              "out2"),
     file = paste0("output/",sp,"slope_results.RData"))


# gg = ggs(out2$samples)
# 
# ggy = ggs(out2$samples,family = "B")
# bby2 = ggs(out2$samples,family = "beta")
# gga = ggs(out2$samples,family = "alpha")
# ggall = rbind(ggy,gga)
# ggmcmc(ggall,file = paste0("output/mcmc_",sp,".pdf"))
# ggmcmc(bby2,file = paste0("output/mcmc_bby2_",sp,".pdf"))
# 
# ggsd = ggs(out2$samples,family = "sd")
# ggmcmc(ggsd,file = paste0("output/mcmc_ggsd_",sp,".pdf"))




}#end species loops


stopCluster(cl = cluster)





# Plotting ----------------------------------------------------------------

library(tidybayes)

load("data/allShorebirdPrismFallCounts.RData")
source("functions/Utility_functions.R")

for(sp in sps){
  
  
  
  load(paste0("output/",sp,"slope_results.RData"))
  
  strats = unique(dts[,c("strat","Region")])
  strats = rename(strats,s = strat)
  
  
  
  sums = data.frame(out2$summary)
  names(sums) <- c("mean","sd","lci","lqrt","median","uqrt","uci","Rhat","n.eff","overlap0","f")
  sums$Parameter = row.names(sums)
  
  # compiling indices -------------------------------------
  
  n_inds <- extr_inds(param = "n_s")
  N_inds <- extr_inds(param = "N",regions = FALSE)
 
  
  
  
  # calculating trends  -----------------------------------------------------
  
  NSamples <- out2$samples %>% gather_draws(N[y])
  NSamples$year <- NSamples$y + 1973
  

  
  
  n_sSamples <- out2$samples %>% gather_draws(n_s[s,y])
  n_sSamples$year <- n_sSamples$y + 1973
  n_sSamples <- left_join(n_sSamples,strats,by = "s")
  
 
  
  
 
  t_n_s <- ItoT(inds = n_sSamples,regions = TRUE)
  
  t_N <- ItoT(inds = NSamples,regions = FALSE)
  
 
  
  # plotting indices --------------------------------------------------------
  
  
  plot_by_st <- plot_ind(inds = n_inds,
                         #smooth_inds = n_sm_inds,
                         raw = dts,
                         add_observed = TRUE,
                         add_samplesize = TRUE,
                         species = sp,
                         regions = TRUE,
                         title_size = 20,
                         axis_title_size = 18,
                         axis_text_size = 16)  
  
  
}

