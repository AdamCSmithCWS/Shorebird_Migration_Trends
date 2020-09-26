# SPECIES MCMC data-prep -------------------------------------------------------

library(jagsUI)
library(ggmcmc)


library(doParallel)
library(foreach)


 load("data/allShorebirdPrismFallCounts.RData")
# source("functions/GAM_basis_function.R")

n_cores <- 4
cluster <- makeCluster(n_cores, type = "PSOCK")
registerDoParallel(cluster)



fullrun <- foreach(sp = sps[c(11,25,3,8)],
                   .packages = c("jagsUI","tidyverse","ggmcmc"),
                   .inorder = FALSE,
                   .errorhandling = "pass") %dopar%
  {
 # sp = sps[25]  
    load("data/allShorebirdPrismFallCounts.RData")
    
    source("functions/GAM_basis_function.R")
    

#for(sp in sps){
#sp = sps[25]

dts <- filter(ssData,CommonName == sp,
              YearCollected > 1977)
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

#number of years with non-zero observations, by region
nyrs_region <- dts %>% 
  filter(present == TRUE) %>%  
  group_by(Region,YearCollected) %>% 
  summarise(nobs = n()) %>% 
  group_by(Region) %>% 
  summarise(nyears = n(),
            span_years = yspn(YearCollected))

#strats with 30 or more years of non-zero, observations - species has to be observed in a region in at least 2/3 of the years in the time-series
regions_keep <- nyrs_region[which(nyrs_region$nyears > 29),"Region"]

dts <- filter(dts,Region %in% regions_keep$Region) 


# prepare jags data -------------------------------------------------------

fday = min(dts$doy)-1

syear = min(dts$YearCollected)
dts <- dts %>% mutate(count = as.integer(ObservationCount),
                      year = as.integer(YearCollected),
                      yr = as.integer(year-syear),
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

ndays <- basis_season$npredpoints_season






jags_data <- list(count = as.integer(unlist(dts$count)),
                  yr = as.integer(unlist(dts$yr)),
                  site = as.integer(unlist(dts$site)),
                  strat = as.integer(unlist(dts$strat)),
                  
                  nyears = nyears,
                  nstrata = nstrata,
                  nsites = nsites,
                  ncounts = ncounts,
                  ndays = ndays,
                  
                  season_basis = basis_season$season_basis,
                  season_basispred = basis_season$season_basispred,
                  nknots_season = basis_season$nknots_season,
                  
                  midyear = midyear)



mod.file = "models/AnnualSlopeSeasonalGAM_ZIP.R"




parms = c("sdnoise",
           "nu", #if optional heavy-tailed noise
          "nu_site",
          "sdgam_season",
          "B",
          "b",
          "sdsite",
          "N",
          "n_s",
          "N_comp",
          # "n_s_a1",
          # "n_s_a2",
          # "N_sc2",
          # "N_sc",
          # "n_s_scaled",
          # "n_s_scaled2",
          "alpha",
          "vis.sm_season",
          "psi")


#adaptSteps = 200              # Number of steps to "tune" the samplers.
burnInSteps = 100            # Number of steps to "burn-in" the samplers.
nChains = 3                   # Number of chains to run.
numSavedSteps=400          # Total number of steps in each chain to save.
thinSteps=5                   # Number of steps to "thin" (1=keep every step).
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
     file = paste0("output/",sp,"slope_ZIP_results.RData"))


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
  
  if(file.exists(paste0("output/",sp,"slope_ZIP_results.RData"))){
  
  load(paste0("output/",sp,"slope_ZIP_results.RData"))
    fyear = (min(dts$YearCollected))
    
  strats = unique(dts[,c("strat","Region")])
  strats = rename(strats,s = strat)
  
  
  
  sums = data.frame(out2$summary)
  names(sums) <- c("mean","sd","lci","lqrt","median","uqrt","uci","Rhat","n.eff","overlap0","f")
  sums$Parameter = row.names(sums)
  
  # compiling indices -------------------------------------
  
  n_inds <- extr_inds(param = "n_s")
  N_inds <- extr_inds(param = "N",regions = FALSE)
 
 
  
  sdnoise_st = extr_sum(param = "sdnoise",
                    index = c("s"),
                    log_retrans = F) 
  sdnoise_st <- left_join(sdnoise_st,strats,by = "s")
  
  nu_st = extr_sum(param = "nu",
                        index = c("s"),
                        log_retrans = F) 
  nu_st <- left_join(nu_st,strats,by = "s")
  
  
  
  psi_st = extr_sum(param = "psi",
                       index = c("s"),
                       log_retrans = F) 
 
  psiSamples <- out2$samples %>% gather_draws(psi)
  sdnoiseSamples <- out2$samples %>% gather_draws(sdnoise[s])
  
# extracting the seasonal smooth ------------------------------------------

  season_sm = extr_sum(param = "vis.sm_season",
                       index = c("day","s"),
                       log_retrans = TRUE) 
  
  season_sm <- left_join(season_sm,strats)
  pp <- ggplot()+
    geom_line(data = season_sm,aes(x = day,y = mean,colour = Region,group = Region))+
    geom_ribbon(data = season_sm,aes(x = day,ymax = uci,ymin = lci),alpha = 0.2)+
    ylab("")+
    xlab("Days since July 1")+
    facet_wrap(facets = ~Region,ncol = 3,scales = "free")
  
  
  
  pdf(file = paste0("Figures/",sp,"_Season_ZIP_slope.pdf"),
      width = 8.5,
      height = 11)
  print(pp)
  dev.off()
  
  # calculating trends  -----------------------------------------------------
  

 
  NSamples <- out2$samples %>% gather_draws(N[y])
  NSamples$year <- NSamples$y + fyear-1
  
  
  
  N_compSamples <- out2$samples %>% gather_draws(N_comp[y])
  N_compSamples$year <- N_compSamples$y + fyear-1
  
  n_sSamples <- out2$samples %>% gather_draws(n_s[s,y])
  n_sSamples$year <- n_sSamples$y + fyear-1
  n_sSamples <- left_join(n_sSamples,strats,by = "s")
  
  
  
  t_n_s <- ItoT(inds = n_sSamples,regions = TRUE)
  
  
  
  
  t_n_sS <- ItoT_slope(inds = n_sSamples,regions = TRUE)
  
  
  # t_n_s_a1 <- ItoT(inds = n_s_a1Samples,regions = TRUE,retransformation_type = "lognormal_only")
  # 
  # 
  # t_n_s_a2 <- ItoT(inds = n_s_a2Samples,regions = TRUE,retransformation_type = "none")
  # 
  # 
  # 
  t_N <- ItoT(inds = NSamples,regions = FALSE)
  
  
  t_n_s_15 <- ItoT(inds = n_sSamples,regions = TRUE,start= 2004)
  
  
  t_n_sS_15 <- ItoT_slope(inds = n_sSamples,regions = TRUE,start= 2004)
  
  
  t_N_15 <- ItoT(inds = NSamples,regions = FALSE,start= 2004)
  
  t_NS <- ItoT_slope(inds = NSamples,regions = FALSE)
  t_NS_15 <- ItoT_slope(inds = NSamples,regions = FALSE,start= 2004)
  
  
  
  t_N_comp <- ItoT(inds = N_compSamples,regions = FALSE,start = 1978)
  
  t_N_comp_15 <- ItoT(inds = N_compSamples,regions = FALSE,start= 2004)
  
  t_N_compS <- ItoT_slope(inds = N_compSamples,regions = FALSE)
  t_N_compS_15 <- ItoT_slope(inds = N_compSamples,regions = FALSE,start= 2004)
  
  
  
  trend_out <- bind_rows(t_N,
                         t_N_15,
                         t_NS,
                         t_NS_15,
                         t_N_comp,
                         t_N_comp_15,
                         t_N_compS,
                         t_N_compS_15,
                         t_n_s,
                         t_n_sS,
                         t_n_s_15,
                         t_n_sS_15)
  
  write.csv(trend_out,file = paste0("Trends/trends_slope_ZIP_",sp,".csv"),row.names = F)
  
  # plotting indices --------------------------------------------------------
  
  
  # plot_Hyper <- plot_ind(inds = N_inds,
  #                        #smooth_inds = ,
  #                        raw = dts,
  #                        add_observed = TRUE,
  #                        add_samplesize = TRUE,
  #                        species = sp,
  #                        regions = FALSE,
  #                        title_size = 20,
  #                        axis_title_size = 18,
  #                        axis_text_size = 16)  
  # 
  # pdf(file = paste0("Figures/",sp,"_Hyperparameter_slope.pdf"),
  #     width = 8.5,
  #     height = 11)
  # print(plot_Hyper)
  # dev.off()
  
  
  
  



plot_by_st <- plot_ind(inds = n_inds,
                       smooth_inds = NULL,
                       raw = dts,
                       add_observed = TRUE,
                       add_samplesize = TRUE,
                       species = sp,
                       regions = TRUE,
                       title_size = 20,
                       axis_title_size = 18,
                       axis_text_size = 16)  

pdf(file = paste0("Figures/",sp,"_slope_ZIP.pdf"),
    width = 8.5,
    height = 11)
print(plot_by_st)
dev.off()




}# end if output exists

}

