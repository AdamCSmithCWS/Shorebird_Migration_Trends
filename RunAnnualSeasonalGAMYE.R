# SPECIES MCMC data-prep -------------------------------------------------------

library(jagsUI)
library(ggmcmc)

library(doParallel)
library(foreach)


 load("data/allShorebirdPrismFallCounts.RData")
# source("functions/GAM_basis_function.R")

n_cores <- 6
cluster <- makeCluster(n_cores, type = "PSOCK")
registerDoParallel(cluster)



fullrun <- foreach(sp = sps[c(3,4,6,8,21,22)],
                   .packages = c("jagsUI","tidyverse","ggmcmc"),
                   .inorder = FALSE,
                   .errorhandling = "pass") %dopar%
  {
    
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



mod.file = "models/AnnualSeasonalGAMYE.R"




parms = c("sdnoise",
          # "nu", #if optional heavy-tailed noise
          "sdgam_season",
          "sdgam_year",
          "sdgam_year_b",
          "beta_season",
          "b_year",
          "B_year",
          "sd_year",
          "sd_ye",
          "YE",
          "beta_year",
          "sdsite",
          "N",
          "n_s",
          "N_sm",
          "n_s_sm",
          "n_s_a1",
          "n_s_sm_a1",
          "n_s_a2",
          "n_s_sm_a2",
          "alpha",
          "vis.sm_season")


#adaptSteps = 200              # Number of steps to "tune" the samplers.
burnInSteps = 10000            # Number of steps to "burn-in" the samplers.
nChains = 3                   # Number of chains to run.
numSavedSteps=2000          # Total number of steps in each chain to save.
thinSteps=20                   # Number of steps to "thin" (1=keep every step).
nIter = ceiling( ( (numSavedSteps * thinSteps )+burnInSteps)) # Steps per chain.

t1 = Sys.time()




# MCMC sampling -----------------------------------------------------------



out2 = jagsUI(data = jags_data,
              parameters.to.save = parms,
              n.chains = nChains,
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
              "basis_year",
              "dts",
              "t2",
              "t1",
              "out2"),
     file = paste0("output/",sp,"GAMYE_results.RData"))

# gg = ggs(out2$samples)
# 
# ggy = ggs(out2$samples,family = "B_year")
# bby2 = ggs(out2$samples,family = "beta_year")
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

  
  if(file.exists(paste0("output/",sp,"GAMYE_results.RData"))){
    
  load(paste0("output/",sp,"GAMYE_results.RData"))

  strats = unique(dts[,c("strat","Region")])
  strats = rename(strats,s = strat)
  
  

sums = data.frame(out2$summary)
names(sums) <- c("mean","sd","lci","lqrt","median","uqrt","uci","Rhat","n.eff","overlap0","f")
sums$Parameter = row.names(sums)



# extracting the seasonal smooth ------------------------------------------

season_sm = extr_sum(param = "vis.sm_season",
                     index = "day",
                     log_retrans = TRUE) 


pp <- ggplot()+
  geom_line(data = season_sm,aes(x = day,y = mean))+
  geom_ribbon(data = season_sm,aes(x = day,ymax = uci,ymin = lci),alpha = 0.2)+
  ylab("")+
  xlab("Days since July 1")

pdf(file = paste0("Figures/",sp,"_Season_GAMYE.pdf"),
    width = 8.5,
    height = 11)
print(pp)
dev.off()


# compiling indices -------------------------------------

n_inds <- extr_inds(param = "n_s")
n_sm_inds <- extr_inds(param = "n_s_sm")
N_inds <- extr_inds(param = "N",regions = FALSE)
N_sm_inds <- extr_inds(param = "N_sm",regions = FALSE)




n_inds_a1 <- extr_inds(param = "n_s_a1")
n_sm_inds_a1 <- extr_inds(param = "n_s_sm_a1")






n_inds_a2 <- extr_inds(param = "n_s_a2")
n_sm_inds_a2 <- extr_inds(param = "n_s_sm_a2")





# Trends  -----------------------------------------------------

  NSamples <- out2$samples %>% gather_draws(N[y])
  NSamples$year <- NSamples$y + 1973
 
  N_smSamples <- out2$samples %>% gather_draws(N_sm[y])
  N_smSamples$year <- N_smSamples$y + 1973
 
   
  
  n_sSamples <- out2$samples %>% gather_draws(n_s[s,y])
  n_sSamples$year <- n_sSamples$y + 1973
 n_sSamples <- left_join(n_sSamples,strats,by = "s")
 
 n_s_smSamples <- out2$samples %>% gather_draws(n_s_sm[s,y])
 n_s_smSamples$year <- n_s_smSamples$y + 1973
 n_s_smSamples <- left_join(n_s_smSamples,strats,by = "s")
 
 
 n_s_a1Samples <- out2$samples %>% gather_draws(n_s_a1[s,y])
 n_s_a1Samples$year <- n_s_a1Samples$y + 1973
 n_s_a1Samples <- left_join(n_s_a1Samples,strats,by = "s")
 
 n_s_sm_a1Samples <- out2$samples %>% gather_draws(n_s_sm_a1[s,y])
 n_s_sm_a1Samples$year <- n_s_sm_a1Samples$y + 1973
 n_s_sm_a1Samples <- left_join(n_s_sm_a1Samples,strats,by = "s")
 
 
 
 n_s_a2Samples <- out2$samples %>% gather_draws(n_s_a2[s,y])
 n_s_a2Samples$year <- n_s_a2Samples$y + 1973
 n_s_a2Samples <- left_join(n_s_a2Samples,strats,by = "s")
 
 n_s_sm_a2Samples <- out2$samples %>% gather_draws(n_s_sm_a2[s,y])
 n_s_sm_a2Samples$year <- n_s_sm_a2Samples$y + 1973
 n_s_sm_a2Samples <- left_join(n_s_sm_a2Samples,strats,by = "s")
 
  
t_n_s_sm_a1 <- ItoT(inds = n_s_sm_a1Samples,regions = TRUE,index_type = "smoothed",retransformation_type = "lognormal_only")
t_n_s_sm_15_a1 <- ItoT(inds = n_s_sm_a1Samples,regions = TRUE,start = 2004,index_type = "smoothed",retransformation_type = "lognormal_only")


t_n_s_sm_a2 <- ItoT(inds = n_s_sm_a2Samples,regions = TRUE,index_type = "smoothed",retransformation_type = "none")
t_n_s_sm_15_a2 <- ItoT(inds = n_s_sm_a2Samples,regions = TRUE,start = 2004,index_type = "smoothed",retransformation_type = "none")


t_n_s_sm <- ItoT(inds = n_s_smSamples,regions = TRUE,index_type = "smoothed")


t_n_s_sm_15 <- ItoT(inds = n_s_smSamples,regions = TRUE,start = 2004,index_type = "smoothed")

t_n_s <- ItoT(inds = n_sSamples,regions = TRUE)
t_n_s_15 <- ItoT(inds = n_sSamples,regions = TRUE,start = 2004)

t_N <- ItoT(inds = NSamples,regions = FALSE)


t_N_sm <- ItoT(inds = N_smSamples,regions = FALSE,index_type = "smoothed")
t_N_sm_15 <- ItoT(inds = N_smSamples,regions = FALSE,start = 2004,index_type = "smoothed")



# Slope Trends ------------------------------------------------------------


t_n_s_slope <- ItoT_slope(inds = n_sSamples,regions = TRUE)

t_N_slope <- ItoT_slope(inds = NSamples,regions = FALSE)


t_n_s_slope_15 <- ItoT_slope(inds = n_sSamples,regions = TRUE,start = 2004)

t_N_slope_15 <- ItoT_slope(inds = NSamples,regions = FALSE,start = 2004)

trend_out <- bind_rows(t_N,
                       t_N_slope,
                       t_N_slope_15,
                       t_N_sm,
                       t_N_sm_15,
                       t_n_s_sm,
                       t_n_s_sm_a1,
                       t_n_s_sm_a2,
                       t_n_s_sm_15,
                       t_n_s_sm_15_a1,
                       t_n_s_sm_15_a2,
                       t_n_s,
                       t_n_s_slope,
                       t_n_s_15,
                       t_n_s_slope_15)

write.csv(trend_out,file = paste0("Trends/trends_GAMYE_",sp,".csv"),row.names = F)

# plotting indices --------------------------------------------------------



plot_by_st <- plot_ind(inds = n_inds_a2,
                       #smooth_inds = ,
                       raw = dts,
                       add_observed = TRUE,
                       add_samplesize = TRUE,
                       species = sp,
                       regions = TRUE,
                       title_size = 20,
                       axis_title_size = 18,
                       axis_text_size = 16)  

pdf(file = paste0("Figures/",sp,"GAMYE_A2.pdf"),
    width = 8.5,
    height = 11)
print(plot_by_st)
dev.off()


plot_by_st <- plot_ind(inds = n_inds_a1,
                       smooth_inds = n_inds,
                       raw = dts,
                       add_observed = TRUE,
                       add_samplesize = TRUE,
                       species = sp,
                       regions = TRUE,
                       title_size = 20,
                       axis_title_size = 18,
                       axis_text_size = 16)  

pdf(file = paste0("Figures/",sp,"GAMYE_A1.pdf"),
    width = 8.5,
    height = 11)
print(plot_by_st)
dev.off()



plot_by_st <- plot_ind(inds = n_inds,
                       smooth_inds = n_sm_inds,
                       raw = dts,
                       add_observed = TRUE,
                       add_samplesize = TRUE,
                       species = sp,
                       regions = TRUE,
                       title_size = 20,
                       axis_title_size = 18,
                       axis_text_size = 16)  
 
pdf(file = paste0("Figures/",sp,"GAMYE.pdf"),
    width = 8.5,
    height = 11)
print(plot_by_st)
dev.off()


}#end if output exists

}


