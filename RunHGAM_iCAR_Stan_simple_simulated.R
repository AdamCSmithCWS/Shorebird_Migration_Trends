# SPECIES MCMC data-prep -------------------------------------------------------
library(tidyverse)
library(rstan)
rstan_options(auto_write = TRUE, javascript = FALSE)
library(shinystan)
library(sf)
library(spdep)
library(ggforce)
library(tidybayes)


# library(doParallel)
# library(foreach)

#load data
load("data/allShorebirdPrismFallCounts.RData")
grid_spacing <- 300000  # size of squares, in units of the CRS (i.e. meters for lae)

 
FYYYY = 1980
# 
# w_cosewic = sps[c(2:4,7,10,12:20,22,11,25)]
# 
# # for(sp in w_cosewic){
# #   
# #   if(file.exists(paste0("output/",sp,"_GAMYE_strat_simple",grid_spacing/1000,".RData"))){w_cosewic <- w_cosewic[-which(w_cosewic == sp)]}
# # 
# # }  
# 
# sps_remain = sps[-which(sps %in% w_cosewic)]


 
   
    #load(paste0("data/data_simulated_stable_upturn_GAMYE_strat_simple",".RData"))
#load(paste0("data/data_simulated_stable_Balanced_upturn_GAMYE_strat_simple",".RData"))

#load(paste0("data/data_simulated_stable_downturn2_GAMYE_strat_simple.RData"))
    #load(paste0("data/data_simulated_stable_Balanced_downturn2_GAMYE_strat_simple.RData"))
    load(paste0("data/data_simulated_stable_decline_GAMYE_strat_simple.RData"))
    #load(paste0("data/data_simulated_balanced_stable_decline_GAMYE_strat_simple.RData"))
    
    
mod.file = c("models/GAMYE_strata_two_season_normal_tail_noCon.stan")

## compile model
slope_icar_model = stan_model(file=mod.file)

## run sampler on model, data
slope_icar_stanfit <- sampling(slope_icar_model,
                               data=stan_data,
                               verbose=TRUE, refresh=100,
                               chains=4, iter=1800,
                               warmup=1200,
                               cores = 4,
                               pars = c(parms),
                               control = list(adapt_delta = 0.9,
                                              max_treedepth = 14))



save(list = c("slope_icar_stanfit",
              "stan_data",
              "dts",
              "real_grid",
              "strats_dts",
              "strat_regions",
              "mod.file",
              "parms"),
     file = paste0("output/simulated_stable_decline_noCon_GAMYE_strat_simple.RData"))



save(list = c("slope_icar_stanfit",
              "stan_data",
              "dts",
              "real_grid",
              "strats_dts",
              "strat_regions",
              "mod.file",
              "parms"),
     file = paste0("output/simulated_stable_decline_alt_spati_prior_GAMYE_strat_simple.RData"))


save(list = c("slope_icar_stanfit",
              "stan_data",
              "dts",
              "real_grid",
              "strats_dts",
              "strat_regions",
              "mod.file",
              "parms"),
     file = paste0("output/simulated_balanced_stable_decline_upturn_GAMYE_strat_simple.RData"))


save(list = c("slope_icar_stanfit",
              "stan_data",
              "dts",
              "real_grid",
              "strats_dts",
              "strat_regions",
              "mod.file",
              "parms"),
     file = paste0("output/simulated_stable_decline_upturn_GAMYE_strat_simple.RData"))



save(list = c("slope_icar_stanfit",
              "stan_data",
              "dts",
              "real_grid",
              "strats_dts",
              "strat_regions",
              "mod.file",
              "parms"),
     file = paste0("output/simulated_stable_Balanced_upturn_GAMYE_strat_simple.RData"))

save(list = c("slope_icar_stanfit",
              "stan_data",
              "dts",
              "real_grid",
              "strats_dts",
              "strat_regions",
              "mod.file",
              "parms"),
     file = paste0("output/simulated_stable_upturn_GAMYE_strat_simple.RData"))




save(list = c("slope_icar_stanfit",
              "stan_data",
              "dts",
              "real_grid",
              "strats_dts",
              "strat_regions",
              "mod.file",
              "parms"),
     file = paste0("output/simulated_stable_downturn2altprior_GAMYE_strat_simple.RData"))



save(list = c("slope_icar_stanfit",
              "stan_data",
              "dts",
              "real_grid",
              "strats_dts",
              "strat_regions",
              "mod.file",
              "parms"),
     file = paste0("output/simulated_stable_downturn2altsdyeargam_GAMYE_strat_simple.RData"))

save(list = c("slope_icar_stanfit",
              "stan_data",
              "dts",
              "real_grid",
              "strats_dts",
              "strat_regions",
              "mod.file",
              "parms"),
     file = paste0("output/simulated_stable_uptick_GAMYE_strat_simple.RData"))


save(list = c("slope_icar_stanfit",
              "stan_data",
              "dts",
              "real_grid",
              "strats_dts",
              "strat_regions",
              "mod.file",
              "parms"),
     file = paste0("output/simulated_stable_Balanced_downturn2_GAMYE_strat_simple.RData"))

















loo_out = loo(slope_icar_stanfit)



launch_shinystan(slope_icar_stanfit) 

source("functions/utility_functions.R")
source("functions/palettes.R")


sp = "Balanced stable decline"
three_gen = 20
y3g <- 2019-three_gen

syear = min(dts$YearCollected)


Nsamples <- slope_icar_stanfit %>% gather_draws(N[y])
Nsamples$year <- Nsamples$y + (syear-1)

NSmoothsamples <- slope_icar_stanfit %>% gather_draws(NSmooth[y])
NSmoothsamples$year <- NSmoothsamples$y + (syear-1)

# Alphas by site ----------------------------------------------------------
alpha_samples <- slope_icar_stanfit %>% gather_draws(alpha[s])

sites_strat = (stan_data$sites)
nstrata = stan_data$nstrata
nsites_strat = stan_data$nsites_strat
#offsets = stan_data$site_size

sitesbystrat = NULL
for(st in 1:stan_data$nstrata){
   tmp = data.frame(strat = st,
                    site = sites_strat[1:nsites_strat[st],st])
   sitesbystrat <- bind_rows(sitesbystrat,tmp)
}
sitesbystrat <- arrange(sitesbystrat,site)

alphas = alpha_samples %>% group_by(s) %>% 
   summarise(mean = mean(exp(.value)),
             lci = quantile(exp(.value),0.025),
             uci = quantile(exp(.value),0.975)) %>% 
   mutate(site = s) %>% left_join(.,sitesbystrat)

## site-effects against the observed counts at each site
nr = ceiling(sqrt(nstrata))
obs_by_alpha = ggplot(data = alphas,aes(x = site,y = mean))+
   geom_pointrange(aes(ymin = lci,ymax = uci),colour = "red")+
   geom_point(data = dts, inherit.aes = FALSE,aes(x = site,y = count),
              fill = grDevices::grey(0.6),colour = grDevices::grey(0),alpha = 0.1,
              position = position_jitter(width = 0.2,height = 0))+
   labs(title = sp)+
   facet_wrap(~strat,nrow = nr,ncol = nr,scales = "free")+
   theme_minimal()

pdf(paste0("Figures/",sp,"_obs_by_alpha.pdf"),
    width = 11,
    height = 11)
print(obs_by_alpha)
dev.off()

scale_adj_means <- summary(slope_icar_stanfit,pars = c("sdnoise","ALPHA1"))$summary[,1]

scale_adj <- ((scale_adj_means[["sdnoise"]]^2)*0.5 + scale_adj_means[["ALPHA1"]])


# n_by_strat <- dts %>% group_by(strat) %>% 
#   summarise(n_cst = n())
# n_countsby_site <- dts %>% group_by(site,strat) %>% 
#   summarise(n_c = n()) %>% 
#   left_join(.,n_by_strat,by = "strat") %>% 
#   mutate(n_c = n_c/n_cst) %>% 
#   select(site,strat,n_c)


n_countsby_site <- dts %>% group_by(site) %>% 
   summarise(n_c = n()) 

if(grepl(x = mod.file,pattern = "two_season")){
   
   strat_season_strat <- dts %>% distinct(seas_strat,strat,hex_name)
   
   alphas <- left_join(alphas,strat_season_strat,by = "strat")
   
   
   strat_offs <- alphas %>% left_join(.,n_countsby_site,by = "site") %>% 
      group_by(strat,seas_strat) %>% 
      summarise(adjs = mean(mean))
   
   season_samples <- slope_icar_stanfit %>% gather_draws(season_pred[d,s])
   seasonEffectT = season_samples %>% group_by(d,s) %>% 
      summarise(mean = mean(exp(.value+scale_adj)),
                lci = quantile(exp(.value+scale_adj),0.025),
                uci = quantile(exp(.value+scale_adj),0.975)) %>% 
      mutate(day = d,
             seas_strat = s) 
   
   obs_season <- dts %>% group_by(date,seas_strat) %>% 
      summarise(mean = mean(count),
                median = median(count),
                lqrt = quantile(count,0.05),
                uqrt = quantile(count,0.95)) %>% 
      mutate(day = date)
   
   sse <- seasonEffectT %>% group_by(seas_strat) %>% 
      group_split()
   
   #seasonEffect <- expand_grid(seasonEffect,strat = c(1:nstrata))
   tout <- NULL
   for(j in 1:nstrata){
      wg = as.integer(strat_offs[which(strat_offs$strat == j),"seas_strat"])
      tmp <- sse[[wg]]
      tmp$strat <- j
      tout <- bind_rows(tout,tmp)
   }
   
   seasonEffect <- left_join(tout,strat_offs,by = "strat") %>% 
      mutate(mean = mean*adjs,
             lci = lci*adjs,
             uci = uci*adjs)
   
   
   
   ncl = 3
   nrr = 3
   ppag = ncl*nrr
   rem = nstrata-(floor(nstrata/ppag)*ppag)
   if(rem < ncl){
      nrr = 4
      ppag = ncl*nrr
      rem = nstrata-(floor(nstrata/ppag)*ppag)
   }
   if(rem < ncl){
      nrr = 3
      ncl = 2
      ppag = ncl*nrr
      rem = nstrata-(floor(nstrata/ppag)*ppag)
   }
   if(rem < ncl){
      nrr = 5
      ncl = 3
      ppag = ncl*nrr
      rem = nstrata-(floor(nstrata/ppag)*ppag)
   }
   
   tmp_season_graphs <- vector(mode = "list",length = ceiling(nstrata/ppag))
   
   pdf(file = paste0("Figures/",sp,"simple_Season.pdf"),
       width = 8.5,
       height = 8.5)
   
   
   for(jj in 1:ceiling(nstrata/ppag)){
      #yup <- quantile(dts$count,0.99)
      pp <- ggplot()+
         # geom_pointrange(data = obs_season,inherit.aes = FALSE,
         #                 aes(x = day,y = mean,ymin = lqrt,ymax = uqrt),alpha = 0.1)+
         geom_point(data = dts,aes(x = date,y = count,colour = year),alpha = 0.5,size = 1,position = position_jitter(width = 0.7,height = 0))+
         scale_colour_viridis_c()+
         geom_line(data = seasonEffect,aes(x = day,y = mean),inherit.aes = FALSE)+
         #coord_cartesian(ylim = c(0,yup))+
         geom_ribbon(data = seasonEffect,aes(x = day,y = mean,ymax = uci,ymin = lci),alpha = 0.2,inherit.aes = FALSE)+
         ylab("")+
         xlab("Days since July 1")+
         facet_wrap_paginate(facets = ~strat,page = jj,nrow = nrr, ncol = ncl,scales = "free")+
         labs(title = sp)
      tmp_season_graphs[[jj]] <- pp
      print(pp)
   }
   dev.off()
   
   # pp <- ggplot(data = seasonEffect,aes(x = day,y = mean))+
   #   # geom_pointrange(data = obs_season,inherit.aes = FALSE,
   #   #                 aes(x = day,y = mean,ymin = lqrt,ymax = uqrt),alpha = 0.1)+
   #   geom_point(data = dts,inherit.aes = FALSE,aes(x = day,y = count),alpha = 0.1,size = 1)+
   #   geom_line()+
   #   geom_ribbon(aes(x = day,y = mean,ymax = uci,ymin = lci),alpha = 0.2)+
   #   ylab("")+
   #   xlab("Days since July 1")+
   #   labs(title = sp)+
   #   facet_wrap(facets = ~seas_strat,ncol = 3,scales = "free")
   # 
   
}else{
   strat_offs <- alphas %>% group_by(strat) %>% 
      summarise(adjs = mean(mean))
   
   season_samples <- slope_icar_stanfit %>% gather_draws(season_pred[d])
   seasonEffect = season_samples %>% group_by(d) %>% 
      summarise(mean = mean(exp(.value+scale_adj)),
                lci = quantile(exp(.value+scale_adj),0.025),
                uci = quantile(exp(.value+scale_adj),0.975)) %>% 
      mutate(day = d) 
   
   obs_season <- dts %>% group_by(date) %>% 
      summarise(mean = mean(count),
                median = median(count),
                lqrt = quantile(count,0.05),
                uqrt = quantile(count,0.95)) %>% 
      mutate(day = date)
   
   seasonEffect <- expand_grid(seasonEffect,strat = c(1:nstrata))
   seasonEffect <- left_join(seasonEffect,strat_offs) %>% 
      mutate(mean = mean*adjs,
             lci = lci*adjs,
             uci = uci*adjs)
   
   
   ncl = 3
   nrr = 3
   ppag = ncl*nrr
   rem = nstrata-(floor(nstrata/ppag)*ppag)
   if(rem < ncl){
      nrr = 4
      ppag = ncl*nrr
      rem = nstrata-(floor(nstrata/ppag)*ppag)
   }
   if(rem < ncl){
      nrr = 3
      ncl = 2
      ppag = ncl*nrr
      rem = nstrata-(floor(nstrata/ppag)*ppag)
   }
   if(rem < ncl){
      nrr = 5
      ncl = 3
      ppag = ncl*nrr
      rem = nstrata-(floor(nstrata/ppag)*ppag)
   }
   tmp_season_graphs <- vector(mode = "list",length = ceiling(nstrata/ppag))
   
   pdf(file = paste0("Figures/",sp,"simple_Season.pdf"),
       width = 8.5,
       height = 8.5)
   
   for(jj in 1:ceiling(nstrata/ppag)){
      #yup <- quantile(dts$count,0.99)
      pp <- ggplot()+
         # geom_pointrange(data = obs_season,inherit.aes = FALSE,
         #                 aes(x = day,y = mean,ymin = lqrt,ymax = uqrt),alpha = 0.1)+
         geom_point(data = dts,aes(x = date,y = count,colour = year),alpha = 0.5,size = 1,position = position_jitter(width = 0.7,height = 0))+
         scale_colour_viridis_c()+
         geom_line(data = seasonEffect,aes(x = day,y = mean),inherit.aes = FALSE)+
         #coord_cartesian(ylim = c(0,yup))+
         geom_ribbon(data = seasonEffect,aes(x = day,y = mean,ymax = uci,ymin = lci),alpha = 0.2,inherit.aes = FALSE)+
         ylab("")+
         xlab("Days since July 1")+
         facet_wrap_paginate(facets = ~strat,page = jj,nrow = nrr, ncol = ncl,scales = "free")+
         labs(title = sp)
      tmp_season_graphs[[jj]] <- pp
      print(pp)
   }
   dev.off()
   
}

season_graphs[[sp]] <- tmp_season_graphs



TRENDSout <- NULL
# calculate trends continent --------------------------------------------------------

t_NSmooth_80 <- ItoT(inds = NSmoothsamples,
                     start = 1980,
                     end = 2019,
                     regions = NULL,
                     qs = 95,
                     sp = sp,
                     type = "Long-term")
TRENDSout <- bind_rows(TRENDSout,t_NSmooth_80)

t_NSmooth_04 <- ItoT(inds = NSmoothsamples,
                     start = 2004,
                     end = 2019,
                     regions = NULL,
                     qs = 95,
                     sp = sp,
                     type = "15-year")
TRENDSout <- bind_rows(TRENDSout,t_NSmooth_04)

t_NSmooth_3g <- ItoT(inds = NSmoothsamples,
                     start = y3g,
                     end = 2019,
                     regions = NULL,
                     qs = 95,
                     sp = sp,
                     type = "Recent-three-generation")
TRENDSout <- bind_rows(TRENDSout,t_NSmooth_3g)


syL3g = y3g-(2019-y3g)
syL3g = max(1980,syL3g)
t_NSmooth_L3g <- ItoT(inds = NSmoothsamples,
                      start = syL3g,
                      end = y3g,
                      regions = NULL,
                      qs = 95,
                      sp = sp,
                      type = "Previous-three-generation")

TRENDSout <- bind_rows(TRENDSout,t_NSmooth_L3g)



anot_funct <- function(x){
   ant = paste(signif(x$percent_change,3),
               "% ",
               x$start_year,
               ":",
               x$end_year,
               "[",
               signif(x$p_ch_lci,3),
               " : ",
               signif(x$p_ch_uci,3),
               "]")
}

anot_90 = anot_funct(t_NSmooth_L3g)
anot_80 = anot_funct(t_NSmooth_80)
anot_07 = anot_funct(t_NSmooth_3g)


# Extract Annual Indices ------------------------------------------------


indicesN <- index_summary(parm = "N",
                          dims = "year",
                          site_scale = FALSE,
                          season_scale = FALSE,
                          strat_offsets = strat_offs)


indicesNSmooth <- index_summary(parm = "NSmooth",
                                dims = "year",
                                site_scale = FALSE,
                                season_scale = FALSE,
                                strat_offsets = strat_offs)






indices = bind_rows(indicesN,indicesNSmooth)
indices$year = indices$year + (syear-1)
yup = max(max(indices$PI97_5)*1.5,quantile(indices$obsmean,0.7))

indices$parm <- factor(indices$parm,ordered = T,levels = c("NSmooth","N"))

N_gg = ggplot(data = indices,aes(x = year, y = PI50,fill = parm))+
   geom_ribbon(aes(ymin = PI2_5,ymax = PI97_5),alpha = 0.2)+
   geom_line(aes(colour = parm))+
   labs(title = paste(sp,"Survey-wide trajectory (full and smooth) with obs means"))+
   annotate("text", x = 1997, y = yup*0.9, label = anot_80)+
   annotate("text", x = 1997, y = yup*0.8, label = anot_90)+
   annotate("text", x = 1997, y = yup*0.7, label = anot_07)+
   my_col2_traj+
   coord_cartesian(ylim = c(0,yup))+
   geom_point(aes(y = obsmean,size = mean_counts_incl_strata),colour = grey(0.5),alpha = 0.3)+
   theme_classic()+
   theme(legend.position = "none")#+
#scale_size_area()


pdf(paste0("figures/",sp,FYYYY,"_GAMYE_survey_wide_trajectory_simple",grid_spacing/1000,".pdf"))
print(N_gg)
dev.off()

#sp_ind_plots_diagnostic[[sp]] <- N_gg




N_gg_simple = ggplot(data = indices,aes(x = year, y = PI50,fill = parm))+
   geom_ribbon(aes(ymin = PI2_5,ymax = PI97_5),alpha = 0.2)+
   geom_line(aes(colour = parm))+
   labs(title = paste(sp,"Survey-wide trajectory (full and smooth)"))+
   my_col2_traj+
   xlab("")+
   ylab("Modeled mean count")+
   theme_classic()+
   theme(legend.position = "none")#+
#scale_size_area()

#print(N_gg_simple)

#sp_ind_plots[[sp]] <- N_gg_simple


# Continental Smooth spaghetti plot ---------------------------------------
indicesNSmooth$year = indicesNSmooth$year + (syear-1)

set.seed(2019)
r_draws <- sample(size = 50,x = 1:max(NSmoothsamples$.draw))
spg_trajs <- NSmoothsamples %>% filter(.draw %in% r_draws) %>% 
   mutate(draw_f = factor(.draw))

lev = spg_trajs %>% 
   ungroup() %>% 
   filter(year == 2019) %>% 
   select(draw_f,.value) %>% 
   mutate(midp = log(.value)) %>% 
   select(draw_f,midp)

spg_trajs <- left_join(spg_trajs,lev,by = "draw_f")


n_gg_spag = ggplot(data = indicesNSmooth,aes(x = year, y = PI50))+
   geom_ribbon(aes(ymin = PI2_5,ymax = PI97_5),alpha = 0.25)+
   geom_line(size =2)+
   labs(title = paste(sp,"Random selection of posterior draws of survey-wide trajectories"),
        subtitle = "Colour of each posterior draw reflects the value in 1980, demonstrating similar smooths across draws")+
   xlab("")+
   ylab("Modeled mean count")+
   theme_classic()+
   theme(legend.position = "none") + 
   geom_line(data = spg_trajs,aes(x = year,y = .value,group = draw_f,colour = midp),alpha = 0.8,inherit.aes = FALSE,size = 1)+
   scale_color_viridis_c(aesthetics = "colour")+
   scale_y_log10()

sp_spag_plots_diagnostic[[sp]] <- n_gg_spag


print(n_gg_spag)

pdf(paste0("figures/",sp,FYYYY,"_GAMYE_survey_wide_trajectory_simple",grid_spacing/1000,".pdf"))
print(N_gg_simple)
dev.off()


alpha_adjs <- alphas %>% mutate(adjs = mean) %>% select(site,adjs)


nstrata = stan_data$nstrata
indicesnsmooth <- index_summary(parm = "nsmooth",
                                dims = c("stratn","year"),
                                season_scale = FALSE,
                                site_scale = FALSE,
                                site_offsets = alpha_adjs)

indicesn <- index_summary(parm = "n",
                          dims = c("stratn","year"),
                          season_scale = FALSE,
                          site_scale = FALSE,
                          site_offsets = alpha_adjs)

indices_strat = bind_rows(indicesn,indicesnsmooth)
indices_strat$year = indices_strat$year + (syear-1)
indices_strat <- left_join(indices_strat,strats_dts, by = "stratn")
indices_strat$parm <- factor(indices_strat$parm,ordered = T,levels = c("nsmooth","n"))



pdf(file = paste0("figures/", sp,FYYYY,"_GAMYE_Strata_trajectories_simple",grid_spacing/1000,".pdf"),
    width = 8.5,
    height = 11)
print(N_gg)
# ncl = 3
# nrr = 3
# ppag = ncl*nrr
# rem = nstrata-(floor(nstrata/ppag)*ppag)
# if(rem < ncl){
#   nrr = 4
#   ppag = ncl*nrr
#   rem = nstrata-(floor(nstrata/ppag)*ppag)
# }
# if(rem < ncl){
#   nrr = 3
#   ncl = 2
#   ppag = ncl*nrr
#   rem = nstrata-(floor(nstrata/ppag)*ppag)
# }
# if(rem < ncl){
#   nrr = 5
#   ncl = 3
#   ppag = ncl*nrr
#   rem = nstrata-(floor(nstrata/ppag)*ppag)
# }
tmp_sp_ind_plots <- vector(mode = "list",length = ceiling(nstrata/ppag))
tmp_sp_ind_plots_diagnostic <- tmp_sp_ind_plots
for(jj in 1:ceiling(nstrata/ppag)){
   n_gg_simple = ggplot(data = indices_strat,aes(x = year, y = PI50,fill = parm))+
      geom_ribbon(aes(ymin = PI2_5,ymax = PI97_5),alpha = 0.2)+
      geom_line(aes(colour = parm))+
      my_col2_traj+
      xlab("")+
      ylab("Modeled mean count")+
      theme_classic()+
      theme(legend.position = "none")+
      labs(title = sp)+
      facet_wrap_paginate(facets = ~stratn,page = jj,nrow = nrr, ncol = ncl,scales = "free")
   #tmp_sp_ind_plots[[jj]] <- n_gg_simple
   
   
   n_gg = ggplot(data = indices_strat,aes(x = year, y = PI50,fill = parm))+
      geom_ribbon(aes(ymin = PI2_5,ymax = PI97_5),alpha = 0.2)+
      geom_line(aes(colour = parm))+
      my_col2_traj+
      geom_point(aes(y = obsmean,size = mean_counts_incl_sites),colour = grey(0.5),alpha = 0.3)+
      theme_classic()+
      labs(title = sp)+
      facet_wrap_paginate(facets = ~stratn,page = jj,nrow = nrr, ncol = ncl,scales = "free")+
      scale_size_area()
   #tmp_sp_ind_plots_diagnostic[[jj]] <- n_gg
   
   
   
   print(n_gg)
}
print(n_gg_spag)
dev.off()

# sp_ind_plots_strat[[sp]] <- tmp_sp_ind_plots
# sp_ind_plots_strat_diagnostic[[sp]] <- tmp_sp_ind_plots_diagnostic


# Trajectory Overplot log-scale -------------------------------------------

indices_strat_smooth <- indices_strat %>% filter(parm == "nsmooth")



n_gg_over = ggplot(data = indices_strat_smooth,aes(x = year, y = PI50,colour = hex_name))+
   #geom_ribbon(aes(ymin = PI2_5,ymax = PI97_5),alpha = 0.2)+
   geom_line(alpha = 0.8)+
   scale_colour_viridis_d()+
   geom_line(data = indicesNSmooth,aes(x = year,y = PI50),inherit.aes = FALSE)+
   theme_classic()+
   labs(title = sp)+
   theme(legend.position = "none")+
   scale_y_log10()

#traj_overplots[[jj]] <- n_gg_over



# Beta plots --------------------------------------------------------------

B_samples <- slope_icar_stanfit %>% gather_draws(B[k])    
b_samples <- slope_icar_stanfit %>% gather_draws(b[s,k])    
b_samples <- left_join(b_samples,strats_dts,by = c("s" = "stratn"))

B <- B_samples %>% group_by(k) %>% 
   summarise(mean = mean(.value),
             lci = quantile(.value,0.05),
             uci = quantile(.value,0.95))
b <- b_samples %>% group_by(k,hex_name) %>% 
   summarise(mean = mean(.value),
             lci = quantile(.value,0.025),
             uci = quantile(.value,0.975))

sdbeta_samples <- slope_icar_stanfit %>% gather_draws(sdyear_gam_strat[k])
sdbeta <- sdbeta_samples %>% group_by(k) %>% 
   summarise(mean = mean(.value),
             lci = quantile(.value,0.025),
             uci = quantile(.value,0.975)) %>% 
   mutate(k = k+0.1)
#B_b <- bind_rows(B,b)

b_plot <- ggplot(data = b,aes(x = k,y = mean,colour = hex_name))+
   geom_point(data = B, inherit.aes = FALSE, 
              aes(x = k,y = mean),size = 1)+
   geom_errorbar(data = B, inherit.aes = FALSE, 
                 aes(x = k,y = mean,ymin = lci,ymax = uci),width = 0,alpha = 0.5)+
   geom_point()+
   theme_classic()+
   labs(title = sp)+
   geom_pointrange(data = sdbeta,inherit.aes = FALSE,aes(x = k,y = mean,ymin = lci,ymax = uci),colour = "red")+
   scale_colour_viridis_d()+
   theme(legend.position = "none")
print(b_plot)


#beta_overplots[[jj]] <- b_plot

# combine the trajectories within the original regional strata ------------




# map the trend estimates -------------------------------------------------

# real_grid


# reg_strats <- data.frame(hex_name = real_grid_regs$hex_name,
#                          Region = real_grid_regs$Region)
# 
strats_dts <- left_join(strats_dts,strat_regions,by = c("stratn" = "strat"))


nsmoothsamples <- slope_icar_stanfit %>% gather_draws(nsmooth[s,y])
nsmoothsamples$year <- nsmoothsamples$y + (syear-1)
nsmoothsamples <- left_join(nsmoothsamples,strats_dts,by = c("s" = "stratn"))
# st_drop <- "2305040_1068494"
#nsmoothsamples <- nsmoothsamples %>% filter(hex_name != st_drop)


t_nsmooth_strat_04 <- ItoT(inds = nsmoothsamples,
                           start = 2004,
                           end = 2019,
                           regions = "hex_name",
                           qs = 95,
                           sp = sp,
                           type = "15-year")

#trendsout <- bind_rows(trendsout,
#                       t_nsmooth_strat_04)


t_nsmooth_strat_80 <- ItoT(inds = nsmoothsamples,
                           start = 1980,
                           end = 2019,
                           regions = "hex_name",
                           qs = 95,
                           sp = sp,
                           type = "Long-term")
#trendsout <- bind_rows(trendsout,
#                       t_nsmooth_strat_80)



t_nsmooth_strat_3g <- ItoT(inds = nsmoothsamples,
                           start = y3g,
                           end = 2019,
                           regions = "hex_name",
                           qs = 95,
                           sp = sp,
                           type = "Recent-three-generation")
# trendsout <- bind_rows(trendsout,
#                       t_nsmooth_strat_3g)

t_nsmooth_strat_L3g <- ItoT(inds = nsmoothsamples,
                            start = syL3g,
                            end = y3g,
                            regions = "hex_name",
                            qs = 95,
                            sp = sp,
                            type = "Previous-three-generation")
#trendsout <- bind_rows(trendsout,
#                       t_nsmooth_strat_L3g)




t_nsmooth_reg_04 <- ItoT(inds = nsmoothsamples,
                         start = 2004,
                         end = 2019,
                         regions = "Region",
                         qs = 95,
                         sp = sp,
                         type = "15-year")
#trendsout <- bind_rows(trendsout,
#                       t_nsmooth_reg_04)

t_nsmooth_reg_80 <- ItoT(inds = nsmoothsamples,
                         start = 1980,
                         end = 2019,
                         regions = "Region",
                         qs = 95,
                         sp = sp,
                         type = "Long-term")
#trendsout <- bind_rows(trendsout,
#                       t_nsmooth_reg_80)

t_nsmooth_reg_3g <- ItoT(inds = nsmoothsamples,
                         start = y3g,
                         end = 2019,
                         regions = "Region",
                         qs = 95,
                         sp = sp,
                         type = "Recent-three-generation")
#trendsout <- bind_rows(trendsout,
#                       t_nsmooth_reg_3g)

t_nsmooth_reg_L3g <- ItoT(inds = nsmoothsamples,
                          start = syL3g,
                          end = y3g,
                          regions = "Region",
                          qs = 95,
                          sp = sp,
                          type = "Previous-three-generation")
#trendsout <- bind_rows(trendsout,
#                       t_nsmooth_reg_L3g)



nsamples <- slope_icar_stanfit %>% gather_draws(n[s,y])
nsamples$year <- nsamples$y + (syear-1)
nsamples <- left_join(nsamples,strats_dts,by = c("s" = "stratn"))


# st_drop <- "2305040_1068494"
# nsamples <- nsamples %>% filter(hex_name != st_drop)
# 
indsnsmooth_region = ItoI(inds = nsmoothsamples,
                          regions = "Region")
indsnsmooth_region$type = "Smooth"

indsn_region = ItoI(inds = nsamples,
                    regions = "Region")
indsn_region$type = "Full"
indsn_region <- bind_rows(indsnsmooth_region,indsn_region)

ind_fc = ggplot(data = indsn_region,aes(x = year,y = median,group = type))+
   geom_ribbon(aes(ymin = lci,ymax = uci,fill = type),alpha = 0.2)+
   geom_line(aes(colour = type))+
   coord_cartesian(ylim = c(0,NA))+
   #scale_colour_viridis_d(aesthetics = c("fill","colour"))+
   my_col2_traj+
   theme(legend.position = "none")+
   theme_classic()+
   labs(title = sp)+
   facet_wrap(~region,scales = "free")
print(ind_fc)

#composite_trajectories[[sp]] <- ind_fc  


t_comparison_1 <- inner_join(t_nsmooth_strat_L3g,
                             strats_xy,
                             by = c("region" = "hex_name")) %>% 
   mutate(true_trend = (exp(b1)-1)*100)
plot(t_comparison_1$trend,t_comparison_1$true_trend)#, xlim = c(-2.5,2),ylim = c(-2.5,2))
abline(0,1)

t_comparison_2 <- inner_join(t_nsmooth_strat_3g,
                             strats_xy,
                             by = c("region" = "hex_name")) %>% 
   mutate(true_trend = (exp(b2)-1)*100)
plot(t_comparison_2$trend,t_comparison_2$true_trend)#, xlim = c(-9,-2),ylim = c(-9,-2))
#plot(t_comparison_2$trend,t_comparison_2$true_trend, xlim = rev(-1*c(-9,-2)),ylim = rev(-1*c(-9,-2)))
abline(0,1)

t_comparison_all <- inner_join(t_nsmooth_strat_80,
                               strats_xy,
                               by = c("region" = "hex_name")) %>% 
   mutate(true_trend = (exp(b1+b2)-1)*100)
plot(t_comparison_all$trend,t_comparison_all$true_trend, xlim = c(-8.5,-2),ylim = c(-8.5,-2))
plot(t_comparison_all$trend,t_comparison_all$true_trend, xlim = rev(-1*c(-9,-2)),ylim = rev(-1*c(-9,-2)))
abline(0,1)


# Trend heatmaps ----------------------------------------------------------
## consider adding stratum labels, etc.

pdf(paste0("Figures/Trend_Heat_maps_",sp,".pdf"),
    width = 11,
    height = 8.5)
t_80 = trend_map(t_nsmooth_strat_80,
                 size_value = "Mean Observed Count")
print(t_80)
t_04 = trend_map(t_nsmooth_strat_04,
                 size_value = "Mean Observed Count")
print(t_04)

t_3g = trend_map(t_nsmooth_strat_3g,
                 size_value = "Mean Observed Count")
print(t_3g)

t_L3g = trend_map(t_nsmooth_strat_L3g,
                  size_value = "Mean Observed Count")
print(t_L3g)

dev.off()

TRENDSout <- TRENDSout %>% 
   mutate(span = end_year-start_year,
          yrs = paste(start_year,end_year,sep = "-"))

tplot <- ggplot(data = TRENDSout,aes(x = yrs,y = trend,colour = end_year))+
   geom_point()+
   geom_errorbar(aes(ymin = lci,ymax = uci),alpha = 0.3,width = 0)+
   labs(title = sp)
print(tplot)


pdf(file = paste0("Figures/Trends_",sp,".pdf"),
    width = 8.5,
    height = 11)
print(tplot)
dev.off()






 

