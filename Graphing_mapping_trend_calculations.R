library(tidyverse)
library(sf)
library(spdep)
library(ggforce)
library(tidybayes)
source("functions/utility_functions.R")

library(loo)
load("data/allShorebirdPrismFallCounts.RData")


#lists for stored figures

blank_list <- vector(mode = "list",length = length(sps))
names(blank_list) <- sps

trend_maps_1980 <- blank_list
trend_maps_2004 <- blank_list
trend_maps_3gen <- blank_list
trend_maps_L3gen <- blank_list

composite_trajectories <- blank_list

sp_ind_plots_strat <- blank_list
sp_ind_plots <- blank_list

season_graphs <- blank_list
loo_ic <- blank_list

trendsout <- NULL
TRENDSout <- NULL


# Loading gen times from Bird et al 2020 --------
gens = read.csv("data/cobi13486-sup-0004-tables4.csv")
fullgensnames = read.csv("data/cobi13486-sup-0001-tables1.csv")

fullgensnames <- fullgensnames %>% select(Scientific_name,Common_name)
gens <- gens %>% select(Scientific_name,
                        GenLength) %>% 
  left_join(.,fullgensnames)

gens[which(gens$Common_name == "American Golden Plover"),"Common_name"] <- "American Golden-Plover"
gens[which(gens$Common_name == "Grey Plover"),"Common_name"] <- "Black-bellied Plover"

sps[-which(sps %in% gens$Common_name)]
 
  gens <- gens %>% filter(Common_name %in% sps)



spslist <- left_join(spslist,gens,by = c("sci_name" = "Scientific_name"))




grid_spacing <- 300000  # size of squares, in units of the CRS (i.e. meters for lae)


FYYYY = 1980


for(sp in sps){
  if(file.exists(paste0("output/",sp,"_GAMYE_strat_simple",grid_spacing/1000,".RData"))){
    load(paste0("output/",sp,"_GAMYE_strat_simple",grid_spacing/1000,".RData"))
    
    three_gen <- ceiling(gens[which(gens$Common_name == sp),"GenLength"]*3)
    #Three generation assessment time in COSEWIC report
    y3g <- 2019-three_gen
      
      # CONSIDERATIONS ----------------------------------------------------------
    
    ## map the looic values to see if some areas are being poorly predicted
    
    ## add some variation to the plotting of the observed means (box and whisker info)
    
    ## assess the fit of the site adjustment model and the non-adjusted model
    
    ## export a table with an n-surveys by site by year matrix for the data holders to confirm
    ## send the model summary to Paul - seasonal-split, maps of regions, etc.
    ## rationalize the 1980 start date (see above table), plot number of surveys and sites by time
    
    # Calculate Annual indices using samples ----------------------------------
    
    syear = min(dts$YearCollected)
    
    
    Nsamples <- slope_icar_stanfit %>% gather_draws(N[y])
    Nsamples$year <- Nsamples$y + (syear-1)
    
    NSmoothsamples <- slope_icar_stanfit %>% gather_draws(NSmooth[y])
    NSmoothsamples$year <- NSmoothsamples$y + (syear-1)
    
    
    
    # looic -------------------------------------------------------------------
    
    
    loo_ic[[sp]] = loo(slope_icar_stanfit)
    
    
    
    
    # Alphas by year ----------------------------------------------------------
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
    
    pdf(paste0("Figures/",sp,"obs_by_alpha.pdf"),
        width = 11,
        height = 11)
    print(obs_by_alpha)
    dev.off()
    
    # visualize the seasonal corrections --------------------------------------
    
    # # extracting the seasonal smooth ------------------------------------------
    # 
    
    scale_adj_means <- summary(slope_icar_stanfit,pars = c("sdnoise","ALPHA1"))$summary[,1]
    
    scale_adj <- ((scale_adj_means[["sdnoise"]]^2)*0.5 + scale_adj_means[["ALPHA1"]])
    
    if(mod.file == "models/GAMYE_strata_two_season_simple.stan"){
    season_samples <- slope_icar_stanfit %>% gather_draws(season_pred[d,s])
    seasonEffect = season_samples %>% group_by(d,s) %>% 
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
    
    pp <- ggplot(data = seasonEffect,aes(x = day,y = mean))+
      geom_pointrange(data = obs_season,inherit.aes = FALSE,
                      aes(x = day,y = mean,ymin = lqrt,ymax = uqrt),alpha = 0.1)+
      geom_line()+
      geom_ribbon(aes(x = day,y = mean,ymax = uci,ymin = lci),alpha = 0.2)+
      ylab("")+
      xlab("Days since July 1")+
      labs(title = sp)+
      facet_wrap(facets = ~seas_strat,ncol = 3,scales = "free")
    
    
    }else{
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
      
      pp <- ggplot(data = seasonEffect,aes(x = day,y = mean))+
        geom_pointrange(data = obs_season,inherit.aes = FALSE,
                        aes(x = day,y = mean,ymin = lqrt,ymax = uqrt),alpha = 0.1)+
        geom_line()+
        geom_ribbon(aes(x = day,y = mean,ymax = uci,ymin = lci),alpha = 0.2)+
        ylab("")+
        xlab("Days since July 1")+
        labs(title = sp)
      
    }
    

   
    
    
    pdf(file = paste0("Figures/",sp,"simple_Season.pdf"),
        width = 8.5,
        height = 8.5)
    print(pp)
    dev.off()
    
    season_graphs[[sp]] <- pp
    
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
                         type = "Three-generation")
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
                  "% since",
                  x$start_year,
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
                              season_scale = FALSE)
    
    
    indicesNSmooth <- index_summary(parm = "NSmooth",
                                    dims = "year",
                                    site_scale = FALSE,
                                    season_scale = FALSE)
    
    
    
    
    
    
    indices = bind_rows(indicesN,indicesNSmooth)
    indices$year = indices$year + (syear-1)
    yup = max(max(indices$PI97_5)*1.5,quantile(indices$obsmean,0.7))
    
    N_gg = ggplot(data = indices,aes(x = year, y = PI50,fill = parm))+
      geom_ribbon(aes(ymin = PI2_5,ymax = PI97_5),alpha = 0.2)+
      geom_line(aes(colour = parm))+
      labs(title = paste(sp,"Survey-wide trajectory (full and smooth) with obs means"))+
      theme(legend.position = "none")+
      annotate("text", x = 1997, y = yup*0.9, label = anot_80)+
      annotate("text", x = 1997, y = yup*0.8, label = anot_90)+
      annotate("text", x = 1997, y = yup*0.7, label = anot_07)+
      coord_cartesian(ylim = c(0,yup))+
      geom_point(aes(y = obsmean,size = nsurveys),colour = grey(0.5),alpha = 0.3)+
      theme_classic()+
      scale_size_area()
    
    
    pdf(paste0("figures/",sp,FYYYY,"_GAMYE_survey_wide_trajectory_simple",grid_spacing/1000,".pdf"))
    print(N_gg)
    dev.off()
    sp_ind_plots[[sp]] <- N_gg
    
    
    nstrata = stan_data$nstrata
    indicesnsmooth <- index_summary(parm = "nsmooth",
                                    dims = c("stratn","year"),
                                    season_scale = FALSE,
                                    site_scale = FALSE)
    
    indicesn <- index_summary(parm = "n",
                              dims = c("stratn","year"),
                              season_scale = FALSE,
                              site_scale = FALSE)
    
    indices_strat = bind_rows(indicesn,indicesnsmooth)
    indices_strat$year = indices_strat$year + (syear-1)
    indices_strat <- left_join(indices_strat,strats_dts, by = "stratn")
    
    pdf(file = paste0("figures/", sp,FYYYY,"_GAMYE_Strata_trajectories_simple",grid_spacing/1000,".pdf"),
        width = 8.5,
        height = 11)
    print(N_gg)
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
    
    tmp_sp_ind_plots <- vector(mode = "list",length = ceiling(nstrata/ppag))
    
    for(jj in 1:ceiling(nstrata/ppag)){
      n_gg = ggplot(data = indices_strat,aes(x = year, y = PI50,fill = parm))+
        geom_ribbon(aes(ymin = PI2_5,ymax = PI97_5),alpha = 0.2)+
        geom_line(aes(colour = parm))+
        geom_point(aes(y = obsmean,size = nsurveys),colour = grey(0.5),alpha = 0.3)+
        theme_classic()+
        facet_wrap_paginate(facets = ~stratn,page = jj,nrow = nrr, ncol = ncl,scales = "free")+
        scale_size_area()
      tmp_sp_ind_plots[[jj]] <- n_gg
      print(n_gg)
    }
    dev.off()
    
    sp_ind_plots_strat[[sp]] <- tmp_sp_ind_plots
    
    
    
    
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
    
    trendsout <- bind_rows(trendsout,
                           t_nsmooth_strat_04)
    
    
    t_nsmooth_strat_80 <- ItoT(inds = nsmoothsamples,
                               start = 1980,
                               end = 2019,
                               regions = "hex_name",
                               qs = 95,
                               sp = sp,
                               type = "Long-term")
    trendsout <- bind_rows(trendsout,
                           t_nsmooth_strat_80)
    
    
    
    t_nsmooth_strat_3g <- ItoT(inds = nsmoothsamples,
                               start = y3g,
                               end = 2019,
                               regions = "hex_name",
                               qs = 95,
                               sp = sp,
                               type = "Three-generation")
    trendsout <- bind_rows(trendsout,
                           t_nsmooth_strat_3g)
    
    t_nsmooth_strat_L3g <- ItoT(inds = nsmoothsamples,
                                start = syL3g,
                                end = y3g,
                                regions = "hex_name",
                                qs = 95,
                                sp = sp,
                                type = "Previous-three-generation")
    trendsout <- bind_rows(trendsout,
                           t_nsmooth_strat_L3g)
    
    
    
    
    t_nsmooth_reg_04 <- ItoT(inds = nsmoothsamples,
                             start = 2004,
                             end = 2019,
                             regions = "Region",
                             qs = 95,
                             sp = sp,
                             type = "15-year")
    trendsout <- bind_rows(trendsout,
                           t_nsmooth_reg_04)
    
    t_nsmooth_reg_80 <- ItoT(inds = nsmoothsamples,
                             start = 1980,
                             end = 2019,
                             regions = "Region",
                             qs = 95,
                             sp = sp,
                             type = "Long-term")
    trendsout <- bind_rows(trendsout,
                           t_nsmooth_reg_80)
    
    t_nsmooth_reg_3g <- ItoT(inds = nsmoothsamples,
                               start = y3g,
                               end = 2019,
                               regions = "Region",
                               qs = 95,
                               sp = sp,
                               type = "Three-generation")
    trendsout <- bind_rows(trendsout,
                           t_nsmooth_reg_3g)
    
    t_nsmooth_reg_L3g <- ItoT(inds = nsmoothsamples,
                                start = syL3g,
                                end = y3g,
                                regions = "Region",
                                qs = 95,
                                sp = sp,
                                type = "Previous-three-generation")
    trendsout <- bind_rows(trendsout,
                           t_nsmooth_reg_L3g)
    
    
    
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
    
    ind_fc = ggplot(data = indsn_region,aes(x = year,y = median,group = region))+
      geom_ribbon(aes(ymin = lci,ymax = uci,fill = type),alpha = 0.2)+
      geom_line(aes(colour = type))+
      coord_cartesian(ylim = c(0,NA))+
      scale_colour_viridis_d()+
      theme(legend.position = "none")+
      facet_wrap(~region,scales = "free")
    
    composite_trajectories[[sp]] <- ind_fc  
    
    
    
    
    
    
    # Trend heatmaps ----------------------------------------------------------
    ## consider adding site-locations, stratum labels, etc.
    
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
    
    trend_maps_1980[[sp]] <- t_80
    trend_maps_2004[[sp]] <- t_04
    trend_maps_3gen[[sp]] <- t_3g
    trend_maps_L3gen[[sp]] <- t_L3g
    
    composite_trajectories
    
  }# end if species output data exists
}#end species loop


write.csv(trendsout,"trends/All_regional_trends.csv",row.names = FALSE)
write.csv(TRENDSout,"trends/All_survey_wide_trends.csv",row.names = FALSE)

save(list = c("trend_maps_1980",
              "trend_maps_2004",
              "trend_maps_3gen",
              "trend_maps_L3gen",
              
              "composite_trajectories",
              
              "sp_ind_plots_strat",
              "sp_ind_plots",
              
              "season_graphs",
              "loo_ic"),
     file = "Figures/All_stored_maps.RData")

