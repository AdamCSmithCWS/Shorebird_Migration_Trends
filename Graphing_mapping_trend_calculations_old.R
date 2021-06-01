library(tidyverse)
library(sf)
library(spdep)
library(ggforce)
library(tidybayes)
library(GGally)
source("functions/utility_functions.R")

library(loo)
load("data/allShorebirdPrismFallCounts.RData")
source("Functions/palettes.R")


#lists for stored figures

blank_list <- vector(mode = "list",length = length(sps))
names(blank_list) <- sps

trend_maps_1980 <- blank_list
trend_maps_2004 <- blank_list
trend_maps_3gen <- blank_list
trend_maps_L3gen <- blank_list

composite_trajectories <- blank_list

sp_ind_plots_strat_diagnostic <- blank_list
sp_ind_plots_diagnostic <- blank_list

sp_ind_plots_strat <- blank_list
sp_ind_plots <- blank_list

season_graphs <- blank_list
loo_ic <- blank_list
sp_spag_plots_diagnostic <- blank_list
beta_overplots <- blank_list
traj_overplots <- blank_list
alternate_traj_overplots <- blank_list
out_simple_season_graphs <- blank_list
out_alphas_by_yr <- blank_list

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



grid_spacing <- 300000  # size of squares, in units of the CRS (i.e. meters for lae)


FYYYY = 1980

t1 = Sys.time()


w_cosewic = sps[c(2:4,7,10,12:20,22,11,25)]

# Species loop ------------------------------------------------------------


for(sp in sps){
  #if(sp == "Semipalmated Sandpiper"){next}
  if(file.exists(paste0("output/",sp,"_GAMYE_strat_simple",grid_spacing/1000,".RData"))){
    load(paste0("output/",sp,"_GAMYE_strat_simple",grid_spacing/1000,".RData"))
    
    three_gen <- max(10,ceiling(gens[which(gens$Common_name == sp),"GenLength"]*3))
    #Three generation assessment time in COSEWIC report
    y3g <- 2019-three_gen
      
      # CONSIDERATIONS ----------------------------------------------------------
    
    ## map the looic values to see if some areas are being poorly predicted
    
    ## add some variation to the plotting of the observed means (box and whisker info)
    
    #################################
    ## plot the changing mean alpha through time - are smaller sites being surveyed
    ## more in recent years?
    #################################
    

    ## generate some fake data to ensure there isn't a negative bias
    
    
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
    
    
    pw_loo <- loo_ic[[sp]]
    loopoint = as.data.frame(pw_loo$pointwise)
    
    dts_loo <- bind_cols(dts,loopoint) %>% 
      mutate(log_countp1 <- log(count+1))
    
    wcl <- which(names(dts_loo) %in% c("stratn","log_countp1","year","date","influence_pareto_k","looic"))
    prs_plot <- ggpairs(data = dts_loo,columns = wcl)
    
    
    
    
    loo_by_strat <- dts_loo %>% group_by(hex_name) %>% 
      summarise(mean_looic = mean(looic),
                median_looic = median(looic),
                q75_looic = as.numeric(quantile(looic,0.75)),
                mean_k = mean(influence_pareto_k),
                median_k = median(influence_pareto_k),
                q75_k = as.numeric(quantile(influence_pareto_k,0.75)),
                n_counts = n(),
                median_count = median(count),
                mean_count = mean(count))
    
    strat_pairs <- ggpairs(data = loo_by_strat,columns = 2:ncol(loo_by_strat))
    pdf(paste0("Figures/",sp,"_loo_pairs.pdf"),
        width = 11,
        height = 11)
    print(prs_plot)
    print(strat_pairs)
    dev.off()
    
    
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
    
    
 

    # visualize the seasonal corrections --------------------------------------
    
    
 
    # # extracting the seasonal smooth ------------------------------------------
    # 
    
    
    
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
    
    
    seas_strat_offs <- alphas %>% left_join(.,n_countsby_site,by = "site") %>% 
      group_by(seas_strat) %>% 
      summarise(adjs = mean(mean),
                adjs2 = median(mean))
    
    seasonEffect_plot <- left_join(seasonEffectT,seas_strat_offs,by = "seas_strat") %>% 
      mutate(mean = mean*adjs + 1,
             lci = lci*adjs + 1,
             uci = uci*adjs + 1)
    
    
    pp_simple <- ggplot()+
      geom_point(data = dts,aes(x = date,y = count+1,colour = year),alpha = 0.5,size = 1,position = position_jitter(width = 0.7,height = 0))+
      #geom_smooth(data = dts,aes(x = date,y = count+1))+
      scale_colour_viridis_c()+
      geom_line(data = seasonEffect_plot,aes(x = day,y = mean),inherit.aes = FALSE)+
      #coord_cartesian(ylim = c(0,yup))+
      geom_ribbon(data = seasonEffect_plot,aes(x = day,y = mean,ymax = uci,ymin = lci),alpha = 0.2,inherit.aes = FALSE)+
      ylab("")+
      xlab("Days since July 1")+
      scale_y_log10()+
      facet_wrap(facets = ~seas_strat,nrow = 2, ncol = 1,scales = "free")+
      labs(title = sp)
    #print(pp_simple)
    out_simple_season_graphs[[sp]] <- pp_simple
    
    
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
      seasonEffectT = season_samples %>% group_by(d) %>% 
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
      
      seasonEffect <- expand_grid(seasonEffectT,strat = c(1:nstrata))
      seasonEffect <- left_join(seasonEffect,strat_offs) %>% 
        mutate(mean = mean*adjs,
               lci = lci*adjs,
               uci = uci*adjs)
      
      
      
      # seas_strat_offs <- alphas %>% left_join(.,n_countsby_site,by = "site") %>% 
      #   group_by(seas_strat) %>% 
      #   summarise(adjs = mean(mean),
      #             adjs2 = median(mean))
      
      seasonEffect_plot <- seasonEffectT %>% 
        mutate(mean = mean*mean(alphas$mean)+1,
               lci = lci*mean(alphas$mean)+1,
               uci = uci*mean(alphas$mean)+1)
      
      
      pp_simple <- ggplot()+
        geom_point(data = dts,aes(x = date,y = count+1,colour = year),alpha = 0.5,size = 1,position = position_jitter(width = 0.7,height = 0))+
        #geom_smooth(data = dts,aes(x = date,y = count+1))+
        scale_colour_viridis_c()+
        geom_line(data = seasonEffect_plot,aes(x = day,y = mean),inherit.aes = FALSE)+
        #coord_cartesian(ylim = c(0,yup))+
        geom_ribbon(data = seasonEffect_plot,aes(x = day,y = mean,ymax = uci,ymin = lci),alpha = 0.2,inherit.aes = FALSE)+
        ylab("")+
        xlab("Days since July 1")+
        scale_y_log10()+
        labs(title = sp)
      #print(pp_simple)
      out_simple_season_graphs[[sp]] <- pp_simple
      
      
      
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
    
    

    
    # plot the changing mean site-effect over years ---------------------------
    
    
    # overall -----------------------------------------------------------------
    
    alphasE = alpha_samples %>% group_by(s) %>% 
      summarise(mean = mean((.value)),
                lci = quantile((.value),0.025),
                uci = quantile((.value),0.975)) %>% 
      mutate(site = s) %>% left_join(.,sitesbystrat)
    
    dts_alpha <- left_join(dts,alphasE,by = c("site","strat"))
    
    alphas_yr <- dts_alpha %>% group_by(YearCollected,site,strat) %>% 
      summarise(mean_alpha = mean(mean)) %>% ungroup() %>% 
      group_by(YearCollected) %>% 
      summarise(mean_alpha = mean(mean_alpha))
    
    AA_y_p <- ggplot(data = alphas_yr,aes(x = YearCollected,y = mean_alpha))+
      geom_point()+
      geom_smooth()+
      geom_abline(slope = 0,intercept = 0,colour = grey(0.3))+
      xlab("")+
      labs(title = paste(sp,"mean site-effect of included sites by year"))
    
    pdf(file = paste0("Figures/",sp,"alphas_by_yr.pdf"),
        width = 8.5,
        height = 8.5)
    print(AA_y_p)
    
    tmp_out_alphas_by_yr <- vector(mode = "list",length = 1+(ceiling(nstrata/ppag)))
    tmp_out_alphas_by_yr[[1]] <- AA_y_p
    # by strata ---------------------------------------------------------------
    
    alphas_yr_s <- dts_alpha %>% group_by(YearCollected,site,strat) %>% 
      summarise(mean_alpha = mean(mean)) %>% ungroup() %>% 
      group_by(YearCollected,strat) %>% 
      summarise(mean_alpha = mean(mean_alpha))
    
    for(jj in 1:ceiling(nstrata/ppag)){
      
      a_y_p <- ggplot(data = alphas_yr_s,aes(x = YearCollected,y = mean_alpha))+
        geom_point()+
        geom_smooth()+
        geom_abline(slope = 0,intercept = 0,colour = grey(0.3))+
        xlab("")+
        labs(title = paste(sp,"mean site-effect of included sites by year"))+
        facet_wrap_paginate(facets = ~strat,page = jj,nrow = nrr, ncol = ncl,scales = "free")
      
      print(a_y_p)
      tmp_out_alphas_by_yr[[jj+1]] <- a_y_p
    }
    
   dev.off()
    
   out_alphas_by_yr[[sp]] <- tmp_out_alphas_by_yr
    

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
    
    t_NSmooth_F15 <- ItoT(inds = NSmoothsamples,
                         start = 1980,
                         end = 1995,
                         regions = NULL,
                         qs = 95,
                         sp = sp,
                         type = "First-15-year")
    TRENDSout <- bind_rows(TRENDSout,t_NSmooth_F15)
    
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
      scale_y_continuous(limits = c(0,NA))+
      theme(legend.position = "none")#+
      #scale_size_area()
    
    
    pdf(paste0("figures/",sp,FYYYY,"_GAMYE_survey_wide_trajectory_simple",grid_spacing/1000,".pdf"))
    print(N_gg)
    dev.off()
    
    sp_ind_plots_diagnostic[[sp]] <- N_gg
    
    
    
    
    N_gg_simple = ggplot(data = indices,aes(x = year, y = PI50,fill = parm))+
      geom_ribbon(aes(ymin = PI2_5,ymax = PI97_5),alpha = 0.2)+
      geom_line(aes(colour = parm))+
      labs(title = paste(sp,"Survey-wide trajectory (full and smooth)"))+
      my_col2_traj+
      xlab("")+
      ylab("Modeled mean count")+
      theme_classic()+
      scale_y_continuous(limits = c(0,NA))+
      theme(legend.position = "none")#+
    #scale_size_area()
    
    #print(N_gg_simple)
    
    sp_ind_plots[[sp]] <- N_gg_simple
    

# Continental Smooth spaghetti plot ---------------------------------------
    indicesNSmooth$year = indicesNSmooth$year + (syear-1)
    
    set.seed(2019)
    r_draws <- sample(size = 100,x = 1:max(NSmoothsamples$.draw))
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
      labs(title = paste(sp,"Random selection of 100 posterior draws of survey-wide trajectories"),
           subtitle = "Colour of each posterior draw reflects the value in 2019, demonstrating similar smooths across draws")+
      xlab("")+
      ylab("Modeled mean count")+
      theme_classic()+
      theme(legend.position = "none") + 
      geom_line(data = spg_trajs,aes(x = year,y = .value,group = draw_f,colour = midp),alpha = 0.8,inherit.aes = FALSE,size = 1)+
      scale_color_viridis_c(aesthetics = "colour")+
      scale_y_log10()
    
    sp_spag_plots_diagnostic[[sp]] <- n_gg_spag
    
        
    #print(n_gg_spag)
    
    # pdf(paste0("figures/",sp,FYYYY,"_GAMYE_survey_wide_trajectory_simple",grid_spacing/1000,".pdf"))
    # print(N_gg_simple)
    # dev.off()
    

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
      tmp_sp_ind_plots[[jj]] <- n_gg_simple
      
      
      n_gg = ggplot(data = indices_strat,aes(x = year, y = PI50,fill = parm))+
        geom_ribbon(aes(ymin = PI2_5,ymax = PI97_5),alpha = 0.2)+
        geom_line(aes(colour = parm))+
        my_col2_traj+
        geom_point(aes(y = obsmean,size = mean_counts_incl_sites),colour = grey(0.5),alpha = 0.3)+
        theme_classic()+
        labs(title = sp)+
        facet_wrap_paginate(facets = ~stratn,page = jj,nrow = nrr, ncol = ncl,scales = "free")+
        scale_size_area()
      tmp_sp_ind_plots_diagnostic[[jj]] <- n_gg
      
      
      
      print(n_gg)
    }
    print(n_gg_spag)
    dev.off()
    
    sp_ind_plots_strat[[sp]] <- tmp_sp_ind_plots
    sp_ind_plots_strat_diagnostic[[sp]] <- tmp_sp_ind_plots_diagnostic
    

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
    
    traj_overplots[[sp]] <- n_gg_over
    


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
    #print(b_plot)
    
    
    beta_overplots[[sp]] <- b_plot
    
    

# Trajectories without any intercepts -------------------------------------

    traj_funct <- function(basis,bs){
      out_vec = as.numeric(basis %*% bs)
      return(out_vec)
    }
    yr_funct <- function(basis){
      out_vec = 1:nrow(basis)+(FYYYY-1)
      return(out_vec)
    }

    btraj <- b_samples %>% group_by(.draw,hex_name) %>% 
      summarise(traj = traj_funct(stan_data$year_basispred,.value),
                year = yr_funct(stan_data$year_basispred),
                .groups = "drop") %>%
      group_by(hex_name,year) %>% 
      summarise(mean = mean((traj)),
                lci = quantile((traj),0.025),
                uci = quantile((traj),0.975))
      
    Btraj <- B_samples %>% group_by(.draw) %>% 
      summarise(traj = traj_funct(stan_data$year_basispred,.value),
                year = yr_funct(stan_data$year_basispred),
                .groups = "drop") %>%
      group_by(year) %>% 
      summarise(mean = mean((traj)),
                lci = quantile((traj),0.025),
                uci = quantile((traj),0.975))
    
    

    
    b_plot_alt = ggplot(data = btraj,aes(x = year, y = mean))+
      geom_line(data = Btraj,aes(x = year,y = mean),inherit.aes = FALSE,size = 2)+
      geom_ribbon(data = Btraj,aes(ymin = uci,ymax = lci),alpha = 0.1)+
      geom_line(alpha = 0.8,aes(colour = hex_name))+
      scale_colour_viridis_d(aesthetics = c("colour","fill"))+
      theme_classic()+
      labs(title = sp)+
      ylab("Centered log-scale smooths")+
      #scale_y_log10()+
      theme(legend.position = "none")
    #print(b_plot_alt)
    alternate_traj_overplots[[sp]] <- b_plot_alt
    
    
    #(stan_data$year_basispred * transpose(b[s,]))
    
        
    
    # combine the trajectories within the original regional strata ------------
    
    

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
                               type = "15-year",
                               centered_trends = TRUE)
    
    trendsout <- bind_rows(trendsout,
                           t_nsmooth_strat_04)
    
    
    t_nsmooth_strat_80 <- ItoT(inds = nsmoothsamples,
                               start = 1980,
                               end = 2019,
                               regions = "hex_name",
                               qs = 95,
                               sp = sp,
                               type = "Long-term",
                               centered_trends = TRUE)
    trendsout <- bind_rows(trendsout,
                           t_nsmooth_strat_80)
    
    
    
    t_nsmooth_strat_3g <- ItoT(inds = nsmoothsamples,
                               start = y3g,
                               end = 2019,
                               regions = "hex_name",
                               qs = 95,
                               sp = sp,
                               type = "Recent-three-generation",
                               centered_trends = TRUE)
    trendsout <- bind_rows(trendsout,
                           t_nsmooth_strat_3g)
    
    t_nsmooth_strat_L3g <- ItoT(inds = nsmoothsamples,
                                start = syL3g,
                                end = y3g,
                                regions = "hex_name",
                                qs = 95,
                                sp = sp,
                                type = "Previous-three-generation",
                                centered_trends = TRUE)
    trendsout <- bind_rows(trendsout,
                           t_nsmooth_strat_L3g)
    
    
    
    
    t_nsmooth_reg_04 <- ItoT(inds = nsmoothsamples,
                             start = 2004,
                             end = 2019,
                             regions = "Region",
                             qs = 95,
                             sp = sp,
                             type = "15-year",
                             centered_trends = TRUE)
    trendsout <- bind_rows(trendsout,
                           t_nsmooth_reg_04)
    
    t_nsmooth_reg_80 <- ItoT(inds = nsmoothsamples,
                             start = 1980,
                             end = 2019,
                             regions = "Region",
                             qs = 95,
                             sp = sp,
                             type = "Long-term",
                             centered_trends = TRUE)
    trendsout <- bind_rows(trendsout,
                           t_nsmooth_reg_80)
    
    t_nsmooth_reg_3g <- ItoT(inds = nsmoothsamples,
                               start = y3g,
                               end = 2019,
                               regions = "Region",
                               qs = 95,
                               sp = sp,
                               type = "Recent-three-generation",
                             centered_trends = TRUE)
    trendsout <- bind_rows(trendsout,
                           t_nsmooth_reg_3g)
    
    t_nsmooth_reg_L3g <- ItoT(inds = nsmoothsamples,
                                start = syL3g,
                                end = y3g,
                                regions = "Region",
                                qs = 95,
                                sp = sp,
                                type = "Previous-three-generation",
                              centered_trends = TRUE)
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
    #print(ind_fc)
    
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
    
    
    
  }# end if species output data exists
}#end species loop

t2 = Sys.time()
t2-t1

write.csv(trendsout,"trends/All_region_strata_composite_trends.csv",row.names = FALSE)

TRENDSout <- TRENDSout %>% relocate(species,start_year,end_year,trend_type) %>% 
  select(-parameter)

trendsoutsplit <- trendsout %>% relocate(species,start_year,end_year,trend_type,region) %>% 
  select(-c(parameter)) %>% 
  filter(region != "Composite") %>% 
           group_by(region_type) %>% 
           group_split()


write.csv(trendsoutsplit[[1]],"trends/All_strata_level_trends.csv",row.names = FALSE)
write.csv(trendsoutsplit[[2]],"trends/All_region_level_trends.csv",row.names = FALSE)

write.csv(TRENDSout,"trends/All_survey_wide_trends.csv",row.names = FALSE)

save(list = c("trend_maps_1980",
              "trend_maps_2004",
              "trend_maps_3gen",
              "trend_maps_L3gen",
              
              "composite_trajectories",
              
              "sp_ind_plots_strat",
              "sp_ind_plots",
              "sp_ind_plots_strat_diagnostic",
              "sp_ind_plots_diagnostic",              
              "season_graphs",
              "loo_ic"),
     file = "Figures/All_stored_maps.RData")



# Printing the maps -------------------------------------------------------


# Trend maps --------------------------------------------------------------
pdf(file = "Figures/All_trend_maps.pdf",
    width = 9, height = 6.5)
for(sp in sps){
  if(!is.null(trend_maps_1980[[sp]])){
  print(trend_maps_1980[[sp]])
    print(trend_maps_2004[[sp]])
    print(trend_maps_L3gen[[sp]])
    print(trend_maps_3gen[[sp]])
      }
}
dev.off()




# Season Graphs --------------------------------------------------------------
pdf(file = "Figures/All_Season_graphs.pdf",
    width = 9, height = 6.5)
for(sp in sps){
  if(!is.null(season_graphs[[sp]])){
    print(season_graphs[[sp]])

  }
}
dev.off()


# Simple Season Graphs --------------------------------------------------------------
pdf(file = "Figures/All_Simple_Season_graphs.pdf",
    width = 9, height = 6.5)
for(sp in sps){
  if(!is.null(out_simple_season_graphs[[sp]])){
    print(out_simple_season_graphs[[sp]])
    
  }
}
dev.off()


# Continental Trajectories ------------------------------------------------

pdf(file = "Figures/All_Continental_Trajectories.pdf",
    width = 9, height = 6.5)
for(sp in sps){
  if(!is.null(sp_ind_plots[[sp]])){
    print(sp_ind_plots[[sp]])
    
  }
}
dev.off()


# Stratum-level Trajectories ------------------------------------------------

pdf(file = "Figures/All_Strata_Trajectories.pdf",
    width = 9, height = 6.5)
for(sp in sps){
  if(!is.null(sp_ind_plots_strat[[sp]])){
    tmp = sp_ind_plots_strat[[sp]]
    for(j in 1:length(tmp))
    print(tmp[[j]])
  }
}
dev.off()



# Composite Regional Trajectories --------------------------------------------------------------
pdf(file = "Figures/All_Composite_Trajectories.pdf",
    width = 9, height = 6.5)
for(sp in sps){
  if(!is.null(composite_trajectories[[sp]])){
    print(composite_trajectories[[sp]])
    
  }
}
dev.off()


# Continental Diagnostic Trajectories ------------------------------------------------

pdf(file = "Figures/All_Diagnostic_Continental_Trajectories.pdf",
    width = 9, height = 6.5)
for(sp in sps){
  if(!is.null(sp_ind_plots_diagnostic[[sp]])){
    print(sp_ind_plots_diagnostic[[sp]])
    print(sp_spag_plots_diagnostic[[sp]])
  }
}
dev.off()


# Stratum-level Diagnostic Trajectories ------------------------------------------------

pdf(file = "Figures/All_Diagnostic_Strata_Trajectories.pdf",
    width = 9, height = 6.5)
for(sp in sps){
  if(!is.null(sp_ind_plots_strat_diagnostic[[sp]])){
    tmp = sp_ind_plots_strat_diagnostic[[sp]]
    for(j in 1:length(tmp))
      print(tmp[[j]])
  }
}
dev.off()

# alphas by year ------------------------------------------------

pdf(file = "Figures/All_alphas_by_yr.pdf",
    width = 9, height = 6.5)
for(sp in sps){
  if(!is.null(out_alphas_by_yr[[sp]])){
    tmp = out_alphas_by_yr[[sp]]
    for(j in 1:length(tmp))
      print(tmp[[j]])
  }
}
dev.off()


# Trend plots -------------------------------------------------------------
TRENDSout$trend_type <- factor(TRENDSout$trend_type,ordered = TRUE,levels = c("Long-term","First-15-year",
                                                                              "15-year",
                                                                              "Previous-three-generation",
                                                                              "Recent-three-generation"))

LT_trends <- TRENDSout %>% filter(trend_type %in% c("Long-term","15-year"))

lt_tplot <- ggplot(data = LT_trends,aes(x = species,y = trend,colour = trend_type))+
  geom_pointrange(aes(ymax = uci,ymin = lci),position = position_dodge(width = 0.2))+
  geom_abline(slope = 0,intercept = 0,alpha = 0.7)+
  ylab("Trend (%/year)")+
  xlab("")+
  my_col_sim+
  theme_classic()+
  theme(legend.position = "bottom")+
  coord_flip()

pdf(file = "Figures/All_long_short_term_trends.pdf",
    height = 9,
    width = 6.5)
print(lt_tplot)
dev.off()


FL15_trends <- TRENDSout %>% filter(trend_type %in% c("First-15-year","15-year"))

FL15_tplot <- ggplot(data = FL15_trends,aes(x = species,y = trend,colour = trend_type))+
  geom_pointrange(aes(ymax = uci,ymin = lci),position = position_dodge(width = 0.2))+
  geom_abline(slope = 0,intercept = 0,alpha = 0.7)+
  ylab("Trend (%/year)")+
  xlab("")+
  my_col_sim+
  theme_classic()+
  theme(legend.position = "bottom")+
  coord_flip()

pdf(file = "Figures/First_last_15Year_trends.pdf",
    height = 9,
    width = 6.5)
print(FL15_tplot)
dev.off()



EL_trends <- TRENDSout %>% filter(trend_type %in% c("Previous-three-generation","Recent-three-generation"))

el_tplot <- ggplot(data = EL_trends,aes(x = species,y = trend,colour = trend_type))+
  geom_pointrange(aes(ymax = uci,ymin = lci),position = position_dodge(width = 0.2))+
  geom_abline(slope = 0,intercept = 0,alpha = 0.7)+
  ylab("Trend (%/year)")+
  xlab("")+
  my_col_sim+
  theme_classic()+
  theme(legend.position = "bottom")+
  coord_flip()


pdf(file = "Figures/All_early_recent_3generation_trends.pdf",
    height = 9,
    width = 6.5)
print(el_tplot)
dev.off()



LT3_trends <- TRENDSout %>% filter(trend_type %in% c("Long-term","Recent-three-generation"))

lt3_tplot <- ggplot(data = LT3_trends,aes(x = species,y = trend,colour = trend_type))+
  geom_pointrange(aes(ymax = uci,ymin = lci),position = position_dodge(width = 0.2))+
  geom_abline(slope = 0,intercept = 0,alpha = 0.7)+
  ylab("Trend (%/year)")+
  xlab("")+
  my_col_sim+
  theme_classic()+
  theme(legend.position = "bottom")+
  coord_flip()

pdf(file = "Figures/All_long_term_3Gen_trends.pdf",
    height = 9,
    width = 6.5)
print(lt3_tplot)
dev.off()



pdf(file = "Figures/All_beta_plots.pdf",
    width = 9, height = 6.5)
for(sp in sps){
  if(!is.null(sp_ind_plots_strat_diagnostic[[sp]])){
    print(beta_overplots[[sp]])
}
}
dev.off()


pdf(file = "Figures/All_trajectory_plots.pdf",
    width = 9, height = 6.5)
for(sp in sps){
  if(!is.null(sp_ind_plots_strat_diagnostic[[sp]])){
    print(traj_overplots[[sp]])
  }
}
dev.off()



