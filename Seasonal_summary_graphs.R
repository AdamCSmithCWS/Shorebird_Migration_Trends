library(tidyverse)

library(posterior)
library(cmdstanr)
source("functions/utility_functions.R")
source("functions/posterior_summary_functions.R")
source("Functions/palettes.R")

load("data/allShorebirdPrismFallCounts.RData")


#lists for stored figures

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

season_out<- NULL
seasonEffect_out <- NULL


# Species loop ------------------------------------------------------------
output_dir <- "g:/Shorebird_Migration_Trends/output"


for(sp in sps){
  #if(sp == "Semipalmated Sandpiper"){next}
  spf = gsub(sp,pattern = " ",replacement = "_")
  spf = gsub(pattern = "\'",replacement = "",
             x = spf)
  
  load(paste0("data/data",sp,"_cmdstanr_data.RData"))
  
  noise_dist_sel <- noise_dist2
  
  sp_file_name <- paste0(spf,"-",prior,"-",noise_dist_sel)

  if(file.exists(paste0(output_dir,"/",sp_file_name,"_fit_add.RData"))){
    load(paste0(output_dir,"/",sp_file_name,"_fit_add.RData"))
    if(length(csvfl) > 4){csvfl <- csvfl[1:4]}
    cmdstanfit <- as_cmdstan_fit(csvfl) 
    #<- readRDS(paste0(output_dir,"/",sp_file_name,".RDS"))
    #load(paste0(output_dir,"/",spf,"_fit_add.RData"))
    three_gen <- max(10,ceiling(gens[which(gens$Common_name == sp),"GenLength"]*3))
    #Three generation assessment time in COSEWIC report
    y3g <- 2019-three_gen
    
   
    syear = min(dts$YearCollected)
    
 
    alpha_samples <- posterior_samples(fit = cmdstanfit,
                                       parm = "alpha",
                                       dims = c("s"))

    sdalpha_samples <- posterior_samples(fit = cmdstanfit,
                                         parm = "sdalpha")
    sdalpha = sdalpha_samples %>%
      summarise(mean = mean((.value)),
                lci = quantile((.value),0.025),
                uci = quantile((.value),0.975))


    sites_strat = (stan_data$sites)
    nstrata = stan_data$nstrata
    nsites_strat = stan_data$nsites_strat
    # #offsets = stan_data$site_size
    #
    sitesbystrat = NULL
    for(st in 1:stan_data$nstrata){
      tmp = data.frame(strat = st,
                       site = sites_strat[1:nsites_strat[st],st])
      sitesbystrat <- bind_rows(sitesbystrat,tmp)
    }
    sitesbystrat <- arrange(sitesbystrat,site)
    #
    alphas = alpha_samples %>% group_by(s) %>%
      summarise(mean = mean(exp(.value)),
                lci = quantile(exp(.value),0.025),
                uci = quantile(exp(.value),0.975)) %>%
      mutate(site = s) %>% left_join(.,sitesbystrat)
    #
    # # ## site-effects against the observed counts at each site
    # nr = ceiling(sqrt(nstrata))
    # obs_by_alpha = ggplot(data = alphas,aes(x = site,y = mean))+
    #   geom_pointrange(aes(ymin = lci,ymax = uci),colour = "red")+
    #   geom_point(data = dts, inherit.aes = FALSE,aes(x = site,y = count),
    #              fill = grDevices::grey(0.6),colour = grDevices::grey(0),alpha = 0.1,
    #              position = position_jitter(width = 0.2,height = 0))+
    #   labs(title = sp)+
    #   facet_wrap(~strat,nrow = nr,ncol = nr,scales = "free")+
    #   theme_minimal()
    # 
    # # pdf(paste0("Figures/",sp,prior,"_",noise_dist_sel,"_obs_by_alpha.pdf"),
    # #     width = 11,
    # #     height = 11)
    # # print(obs_by_alpha)
    # # dev.off()
    # 
    # # 
    # # 
    # # 
    # # # visualize the seasonal corrections --------------------------------------
    # # 
    # # 
    # # 
    # # # # extracting the seasonal smooth ------------------------------------------
    # # # 
    # # 
    # # 
    # 
    ALm <- posterior_samples(fit = cmdstanfit,
                             parm = "ALPHA1",
                             dims = NULL) %>%
      summarise(mean = mean(.value)) %>%
      as.numeric()


    sdnm <- posterior_samples(fit = cmdstanfit,
                              parm = "sdnoise",
                              dims = NULL) %>%
      summarise(mean = mean(.value)) %>%
      as.numeric()

    scale_adj <- ((sdnm^2)*0.5 + ALm)
    # 
    # 
    # # n_by_strat <- dts %>% group_by(strat) %>% 
    # #   summarise(n_cst = n())
    # # n_countsby_site <- dts %>% group_by(site,strat) %>% 
    # #   summarise(n_c = n()) %>% 
    # #   left_join(.,n_by_strat,by = "strat") %>% 
    # #   mutate(n_c = n_c/n_cst) %>% 
    # #   select(site,strat,n_c)
    # 
    # 
    n_countsby_site <- dts %>% group_by(site) %>%
      summarise(n_c = n())
    # 
    if(grepl(x = mod.file1,pattern = "two_season")){
      
      strat_season_strat <- dts %>% distinct(seas_strat,strat,hex_name)
      
      alphas <- left_join(alphas,strat_season_strat,by = "strat")
      
      
      strat_offs <- alphas %>% left_join(.,n_countsby_site,by = "site") %>% 
        group_by(strat,seas_strat) %>% 
        summarise(adjs = mean(mean))
      
      
      
      
      season_samples <- posterior_samples(fit = cmdstanfit,
                                          parm = "season_pred",
                                          dims = c("d","s"))
      
      seasonEffect2 = season_samples %>% group_by(d,s) %>% 
        summarise(mean = mean(exp(.value)),
                  lci = quantile(exp(.value),0.025),
                  uci = quantile(exp(.value),0.975)) %>% 
        mutate(day = d,
               seas_strat = s,
               species = sp) 
      
      seasonEffect_out <- bind_rows(seasonEffect_out,seasonEffect2)
      
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
               uci = uci*adjs + 1,
               species = sp)
      
      season_out <- bind_rows(season_out,seasonEffect_plot)
      
      # pp_simple <- ggplot()+
      #   geom_point(data = dts,aes(x = date,y = count+1,colour = year),alpha = 0.5,size = 1,position = position_jitter(width = 0.7,height = 0))+
      #   #geom_smooth(data = dts,aes(x = date,y = count+1))+
      #   scale_colour_viridis_c()+
      #   geom_line(data = seasonEffect_plot,aes(x = day,y = mean),inherit.aes = FALSE)+
      #   #coord_cartesian(ylim = c(0,yup))+
      #   geom_ribbon(data = seasonEffect_plot,aes(x = day,y = mean,ymax = uci,ymin = lci),alpha = 0.2,inherit.aes = FALSE)+
      #   ylab("")+
      #   xlab("Days since July 1")+
      #   scale_y_log10()+
      #   facet_wrap(facets = ~seas_strat,nrow = 2, ncol = 1,scales = "free")+
      #   labs(title = sp)
      # 
      # print(pp_simple)
      # 
      # 
      
      # out_simple_season_graphs[[sp]] <- pp_simple
      # 
      # 
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
      # 
      # tmp_season_graphs <- vector(mode = "list",length = ceiling(nstrata/ppag))
      # 
      # pdf(file = paste0("Figures/",sp,prior,"_",noise_dist_sel,"_cmd_simple_Season.pdf"),
      #     width = 8.5,
      #     height = 8.5)
      # 
      # 
      # for(jj in 1:ceiling(nstrata/ppag)){
      #   #yup <- quantile(dts$count,0.99)
      #   pp <- ggplot()+
      #     # geom_pointrange(data = obs_season,inherit.aes = FALSE,
      #     #                 aes(x = day,y = mean,ymin = lqrt,ymax = uqrt),alpha = 0.1)+
      #     geom_point(data = dts,aes(x = date,y = count,colour = year),alpha = 0.5,size = 1,position = position_jitter(width = 0.7,height = 0))+
      #     scale_colour_viridis_c()+
      #     geom_line(data = seasonEffect,aes(x = day,y = mean),inherit.aes = FALSE)+
      #     #coord_cartesian(ylim = c(0,yup))+
      #     geom_ribbon(data = seasonEffect,aes(x = day,y = mean,ymax = uci,ymin = lci),alpha = 0.2,inherit.aes = FALSE)+
      #     ylab("")+
      #     xlab("Days since July 1")+
      #     facet_wrap_paginate(facets = ~strat,page = jj,nrow = nrr, ncol = ncl,scales = "free")+
      #     labs(title = sp)
      #   tmp_season_graphs[[jj]] <- pp
      #   print(pp)
      # }
      # dev.off()
      # 
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
      
      
      season_samples <- posterior_samples(fit = cmdstanfit,
                                          parm = "season_pred",
                                          dims = c("d"))
      
      
      #season_samples <- slope_icar_stanfit %>% gather_draws(season_pred[d])
      seasonEffectT = season_samples %>% group_by(d) %>% 
        summarise(mean = mean(exp(.value+scale_adj)),
                  lci = quantile(exp(.value+scale_adj),0.025),
                  uci = quantile(exp(.value+scale_adj),0.975)) %>% 
        mutate(day = d) 
      
      
      seasonEffect2 = season_samples %>% group_by(d) %>% 
        summarise(mean = mean(exp(.value)),
                  lci = quantile(exp(.value),0.025),
                  uci = quantile(exp(.value),0.975)) %>% 
        mutate(day = d,
               species = sp) 
      
      seasonEffect_out <- bind_rows(seasonEffect_out,seasonEffect2)
      
      
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
               uci = uci*mean(alphas$mean)+1,
               species = sp)
      
      season_out <- bind_rows(season_out,seasonEffect_plot)
      # 
      # pp_simple <- ggplot()+
      #   geom_point(data = dts,aes(x = date,y = count+1,colour = year),alpha = 0.5,size = 1,position = position_jitter(width = 0.7,height = 0))+
      #   #geom_smooth(data = dts,aes(x = date,y = count+1))+
      #   scale_colour_viridis_c()+
      #   geom_line(data = seasonEffect_plot,aes(x = day,y = mean),inherit.aes = FALSE)+
      #   #coord_cartesian(ylim = c(0,yup))+
      #   geom_ribbon(data = seasonEffect_plot,aes(x = day,y = mean,ymax = uci,ymin = lci),alpha = 0.2,inherit.aes = FALSE)+
      #   ylab("")+
      #   xlab("Days since July 1")+
      #   scale_y_log10()+
      #   labs(title = sp)
      # #print(pp_simple)
      # out_simple_season_graphs[[sp]] <- pp_simple
      # 
      # pdf(file = paste0("Figures/",sp,prior,"_",noise_dist_sel,"_cmd_simple_Season_simplified.pdf"),
      #     width = 8.5,
      #     height = 8.5)
      # print(pp_simple)
      # dev.off()
      # 
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
      # tmp_season_graphs <- vector(mode = "list",length = ceiling(nstrata/ppag))
      # 
      # pdf(file = paste0("Figures/",sp,prior,"_",noise_dist_sel,"_cmd_simple_Season.pdf"),
      #     width = 8.5,
      #     height = 8.5)
      # 
      # for(jj in 1:ceiling(nstrata/ppag)){
      #   #yup <- quantile(dts$count,0.99)
      #   pp <- ggplot()+
      #     # geom_pointrange(data = obs_season,inherit.aes = FALSE,
      #     #                 aes(x = day,y = mean,ymin = lqrt,ymax = uqrt),alpha = 0.1)+
      #     geom_point(data = dts,aes(x = date,y = count,colour = year),alpha = 0.5,size = 1,position = position_jitter(width = 0.7,height = 0))+
      #     scale_colour_viridis_c()+
      #     geom_line(data = seasonEffect,aes(x = day,y = mean),inherit.aes = FALSE)+
      #     #coord_cartesian(ylim = c(0,yup))+
      #     geom_ribbon(data = seasonEffect,aes(x = day,y = mean,ymax = uci,ymin = lci),alpha = 0.2,inherit.aes = FALSE)+
      #     ylab("")+
      #     xlab("Days since July 1")+
      #     facet_wrap_paginate(facets = ~strat,page = jj,nrow = nrr, ncol = ncl,scales = "free")+
      #     labs(title = sp)
      #   tmp_season_graphs[[jj]] <- pp
      #   print(pp)
      # }
      # dev.off()
      # 
    }
    
save(list = c("seasonEffect_out",
              "season_out"),
     file = "Data/Season_effects.RData")
    
print(sp)
  }
}





# Plotting ----------------------------------------------------------------



load("Data/Season_effects.RData")

d1 <- yday(as.Date.character("2019-6-30","%Y-%m-%d"))
od <- function(x){
  y <- as.Date(paste0("2019-",x + d1),"%Y-%j")
}
season_out <- season_out %>% 
  mutate(n_season = factor(ifelse(is.na(seas_strat),"Single","Split"),
                           levels = c("Single","Split"),
                           ordered = TRUE),
         Season_Region = factor(ifelse(is.na(seas_strat)|seas_strat == 1,
                                       "Most or all of survey area",
                                       "Southern portion")),
         date = od(day)) %>% 
  arrange(n_season)
spl <- unique(season_out$species)
season_out <- season_out %>% 
  mutate(species = factor(species,levels = spl,
                          ordered = TRUE))

season_spag2 <- ggplot(data = season_out,aes(x = date,y = mean,
                                             group = Season_Region,
                                             colour = Season_Region))+
  geom_ribbon(aes(x = date,y = mean,
                  ymin = lci,ymax = uci,fill = Season_Region),
              inherit.aes = FALSE,alpha = 0.2)+
  geom_line(alpha = 0.5)+
  my_col2_traj+
  guides(fill = guide_legend(title = "Region of survey area"),
         colour = guide_legend(title = "Region of survey area"))+
  scale_y_continuous(trans = "log10")+
  scale_x_date(date_breaks = "1 month",
               #date_minor_breaks = "1 week",
               date_labels =  "%b %d")+
  theme_bw()+
  ylab("Mean count across all sites")+
  #xlab("day")+
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 7.5),
        axis.text.x=element_text(angle=45, hjust=1),
        legend.position = "top")+
  facet_wrap(vars(species),
             nrow = 5,
             ncol = 6,
             scales = "free_y")

print(season_spag2)


pdf("Figures/Seasonal_means.pdf",
    width = 10,
    height = 7)
print(season_spag2)
dev.off()





seasonEffect_out <- seasonEffect_out %>% 
  mutate(Season_Region = factor(ifelse(is.na(seas_strat)|seas_strat == 1,"Northern",
                                "Southern")))


season_spag1 <- ggplot(data = seasonEffect_out,aes(x = day,y = mean,
                                            group = Season_Region,
                                            colour = Season_Region))+
  geom_ribbon(aes(x = day,y = mean,
                  ymin = lci,ymax = uci,fill = Season_Region),
              inherit.aes = FALSE,alpha = 0.2)+
  geom_line(alpha = 0.5)+
  my_col_sim+
  scale_y_continuous(trans = "log10")+
  facet_wrap(vars(species),
             nrow = 5,
             ncol = 6,
             scales = "free_y")

print(season_spag1)







