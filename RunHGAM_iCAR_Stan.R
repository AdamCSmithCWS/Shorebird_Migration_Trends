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

 
 for(sp in sps){
  
   
    load(paste0("data/data",sp,"_GAMYE_strat_offset_",grid_spacing/1000,".RData"))



## compile model
slope_icar_model = stan_model(file=mod.file)

print(sp)
## run sampler on model, data
slope_icar_stanfit <- sampling(slope_icar_model,
                               data=stan_data,
                               verbose=TRUE, refresh=100,
                               chains=4, iter=1500,
                               warmup=1000,
                               cores = 4,
                               pars = parms,
                               control = list(adapt_delta = 0.85,
                                              max_treedepth = 14))





save(list = c("slope_icar_stanfit",
              "stan_data",
              "dts",
              "real_grid",
              "strats_dts",
              "strat_regions",
              "mod.file",
              "parms"),
     file = paste0("output/",sp,"_GAMYE_strat_offset",grid_spacing/1000,".RData"))


}#end modeling loop



















launch_shinystan(slope_icar_stanfit) 

load(paste0(sp,"_GAMYE_strat_offset",grid_spacing/1000,".RData"))


# CONSIDERATIONS ----------------------------------------------------------

## map the looic values to see if some areas are being poorly predicted

## add some variation to the plotting of the observed means (box and whisker info)

## assess the fit of the site adjuestment model and the non-adjusted model



# Calculate Annual indices using samples ----------------------------------


source("functions/utility_functions.R")


Nsamples <- slope_icar_stanfit %>% gather_draws(N[y])
Nsamples$year <- Nsamples$y + (syear-1)

NSmoothsamples <- slope_icar_stanfit %>% gather_draws(NSmooth[y])
NSmoothsamples$year <- NSmoothsamples$y + (syear-1)



# calculate trends continent --------------------------------------------------------


t_NSmooth_90 <- ItoT(inds = NSmoothsamples,
                  start = 1990,
                  end = 2019,
                  regions = NULL,
                  qs = 95,
                  trend_type = "endpoint",
                  index_type = "smooth",
                  retransformation_type = "standard")
t_NSmooth_80 <- ItoT(inds = NSmoothsamples,
                  start = 1980,
                  end = 2019,
                  regions = NULL,
                  qs = 95,
                  trend_type = "endpoint",
                  index_type = "smooth",
                  retransformation_type = "standard")

t_NSmooth_07 <- ItoT(inds = NSmoothsamples,
                        start = 2007,
                        end = 2019,
                        regions = NULL,
                        qs = 95,
                        index_type = "smooth")

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

anot_90 = anot_funct(t_NSmooth_90)
anot_80 = anot_funct(t_NSmooth_80)
anot_07 = anot_funct(t_NSmooth_07)


# Extract Annual Indices ------------------------------------------------


indicesN <- index_summary(parm = "N",
                          dims = "year",
                          site_scale = FALSE,
                          season_scale = TRUE)


indicesNSmooth <- index_summary(parm = "NSmooth",
                                dims = "year",
                                site_scale = FALSE,
                                season_scale = TRUE)





  
  indices = bind_rows(indicesN,indicesNSmooth)
indices$year = indices$year + (syear-1)
yup = max(max(indices$PI97_5)*1.5,quantile(indices$obsmean,0.95))

N_gg = ggplot(data = indices,aes(x = year, y = PI50,fill = parm))+
  geom_ribbon(aes(ymin = PI2_5,ymax = PI97_5),alpha = 0.2)+
  geom_line(aes(colour = parm))+
  labs(title = paste(sp,"GAMYE continental traj (full and smooth) with obs means"))+
  theme(legend.position = "none")+
  annotate("text", x = 1997, y = yup*0.9, label = anot_80)+
  annotate("text", x = 1997, y = yup*0.8, label = anot_90)+
  annotate("text", x = 1997, y = yup*0.7, label = anot_07)+
  coord_cartesian(ylim = c(0,yup))+
  geom_point(aes(y = obsmean),colour = grey(0.5),alpha = 0.3)


pdf(paste0("figures/",sp,FYYYY,"_GAMYE_survey_wide_trajectory",grid_spacing/1000,".pdf"))
print(N_gg)
dev.off()

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

pdf(file = paste0("figures/", sp,FYYYY,"_GAMYE_Strata_trajectories",grid_spacing/1000,".pdf"),
    width = 8.5,
    height = 11)
print(N_gg)
for(jj in 1:ceiling(nstrata/12)){
n_gg = ggplot(data = indices_strat,aes(x = year, y = PI50,fill = parm))+
  geom_ribbon(aes(ymin = PI2_5,ymax = PI97_5),alpha = 0.2)+
  geom_line(aes(colour = parm))+
  geom_point(aes(y = obsmean),colour = grey(0.5),alpha = 0.1)+
  facet_wrap_paginate(facets = ~hex_name,page = jj,nrow = 4, ncol = 3,scales = "free")
print(n_gg)
}
  dev.off()





# combine the trajectories within the original regional strata ------------




# map the trend estimates -------------------------------------------------

  # real_grid
  

  reg_strats <- data.frame(hex_name = real_grid_regs$hex_name,
                           Region = real_grid_regs$Region)
  
 strats_dts <- left_join(strats_dts,reg_strats,by = "hex_name")
  
 
  nsmoothsamples <- slope_icar_stanfit %>% gather_draws(nsmooth[s,y])
  nsmoothsamples$year <- nsmoothsamples$y + (syear-1)
  nsmoothsamples <- left_join(nsmoothsamples,strats_dts,by = c("s" = "stratn"))
 # st_drop <- "2305040_1068494"
  #nsmoothsamples <- nsmoothsamples %>% filter(hex_name != st_drop)
  
  
  t_nsmooth_strat_07 <- ItoT(inds = nsmoothsamples,
                       start = 2007,
                       end = 2019,
                       regions = "hex_name",
                       qs = 95,
                       index_type = "smooth")
  
  t_nsmooth_reg_07 <- ItoT(inds = nsmoothsamples,
                       start = 2007,
                       end = 2019,
                       regions = "Region",
                       qs = 95,
                       index_type = "smooth")
  
  t_nsmooth_reg_80 <- ItoT(inds = nsmoothsamples,
                           start = 1980,
                           end = 2019,
                           regions = "Region",
                           qs = 95,
                           index_type = "smooth")
  
  
  t_nsmooth_reg_00 <- ItoT(inds = nsmoothsamples,
                           start = 2000,
                           end = 2019,
                           regions = "Region",
                           qs = 95,
                           index_type = "smooth")
  
  
  
  
  
  nsamples <- slope_icar_stanfit %>% gather_draws(n[s,y])
  nsamples$year <- nsamples$y + (syear-1)
  nsamples <- left_join(nsamples,strats_dts,by = c("s" = "stratn"))
  
  st_drop <- "2305040_1068494"
  nsamples <- nsamples %>% filter(hex_name != st_drop)
  
  indsn_region = ItoI(inds = nsmoothsamples,
               regions = "Region")
  
  ind_fc = ggplot(data = indsn_region,aes(x = year,y = median,group = region))+
    geom_ribbon(aes(ymin = lci,ymax = uci),alpha = 0.2)+
    geom_line()+
    coord_cartesian(ylim = c(0,NA))+
    facet_wrap(~region,scales = "free")
  
print(ind_fc)  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
slopes = as.data.frame(summary(slope_icar_stanfit,
                 pars = c("b"),
                 probs = c(0.025,0.5,0.975))$summary)

myrename = function(fit){
  rename_with(fit,~ paste0("PI",gsub(".","_",gsub("%", "", .x, fixed = TRUE), fixed = TRUE)),ends_with("%"))

}
slopes = myrename(slopes)
slopes$stratn = 1:nrow(slopes)
slopes$trend = (exp(slopes$mean)-1)*100

breaks <- c(-7, -4, -2, -1, -0.5, 0.5, 1, 2, 4, 7)
labls = c(paste0("< ",breaks[1]),paste0(breaks[-c(length(breaks))],":", breaks[-c(1)]),paste0("> ",breaks[length(breaks)]))
labls = paste0(labls, " %")

slopes$Tplot <- cut(slopes$trend,breaks = c(-Inf, breaks, Inf),labels = labls)

map_palette <- c("#a50026", "#d73027", "#f46d43", "#fdae61", "#fee090", "#ffffbf",
                 "#e0f3f8", "#abd9e9", "#74add1", "#4575b4", "#313695")
names(map_palette) <- labls

hex_map <- left_join(real_grid,strats_dts)

hex_map <- left_join(hex_map,slopes,by = "stratn")

trend_map <- ggplot()+
  geom_sf(data = prov_state,colour = grey(0.8))+
  geom_sf(data = hex_map,aes(fill = Tplot,colour = Tplot))+
  coord_sf(ylim = c(-1509123,1890567),expand = TRUE)+
  scale_colour_manual(values = map_palette, aesthetics = c("fill","colour"),
                      guide = ggplot2::guide_legend(reverse=TRUE),
                      name = "Trends")#paste0("Trend\n",fyr,"-",lyr))
print(trend_map)





save(list = c("stan_data",
              #"basis_season",
              "dts",
              "slope_icar_model",
              "hex_map",
              "strats_dts"),
     file = paste0("output/",sp,"slope_iCAR_results.RData"))





}#end species loops


stopCluster(cl = cluster)


# 
# 
# 
# # Plotting ----------------------------------------------------------------
# 
# library(tidybayes)
# 
# load("data/allShorebirdPrismFallCounts.RData")
# source("functions/Utility_functions.R")
# 
# for(sp in sps){
#   
#   if(file.exists(paste0("output/",sp,"slope_results.RData"))){
#   
#   load(paste0("output/",sp,"slope_ZIP_results.RData"))
#   
#     fyear = (min(dts$YearCollected))
#     
#   strats = unique(dts[,c("strat","Region")])
#   strats = rename(strats,s = strat)
#   
#   
#   
#   sums = data.frame(out2$summary)
#   names(sums) <- c("mean","sd","lci","lqrt","median","uqrt","uci","Rhat","n.eff","overlap0","f")
#   sums$Parameter = row.names(sums)
#   
#   # compiling indices -------------------------------------
#   
#   n_inds <- extr_inds(param = "n_s")
#   N_inds <- extr_inds(param = "N",regions = FALSE)
#  
#  
#   
#   sdnoise_st = extr_sum(param = "sdnoise",
#                     index = c("s"),
#                     log_retrans = F) 
#   sdnoise_st <- left_join(sdnoise_st,strats,by = "s")
#   
#   nu_st = extr_sum(param = "nu",
#                         index = c("s"),
#                         log_retrans = F) 
#   nu_st <- left_join(nu_st,strats,by = "s")
#   
#   
#   
#   psi_st = extr_sum(param = "psi",
#                        index = c("j","s"),
#                        log_retrans = F) 
#  
#   psiSamples <- out2$samples %>% gather_draws(psi)
#   sdnoiseSamples <- out2$samples %>% gather_draws(sdnoise[s])
#   
# # extracting the seasonal smooth ------------------------------------------
# 
#   season_sm = extr_sum(param = "vis.sm_season",
#                        index = c("day","s"),
#                        log_retrans = TRUE) 
#   
#   season_sm <- left_join(season_sm,strats)
#   pp <- ggplot()+
#     geom_line(data = season_sm,aes(x = day,y = mean,colour = Region,group = Region))+
#     geom_ribbon(data = season_sm,aes(x = day,ymax = uci,ymin = lci),alpha = 0.2)+
#     ylab("")+
#     xlab("Days since July 1")+
#     facet_wrap(facets = ~Region,ncol = 3,scales = "free")
#   
#   
#   
#   pdf(file = paste0("Figures/",sp,"_Season_slope.pdf"),
#       width = 8.5,
#       height = 11)
#   print(pp)
#   dev.off()
#   
#   # calculating trends  -----------------------------------------------------
#   
#   NSamples <- out2$samples %>% gather_draws(N[y])
#   NSamples$year <- NSamples$y + fyear-1
#   
# 
#   
#   N_compSamples <- out2$samples %>% gather_draws(N_comp[y])
#   N_compSamples$year <- N_compSamples$y + fyear-1
#   
#   n_sSamples <- out2$samples %>% gather_draws(n_s[s,y])
#   n_sSamples$year <- n_sSamples$y + fyear-1
#   n_sSamples <- left_join(n_sSamples,strats,by = "s")
#   
#  
# 
#   t_n_s <- ItoT(inds = n_sSamples,regions = TRUE)
# 
#   
#   
#   
#   t_n_sS <- ItoT_slope(inds = n_sSamples,regions = TRUE)
# 
#   
#   # t_n_s_a1 <- ItoT(inds = n_s_a1Samples,regions = TRUE,retransformation_type = "lognormal_only")
#   # 
#   # 
#   # t_n_s_a2 <- ItoT(inds = n_s_a2Samples,regions = TRUE,retransformation_type = "none")
#   # 
#   # 
#   # 
#   t_N <- ItoT(inds = NSamples,regions = FALSE)
# 
#   
#   t_n_s_15 <- ItoT(inds = n_sSamples,regions = TRUE,start= 2004)
# 
#   
#   t_n_sS_15 <- ItoT_slope(inds = n_sSamples,regions = TRUE,start= 2004)
# 
#   
#   t_N_15 <- ItoT(inds = NSamples,regions = FALSE,start= 2004)
# 
#   t_NS <- ItoT_slope(inds = NSamples,regions = FALSE)
#   t_NS_15 <- ItoT_slope(inds = NSamples,regions = FALSE,start= 2004)
#   
#   
#   
#   t_N_comp <- ItoT(inds = N_compSamples,regions = FALSE,start = 1978)
#   
#   t_N_comp_15 <- ItoT(inds = N_compSamples,regions = FALSE,start= 2004)
#   
#   t_N_compS <- ItoT_slope(inds = N_compSamples,regions = FALSE)
#   t_N_compS_15 <- ItoT_slope(inds = N_compSamples,regions = FALSE,start= 2004)
#   
#   
# 
#   trend_out <- bind_rows(t_N,
#                          t_N_15,
#                          t_NS,
#                          t_NS_15,
#                          t_N_comp,
#                          t_N_comp_15,
#                          t_N_compS,
#                          t_N_compS_15,
#                          t_n_s,
#                          t_n_sS,
#                          t_n_s_15,
#                          t_n_sS_15)
#   
#   write.csv(trend_out,file = paste0("Trends/trends_slope_",sp,".csv"),row.names = F)
#   
#   # plotting indices --------------------------------------------------------
#   
#   
#   # plot_Hyper <- plot_ind(inds = N_inds,
#   #                        #smooth_inds = ,
#   #                        raw = dts,
#   #                        add_observed = TRUE,
#   #                        add_samplesize = TRUE,
#   #                        species = sp,
#   #                        regions = FALSE,
#   #                        title_size = 20,
#   #                        axis_title_size = 18,
#   #                        axis_text_size = 16)  
#   # 
#   # pdf(file = paste0("Figures/",sp,"_Hyperparameter_slope.pdf"),
#   #     width = 8.5,
#   #     height = 11)
#   # print(plot_Hyper)
#   # dev.off()
#   
#   
#   
#   
# #   plot_by_st <- plot_ind(inds = n_inds_a2,
# #                          #smooth_inds = ,
# #                          raw = dts,
# #                          add_observed = TRUE,
# #                          add_samplesize = TRUE,
# #                          species = sp,
# #                          regions = TRUE,
# #                          title_size = 20,
# #                          axis_title_size = 18,
# #                          axis_text_size = 16)  
# #   
# #   pdf(file = paste0("Figures/",sp,"_A2_slope.pdf"),
# #       width = 8.5,
# #       height = 11)
# # print(plot_by_st)
# # dev.off()
# 
# 
# plot_by_st <- plot_ind(inds = n_inds,
#                        smooth_inds = NULL,
#                        raw = dts,
#                        add_observed = TRUE,
#                        add_samplesize = TRUE,
#                        species = sp,
#                        regions = TRUE,
#                        title_size = 20,
#                        axis_title_size = 18,
#                        axis_text_size = 16)  
# 
# pdf(file = paste0("Figures/",sp,"_A1_slope.pdf"),
#     width = 8.5,
#     height = 11)
# print(plot_by_st)
# dev.off()
# 
# 
# 
# 
# }# end if output exists
# 
# }

