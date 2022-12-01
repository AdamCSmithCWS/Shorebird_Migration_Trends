###
# this script generates simulated data for each species,
# uses the real survey events for each species (strata, sites, dates, etc.)
# it also uses the estimated long-term trend from the real data to generated the simulated data
# it assumes that the true data-generating process includes a constant long-term trend
# and that the regional trends vary by latitude, with more negative trends in teh northeast
# and more positive trends in the south and west. 
# this particular simulation of a constant long-term trend and the spatial variation
# was designed specifically to test the rigor of the conclusion that shorebird population declines
# are accelerating.
# after generating the data, the script fits the model to the simulated data
# then compares the estimated trends to the true trends 

##### this script generates information that is then used to produce one of the 
##### supplements for the Smith et al. paper.

# SPECIES MCMC data-prep -------------------------------------------------------
library(tidyverse)
library(cmdstanr)
library(shinystan)
library(sf)
library(spdep)
library(ggforce)
library(tidybayes)

source("functions/utility_functions.R")
source("functions/posterior_summary_functions.R")
source("Functions/palettes.R")
source("functions/GAM_basis_function_mgcv.R")




sps <- readRDS("data/species_vector.rds")
grid_spacing <- 300000  # size of squares, in units of the CRS (i.e. meters for lae)


FYYYY = 1980

output_dir <- "output"
#output_dir2 <- "g:/Shorebird_Migration_Trends/output"

# Loading gen times from Bird et al 2020 supplementary material --------
gens = read.csv("data/cobi13486-sup-0004-tables4.csv")
fullgensnames = read.csv("data/cobi13486-sup-0001-tables1.csv")

fullgensnames <- fullgensnames %>% select(Scientific_name,Common_name)
gens <- gens %>% select(Scientific_name,
                        GenLength) %>% 
   left_join(.,fullgensnames)

gens[which(gens$Common_name == "American Golden Plover"),"Common_name"] <- "American Golden-Plover"
gens[which(gens$Common_name == "Grey Plover"),"Common_name"] <- "Black-bellied Plover"

#looking for missing species
sps[-which(sps %in% gens$Common_name)]

gens <- gens %>% filter(Common_name %in% sps)


balanced_sim = FALSE # script includes an option to generate a fake dataset that does not have temporal bias in the sampling at each site (which we have not highlighted because it is not a sufficiently challenging test of the model). The results from this balanced version are the same as for the version that uses the real survey distributions for each species. 

 for(sp1 in sps){

 
    spf1 = gsub(sp1,pattern = " ",replacement = "_")
    spf1 = gsub(pattern = "\'",replacement = "",
               x = spf1)
    sp_file_name1 <- paste0(spf1,"-gamma-t")
    
    load(paste0("data/data",sp1,"_cmdstanr_data.RData"))
    load(paste0(output_dir,"/",sp_file_name1,"_fit_add.RData"))
    
    cmdstanfit <- rstan::read_stan_csv(csvfl) 


# Gather estimates from fitted model --------------------------------------

    alpha <- posterior_samples(fit = cmdstanfit,
                                       parm = "alpha",
                                       dims = c("s")) %>% 
       group_by(s) %>%
       summarise(mean = mean((.value))) %>% 
       ungroup() %>% 
       select(mean) %>% 
       unlist() %>% 
       as.numeric()
    
    ALPHA1 <- posterior_samples(fit = cmdstanfit,
                             parm = "ALPHA1",
                             dims = NULL) %>% 
       summarise(mean = mean(.value)) %>% 
       as.numeric()
    

    
    sdnoise <- posterior_samples(fit = cmdstanfit,
                                parm = "sdnoise",
                                dims = NULL) %>% 
       summarise(mean = mean(.value)) %>% 
       as.numeric()
    
    noise = rnorm(stan_data$ncounts,0,sdnoise) 
    # noise2 = rt(stan_data$ncounts,20)*sdnoise
    # 
# Seasonal effects --------------------------------------------------------

    
    if(grepl(x = mod.file2,pattern = "two_season")){
       season_pred_df <- posterior_samples(fit = cmdstanfit,
                                           parm = "season_pred",
                                           dims = c("d","s"))%>% 
          group_by(d,s) %>%
          summarise(mean = mean((.value)))
       
       tmps1 <- season_pred_df %>% filter(s == 1) %>% 
          ungroup() %>% 
          select(mean) %>% 
          unlist() %>% 
          as.numeric()
       tmps2 <- season_pred_df %>% filter(s == 2) %>% 
          ungroup() %>% 
          select(mean) %>% 
          unlist() %>% 
          as.numeric()
       
       season_pred <- matrix(c(tmps1,tmps2),ncol = 2)
       
    }else{
       season_pred_df <- posterior_samples(fit = cmdstanfit,
                                        parm = "season_pred",
                                        dims = c("d"))%>% 
          group_by(d) %>%
          summarise(mean = mean((.value)))
       
       season_pred <- season_pred_df %>% 
          ungroup() %>% 
          select(mean) %>% 
          unlist() %>% 
          as.numeric()
       
    }
    
    

# fake linear trend to contrast with increasing rate of decline -----------

    year_pred <- matrix(NA,nrow = stan_data$nyears,ncol = stan_data$nstrata)
    
    

# use the mean estimated survey wide trend as the base for this sp --------
## then build a trend model that assumes a constant, log-linear slope at that rate
    NSmoothsamples <- posterior_samples(fit = cmdstanfit,
                                        parm = "NSmooth",
                                        dims = c("y"))
    NSmoothsamples$year <- NSmoothsamples$y + (1980-1)
    
    
    t_NSmooth_80 <- ItoT(inds = NSmoothsamples,
                         start = 1980,
                         end = 2019,
                         regions = NULL,
                         qs = 95,
                         sp = sp,
                         type = "Long-term")
    
    
    B1 = 0.5*(t_NSmooth_80$trend/100)  #linear trend at center of range = 2*B1
    bx = -0.005 #implies a 1-2%/year range in trends by longitude
    by = -0.02 #implies a 2-4%/year range in trends by latitude
    
    strats_xy <- strats_dts %>% 
       mutate(x = as.numeric(str_split(hex_name,pattern = "_",simplify = TRUE)[,1]),
              y = as.numeric(str_split(hex_name,pattern = "_",simplify = TRUE)[,2]),
              x_scale = (x-mean(x))/(0.5*diff(range(x))),#scaled x coordinate to create a longitudinal gradient in trends
              y_scale = (y-mean(y))/(0.5*diff(range(y))),#scaled y coordinate to create a latitudinal gradient in trends
              b1 = (x_scale*bx+B1) + (y_scale*by+B1)) #initial log-linear slopes for each stratum with spatial gradients
    
    midyear = 20
    b1 = strats_xy$b1 #initial slopes for each stratum

    #generating a smooth log-linear trajectory for each stratum
    for(s in 1:stan_data$nstrata){
       year_pred[1:midyear,s] <- b1[s]*((1:midyear)-midyear)
       year_pred[(midyear+1):stan_data$nyears,s] <- b1[s]*((midyear+1):stan_data$nyears-midyear)
    }
    
    
    
# fake - evenly distributed and oversimplified year-effects ----------------------------------
      # the even distribution and simple, systematic year-effects
      # ensure that there's no accidental pattern in the year-effects that might
      # induce bias in estimated smooth component that tracks the species trend
    
    sdyear <- posterior_samples(fit = cmdstanfit,
                                parm = "sdyear",
                                dims = NULL) %>% 
       summarise(mean = mean(.value)) %>% 
       as.numeric()
    
    year_effect <- rep(c(-1,1),length = stan_data$nyears)*(sdyear/2)
    year_effect[c(1,stan_data$nyears)] <- 0
    
    
    
    

# Generate fake counts using above estimates ------------------------------


    stan_data_sim = stan_data
    
    if(grepl(x = mod.file2,pattern = "two_season")){
    for(i in 1:stan_data_sim$ncounts){
       stan_data_sim$count[i] = rpois(n = 1,exp(ALPHA1 + year_pred[stan_data_sim$year_raw[i],stan_data_sim$strat[i]] + 
                                                   alpha[stan_data_sim$site[i]] + year_effect[stan_data_sim$year_raw[i]] + 
                                                   season_pred[stan_data_sim$date[i],stan_data_sim$seas_strat[i]] + noise[i]));
       
    }
    }else{
       for(i in 1:stan_data_sim$ncounts){
          stan_data_sim$count[i] = rpois(n = 1,exp(ALPHA1 + year_pred[stan_data_sim$year_raw[i],stan_data_sim$strat[i]] + 
                                                   alpha[stan_data_sim$site[i]] + year_effect[stan_data_sim$year_raw[i]] + 
                                                   season_pred[stan_data_sim$date[i]] + noise[i]));
       }
    }
    

    
    # catch a few extreme counts -----------------------------------
    
    tt = (which(stan_data_sim$count > max(stan_data$count)*2))
    if(length(tt)/stan_data$ncounts > 0.002){
       warning(paste(sp,"simulated data has many extreme values"))
    }
    stan_data_sim$count[tt] <- max(stan_data$count)*2
    
# compare fake to observed data -------------------------------------------

    

   
    
    comp_dat1 = data.frame(year_raw = stan_data$year_raw,
                          site = stan_data$site,
                          strat = as.factor(stan_data$strat),
                          date = stan_data$date,
                          count = stan_data$count,
                          sim = "observed")
    
    comp_dat2 = data.frame(year_raw = stan_data$year_raw,
                           site = stan_data$site,
                           strat = as.factor(stan_data$strat),
                           date = stan_data$date,
                           count = stan_data_sim$count,
                           sim = "simulated")
    comp_dat <- bind_rows(comp_dat1,comp_dat2)
   
    #remove the fitted object from the real-data analysis
rm(list = "cmdstanfit")
    
# fit simulated data ------------------------------------------------------

    
sp = paste("Simulated",sp1,sep = "_")    

      
   spf = gsub(sp,pattern = " ",replacement = "_")
   spf = gsub(pattern = "\'",replacement = "",
           x = spf)
   sp_file_name <- paste0(spf,"-gamma-t")

   pdim = ceiling(sqrt(stan_data$nstrata)) # just for plotting
   
   pdf(paste0("figures/",sp_file_name,"_sim_vs_obs.pdf"),
       width = 11,
       height = 8)
   # comp_box = ggplot(data = comp_dat,aes(y = count,colour = sim))+
   #    geom_boxplot()+
   #    facet_wrap(~strat,nrow = pdim, ncol = pdim,scales = "free")
   # 
   # print(comp_box)
   
   comp_point = ggplot(data = comp_dat,aes(y = count+1,x = sim,colour = sim))+
      geom_point(position = position_jitter(width = 0.4))+
      scale_y_log10()+
      facet_wrap(~strat,nrow = pdim, ncol = pdim,scales = "free")
   
   print(comp_point)
   
   dev.off()
   
   
   
   #if(!file.exists(paste0(output_dir,"/",sp_file_name,"-",1,".csv"))){
   
   
 

    
    ## compile model
    modl = cmdstan_model(stan_file=mod.file2)

print(sp)

if(grepl(x = mod.file2,pattern = "two_season")){
init_def <- function(){ list(noise_raw = rnorm(stan_data$ncounts,0,0.1),
                             alpha_raw = rnorm(stan_data$nsites,0,0.1),
                             ALPHA1 = 0,
                             year_effect_raw = rnorm(stan_data$nyears,0,0.1),
                             B_season_raw1 = rnorm(stan_data$ndays,0,0.1),
                             B_season_raw2 = rnorm(stan_data$ndays,0,0.1),
                             sdnoise = 0.2,
                             sdalpha = 0.1,
                             sdyear_gam = 1,
                             sdyear_gam_strat = runif(stan_data$nknots_year,0.1,0.2),
                             sdseason = c(0.1,0.1),
                             sdyear = 0.1,
                             B_raw = rnorm(stan_data$nknots_year,0,0.1),
                             b_raw = matrix(rnorm(stan_data$nknots_year*stan_data$nstrata,0,0.01),
                                            nrow = stan_data$nstrata,ncol = stan_data$nknots_year))}

}else{
   init_def <- function(){ list(noise_raw = rnorm(stan_data$ncounts,0,0.1),
                                alpha_raw = rnorm(stan_data$nsites,0,0.1),
                                ALPHA1 = 0,
                                year_effect_raw = rnorm(stan_data$nyears,0,0.1),
                                B_season_raw = rnorm(stan_data$ndays,0,0.1),
                                #B_season_raw2 = rnorm(stan_data$ndays,0,0.1),
                                sdnoise = 0.2,
                                sdalpha = 0.1,
                                sdyear_gam = 1,
                                sdyear_gam_strat = runif(stan_data$nknots_year,0.1,0.2),
                                sdseason = c(0.1),
                                sdyear = 0.1,
                                B_raw = rnorm(stan_data$nknots_year,0,0.1),
                                b_raw = matrix(rnorm(stan_data$nknots_year*stan_data$nstrata,0,0.01),
                                               nrow = stan_data$nstrata,ncol = stan_data$nknots_year))}
}

cmdstanfit<- modl$sample(data=stan_data_sim,
               refresh=100,
               chains=4, iter_sampling =800,
               iter_warmup=1000,
               parallel_chains = 4,
               max_treedepth = 15,
               adapt_delta = 0.8,
               init = init_def)


cmdstanfit$save_output_files(dir = output_dir,
                             basename = sp_file_name,
                             timestamp = FALSE,
                             random = FALSE)

csvfl = paste0(output_dir,"/",sp_file_name,"-",1:4,".csv")
#cmdstanfit$save_object(file = paste0(output_dir,"/",sp_file_name,".RDS"))

shiny_explore <- FALSE
if(shiny_explore){
   #load(paste0(output_dir,"/",sp_file_name,"_fit_add.RData"))
  sl_rstan <- rstan::read_stan_csv(csvfl)
  launch_shinystan(as.shinystan(sl_rstan))
}




save(list = c("cmdstanfit",
              "stan_data_sim",
              "dts",
              "real_grid",
              "strats_dts",
              "strat_regions",
              "parms",
              "sp",
              "spf",
              "sp_file_name",
              "output_dir",
              "prior",
              "mod.file2",
              "noise_dist2",
              "csvfl",
              "strats_xy"),
     file = paste0(output_dir,"/",sp_file_name,"_fit_add.RData"))






if(balanced_sim){
# Fit simulated data with balanced sampling -------------------------------

   
   
   
stan_data_bal <- stan_data

#table of the sites to sample
sites_strats <- unique(data.frame(strat = stan_data$strat,
                                  site = stan_data$site,
                                  seas_strat = stan_data$seas_strat))

# balanced sequence of sampling within season - assuming each site sample 12 times each year
sample_days = floor(seq(6,stan_data$ndays-6,length = 12))

bal_days <- NULL
for(j in 1:nrow(sites_strats)){ #adding some random variation to which days each site is sampled
   tmp = sites_strats[j,]
   tmp2 = sample_days + round(runif(length(sample_days),-5,5))
   tmp = expand_grid(tmp,date = as.integer(tmp2))
   
   bal_days <- bind_rows(bal_days,tmp)
}

bal_dat <- expand_grid(bal_days,year_raw = c(1:stan_data$nyears))


stan_data_bal$ncounts <- nrow(bal_dat)
stan_data_bal$year_raw <- bal_dat$year_raw
stan_data_bal$site <- bal_dat$site
stan_data_bal$strat <- bal_dat$strat
stan_data_bal$date <- bal_dat$date
stan_data_bal$seas_strat <- bal_dat$seas_strat


noise = rnorm(stan_data_bal$ncounts,0,sdnoise)

stan_data_bal$count <- vector(mode = "integer",length = stan_data_bal$ncounts)


# Generate fake counts using above estimates ------------------------------



if(grepl(x = mod.file2,pattern = "two_season")){
   for(i in 1:stan_data_sim$ncounts){
      stan_data_bal$count[i] = rpois(n = 1,exp(ALPHA1 + year_pred[stan_data_bal$year_raw[i],stan_data_bal$strat[i]] + 
                                                  alpha[stan_data_bal$site[i]] + year_effect[stan_data_bal$year_raw[i]] + 
                                                  season_pred[stan_data_bal$date[i],stan_data_bal$seas_strat[i]] + noise[i]));
      
   }
}else{
   stan_data_bal$count[i] = rpois(n = 1,exp(ALPHA1 + year_pred[stan_data_bal$year_raw[i],stan_data_bal$strat[i]] + 
                                               alpha[stan_data_bal$site[i]] + year_effect[stan_data_bal$year_raw[i]] + 
                                               season_pred[stan_data_bal$date[i]] + noise[i]));
   
}



# catch a few extremely unlikely counts -----------------------------------

tt = (which(stan_data_bal$count > max(stan_data$count)*2))
if(length(tt)/stan_data$ncounts > 0.002){
   warning(paste(sp,"simulated data have many extreme values"))
}
stan_data_bal$count[tt] <- max(stan_data$count)*2



# fit Balanced simulated data ------------------------------------------------------


#sp = paste("Simulated",sp1,sep = "_")    
sp = paste0("Balanced_",sp) 


spf = gsub(sp,pattern = " ",replacement = "_")
spf = gsub(pattern = "\'",replacement = "",
           x = spf)
sp_file_name <- paste0(spf,"-gamma-t")


#if(!file.exists(paste0(output_dir,"/",sp_file_name,"-",1,".csv"))){





## compile model
modl = cmdstan_model(stan_file=mod.file2)

print(sp)

if(grepl(x = mod.file2,pattern = "two_season")){
   
init_def <- function(){ list(noise_raw = rnorm(stan_data_bal$ncounts,0,0.1),
                             alpha_raw = rnorm(stan_data_bal$nsites,0,0.1),
                             ALPHA1 = 0,
                             year_effect_raw = rnorm(stan_data_bal$nyears,0,0.1),
                             B_season_raw1 = rnorm(stan_data_bal$ndays,0,0.1),
                             B_season_raw2 = rnorm(stan_data_bal$ndays,0,0.1),
                             sdnoise = 0.2,
                             sdalpha = 0.1,
                             sdyear_gam = 1,
                             sdyear_gam_strat = runif(stan_data_bal$nknots_year,0.1,0.2),
                             sdseason = c(0.1,0.1),
                             sdyear = 0.1,
                             B_raw = rnorm(stan_data_bal$nknots_year,0,0.1),
                             b_raw = matrix(rnorm(stan_data_bal$nknots_year*stan_data_bal$nstrata,0,0.01),
                                            nrow = stan_data_bal$nstrata,ncol = stan_data_bal$nknots_year))}


}else{
   init_def <- function(){ list(noise_raw = rnorm(stan_data_bal$ncounts,0,0.1),
                                alpha_raw = rnorm(stan_data_bal$nsites,0,0.1),
                                ALPHA1 = 0,
                                year_effect_raw = rnorm(stan_data_bal$nyears,0,0.1),
                                B_season_raw = rnorm(stan_data_bal$ndays,0,0.1),
                                #B_season_raw2 = rnorm(stan_data_bal$ndays,0,0.1),
                                sdnoise = 0.2,
                                sdalpha = 0.1,
                                sdyear_gam = 1,
                                sdyear_gam_strat = runif(stan_data_bal$nknots_year,0.1,0.2),
                                sdseason = c(0.1),
                                sdyear = 0.1,
                                B_raw = rnorm(stan_data_bal$nknots_year,0,0.1),
                                b_raw = matrix(rnorm(stan_data_bal$nknots_year*stan_data_bal$nstrata,0,0.01),
                                               nrow = stan_data_bal$nstrata,ncol = stan_data_bal$nknots_year))} 
   
}

cmdstanfit<- modl$sample(data=stan_data_bal,
                         refresh=100,
                         chains=4, iter_sampling =800,
                         iter_warmup=1000,
                         parallel_chains = 4,
                         max_treedepth = 15,
                         adapt_delta = 0.8,
                         init = init_def)


cmdstanfit$save_output_files(dir = output_dir,
                             basename = sp_file_name,
                             timestamp = FALSE,
                             random = FALSE)

csvfl = paste0(output_dir,"/",sp_file_name,"-",1:4,".csv")
#cmdstanfit$save_object(file = paste0(output_dir,"/",sp_file_name,".RDS"))

shiny_explore <- FALSE
if(shiny_explore){
   #load(paste0(output_dir,"/",sp_file_name,"_fit_add.RData"))
   sl_rstan <- rstan::read_stan_csv(csvfl)
   launch_shinystan(as.shinystan(sl_rstan))
}




save(list = c("cmdstanfit",
              "stan_data_bal",
              "dts",
              "real_grid",
              "strats_dts",
              "strat_regions",
              "parms",
              "sp",
              "spf",
              "sp_file_name",
              "output_dir",
              "prior",
              "mod.file2",
              "noise_dist2",
              "csvfl",
              "strats_xy"),
     file = paste0(output_dir,"/",sp_file_name,"_fit_add.RData"))






}#end if balanced  
   

 }### end of species loop









# COMPARE FITTED WITH KNOWN ------------------------------------------



   blank_list <- vector(mode = "list",length = length(sps))
names(blank_list) <- sps

trend_maps_1980 <- blank_list
trend_maps_3gen <- blank_list
trend_maps_L3gen <- blank_list

trend_comparison_1980 <- blank_list
trend_comparison_3gen <- blank_list
trend_comparison_L3gen <- blank_list

trendsout <- NULL
TRENDSout <- NULL
indices_out <- NULL
indices_out_strat <- NULL

saved_trends <- read.csv("trends/All_gamma_t_survey_wide_trends.csv")
saved_trends_normal <- read.csv("trends/All_gamma_normal_survey_wide_trends.csv")

# Comarison Species Loop --------------------------------------------------
prior = "gamma"
noise_dist2 = "t"

for(sp1 in sps){
   
   
   
   
   
saved_trend <- saved_trends[which(saved_trends$species == sp1 &
                                     saved_trends$start_year == 1980 &
                                     saved_trends$end_year == 2019), "trend"]
   
    
   sp = paste("Simulated",sp1,sep = "_")    
   
  
   
   spf = gsub(sp,pattern = " ",replacement = "_")
   spf = gsub(pattern = "\'",replacement = "",
              x = spf)
   sp_file_name <- paste0(spf,"-",prior,"-",noise_dist2)
   

   
   output_dir <- "g:/Shorebird_Migration_Trends/output"



       noise_dist_sel <- noise_dist2

      
      if(file.exists(paste0(output_dir,"/",sp_file_name,"_fit_add.RData"))){
         load(paste0(output_dir,"/",sp_file_name,"_fit_add.RData"))
         
         B1 = 0.5*(saved_trend/100)  #linear trend at center of range = 2*B1
         bx = -0.005 #implies a 1-2%/year range in trends by longitude
         by = -0.02 #implies a 2-4%/year range in trends by latitude
         
         strats_xy <- strats_dts %>% 
            mutate(x = as.numeric(str_split(hex_name,pattern = "_",simplify = TRUE)[,1]),
                   y = as.numeric(str_split(hex_name,pattern = "_",simplify = TRUE)[,2]),
                   x_scale = (x-mean(x))/(0.5*diff(range(x))),#scaled x coordinate to create a longitudinal gradient in trends
                   y_scale = (y-mean(y))/(0.5*diff(range(y))),#scaled y coordinate to create a latitudinal gradient in trends
                   b1 = (x_scale*bx+B1) + (y_scale*by+B1)) #initial log-linear slopes for each stratum with spatial gradients
         
         midyear = 20
         b1 = strats_xy$b1 #initial slopes for each stratum
         #b2 = strats_xy$b2 #second slopes for each stratum
         
         
         if(length(csvfl) > 4){csvfl <- csvfl[1:4]}
         cmdstanfit <- rstan::read_stan_csv(csvfl) 
         #<- readRDS(paste0(output_dir,"/",sp_file_name,".RDS"))
         #load(paste0(output_dir,"/",spf,"_fit_add.RData"))
         three_gen <- max(10,ceiling(gens[which(gens$Common_name == sp),"GenLength"]*3)) #20
         #Three generation assessment time in COSEWIC report
         y3g <- 2019-three_gen
         

         
         # Calculate Annual indices using samples ----------------------------------
         
         syear = min(dts$YearCollected)
         
         # gather_draws2 <- funtion(model,
         #                          vars = "",
         #                          dims = 1){
         #   
         #   cmt = paste0("(",paste(vars,sep = "|"),")")
         #   dmt = paste(rep("[:digit:]",times = dims),sep = ",")
         #   nmpt = paste0("(",cmt,"\\[",dmt,"\\]",")")
         #   plong <- draws %>% pivot_longer(
         #     cols = matches(cmt),
         #     names_pattern = regex(nmpt),
         #     names_to = c("variable","group"),
         #     values_to = ".value"
         #   )
         # }
         # 
         
         
         
         Nsamples <- posterior_samples(fit = cmdstanfit,
                                       parm = "N",
                                       dims = c("y"))
         
         
         
         Nsamples$year <- Nsamples$y + (syear-1)
         
         NSmoothsamples <- posterior_samples(fit = cmdstanfit,
                                             parm = "NSmooth",
                                             dims = c("y"))
         NSmoothsamples$year <- NSmoothsamples$y + (syear-1)
         
         #loo_ic[[sp]] <- loo(cmdstanfit)
         
         # # looic -------------------------------------------------------------------
         # 
         # 
         # loo_ic[[sp]] = loo(slope_icar_stanfit)
         # 
         # 
         # library(loo)
         # library(GGally)
         # 
         # pw_loo <- loo(cmdstanfit)
         # loopoint = as.data.frame(pw_loo$pointwise)
         # 
         # dts_loo <- bind_cols(dts,loopoint) %>%
         #    mutate(log_countp1 <- log(count+1))
         # 
         # wcl <- which(names(dts_loo) %in% c("stratn","log_countp1","year","date","influence_pareto_k","looic"))
         # prs_plot <- ggpairs(data = dts_loo,columns = wcl)
         # 
         # dts_loo$species <- sp
         # 
         # #loo_df <- bind_rows(loo_df,dts_loo)
         # 
         # 
         # loo_by_strat <- dts_loo %>% group_by(hex_name) %>%
         #    summarise(mean_looic = mean(looic),
         #              median_looic = median(looic),
         #              q75_looic = as.numeric(quantile(looic,0.75)),
         #              mean_k = mean(influence_pareto_k),
         #              median_k = median(influence_pareto_k),
         #              q75_k = as.numeric(quantile(influence_pareto_k,0.75)),
         #              n_counts = n(),
         #              median_count = median(count),
         #              mean_count = mean(count))
         # 
         # strat_pairs <- ggpairs(data = loo_by_strat,columns = 2:ncol(loo_by_strat))
         # pdf(paste0("Figures/",sp,prior,"_",noise_dist_sel,"_loo_pairs.pdf"),
         #     width = 11,
         #     height = 11)
         # print(prs_plot)
         # print(strat_pairs)
         # dev.off()
         # 
         # 
         # # Alphas by site ----------------------------------------------------------
         #alpha_samples <- slope_icar_stanfit %>% gather_draws(alpha[s])
         
         alpha_samples <- posterior_samples(fit = cmdstanfit,
                                            parm = "alpha",
                                            dims = c("s"))
         
         sdalpha_samples <- posterior_samples(fit = cmdstanfit,
                                              parm = "sdalpha") 
         sdalpha = sdalpha_samples %>% 
            summarise(mean = mean((.value)),
                      lci = quantile((.value),0.025),
                      uci = quantile((.value),0.975))
         
         
         sites_strat = (stan_data_sim$sites)
         nstrata = stan_data_sim$nstrata
         nsites_strat = stan_data_sim$nsites_strat
         # #offsets = stan_data_sim$site_size
         # 
         sitesbystrat = NULL
         for(st in 1:stan_data_sim$nstrata){
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
         # ## site-effects against the observed counts at each site
         nr = ceiling(sqrt(nstrata))
         obs_by_alpha = ggplot(data = alphas,aes(x = site,y = mean))+
            geom_pointrange(aes(ymin = lci,ymax = uci),colour = "red")+
            geom_point(data = dts, inherit.aes = FALSE,aes(x = site,y = count),
                       fill = grDevices::grey(0.6),colour = grDevices::grey(0),alpha = 0.1,
                       position = position_jitter(width = 0.2,height = 0))+
            labs(title = sp)+
            facet_wrap(~strat,nrow = nr,ncol = nr,scales = "free")+
            theme_minimal()
         
         pdf(paste0("Figures/",sp,prior,"_",noise_dist_sel,"_obs_by_alpha.pdf"),
             width = 11,
             height = 11)
         print(obs_by_alpha)
         dev.off()
         
         # 
         # 
         # 
         # # visualize the seasonal corrections --------------------------------------
         # 
         # 
         # 
         # # # extracting the seasonal smooth ------------------------------------------
         # # 
         # 
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
         
         
         # n_by_strat <- dts %>% group_by(strat) %>% 
         #   summarise(n_cst = n())
         # n_countsby_site <- dts %>% group_by(site,strat) %>% 
         #   summarise(n_c = n()) %>% 
         #   left_join(.,n_by_strat,by = "strat") %>% 
         #   mutate(n_c = n_c/n_cst) %>% 
         #   select(site,strat,n_c)
         
         
         n_countsby_site <- dts %>% group_by(site) %>% 
            summarise(n_c = n()) 
         # 
         # if(grepl(x = mod.file1,pattern = "two_season")){
         #    
         #    strat_season_strat <- dts %>% distinct(seas_strat,strat,hex_name)
         #    
         #    alphas <- left_join(alphas,strat_season_strat,by = "strat")
         #    
         #    
         #    strat_offs <- alphas %>% left_join(.,n_countsby_site,by = "site") %>% 
         #       group_by(strat,seas_strat) %>% 
         #       summarise(adjs = mean(mean))
         #    
         #    
         #    
         #    
         #    season_samples <- posterior_samples(fit = cmdstanfit,
         #                                        parm = "season_pred",
         #                                        dims = c("d","s"))
         #    
         #    
         #    seasonEffectT = season_samples %>% group_by(d,s) %>% 
         #       summarise(mean = mean(exp(.value+scale_adj)),
         #                 lci = quantile(exp(.value+scale_adj),0.025),
         #                 uci = quantile(exp(.value+scale_adj),0.975)) %>% 
         #       mutate(day = d,
         #              seas_strat = s) 
         #    
         #    obs_season <- dts %>% group_by(date,seas_strat) %>% 
         #       summarise(mean = mean(count),
         #                 median = median(count),
         #                 lqrt = quantile(count,0.05),
         #                 uqrt = quantile(count,0.95)) %>% 
         #       mutate(day = date)
         #    
         #    sse <- seasonEffectT %>% group_by(seas_strat) %>% 
         #       group_split()
         #    
         #    #seasonEffect <- expand_grid(seasonEffect,strat = c(1:nstrata))
         #    tout <- NULL
         #    for(j in 1:nstrata){
         #       wg = as.integer(strat_offs[which(strat_offs$strat == j),"seas_strat"])
         #       tmp <- sse[[wg]]
         #       tmp$strat <- j
         #       tout <- bind_rows(tout,tmp)
         #    }
         #    
         #    seasonEffect <- left_join(tout,strat_offs,by = "strat") %>% 
         #       mutate(mean = mean*adjs,
         #              lci = lci*adjs,
         #              uci = uci*adjs)
         #    
         #    
         #    seas_strat_offs <- alphas %>% left_join(.,n_countsby_site,by = "site") %>% 
         #       group_by(seas_strat) %>% 
         #       summarise(adjs = mean(mean),
         #                 adjs2 = median(mean))
         #    
         #    seasonEffect_plot <- left_join(seasonEffectT,seas_strat_offs,by = "seas_strat") %>% 
         #       mutate(mean = mean*adjs + 1,
         #              lci = lci*adjs + 1,
         #              uci = uci*adjs + 1)
         #    
         #    
         #    pp_simple <- ggplot()+
         #       geom_point(data = dts,aes(x = date,y = count+1,colour = year),alpha = 0.5,size = 1,position = position_jitter(width = 0.7,height = 0))+
         #       #geom_smooth(data = dts,aes(x = date,y = count+1))+
         #       scale_colour_viridis_c()+
         #       geom_line(data = seasonEffect_plot,aes(x = day,y = mean),inherit.aes = FALSE)+
         #       #coord_cartesian(ylim = c(0,yup))+
         #       geom_ribbon(data = seasonEffect_plot,aes(x = day,y = mean,ymax = uci,ymin = lci),alpha = 0.2,inherit.aes = FALSE)+
         #       ylab("")+
         #       xlab("Days since July 1")+
         #       scale_y_log10()+
         #       facet_wrap(facets = ~seas_strat,nrow = 2, ncol = 1,scales = "free")+
         #       labs(title = sp)
         #    #print(pp_simple)
         #    out_simple_season_graphs[[sp]] <- pp_simple
         #    
         #    
         #    ncl = 3
         #    nrr = 3
         #    ppag = ncl*nrr
         #    rem = nstrata-(floor(nstrata/ppag)*ppag)
         #    if(rem < ncl){
         #       nrr = 4
         #       ppag = ncl*nrr
         #       rem = nstrata-(floor(nstrata/ppag)*ppag)
         #    }
         #    if(rem < ncl){
         #       nrr = 3
         #       ncl = 2
         #       ppag = ncl*nrr
         #       rem = nstrata-(floor(nstrata/ppag)*ppag)
         #    }
         #    if(rem < ncl){
         #       nrr = 5
         #       ncl = 3
         #       ppag = ncl*nrr
         #       rem = nstrata-(floor(nstrata/ppag)*ppag)
         #    }
         #    
         #    tmp_season_graphs <- vector(mode = "list",length = ceiling(nstrata/ppag))
         #    
         #    pdf(file = paste0("Figures/",sp,prior,"_",noise_dist_sel,"_cmd_simple_Season.pdf"),
         #        width = 8.5,
         #        height = 8.5)
         #    
         #    
         #    for(jj in 1:ceiling(nstrata/ppag)){
         #       #yup <- quantile(dts$count,0.99)
         #       pp <- ggplot()+
         #          # geom_pointrange(data = obs_season,inherit.aes = FALSE,
         #          #                 aes(x = day,y = mean,ymin = lqrt,ymax = uqrt),alpha = 0.1)+
         #          geom_point(data = dts,aes(x = date,y = count,colour = year),alpha = 0.5,size = 1,position = position_jitter(width = 0.7,height = 0))+
         #          scale_colour_viridis_c()+
         #          geom_line(data = seasonEffect,aes(x = day,y = mean),inherit.aes = FALSE)+
         #          #coord_cartesian(ylim = c(0,yup))+
         #          geom_ribbon(data = seasonEffect,aes(x = day,y = mean,ymax = uci,ymin = lci),alpha = 0.2,inherit.aes = FALSE)+
         #          ylab("")+
         #          xlab("Days since July 1")+
         #          facet_wrap_paginate(facets = ~strat,page = jj,nrow = nrr, ncol = ncl,scales = "free")+
         #          labs(title = sp)
         #       tmp_season_graphs[[jj]] <- pp
         #       print(pp)
         #    }
         #    dev.off()
         #    
         #    # pp <- ggplot(data = seasonEffect,aes(x = day,y = mean))+
         #    #   # geom_pointrange(data = obs_season,inherit.aes = FALSE,
         #    #   #                 aes(x = day,y = mean,ymin = lqrt,ymax = uqrt),alpha = 0.1)+
         #    #   geom_point(data = dts,inherit.aes = FALSE,aes(x = day,y = count),alpha = 0.1,size = 1)+
         #    #   geom_line()+
         #    #   geom_ribbon(aes(x = day,y = mean,ymax = uci,ymin = lci),alpha = 0.2)+
         #    #   ylab("")+
         #    #   xlab("Days since July 1")+
         #    #   labs(title = sp)+
         #    #   facet_wrap(facets = ~seas_strat,ncol = 3,scales = "free")
         #    # 
         #    
         # }else{
         #    strat_offs <- alphas %>% group_by(strat) %>% 
         #       summarise(adjs = mean(mean))
         #    
         #    
         #    season_samples <- posterior_samples(fit = cmdstanfit,
         #                                        parm = "season_pred",
         #                                        dims = c("d"))
         #    
         #    
         #    #season_samples <- slope_icar_stanfit %>% gather_draws(season_pred[d])
         #    seasonEffectT = season_samples %>% group_by(d) %>% 
         #       summarise(mean = mean(exp(.value+scale_adj)),
         #                 lci = quantile(exp(.value+scale_adj),0.025),
         #                 uci = quantile(exp(.value+scale_adj),0.975)) %>% 
         #       mutate(day = d) 
         #    
         #    obs_season <- dts %>% group_by(date) %>% 
         #       summarise(mean = mean(count),
         #                 median = median(count),
         #                 lqrt = quantile(count,0.05),
         #                 uqrt = quantile(count,0.95)) %>% 
         #       mutate(day = date)
         #    
         #    seasonEffect <- expand_grid(seasonEffectT,strat = c(1:nstrata))
         #    seasonEffect <- left_join(seasonEffect,strat_offs) %>% 
         #       mutate(mean = mean*adjs,
         #              lci = lci*adjs,
         #              uci = uci*adjs)
         #    
         #    
         #    
         #    # seas_strat_offs <- alphas %>% left_join(.,n_countsby_site,by = "site") %>% 
         #    #   group_by(seas_strat) %>% 
         #    #   summarise(adjs = mean(mean),
         #    #             adjs2 = median(mean))
         #    
         #    seasonEffect_plot <- seasonEffectT %>% 
         #       mutate(mean = mean*mean(alphas$mean)+1,
         #              lci = lci*mean(alphas$mean)+1,
         #              uci = uci*mean(alphas$mean)+1)
         #    
         #    
         #    pp_simple <- ggplot()+
         #       geom_point(data = dts,aes(x = date,y = count+1,colour = year),alpha = 0.5,size = 1,position = position_jitter(width = 0.7,height = 0))+
         #       #geom_smooth(data = dts,aes(x = date,y = count+1))+
         #       scale_colour_viridis_c()+
         #       geom_line(data = seasonEffect_plot,aes(x = day,y = mean),inherit.aes = FALSE)+
         #       #coord_cartesian(ylim = c(0,yup))+
         #       geom_ribbon(data = seasonEffect_plot,aes(x = day,y = mean,ymax = uci,ymin = lci),alpha = 0.2,inherit.aes = FALSE)+
         #       ylab("")+
         #       xlab("Days since July 1")+
         #       scale_y_log10()+
         #       labs(title = sp)
         #    #print(pp_simple)
         #    out_simple_season_graphs[[sp]] <- pp_simple
         #    
         #    pdf(file = paste0("Figures/",sp,prior,"_",noise_dist_sel,"_cmd_simple_Season_simplified.pdf"),
         #        width = 8.5,
         #        height = 8.5)
         #    print(pp_simple)
         #    dev.off()
         #    
         #    ncl = 3
         #    nrr = 3
         #    ppag = ncl*nrr
         #    rem = nstrata-(floor(nstrata/ppag)*ppag)
         #    if(rem < ncl){
         #       nrr = 4
         #       ppag = ncl*nrr
         #       rem = nstrata-(floor(nstrata/ppag)*ppag)
         #    }
         #    if(rem < ncl){
         #       nrr = 3
         #       ncl = 2
         #       ppag = ncl*nrr
         #       rem = nstrata-(floor(nstrata/ppag)*ppag)
         #    }
         #    if(rem < ncl){
         #       nrr = 5
         #       ncl = 3
         #       ppag = ncl*nrr
         #       rem = nstrata-(floor(nstrata/ppag)*ppag)
         #    }
         #    tmp_season_graphs <- vector(mode = "list",length = ceiling(nstrata/ppag))
         #    
         #    pdf(file = paste0("Figures/",sp,prior,"_",noise_dist_sel,"_cmd_simple_Season.pdf"),
         #        width = 8.5,
         #        height = 8.5)
         #    
         #    for(jj in 1:ceiling(nstrata/ppag)){
         #       #yup <- quantile(dts$count,0.99)
         #       pp <- ggplot()+
         #          # geom_pointrange(data = obs_season,inherit.aes = FALSE,
         #          #                 aes(x = day,y = mean,ymin = lqrt,ymax = uqrt),alpha = 0.1)+
         #          geom_point(data = dts,aes(x = date,y = count,colour = year),alpha = 0.5,size = 1,position = position_jitter(width = 0.7,height = 0))+
         #          scale_colour_viridis_c()+
         #          geom_line(data = seasonEffect,aes(x = day,y = mean),inherit.aes = FALSE)+
         #          #coord_cartesian(ylim = c(0,yup))+
         #          geom_ribbon(data = seasonEffect,aes(x = day,y = mean,ymax = uci,ymin = lci),alpha = 0.2,inherit.aes = FALSE)+
         #          ylab("")+
         #          xlab("Days since July 1")+
         #          facet_wrap_paginate(facets = ~strat,page = jj,nrow = nrr, ncol = ncl,scales = "free")+
         #          labs(title = sp)
         #       tmp_season_graphs[[jj]] <- pp
         #       print(pp)
         #    }
         #    dev.off()
         #    
         # }
         # 
         # season_graphs[[sp]] <- tmp_season_graphs
         # 
         # 
         
         
         # plot the changing mean site-effect over years ---------------------------
         
         
         # overall -----------------------------------------------------------------
         
         # alphasE = alpha_samples %>% group_by(s) %>% 
         #    summarise(mean = mean((.value)),
         #              lci = quantile((.value),0.025),
         #              uci = quantile((.value),0.975)) %>% 
         #    mutate(site = s) %>% left_join(.,sitesbystrat)
         # 
         # dts_alpha <- left_join(dts,alphasE,by = c("site","strat"))
         # 
         # alphas_yr <- dts_alpha %>% group_by(YearCollected,site,strat) %>% 
         #    summarise(mean_alpha = mean(mean)) %>% ungroup() %>% 
         #    group_by(YearCollected) %>% 
         #    summarise(mean_alpha = mean(mean_alpha))
         # 
         # AA_y_p <- ggplot(data = alphas_yr,aes(x = YearCollected,y = mean_alpha))+
         #    geom_point()+
         #    geom_smooth()+
         #    geom_abline(slope = 0,intercept = 0,colour = grey(0.3))+
         #    xlab("")+
         #    labs(title = paste(sp,prior,"_",noise_dist_sel,"_mean site-effect of included sites by year"))
         # 
         # pdf(file = paste0("Figures/",sp,prior,"_",noise_dist_sel,"_cmd_alphas_by_yr.pdf"),
         #     width = 8.5,
         #     height = 8.5)
         # print(AA_y_p)
         # 
         # tmp_out_alphas_by_yr <- vector(mode = "list",length = 1+(ceiling(nstrata/ppag)))
         # tmp_out_alphas_by_yr[[1]] <- AA_y_p
         # # by strata ---------------------------------------------------------------
         # 
         # alphas_yr_s <- dts_alpha %>% group_by(YearCollected,site,strat) %>% 
         #    summarise(mean_alpha = mean(mean)) %>% ungroup() %>% 
         #    group_by(YearCollected,strat) %>% 
         #    summarise(mean_alpha = mean(mean_alpha))
         # 
         # for(jj in 1:ceiling(nstrata/ppag)){
         #    
         #    a_y_p <- ggplot(data = alphas_yr_s,aes(x = YearCollected,y = mean_alpha))+
         #       geom_point()+
         #       geom_smooth()+
         #       geom_abline(slope = 0,intercept = 0,colour = grey(0.3))+
         #       xlab("")+
         #       labs(title = paste(sp,prior,"_",noise_dist_sel,"_mean site-effect of included sites by year"))+
         #       facet_wrap_paginate(facets = ~strat,page = jj,nrow = nrr, ncol = ncl,scales = "free")
         #    
         #    print(a_y_p)
         #    tmp_out_alphas_by_yr[[jj+1]] <- a_y_p
         # }
         # 
         # dev.off()
         # 
         # out_alphas_by_yr[[sp]] <- tmp_out_alphas_by_yr
         # 
         # 
         # calculate trends continent --------------------------------------------------------
         
         t_NSmooth_80 <- ItoT(inds = NSmoothsamples,
                              start = 1980,
                              end = 2019,
                              regions = NULL,
                              qs = 95,
                              sp = sp,
                              type = "Long-term")
         TRENDSout <- bind_rows(TRENDSout,t_NSmooth_80)
         
         # t_NSmooth_04 <- ItoT(inds = NSmoothsamples,
         #                      start = 2004,
         #                      end = 2019,
         #                      regions = NULL,
         #                      qs = 95,
         #                      sp = sp,
         #                      type = "15-year")
         #TRENDSout <- bind_rows(TRENDSout,t_NSmooth_04)
         
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
         
         
         
         season_samples <- posterior_samples(fit = cmdstanfit,
                                             parm = "season_pred",
                                             dims = c("d","s")) 
         
         
         #drawst = as_draws_df(cmdstanfit$draws(variables = c("N")))
         
         
         indicesN <- index_summary(samples = Nsamples,
                                   parm = "N",
                                   dims = "year",
                                   strat_offsets = NULL)
         
         
         indicesNSmooth <- index_summary(samples = NSmoothsamples,
                                         parm = "NSmooth",
                                         dims = "year",
                                         strat_offsets = NULL)
         
         
         
         
         
         
         indices = bind_rows(indicesN,indicesNSmooth)
         indices$species <- sp
         
         indices_out <- bind_rows(indices_out,indices)
         #indices$year = indices$year + (syear-1)
         yup = max(max(indices$uci),quantile(indices$obsmean,0.7))
         
         indices$parm <- factor(indices$parm,ordered = T,levels = c("NSmooth","N"))
         
         N_gg = ggplot(data = indices,aes(x = year, y = median,fill = parm))+
            geom_ribbon(aes(ymin = lci,ymax = uci),alpha = 0.2)+
            geom_line(aes(colour = parm))+
            labs(title = paste(sp,prior,"_",noise_dist_sel,"_Survey-wide trajectory (full and smooth) with obs means"))+
            annotate("text", x = 1997, y = yup*0.9, label = anot_80)+
            annotate("text", x = 1997, y = yup*0.8, label = anot_90)+
            annotate("text", x = 1997, y = yup*0.7, label = anot_07)+
            my_col2_traj+
            coord_cartesian(ylim = c(0,yup))+
            geom_point(aes(y = obsmean),colour = grey(0.5),alpha = 0.3)+
            theme_classic()+
            scale_y_continuous(limits = c(0,NA))+
            theme(legend.position = "none")#+
         #scale_size_area()
         
         N_gg_simple = ggplot(data = indices,aes(x = year, y = median,fill = parm))+
            geom_ribbon(aes(ymin = lci,ymax = uci),alpha = 0.2)+
            geom_line(aes(colour = parm))+
            labs(title = paste(sp,prior,"_",noise_dist_sel,"_Survey-wide trajectory (full and smooth)"))+
            my_col2_traj+
            xlab("")+
            ylab("Modeled mean count")+
            theme_classic()+
            scale_y_continuous(limits = c(0,NA))+
            theme(legend.position = "none")
         
         pdf(paste0("figures/",sp,prior,FYYYY,"_cmd_GAMYE_survey_wide_trajectory_simple",grid_spacing/1000,".pdf"))
         print(N_gg)
         print(N_gg_simple)
         dev.off()
         
         #sp_ind_plots_diagnostic[[sp]] <- N_gg
         
         
         
         
        
         # 
         # #print(N_gg_simple)
         # 
         # sp_ind_plots[[sp]] <- N_gg_simple
         # 
         
         # Continental Smooth spaghetti plot ---------------------------------------
         #  indicesNSmooth$year = indicesNSmooth$year + (syear-1)
         
         # set.seed(2019)
         # r_draws <- sample(size = 100,x = 1:max(NSmoothsamples$.draw))
         # spg_trajs <- NSmoothsamples %>% filter(.draw %in% r_draws) %>% 
         #    mutate(draw_f = factor(.draw))
         # 
         # lev = spg_trajs %>% 
         #    ungroup() %>% 
         #    filter(year == 2019) %>% 
         #    select(draw_f,.value) %>% 
         #    mutate(midp = log(.value)) %>% 
         #    select(draw_f,midp)
         # 
         # spg_trajs <- left_join(spg_trajs,lev,by = "draw_f")
         # 
         # 
         # n_gg_spag = ggplot(data = indicesNSmooth,aes(x = year, y = median))+
         #    geom_ribbon(aes(ymin = lci,ymax = uci),alpha = 0.25)+
         #    geom_line(size =2)+
         #    labs(title = paste(sp,prior,"_",noise_dist_sel,"_Random selection of 100 posterior draws of survey-wide trajectories"),
         #         subtitle = "Colour of each posterior draw reflects the value in 2019, demonstrating similar smooths across draws")+
         #    xlab("")+
         #    ylab("Modeled mean count")+
         #    theme_classic()+
         #    theme(legend.position = "none") + 
         #    geom_line(data = spg_trajs,aes(x = year,y = .value,group = draw_f,colour = midp),alpha = 0.8,inherit.aes = FALSE,size = 1)+
         #    scale_color_viridis_c(aesthetics = "colour")+
         #    scale_y_log10()
         # 
         # #sp_spag_plots_diagnostic[[sp]] <- n_gg_spag
         # 
         # 
         # #print(n_gg_spag)
         # 
         # pdf(paste0("figures/",sp,prior,FYYYY,"_GAMYE_spaghetti",".pdf"))
         # print(n_gg_spag)
         # dev.off()
         # 
         
         # combine the trajectories within the original regional strata ------------
         
         
         
         strats_dts <- left_join(strats_dts,strat_regions,by = c("stratn" = "strat"))
         
         
         alpha_adjs <- alphas %>% mutate(adjs = mean) %>% select(site,adjs)
         
         
         nstrata = stan_data_sim$nstrata
         
         nsmoothsamples <- posterior_samples(fit = cmdstanfit,
                                             parm = "nsmooth",
                                             dims = c("stratn","y")) 
         
         
         nsmoothsamples$year <- nsmoothsamples$y + (syear-1)
         nsmoothsamples <- left_join(nsmoothsamples,strats_dts,by = c("stratn"))
         
         
         indicesnsmooth <- index_summary(samples = nsmoothsamples,
                                         parm = "nsmooth",
                                         dims = c("stratn","year"),
                                         site_offsets = alpha_adjs)
         
         
         
         
         nsamples <- posterior_samples(fit = cmdstanfit,
                                       parm = "n",
                                       dims = c("stratn","y")) 
         
         
         nsamples$year <- nsamples$y + (syear-1)
         nsamples <- left_join(nsamples,strats_dts,by = c("stratn"))
         
         
         indicesn <- index_summary(samples = nsamples,
                                   parm = "n",
                                   dims = c("stratn","year"),
                                   site_offsets = alpha_adjs)
         
         
         
         
         
         indices_strat = bind_rows(indicesn,indicesnsmooth)
         #indices_strat$year = indices_strat$year + (syear-1)
         indices_strat <- left_join(indices_strat,strats_dts, by = "stratn")
         indices_strat$species <- sp
         
         indices_out_strat <- bind_rows(indices_out_strat,indices_strat)
         
         
         
         indices_strat$parm <- factor(indices_strat$parm,ordered = T,levels = c("nsmooth","n"))
         
         
         
         pdf(file = paste0("figures/", sp,prior,FYYYY,"_cmd_GAMYE_Strata_trajectories_simple",grid_spacing/1000,".pdf"),
             width = 8.5,
             height = 11)
         print(N_gg_simple)
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
         tmp_sp_ind_plots <- vector(mode = "list",length = ceiling(nstrata/ppag))
         tmp_sp_ind_plots_diagnostic <- tmp_sp_ind_plots
         for(jj in 1:ceiling(nstrata/ppag)){
            n_gg_simple = ggplot(data = indices_strat,aes(x = year, y = median,fill = parm))+
               geom_ribbon(aes(ymin = lci,ymax = uci),alpha = 0.2)+
               geom_line(aes(colour = parm))+
               my_col2_traj+
               xlab("")+
               ylab("Modeled mean count")+
               theme_classic()+
               theme(legend.position = "none")+
               labs(title = sp)+
               facet_wrap_paginate(facets = ~stratn,page = jj,nrow = nrr, ncol = ncl,scales = "free")
            #tmp_sp_ind_plots[[jj]] <- n_gg_simple
            
            print(n_gg_simple)
            
            # n_gg = ggplot(data = indices_strat,aes(x = year, y = median,fill = parm))+
            #    geom_ribbon(aes(ymin = lci,ymax = uci),alpha = 0.2)+
            #    geom_line(aes(colour = parm))+
            #    my_col2_traj+
            #    geom_point(aes(y = obsmean,size = mean_counts_incl_sites),colour = grey(0.5),alpha = 0.3)+
            #    theme_classic()+
            #    labs(title = sp)+
            #    facet_wrap_paginate(facets = ~stratn,page = jj,nrow = nrr, ncol = ncl,scales = "free")+
            #    scale_size_area()
            # tmp_sp_ind_plots_diagnostic[[jj]] <- n_gg
            
            
            
           
         }
         # print(n_gg_spag)
         dev.off()
         
         # sp_ind_plots_strat[[sp]] <- tmp_sp_ind_plots
         # sp_ind_plots_strat_diagnostic[[sp]] <- tmp_sp_ind_plots_diagnostic
         # 
         # 
         # # Trajectory Overplot log-scale -------------------------------------------
         # 
         # indices_strat_smooth <- indices_strat %>% filter(parm == "nsmooth")
         # 
         # 
         # 
         # n_gg_over = ggplot(data = indices_strat_smooth,aes(x = year, y = median,colour = hex_name))+
         #    #geom_ribbon(aes(ymin = q2.5,ymax = q97.5),alpha = 0.2)+
         #    geom_line(alpha = 0.8)+
         #    scale_colour_viridis_d()+
         #    geom_line(data = indicesNSmooth,aes(x = year,y = median),inherit.aes = FALSE)+
         #    theme_classic()+
         #    labs(title = sp)+
         #    theme(legend.position = "none")+
         #    scale_y_log10()
         # pdf(paste0("Figures/log_scale_trajectories_",sp,prior,"_",noise_dist_sel,"_cmd_.pdf"),
         #     width = 11,
         #     height = 8.5)
         # print(n_gg_over)
         # dev.off()
         # traj_overplots[[sp]] <- n_gg_over
         # 
         
         
         # Beta plots --------------------------------------------------------------
         
         
         
         # _samples = as_draws_df(fit$draws(variables = c("B")))
         # sums <- summarise_draws(x = _samples,~quantile2(.x,probs = probs))
         # 
         #   sums[,"k"] = jags_dim(dim = 1,
         #                           var = ,
         #                           cl = "variable",
         #                           dat = B)
         
         # B_samples <- posterior_samples(fit = cmdstanfit,
         #                                parm = "B",
         #                                dims = c("k"))
         # b_samples <- posterior_samples(fit = cmdstanfit,
         #                                parm = "b",
         #                                dims = c("s","k")) %>% 
         #    left_join(.,strats_dts,by = c("s" = "stratn"))
         # 
         # sdbeta_samples <- posterior_samples(fit = cmdstanfit,
         #                                     parm = "sdyear_gam_strat",
         #                                     dims = c("k")) 
         # 
         # sdB_samples <- posterior_samples(fit = cmdstanfit,
         #                                  parm = "sdyear_gam",
         #                                  dims = NULL) 
         # 
         # 
         # sdB <- sdB_samples %>% 
         #    summarise(mean = mean(.value),
         #              lci = quantile(.value,0.025),
         #              uci = quantile(.value,0.975)) 
         # 
         # # B_samples <- slope_icar_stanfit %>% gather_draws(B[k])    
         # #   b_samples <- slope_icar_stanfit %>% gather_draws(b[s,k])    
         # #   b_samples <- left_join(b_samples,strats_dts,by = c("s" = "stratn"))
         # # 
         # B <- B_samples %>% group_by(k) %>%
         #    summarise(mean = mean(.value),
         #              lci = quantile(.value,0.05),
         #              uci = quantile(.value,0.95))
         # b <- b_samples %>% group_by(k,hex_name) %>%
         #    summarise(mean = mean(.value),
         #              lci = quantile(.value,0.025),
         #              uci = quantile(.value,0.975))
         # 
         # # sdbeta_samples <- slope_icar_stanfit %>% gather_draws(sdyear_gam_strat[k])
         # sdbeta <- sdbeta_samples %>% group_by(k) %>%
         #    summarise(mean = mean(.value),
         #              lci = quantile(.value,0.025),
         #              uci = quantile(.value,0.975)) %>%
         #    mutate(k = k+0.1)
         # #B_b <- bind_rows(B,b)
         # 
         # b_plot <- ggplot(data = b,aes(x = k,y = mean,colour = hex_name))+
         #    geom_point(data = B, inherit.aes = FALSE, 
         #               aes(x = k,y = mean),size = 1)+
         #    geom_errorbar(data = B, inherit.aes = FALSE, 
         #                  aes(x = k,y = mean,ymin = lci,ymax = uci),width = 0,alpha = 0.5)+
         #    geom_point()+
         #    theme_classic()+
         #    coord_cartesian(ylim = c(-10,10))+
         #    labs(title = sp)+
         #    geom_pointrange(data = sdbeta,inherit.aes = FALSE,aes(x = k,y = mean,ymin = lci,ymax = uci),colour = "red")+
         #    scale_colour_viridis_d()+
         #    geom_hline(yintercept = 0)+
         #    theme(legend.position = "none")
         # 
         # 
         # pdf(paste0("Figures/beta_sdbeta_",sp,prior,"_",noise_dist_sel,"_cmd_.pdf"),
         #     width = 11,
         #     height = 8.5)
         # print(b_plot)
         # dev.off()
         # 
         # 
         # 
         # beta_overplots[[sp]] <- b_plot
         # 
         
         
         # Trajectories without any intercepts -------------------------------------
         # 
         #     traj_funct <- function(basis,bs){
         #       out_vec = as.numeric(basis %*% bs)
         #       return(out_vec)
         #     }
         #     yr_funct <- function(basis){
         #       out_vec = 1:nrow(basis)+(FYYYY-1)
         #       return(out_vec)
         #     }
         # 
         #     btraj <- b_samples %>% group_by(.draw,hex_name) %>% 
         #       summarise(traj = traj_funct(stan_data_sim$year_basispred,.value),
         #                 year = yr_funct(stan_data_sim$year_basispred),
         #                 .groups = "drop") %>%
         #       group_by(hex_name,year) %>% 
         #       summarise(mean = mean((traj)),
         #                 lci = quantile((traj),0.025),
         #                 uci = quantile((traj),0.975))
         #       
         #     Btraj <- B_samples %>% group_by(.draw) %>% 
         #       summarise(traj = traj_funct(stan_data_sim$year_basispred,.value),
         #                 year = yr_funct(stan_data_sim$year_basispred),
         #                 .groups = "drop") %>%
         #       group_by(year) %>% 
         #       summarise(mean = mean((traj)),
         #                 lci = quantile((traj),0.025),
         #                 uci = quantile((traj),0.975))
         #     
         #     
         # 
         #     
         #     b_plot_alt = ggplot(data = btraj,aes(x = year, y = mean))+
         #       geom_line(data = Btraj,aes(x = year,y = mean),inherit.aes = FALSE,size = 2)+
         #       geom_ribbon(data = Btraj,aes(ymin = uci,ymax = lci),alpha = 0.1)+
         #       geom_line(alpha = 0.8,aes(colour = hex_name))+
         #       scale_colour_viridis_d(aesthetics = c("colour","fill"))+
         #       theme_classic()+
         #       labs(title = sp)+
         #       ylab("Centered log-scale smooths")+
         #       #scale_y_log10()+
         #       theme(legend.position = "none")
         #     #print(b_plot_alt)
         #     alternate_traj_overplots[[sp]] <- b_plot_alt
         #     
         
         #(stan_data_sim$year_basispred * transpose(b[s,]))
         
         
         
         # st_drop <- "2305040_1068494"
         #nsmoothsamples <- nsmoothsamples %>% filter(hex_name != st_drop)
         
         
         # t_nsmooth_strat_04 <- ItoT(inds = nsmoothsamples,
         #                            start = 2004,
         #                            end = 2019,
         #                            regions = "hex_name",
         #                            qs = 95,
         #                            sp = sp,
         #                            type = "15-year",
         #                            centered_trends = TRUE)
         
         #trendsout <- bind_rows(trendsout,t_nsmooth_strat_04)
         
         
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
         
         
         
         
         # t_nsmooth_reg_04 <- ItoT(inds = nsmoothsamples,
         #                          start = 2004,
         #                          end = 2019,
         #                          regions = "Region",
         #                          qs = 95,
         #                          sp = sp,
         #                          type = "15-year",
         #                          centered_trends = TRUE)
         # # trendsout <- bind_rows(trendsout,
         # #                        t_nsmooth_reg_04)
         # 
         # t_nsmooth_reg_80 <- ItoT(inds = nsmoothsamples,
         #                          start = 1980,
         #                          end = 2019,
         #                          regions = "Region",
         #                          qs = 95,
         #                          sp = sp,
         #                          type = "Long-term",
         #                          centered_trends = TRUE)
         # # trendsout <- bind_rows(trendsout,
         # #                        t_nsmooth_reg_80)
         # 
         # t_nsmooth_reg_3g <- ItoT(inds = nsmoothsamples,
         #                          start = y3g,
         #                          end = 2019,
         #                          regions = "Region",
         #                          qs = 95,
         #                          sp = sp,
         #                          type = "Recent-three-generation",
         #                          centered_trends = TRUE)
         # # trendsout <- bind_rows(trendsout,
         # #                        t_nsmooth_reg_3g)
         # # 
         # t_nsmooth_reg_L3g <- ItoT(inds = nsmoothsamples,
         #                           start = syL3g,
         #                           end = y3g,
         #                           regions = "Region",
         #                           qs = 95,
         #                           sp = sp,
         #                           type = "Previous-three-generation",
         #                           centered_trends = TRUE)
         # # trendsout <- bind_rows(trendsout,
         # #                        t_nsmooth_reg_L3g)
         # # 
         # 
         # 
         # 
         # 
         # # st_drop <- "2305040_1068494"
         # # nsamples <- nsamples %>% filter(hex_name != st_drop)
         # # 
         # indsnsmooth_region = ItoI(inds = nsmoothsamples,
         #                           regions = "Region")
         # indsnsmooth_region$type = "Smooth"
         # 
         # indsn_region = ItoI(inds = nsamples,
         #                     regions = "Region")
         # indsn_region$type = "Full"
         # indsn_region <- bind_rows(indsnsmooth_region,indsn_region)
         # indsn_region$species <- sp
         # 
         # ind_fc = ggplot(data = indsn_region,aes(x = year,y = median,group = type))+
         #    geom_ribbon(aes(ymin = lci,ymax = uci,fill = type),alpha = 0.2)+
         #    geom_line(aes(colour = type))+
         #    coord_cartesian(ylim = c(0,NA))+
         #    #scale_colour_viridis_d(aesthetics = c("fill","colour"))+
         #    my_col2_traj+
         #    theme(legend.position = "none")+
         #    theme_classic()+
         #    labs(title = sp)+
         #    facet_wrap(~region,scales = "free")
         # #print(ind_fc)
         # 
         # composite_trajectories[[sp]] <- ind_fc  
         # 
         # 
         # indices_out_composite <- bind_rows(indices_out_composite,indsn_region)
         # 
         # 
         # 
         # 
         
         # Trend heatmaps ----------------------------------------------------------
         ## consider adding site-locations, stratum labels, etc.
         
         strats_xy$trend_1 <- 100*(exp(strats_xy$b1)-1)
         #strats_xy$trend_2 <- 100*(exp(strats_xy$b2)-1)
         strats_xy$North_position <- strats_xy$y_scale
         cont_st <- data.frame(hex_name = "Survey Wide",
                               trend_1 = 100*(exp(mean(strats_xy$b1))-1),
                               North_position = 0)
           strats_xy <- bind_rows(strats_xy,cont_st)
           
           t_NSmooth_80$region <- "Survey Wide"
           t_NSmooth_3g$region <- "Survey Wide"
           t_NSmooth_L3g$region <- "Survey Wide"
           
         pdf(paste0("Figures/Trend_Heat_maps_",sp,prior,"_",noise_dist_sel,"_cmd_.pdf"),
             width = 11,
             height = 8.5)
         t_80 = trend_map(t_nsmooth_strat_80,
                          size_value = "Mean Observed Count")
         print(t_80)
         
         tcomp_80 <-  t_nsmooth_strat_80 %>% 
            left_join(.,strats_xy,by = c("region" = "hex_name")) %>% 
            mutate(prec = 1/((uci-lci)/(1.96*2))^2)
         
         
         t_Ncomp_80 <- left_join(t_NSmooth_80,strats_xy,by = c("region" = "hex_name")) 
         
         tcpl_80 <- ggplot(data = tcomp_80,aes(x = trend_1,y = trend))+
            xlab("True trend")+
            ylab("Estimated Full trend")+
            geom_point(aes(size = prec,colour = North_position))+
            geom_errorbar(aes(ymin = lci,ymax = uci),alpha = 0.3,width = 0)+
            geom_point(data = t_Ncomp_80,aes(x = trend_1,y = trend),alpha = 0.7,inherit.aes = FALSE,size = 5,colour = "red")+
            geom_errorbar(data = t_Ncomp_80,aes(x = trend_1,y = trend,ymin = lci,ymax = uci),
                          alpha = 0.7,colour = "red",width = 0,inherit.aes = FALSE)+
            geom_abline(intercept = 0,slope = 1)
         print(tcpl_80)
            

         t_3g = trend_map(t_nsmooth_strat_3g,
                          size_value = "Mean Observed Count")
         print(t_3g)
         
         tcomp_3g <- left_join(t_nsmooth_strat_3g,strats_xy,by = c("region" = "hex_name")) %>% 
            mutate(prec = 1/((uci-lci)/(1.96*2))^2)
         
        
         t_Ncomp_3g <- left_join(t_NSmooth_3g,strats_xy,by = c("region" = "hex_name")) 
         
         tcpl_3g <- ggplot(data = tcomp_3g,aes(x = trend_1,y = trend))+
            xlab("True trend")+
            ylab("Estimated Recent trend")+
            geom_point(aes(size = prec,colour = North_position))+
            geom_errorbar(aes(ymin = lci,ymax = uci),alpha = 0.3,width = 0)+
            geom_point(data = t_Ncomp_3g,aes(x = trend_1,y = trend),alpha = 0.7,inherit.aes = FALSE,size = 5,colour = "red")+
            geom_errorbar(data = t_Ncomp_3g,aes(x = trend_1,y = trend,ymin = lci,ymax = uci),
                          alpha = 0.7,colour = "red",width = 0,inherit.aes = FALSE)+
            geom_abline(intercept = 0,slope = 1)
         print(tcpl_3g)        
         
         t_L3g = trend_map(t_nsmooth_strat_L3g,
                           size_value = "Mean Observed Count")
         print(t_L3g)
         
         tcomp_L3g <- left_join(t_nsmooth_strat_L3g,strats_xy,by = c("region" = "hex_name")) %>% 
            mutate(prec = 1/((uci-lci)/(1.96*2))^2)
         
         t_Ncomp_L3g <- left_join(t_NSmooth_L3g,strats_xy,by = c("region" = "hex_name")) 
         
         
         tcpl_L3g <- ggplot(data = tcomp_L3g,aes(x = trend_1,y = trend))+
            xlab("True trend")+
            ylab("Estimated Early trend")+
            geom_point(aes(size = prec,colour = North_position))+
            geom_errorbar(aes(ymin = lci,ymax = uci),alpha = 0.3,width = 0)+
            geom_point(data = t_Ncomp_L3g,aes(x = trend_1,y = trend),alpha = 0.7,inherit.aes = FALSE,size = 5,colour = "red")+
            geom_errorbar(data = t_Ncomp_L3g,aes(x = trend_1,y = trend,ymin = lci,ymax = uci),
                          alpha = 0.7,colour = "red",width = 0,inherit.aes = FALSE)+
            geom_abline(intercept = 0,slope = 1)
         print(tcpl_L3g)  
         
         
         dev.off()
         
         trend_maps_1980[[sp]] <- t_80
         # trend_maps_2004[[sp]] <- t_04
         trend_maps_3gen[[sp]] <- t_3g
         trend_maps_L3gen[[sp]] <- t_L3g
         trend_comparison_1980[[sp]] <- tcpl_80
         # trend_maps_2004[[sp]] <- t_04
         trend_comparison_3gen[[sp]] <- tcpl_3g
         trend_comparison_L3gen[[sp]] <- tcpl_L3g
         # 
         
         print(sp)
         
      }# end if species output data exists
}   #end species loop
   



write.csv(trendsout,paste0("trends/All_region_Simulated_strata_",prior,"_",noise_dist_sel,"_composite_trends.csv"),row.names = FALSE)

TRENDSout <- TRENDSout %>% relocate(species,start_year,end_year,trend_type) %>% 
   select(-parameter)

trendsoutsplit <- trendsout %>% relocate(species,start_year,end_year,trend_type,region) %>% 
   select(-c(parameter)) %>% 
   filter(region != "Composite") %>% 
   group_by(region_type) %>% 
   group_split()


write.csv(trendsoutsplit[[1]],paste0("trends/All_strata_Simulated_",prior,"_",noise_dist_sel,"_level_trends.csv"),row.names = FALSE)
#write.csv(trendsoutsplit[[2]],paste0("trends/All_region_Simulated_",prior,"_",noise_dist_sel,"_level_trends.csv"),row.names = FALSE)

write.csv(TRENDSout,paste0("trends/All_Simulated_",prior,"_",noise_dist_sel,"_survey_wide_trends.csv"),row.names = FALSE)

save(list = c("trend_maps_1980",
              "trend_maps_3gen",
              "trend_maps_L3gen",
              "trend_comparison_1980",
              "trend_comparison_3gen",
              "trend_comparison_L3gen"),
              
              #"composite_trajectories",
              
              #"sp_ind_plots_strat",
              #"sp_ind_plots",
              #"sp_ind_plots_strat_diagnostic",
              #"sp_ind_plots_diagnostic",              
              #"season_graphs",
              #"loo_ic"),
     file = paste0("Figures/All_Simulated_",prior,"_",noise_dist_sel,"_stored_maps.RData"))




# Trend maps and comparisons --------------------------------------------------------------
pdf(file = paste0("Figures/All_Simulated_",prior,"_",noise_dist_sel,"_trend_maps.pdf"),
    width = 9, height = 6.5)
for(sp1 in sps){
   sp = paste0("Simulated_",sp1)
   if(!is.null(trend_maps_1980[[sp]])){
      print(trend_maps_1980[[sp]])
      print(trend_maps_L3gen[[sp]])
      print(trend_maps_3gen[[sp]])
      print(trend_comparison_1980[[sp]])
      print(trend_comparison_3gen[[sp]])
      print(trend_comparison_L3gen[[sp]])
   }
}
dev.off()




# Survey-wide trend comparison --------------------------------------------

LT3_trends <- TRENDSout %>% filter(trend_type %in% c("Long-term","Recent-three-generation")) %>% 
   mutate(species = gsub(pattern = "Simulated_",replacement = "",x = species),
          version = "Simulated constant trend") 

LT3_saved <- saved_trends %>% filter(trend_type %in% c("Long-term","Recent-three-generation")) %>% 
   select(species,trend_type,trend,lci,uci) %>% 
   mutate(version = "Estimated") #,species != "Semipalmated Sandpiper"

LT3_comb <- bind_rows(LT3_trends,LT3_saved)

my_col_sim <-  scale_color_viridis_d(aesthetics = c("colour","fill"),
                                     option = "viridis", begin = 0.2,end = 0.75,direction = -1)


lt3_tplot <- ggplot(data = LT3_trends,aes(x = species,y = trend,colour = trend_type))+
   # geom_point(data = LT3_saved,aes(x = species,y = trend),
   #            size = 4,shape = 3,colour = "black",show.legend = FALSE, inherit.aes = FALSE)+
   geom_pointrange(aes(ymax = uci,ymin = lci,x = species,y = trend,colour = trend_type),
                   inherit.aes = FALSE,position = position_dodge(width = 0.2))+
   geom_abline(slope = 0,intercept = 0,alpha = 0.7)+
   ylab("Trend (%/year)")+
   xlab("")+
   my_col_sim+
   theme_classic()+
   theme(legend.position = "bottom")+
   coord_flip()


lt3_tplot2panel <- ggplot(data = LT3_comb,aes(x = species,y = trend,colour = trend_type))+
   # geom_point(data = LT3_saved,aes(x = species,y = trend),
   #            size = 4,shape = 3,colour = "black",show.legend = FALSE, inherit.aes = FALSE)+
   geom_pointrange(aes(ymax = uci,ymin = lci,x = species,y = trend,colour = trend_type),
                   inherit.aes = FALSE,position = position_dodge(width = 0.2))+
   geom_abline(slope = 0,intercept = 0,alpha = 0.7)+
   ylab("Trend (%/year)")+
   xlab("")+
   my_col_sim+
   theme_classic()+
   theme(legend.position = "bottom")+
   facet_wrap(facets = ~version) +
   coord_flip()



pdf(file = paste0("Figures/All_simulated",prior,"_",noise_dist_sel,"_long_term_3Gen_trends.pdf"),
    height = 9,
    width = 6.5)
print(lt3_tplot)
dev.off()


pdf(file = paste0("Figures/All_simulated_comparison",prior,"_",noise_dist_sel,"_long_term_3Gen_trends.pdf"),
    height = 9,
    width = 13.5)
print(lt3_tplot2panel)
dev.off()



# COSEWIC only ------------------------------------------------------------
LT3_cosewic <- LT3_comb %>% filter(species %in% c(w_cosewic,"Red Knot","Lesser Yellowlegs"))

lt3_tplot2panel_c <- ggplot(data = LT3_cosewic,aes(x = species,y = trend,colour = trend_type))+
   # geom_point(data = LT3_saved,aes(x = species,y = trend),
   #            size = 4,shape = 3,colour = "black",show.legend = FALSE, inherit.aes = FALSE)+
   geom_pointrange(aes(ymax = uci,ymin = lci,x = species,y = trend,colour = trend_type),
                   inherit.aes = FALSE,position = position_dodge(width = 0.2))+
   geom_abline(slope = 0,intercept = 0,alpha = 0.7)+
   ylab("Trend (%/year)")+
   xlab("")+
   my_col_sim+
   theme_classic()+
   theme(legend.position = "bottom")+
   facet_wrap(facets = ~version) +
   coord_flip()

pdf(file = paste0("Figures/All_simulated_comparison_Cosewic_etc_",prior,"_",noise_dist_sel,"_long_term_3Gen_trends.pdf"),
    height = 9,
    width = 13.5)
print(lt3_tplot2panel_c)
dev.off()




# Normal vs t over-dispersion ---------------------------------------------


LT3_saved_t <- saved_trends %>% filter(trend_type %in% c("Long-term","Recent-three-generation")) %>% 
   select(species,trend_type,trend,lci,uci) %>% 
   mutate(version = "t") #,species != "Semipalmated Sandpiper"

LT3_saved_normal <- saved_trends_normal %>% filter(trend_type %in% c("Long-term","Recent-three-generation")) %>% 
   select(species,trend_type,trend,lci,uci) %>% 
   mutate(version = "normal") #,species != "Semipalmated Sandpiper"

LT3_comb_overd <- bind_rows(LT3_saved_t,LT3_saved_normal)

lt3_tplot_overd <- ggplot(data = LT3_comb_overd,aes(x = species,y = trend,colour = trend_type))+
   # geom_point(data = LT3_saved,aes(x = species,y = trend),
   #            size = 4,shape = 3,colour = "black",show.legend = FALSE, inherit.aes = FALSE)+
   geom_pointrange(aes(ymax = uci,ymin = lci,x = species,y = trend,colour = trend_type),
                   inherit.aes = FALSE,position = position_dodge(width = 0.2))+
   geom_abline(slope = 0,intercept = 0,alpha = 0.7)+
   ylab("Trend (%/year)")+
   xlab("")+
   my_col_sim+
   theme_classic()+
   theme(legend.position = "bottom")+
   facet_wrap(facets = ~version) +
   coord_flip()

pdf(file = paste0("Figures/All_simulated_comparison_overdispersion_",prior,"_",noise_dist_sel,"_long_term_3Gen_trends.pdf"),
    height = 9,
    width = 13.5)
print(lt3_tplot_overd)
dev.off()


