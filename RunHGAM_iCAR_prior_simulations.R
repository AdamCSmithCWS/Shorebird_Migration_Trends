#The posterior predictive distribution can be strongly affected by the prior when there is not much observed data
# and substantial prior mass is concentrated around infeasible values


# SPECIES MCMC data-prep -------------------------------------------------------
library(tidyverse)
library(sf)
library(spdep)
library(ggforce)
library(tidybayes)

source("functions/posterior_summary_functions.R")
source("functions/GAM_basis_function_mgcv.R")


# library(doParallel)
# library(foreach)



#load data
load("data/allShorebirdPrismFallCounts.RData")
grid_spacing <- 300000  # size of squares, in units of the CRS (i.e. meters for lae)

 
FYYYY = 1980

output_dir <- "g:/Shorebird_Migration_Trends/output"
#output_dir2 <- "g:/Shorebird_Migration_Trends/output"

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




sp1 = "Semipalmated Plover"
    
    prior = "gamma"
    noise_dist2 = "t"
    
    spf1 = gsub(sp1,pattern = " ",replacement = "_")
    spf1 = gsub(pattern = "\'",replacement = "",
               x = spf1)
    sp_file_name1 <- paste0(spf1,"-",prior,"-",noise_dist2)
    
    load(paste0("data/data",sp1,"_cmdstanr_data.RData"))
    load(paste0(output_dir,"/",sp_file_name1,"_fit_add.RData"))
    
    cmdstanfit <- rstan::read_stan_csv(csvfl) 


# Gather estimates from fitted model --------------------------------------

    alpha_samples <- posterior_samples(fit = cmdstanfit,
                                       parm = "alpha",
                                       dims = c("s"))
    
    ALPHA1_samples <- posterior_samples(fit = cmdstanfit,
                             parm = "ALPHA1",
                             dims = NULL) 
    

    
    sdnoise <- posterior_samples(fit = cmdstanfit,
                                parm = "sdnoise",
                                dims = NULL) %>% 
       summarise(mean = mean(.value)) %>% 
       as.numeric()
    
    sdyear <- posterior_samples(fit = cmdstanfit,
                                 parm = "sdyear",
                                 dims = NULL) %>% 
      summarise(mean = mean(.value)) %>% 
      as.numeric()
    
    
    
  

    B_samples <- posterior_samples(fit = cmdstanfit,
                                  parm = "B",
                                  dims = c("k"))
    

    
    save(list = c("stan_data",
                  "B_samples",
                  "sdyear",
                  "sdnoise",
                  "alpha_samples",
                  "ALPHA1_samples"),
         file = "data/temp_SEPL_output.RData")
    
    #neighbourhood matrix from the original species data
    nmat = matrix(0,nrow = stan_data$nstrata,ncol = stan_data$nstrata)
    for(j in 1:stan_data$N_edges){
      ss = stan_data$node1[j]
      ss2 = stan_data$node2[j]
      nmat[ss,ss2] <- 1
      nmat[ss2,ss] <- 1
      
    }
    
    
    
    Trends_comp_out <- NULL
    #number of simulations = number of posterior draws from the fitted model
    niter = max(B_samples$.draw)
    rates = c(0.5,1,2,4)
    
    for(prior_rate in rates){
      #samples from the prior
    sdyear_gam_sim <- rgamma(niter,2,prior_rate)
    
    #blank array to hold the simulated population trajectories from which the trends are generated
    nsmooth <- array(NA,dim = c(niter,stan_data$nstrata,stan_data$nyears))
    
    for(i in 1:niter){
      
      #blank matrix to hold the yearly predictions of the smooths
      year_pred <- matrix(NA,nrow = stan_data$nyears,ncol = stan_data$nstrata)
      # blank matrix to hold the simulated stratum-level GAM parameters
      b <- matrix(NA,nrow = stan_data$nstrata,ncol = stan_data$nknots_year)
      #scaling the neighbourhood matrix into a covariance matrix 
      sig_mat = (nmat*(sdyear_gam_sim[i]^2)) 
      diag(sig_mat) <- 1
      
      # minimal modifications to ensure the 
      sig_mat = Matrix::nearPD(sig_mat,keepDiag = FALSE)
      
      B <- B_samples %>% filter(.draw == i) 
      B <- B$.value
      
      ALPHA1 <-  ALPHA1_samples %>% filter(.draw == i) 
      ALPHA1 <- ALPHA1$.value
      
      alpha <-  alpha_samples %>% filter(.draw == i) 
      alpha <- alpha$.value
      
      
      
      #generate the simulated GAM coefficients using the mvrnorm function from MASS package.
      for(k in 1:stan_data$nknots_year){
      b[1:stan_data$nstrata,k] <- MASS::mvrnorm(mu = rep(0,stan_data$nstrata),Sigma = sig_mat$mat) + B[k] 
      
      }
    #generate the smooths using simulated coefficients
    for(s in 1:stan_data$nstrata){
      year_pred[,s] = (stan_data$year_basispred %*% (b[s,]))
    }
    
    for(s in 1:stan_data$nstrata){
      

      for(y in 1:stan_data$nyears){
        #temporary vector to hold the site-level predictions for each year, stratum, and site
        # addin gin the site-level variation is not actually necessary,
        # we've retained it here just to make a clear link to the quantities estimated in the model
        atmp_smo <- vector(mode = "numeric",length = stan_data$nsites_strat[s])
        
       
        for(j in 1:stan_data$nsites_strat[s]){
          atmp_smo[j] = exp(ALPHA1 + year_pred[y,s] + 0.5*(sdyear^2) + 0.5*(sdnoise^2) + alpha[stan_data$sites[j,s]])
        }
        nsmooth[i,s,y] = mean(atmp_smo)
      }
      
    }
    
    
    
    
    }#i end of iterations
    
    #temporary matrix to contain the trend calcualtions
    Trends <- matrix(NA,nrow = niter,ncol = stan_data$nstrata)
    # additional matrix to contain trends scaled as the ratio of end/start values of the trajectories
    # will be used to calculate %/change of the population - a sometimes more intuitive metric of 
    # long-term trend
    p_change <- Trends
    for(i in 1:niter){
        for(s in 1:stan_data$nstrata){
          Trends[i,s] <- (((nsmooth[i,s,stan_data$nyears]/nsmooth[i,s,1])^(1/stan_data$nyears))-1)*100
          p_change[i,s] <- ((nsmooth[i,s,stan_data$nyears]/nsmooth[i,s,1]))
        }
    }

      # function to calculate the span i.e., the difference between two percentiles of a vector
    span = function(x,q = 1){
      quantile(x,q)-quantile(x,1-q)
    }

    
    # creating a dataframe to store the summaries of the simulations for this prior
    Trends_df <- as.data.frame(Trends)
    names(Trends_df) <- gsub(names(Trends_df),pattern = "V",replacement = "Strat_")
    Trends_df$prior_scale = paste0("gamma(2,",prior_rate,")")
    Trends_df$span95 = apply(Trends,1,FUN = span,q = 0.95)
    Trends_df$span = apply(Trends,1,FUN = span,q = 1)
    Trends_df$p_ch_1 = apply(p_change,1,FUN = quantile,p = 0.1)
    Trends_df$p_ch_9 = apply(p_change,1,FUN = quantile,p = 0.9)
    
    # stacking the summaries across prior values
    Trends_comp_out <- bind_rows(Trends_comp_out,Trends_df)
    
    }#prior_rate


 #    br = c(1/100000,1/10000,1/1000,1/100,1/10,1)
 #    lb = paste0(c((br-1)*100),"%")
 # ggh_1 = ggplot(data = Trends_comp_out,aes(x = p_ch_1))+
 #   geom_histogram()+
 #   scale_x_log10(breaks = br,labels = lb)+
 #   xlab("10th percentile of stratum-level population change")+
 #   facet_wrap(~prior_scale)
 # 
 # print(ggh_1)
 #  
 # 
 # 
 # br = c(1,10,100,1000,10000)
 # lb = paste0(c((br-1)*100),"%")
 # ggh_9 = ggplot(data = Trends_comp_out,aes(x = p_ch_9))+
 #   geom_histogram()+
 #   scale_x_log10(breaks = br,labels = lb)+
 #   xlab("90th percentile of stratum-level population change")+
 #   facet_wrap(~prior_scale)
 # 
 # print(ggh_9)
 # 
 # 
 ch_stack <- Trends_comp_out %>% select(prior_scale,p_ch_1,p_ch_9) %>%
   pivot_longer(.,cols = starts_with("p_ch"),
                                              names_to = "quantile_population_change",
                                              values_to = "Percent_change")


 
 br = c(c(1/1000000000,1/1000000,1/10000,1/100),c(1,100,10000,1000000,1000000000))
 lb = paste0(c((br-1)*100),"%")
 ggh_p = ggplot(data = ch_stack,aes(x = Percent_change,fill = quantile_population_change))+
   geom_histogram()+
   scale_x_log10(breaks = br,labels = lb)+
   scale_fill_viridis_d(end = 0.8)+
   xlab("10th and 90th percentile of stratum-level population change")+
   facet_wrap(~prior_scale,scales = "free")
 
 print(ggh_p)
 
 
 br = c(0,5,10,20,30,50,100,150)
 lb = paste0(br,"%")
 ggh_t = ggplot(data = Trends_comp_out,aes(x = span95))+
   geom_histogram(bins = 50)+
   scale_x_continuous(breaks = br,labels = lb,limits = c(0,NA))+
   xlab("Range of stratum-level long-term trend estimates (%/year)")+
   facet_wrap(~prior_scale,scales = "free")

 print(ggh_t)
 
 tmp = Trends_comp_out %>% filter(prior_scale == "gamma(2,0.1)")
 hist(tmp$span)

