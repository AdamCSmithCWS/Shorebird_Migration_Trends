##3 utility functions

jags_dim <- function(dim = 1,
                     var = "",
                     cl = "Parameter",
                     dat = NULL){
  ##3 function to extract the indicator value of a multi-dimension jagsUI summary table
  require(stringr)
  
  pat = paste0("(?<=",var,"\\[")
  
  if(dim > 1){
    for(j in 1:(dim-1)){
      
      pat2 = paste0(pat,")[:digit:]+")
      cl2 = str_extract(unlist(dat[,cl]),pattern = pat2)
      
      d = max(nchar(cl2))
      
      pat = paste0(pat,"[:digit:]{1,",d,"}[:punct:]")
    }
  }
  
  
  pat = paste0(pat,")[:digit:]+")
  dds = as.integer(str_extract(unlist(dat[,cl]),pattern = pat))
  return(dds)
  
}


# dim_ext <- function(dim = 1,
#                     var = "",
#                     cl = "Parameter",
#                     dat = NULL){
#   ##3 function to extract the indicator values from cmdstanr output
#   require(stringr)
#   
#   pat = paste0("(?<=",var,"\\[")
#   
#   if(dim > 1){
#     for(j in 1:(dim-1)){
#       
#       pat2 = paste0(pat,")[:digit:]+")
#       cl2 = str_extract(unlist(dat[,cl]),pattern = pat2)
#       
#       d = max(nchar(cl2))
#       
#       pat = paste0(pat,"[:digit:]{1,",d,"}[:punct:]")
#     }
#   }
#   
#   
#   pat = paste0(pat,")[:digit:]+")
#   dds = as.integer(str_extract(unlist(dat[,cl]),pattern = pat))
#   return(dds)
#   
# }
# 

myrename = function(fit){
  rename_with(fit,~ paste0("PI",gsub(".","_",gsub("%", "", .x, fixed = TRUE), fixed = TRUE)),ends_with("%"))
  
}


# 
# cmd_samples <- function(fit = cmdstanfit,
#                         parm = "nsmooth",
#                         dims = NULL){
#   require(posterior)
#   samples <- as_draws_df(fit$draws(variables = c(parm)))
#   if(length(dims) > 0){
#     parm_ex <- paste0(parm,"\\[")
#   }else{
#     parm_ex <- parm
#   }
#   
#   plong <- suppressWarnings(samples %>% pivot_longer(
#     cols = matches(parm_ex,ignore.case = FALSE),
#     names_to = c(".variable"),
#     values_to = ".value",
#     values_drop_na = TRUE
#   )) 
#   
#   for(dn in 1:length(dims)){
#     dd = dims[dn]
#     plong[,dd] = dim_ext(dim = dn,
#                          var = parm,
#                          cl = ".variable",
#                          dat = plong)
#     
#   }
#   
#   plong <- plong %>% mutate(.variable = parm)
#   return(plong)
#   
# }


# cmd_summary <- function(samples,
#                         parm = "nsmooth",
#                         dims = c("site","year"),
#                         probs = c(0.025,0.5,0.975),
#                         outnames = c("lci","median","uci")){
#                           
#                           
#    sums <- summarise_draws(x = samples,~quantile2(.x,probs = probs))
#   names(sums)[2:4] <- outnames 
#   for(dn in 1:length(dims)){
#     dd = dims[dn]
#   sums[,dd] = jags_dim(dim = dn,
#                         var = parm,
#                         cl = "variable",
#                         dat = sums)
# 
#   }
#   return(sums)
#   
#                         }


index_summary <- function(samples = Nsamples,
                          rawdat = dts,
                          parm = "nsmooth",
                          dims = c("site","year"),
                          probs = c(0.025,0.5,0.975),
                          strat_offsets = NULL,
                          site_offsets = NULL)
{
  
  
  indsout = posterior_sums(samples,
                           dims = dims,
                           quantiles = probs)
  
  # indsout = as.data.frame(summary(fit,
  #                                 pars = parm,
  #                                 probs = probs)$summary)
  
  # indsout$year = 1980:2019
  # 
  #indsout = myrename(indsout)#removes the special characters in the column names
  #indsout$Parameter = indsout$variable
  indsout$parm = parm
  
  # for(dn in 1:length(dims)){
  #   dd = dims[dn]
  #   indsout[,dd] = jags_dim(dim = dn,
  #                           var = parm,
  #                           cl = "Parameter",
  #                           dat = indsout)
  #   
  # }
  

    rawdat$count_scale = rawdat$count
  
  
    if(length(dims) == 1){
      if(!is.null(strat_offsets)){
        rawdat <- left_join(rawdat,strat_offsets)
        obs = rawdat %>% group_by(year) %>% 
          summarise(obsmean = mean(count_scale),
                    obslci = quantile(count_scale,0.05),
                    obsuci = quantile(count_scale,0.95),
                    obsmed = median(count_scale),
                    nsurveys = n(),
                    sqrt_n = sqrt(nsurveys),
                    nstrats = length(unique(strat)),
                    mean_counts_incl_strata = mean(adjs))
      }else{
      obs = rawdat %>% group_by(year) %>% 
        summarise(obsmean = mean(count_scale),
                  obslci = quantile(count_scale,0.05),
                  obsuci = quantile(count_scale,0.95),
                  obsmed = median(count_scale),
                  nsurveys = n(),
                  sqrt_n = sqrt(nsurveys),
                  nstrats = length(unique(strat)))
}
      indsout <- left_join(indsout,obs,by = c("year"))
      
    }else{
      if(dims[1] == "stratn"){
        if(!is.null(site_offsets)){
          rawdat <- left_join(rawdat,site_offsets)
        obs = rawdat %>% group_by(stratn,year) %>% 
          summarise(obsmean = mean(count_scale),
                    obslci = quantile(count_scale,0.05),
                    obsuci = quantile(count_scale,0.95),
                    obsmed = median(count_scale),
                    nsurveys = n(),
                    sqrt_n = sqrt(nsurveys),
                    nsites = length(unique(site)),
                    mean_counts_incl_sites = mean(adjs))
        indsout <- left_join(indsout,obs,by = c("stratn",
                                                "year" ))
        }else{
          obs = rawdat %>% group_by(stratn,year) %>% 
            summarise(obsmean = mean(count_scale),
                      obslci = quantile(count_scale,0.05),
                      obsuci = quantile(count_scale,0.95),
                      obsmed = median(count_scale),
                      nsurveys = n(),
                      sqrt_n = sqrt(nsurveys),
                      nsites = length(unique(site)))
        } 
      }else{
        obs = rawdat %>% group_by(site,year) %>% 
          summarise(obsmean = mean(count_scale),
                    obslci = quantile(count_scale,0.05),
                    obsuci = quantile(count_scale,0.95),
                    obsmed = median(count_scale),
                    nsurveys = n(),
                    sqrt_n = sqrt(nsurveys))
        indsout <- left_join(indsout,obs,by = c("site",
                                                "year"))
      }
    }
    
    
 
  
  
  return(indsout)
}




texp <- function(x,ny = 2019-1974){
  (x^(1/ny)-1)*100
}




chng <- function(x){
  (x-1)*100
}




extr_sum <- function(param = "vis.sm_season",
                      sumtable = sums,
                      index = "day",
                     log_retrans = FALSE){
  
  pat <- paste0(param,"[")
  wn_s <- grep(pattern = pat,
               x = sumtable$Parameter,
               fixed = TRUE)
  inds <- sumtable[wn_s,]
  

    for(i in 1:length(index)){
    inds[,index[i]] <- jags_dim(dim = i,
                           var = param,
                           dat = inds)
    }
 if(log_retrans){
   for(pp in c("mean","sd","lci","lqrt","median","uqrt","uci")){
    inds[,pp] <- exp(inds[,pp])
  }
 }
  
  return(inds)
}


ItoT <- function(inds = Nsamples,
                 start = syear,
                 end = 2019,
                 regions = "hex_name",
                 qs = 95,
                 sp = "species",
                 type = "",
                 raw_data = dts,
                 centered_trends = TRUE){
  
  
  varbl <- unique(inds$.variable)
  
  lq = (1-(qs/100))/2
  uq = ((qs/100))+lq
  nyrs = end-start
  
  if(is.null(regions)){
    indt <- inds %>% filter(year %in% c(start,end)) %>% 
    ungroup() %>% 
    select(-y) %>% 
    pivot_wider(names_from = year,
                values_from = .value)
  
  indt[,"start"] <- indt[,as.character(start)] 
  indt[,"end"] <- indt[,as.character(end)] 
  
  
  tt <- indt %>% group_by(.draw) %>% 
    summarise(t = texp(end/start,ny = nyrs),
              ch = chng(end/start),
              .groups = "keep") %>%
    ungroup() %>% 
    summarise(trend = mean(t),
              lci = quantile(t,lq,names = FALSE),
              uci = quantile(t,uq,names = FALSE),
              percent_change = median(ch),
              p_ch_lci = quantile(ch,lq,names = FALSE),
              p_ch_uci = quantile(ch,uq,names = FALSE))
  
  mn <- inds %>% filter(year %in% c(start:end)) %>% 
    ungroup() %>%
    group_by(year) %>% 
    summarise(ym = mean(.value)) %>% 
    ungroup() %>% 
    summarise(mean_abundance = mean(ym))
  
  obs <- raw_data %>% filter(year %in% c(start:end)) %>% 
    group_by(year) %>% 
    summarise(ym = mean(count),
              n_c = n()) %>% 
    ungroup() %>% 
    summarise(mean_observed = mean(ym),
              mean_n_surveys = mean(n_c))
  
  tt <- bind_cols(tt,mn,obs)
  tt$trend_type <- type
  tt$species <- sp
  }
  
  if(!is.null(regions)){
    inds[,"region"] = inds[,regions]
    indt <- inds %>% filter(year %in% c(start,end)) %>% 
      ungroup %>% 
      select(-y,-stratn) %>% 
      group_by(region) %>% 
      pivot_wider(names_from = year,
                  values_from = .value)
    
    indt[,"start"] <- indt[,as.character(start)] 
    indt[,"end"] <- indt[,as.character(end)] 
    
    
    tt1 <- indt %>% group_by(.draw,region) %>%
      summarise(end = sum(end),
                start = sum(start),.groups = "keep") %>%  
      summarise(t = texp(end/start,ny = nyrs),
                ch = chng(end/start),
                .groups = "keep") %>%
      ungroup() %>% 
      group_by(region) %>% 
      summarise(trend = mean(t),
                lci = quantile(t,lq,names = FALSE),
                uci = quantile(t,uq,names = FALSE),
                percent_change = median(ch),
                p_ch_lci = quantile(ch,lq,names = FALSE),
                p_ch_uci = quantile(ch,uq,names = FALSE),
                .groups = "keep")
    mn <- inds %>% filter(year %in% c(start:end)) %>% 
      ungroup() %>%
      group_by(region,year) %>% 
      summarise(ym = mean(.value)) %>% 
      ungroup() %>% 
      group_by(region) %>% 
      summarise(mean_abundance = mean(ym))
    
    obs <- raw_data %>% mutate(region = hex_name) %>% 
      filter(year %in% c(start:end)) %>% 
      group_by(region,year) %>% 
      summarise(ym = mean(count),
                n_c = n()) %>% 
      ungroup() %>% 
      group_by(region) %>% 
      summarise(mean_observed = mean(ym),
                mean_n_surveys = mean(n_c))
    
    tt1 <- left_join(tt1,mn,by = "region")
    tt1 <- left_join(tt1,obs,by = "region")
    
    
    tt2 <- indt %>% group_by(.draw) %>% 
      summarise(end = sum(end),
                start = sum(start),.groups = "keep") %>%       
      summarise(t = texp(end/start,ny = nyrs),
               ch = chng(end/start),
               .groups = "keep") %>%
      ungroup() %>% 
      summarise(trend = mean(t),
                lci = quantile(t,lq,names = FALSE),
                uci = quantile(t,uq,names = FALSE),
                percent_change = median(ch),
                p_ch_lci = quantile(ch,lq,names = FALSE),
                p_ch_uci = quantile(ch,uq,names = FALSE))
    
    mn <- inds %>% filter(year %in% c(start:end)) %>% 
      ungroup() %>%
      group_by(year) %>% 
      summarise(ym = mean(.value)) %>% 
      ungroup() %>% 
      summarise(mean_abundance = mean(ym))
    
    obs <- raw_data %>% filter(year %in% c(start:end)) %>% 
      group_by(year) %>% 
      summarise(ym = mean(count),
                n_c = n()) %>% 
      ungroup() %>% 
      summarise(mean_observed = mean(ym),
                mean_n_surveys = mean(n_c))
    tt2 <- bind_cols(tt2,mn,obs)
    tt2$region <- "Composite"  
    
    if(centered_trends){
      ttt1 <- indt %>% group_by(.draw,region) %>%
        summarise(end = sum(end),
                  start = sum(start),.groups = "keep") %>%  
        summarise(t_c = log(end/start),
                  .groups = "keep")
      
      ttt2 <- indt %>% group_by(.draw) %>% 
        summarise(end = sum(end),
                  start = sum(start),.groups = "keep") %>%       
        summarise(t = log(end/start),
                  .groups = "keep")
      
      ttt1 <- left_join(ttt1,ttt2,by = ".draw")
      ttt1 <- ttt1 %>% 
        mutate(td = t_c - t) %>% 
        ungroup() %>% 
        group_by(region) %>% 
        summarise(centered_log_trend = mean(td),
                  centered_log_trend_lci = quantile(td,lq,names = FALSE),
                  centered_log_trend_uci = quantile(td,uq,names = FALSE))
      tt1 <- left_join(tt1,ttt1,by = "region")
      
    }
    
    tt <- bind_rows(tt2,tt1)
    
    
  }
  tt$start_year <- start
  tt$end_year <- end
  tt$species <- sp
  tt$region_type <- regions
  tt$trend_type <- type
  
  tt$parameter = varbl
  
  return(tt)
}



slope_trend <- function(x,y){
  x = log(x)
    n = length(y)
    sx = sum(x)
    sy = sum(y)
    ssy = sum(y^2)
    sxx = sum(x*y)
    b = (n*sxx - sy*sx)/(n*ssy - sy^2)
    return(b)
  }

  
  sltexp <- function(x){((exp(x)-1)*100)} 



ItoT_slope <- function(inds = NSamples,
                 start = fyear,
                 end = 2019,
                 regions = FALSE,
                 qs = 95,
                 trend_type = "slope",
                 index_type = "standard",
                 retransformation_type = "standard"){
  
  varbl <- unique(inds$.variable)
  lq = (1-(qs/100))/2
  uq = ((qs/100))+lq
  nyrs = end-start
  
  tt <- inds %>% 
    filter(year %in% c(start:end)) %>% 
    group_by(.draw) %>% 
    summarise(t = sltexp(slope_trend(x = .value,y = year)),.groups = "keep") %>%
    ungroup() %>% 
    summarise(trend = mean(t),
              lci = quantile(t,lq,names = FALSE),
              uci = quantile(t,uq,names = FALSE))
  
  if(regions){
    
    
    tt1 <- inds %>% 
      filter(year %in% c(start:end)) %>% 
      # ungroup %>% 
      # select(-y,-s) %>% 
      group_by(.draw,region) %>% 
      summarise(t = sltexp(slope_trend(x = .value,y = year)),.groups = "keep") %>%
      group_by(region) %>% 
      #ungroup() %>% 
      summarise(trend = mean(t),
                lci = quantile(t,lq,names = FALSE),
                uci = quantile(t,uq,names = FALSE),
                .groups = "keep")
 
    
    tt2 <- inds %>% 
      filter(year %in% c(start:end)) %>% 
      # ungroup %>% 
      # select(-y,-s) %>% 
      group_by(.draw,year) %>% 
      summarise(sm = sum(.value),.groups = "keep") %>%
      group_by(.draw) %>% 
      summarise(t = sltexp(slope_trend(x = sm,y = year)),.groups = "keep") %>%
      ungroup() %>% 
      summarise(trend = mean(t),
                lci = quantile(t,lq,names = FALSE),
                uci = quantile(t,uq,names = FALSE),
                .groups = "drop")
    

    
    tt2$region <- "Composite"     
    tt <- bind_rows(tt2,tt1)

  }
  
  tt$start_year <- start
  tt$end_year <- end
  tt$trend_type <- trend_type
  tt$index_type <- index_type
  tt$retransformation_type <- retransformation_type
  tt$parameter = varbl
  
  return(tt)
}





plot_ind <- function(inds = N_inds,
                     smooth_inds = NULL,
                     #obs = OBS,
                     raw = dts,
                     add_observed = TRUE,
                     add_samplesize = TRUE,
                     species = sp,
                     regions = FALSE,
                     title_size = 20,
                     axis_title_size = 18,
                     axis_text_size = 16){
  
  require(ggplot2)
  require(ggrepel)
  require(ggforce)
  
  
  
  if(regions){
    
      ss <- raw %>% 
      select(region,year) %>% 
      group_by(region,year) %>% 
      slice_sample(prop = 0.2)
    # summarizing observed values ---------------------------------------------
    
    
    obs <- raw %>% 
      group_by(strat,yr) %>% 
      summarise(n_counts = n(),
                mean_counts = mean(count),
                lqrt_counts = quantile(count,probs = 0.25,na.rm = T),
                uqrt_counts = quantile(count,probs = 0.75,na.rm = T))
    obs$s = obs$strat
    
    obs <- left_join(obs,strats,by = "s")
    obs$year <- obs$yr+1973
    
  }else{
    
    
    ss <- raw %>% 
      select(year) %>% 
      group_by(year) %>% 
      slice_sample(prop = 0.1)
    
    obs <- raw %>% 
      group_by(yr) %>% 
      summarise(n_counts = n(),
                mean_counts = mean(count),
                lqrt_counts = quantile(count,probs = 0.25,na.rm = T),
                uqrt_counts = quantile(count,probs = 0.75,na.rm = T))
    
    obs$year <- obs$yr+1973
    annotobs = filter(obs,year == 1980)
    
  
  }

  cls = scales::viridis_pal(begin = 0.2,end = 0.8)(2)
  
  yys = c(1974,seq(1985,2015,by = 10),2019)
  
  if(regions){
    
    p <- ggplot2::ggplot() +
      ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                     panel.grid.minor = ggplot2::element_blank(),
                     panel.background = ggplot2::element_blank(),
                     axis.line = element_line(colour = "black"),
                     plot.title = ggplot2::element_text(size = title_size),
                     axis.title = ggplot2::element_text(size = axis_title_size),
                     axis.text = ggplot2::element_text(size = axis_text_size)) +
      ggplot2::labs(title = paste(species, " Trajectory ",
                                  sep = ""),
                    x = "Year",
                    y = "Annual index of abundance (mean count)") +
      ggplot2::geom_line(data = inds, ggplot2::aes(x = year, y = median,group = region),colour = cls[1]) +
      ggplot2::geom_ribbon(data = inds, ggplot2::aes(x = year, ymin = lci, ymax = uci,group = region),alpha = 0.3,fill = cls[1])+
      
      ggplot2::scale_x_continuous(breaks = yys)+
      ggplot2::scale_y_continuous(limits = c(0,NA))
    
    if(!is.null(smooth_inds)){
      
      
      p <- p+ ggplot2::geom_line(data = smooth_inds, ggplot2::aes(x = year, y = median,group = region),colour = cls[2]) +
        ggplot2::geom_ribbon(data = smooth_inds, ggplot2::aes(x = year, ymin = lci, ymax = uci,group = region),alpha = 0.3,fill = cls[2])
      
      
    }
    if(add_samplesize){
      
      p <- p + ggplot2::geom_dotplot(data = ss,mapping = ggplot2::aes(x = year,group = region),drop = TRUE,binaxis = "x", stackdir = "up",method = "histodot",binwidth = 1,width = 0.2,inherit.aes = FALSE,fill = grDevices::grey(0.6),colour = grDevices::grey(0.6),alpha = 0.2,dotsize = 0.3)
      
    }
    if(add_observed){
      p <- p + 
        ggplot2::geom_pointrange(data = obs,ggplot2::aes(x = year,y = mean_counts,ymin = lqrt_counts,ymax = uqrt_counts,group = region),colour = grDevices::grey(0.6))
        # ggplot2::annotate(geom = "text",x = annotobs$year,y = annotobs$mean_counts,label = "Observed means",colour = grDevices::grey(0.6))
        # 
    }
    
    p <- p + facet_wrap(facets = ~ region,nrow = 3,scales = "free")
    
  }else{
    p <- ggplot2::ggplot() +
      ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                     panel.grid.minor = ggplot2::element_blank(),
                     panel.background = ggplot2::element_blank(),
                     axis.line = element_line(colour = "black"),
                     plot.title = ggplot2::element_text(size = title_size),
                     axis.title = ggplot2::element_text(size = axis_title_size),
                     axis.text = ggplot2::element_text(size = axis_text_size)) +
      ggplot2::labs(title = paste(species, " Trajectory ",
                                  sep = ""),
                    x = "Year",
                    y = "Annual index of abundance (mean count)") +
      ggplot2::geom_line(data = inds, ggplot2::aes(x = year, y = median),colour = cls[1]) +
      ggplot2::geom_ribbon(data = inds, ggplot2::aes(x = year, ymin = lci, ymax = uci),alpha = 0.3,fill = cls[1])+
      
      ggplot2::scale_x_continuous(breaks = yys)+
      ggplot2::scale_y_continuous(limits = c(0,NA))
    
    if(!is.null(smooth_inds)){
      
      
      p <- p+ ggplot2::geom_line(data = smooth_inds, ggplot2::aes(x = year, y = median),colour = cls[2]) +
        ggplot2::geom_ribbon(data = smooth_inds, ggplot2::aes(x = year, ymin = lci, ymax = uci),alpha = 0.3,fill = cls[2])
      
      
    }
    if(add_samplesize){
      
      p <- p + ggplot2::geom_dotplot(data = ss,mapping = ggplot2::aes(x = year),drop = TRUE,binaxis = "x", stackdir = "up",method = "histodot",binwidth = 1,width = 0.2,inherit.aes = FALSE,fill = grDevices::grey(0.6),colour = grDevices::grey(0.6),alpha = 0.2,dotsize = 0.3)
      
    }
    if(add_observed){
      p <- p + ggplot2::geom_point(data = obs,ggplot2::aes(x = year,y = mean_counts),colour = grDevices::grey(0.6))+
        ggplot2::annotate(geom = "text",x = annotobs$year,y = annotobs$mean_counts,label = "Observed means",colour = grDevices::grey(0.6))
      
    }
    
  }
  
  return(p)
  
}#end index plot










ItoI <- function(inds = nsmoothsamples,
                 start = syear,
                 end = 2019,
                 regions = "hex_name",
                 qs = 95){
  
  
  varbl <- unique(inds$.variable)
  
  lq = (1-(qs/100))/2
  uq = ((qs/100))+lq
  nyrs = end-start
  

  
  if(!is.null(regions)){
    inds[,"region"] = inds[,regions]
    
     tt1 <- inds %>% group_by(.draw,region,year) %>%
      summarise(i = sum(.value),.groups = "keep") %>%  
      ungroup() %>% 
      group_by(region,year) %>% 
      summarise(median = median(i),
                mean = mean(i),
                lci = quantile(i,lq,names = FALSE),
                uci = quantile(i,uq,names = FALSE),
                .groups = "keep")
    
    
    tt2 <- inds %>% group_by(.draw,year) %>%
      summarise(i = sum(.value),.groups = "keep") %>%  
      ungroup() %>% 
      group_by(year) %>% 
      summarise(median = median(i),
                mean = mean(i),
                lci = quantile(i,lq,names = FALSE),
                uci = quantile(i,uq,names = FALSE),
                .groups = "keep")
    
    
    
    tt2$region <- "Composite"     
    tt <- bind_rows(tt2,tt1)
    
  }
  
  return(tt)
}





trend_map = function(
  trends = t_nsmooth_strat_80,
  map.file = "BBS_ProvState_strata",
  hex_map = real_grid,
  size_value = "Mean Abundance",
  raw = dts)
{
  laea = st_crs("+proj=laea +lat_0=40 +lon_0=-95") # Lambert equal area coord reference system
  
  locat = system.file("maps",
                      package = "bbsBayes")
  
  strata_map = read_sf(dsn = locat,
                       layer = map.file)
  strata_map = st_transform(strata_map,crs = laea)
  
  site_loc <- raw %>% distinct(SurveyAreaIdentifier,DecimalLatitude,DecimalLongitude,year) %>% 
    group_by(SurveyAreaIdentifier,DecimalLatitude,DecimalLongitude) %>% 
    summarise(n_years = n())
  
  site_loc <- st_as_sf(site_loc,coords = c("DecimalLongitude","DecimalLatitude"), crs = 4326)
  
  site_loc <- st_transform(site_loc,crs = laea)
  #join the hex map with trends
  #bring in the bbsBayes colour ramp
  #map short and long-term versions
  
  centres = suppressWarnings(st_centroid(hex_map))
  t_map <- left_join(centres,trends,by = c("hex_name" = "region"))
  
  
  
  
  breaks <- c(-7, -4, -2, -1, -0.5, 0.5, 1, 2, 4, 7)
  labls = c(paste0("< ",breaks[1]),paste0(breaks[-c(length(breaks))],":", breaks[-c(1)]),paste0("> ",breaks[length(breaks)]))
  labls = paste0(labls, " %")
  map_palette <- c("#a50026", "#d73027", "#f46d43", "#fdae61", "#fee090", "#ffffbf",
                   "#e0f3f8", "#abd9e9", "#74add1", "#4575b4", "#313695")
  names(map_palette) <- labls
  
  
  t_map$Tplot <- as.numeric(as.character(t_map$trend))
  
  t_map$Tplot <- cut(t_map$Tplot,breaks = c(-Inf, breaks, Inf),labels = labls)
  
  if(grepl(pattern = "bund",size_value)){
    t_map = t_map %>% mutate(size_s = mean_abundance)}
  if(grepl(pattern = "umber",size_value)){
    t_map = t_map %>% mutate(size_s = mean_n_surveys) 
  }
  if(grepl(pattern = "bserve",size_value)){
    t_map = t_map %>% mutate(size_s = mean_observed) 
  }
  
  fyr = unique(t_map$start_year)
  lyr = unique(t_map$end_year)
  
  ptit = paste(sp,"trends",fyr,"-",lyr)
  

  bb = st_bbox(hex_map)
  xl = c(bb["xmin"],bb["xmax"])
  yl = c(bb["ymin"],bb["ymax"])
  
  tmap = ggplot()+

    geom_sf(data = strata_map,alpha = 0,colour = grey(0.8))+
    geom_sf(data = hex_map,alpha = 0,colour = grey(0.9))+
    geom_sf(data = t_map,aes(colour = Tplot,size = size_s),alpha = 1)+
    geom_sf(data = site_loc,alpha = 0.3,size = 1)+
    labs(title = paste0(sp," Trends ",fyr,"-",lyr))+
    coord_sf(xlim = xl,ylim = yl)+
    theme_bw()+
    scale_colour_manual(values = map_palette, aesthetics = c("colour"),
                        guide = ggplot2::guide_legend(reverse=TRUE),
                        name = paste0("Trend\n",fyr,"-",lyr))+
    scale_size(name = size_value,range = c(3,9))
  
  
  return(tmap)
  
}#end function








trend_map_composite = function(
  trends = strat_sums,
  map.file = "BBS_ProvState_strata",
  hex_map = poly_grid,
  size_value = "Probability > or < zero",
  tlab = time)
{
  laea = st_crs("+proj=laea +lat_0=40 +lon_0=-95") # Lambert equal area coord reference system
  
  locat = system.file("maps",
                      package = "bbsBayes")
  
  strata_map = read_sf(dsn = locat,
                       layer = map.file)
  strata_map = st_transform(strata_map,crs = laea)
  
 #join the hex map with trends
  #bring in the bbsBayes colour ramp
  #map short and long-term versions
  if(nrow(hex_map) > nrow(trends)){
    hex_map <- filter(hex_map,hex_name %in% trends$region)
  }
  
  centres = suppressWarnings(st_centroid(hex_map))
  t_map <- left_join(centres,trends,by = c("hex_name" = "region"))
  
  
  
  
  #breaks <- c(-7, -4, -2, -1, -0.5, 0.5, 1, 2, 4, 7)
  breaks <- c(-2.5,-1.5, -1, -0.5,-0.25, 0.25, 0.5, 1, 1.5, 2.5)
  
  labls = c(paste0("< ",breaks[1]),paste0(breaks[-c(length(breaks))],":", breaks[-c(1)]),paste0("> ",breaks[length(breaks)]))
  labls = paste0(labls, " %")
  map_palette <- c("#a50026", "#d73027", "#f46d43", "#fdae61", "#fee090", "#ffffbf",
                   "#e0f3f8", "#abd9e9", "#74add1", "#4575b4", "#313695")
  names(map_palette) <- labls
  
  
  t_map$Tplot <- as.numeric(as.character(t_map$trend))
  
  t_map$Tplot <- cut(t_map$Tplot,breaks = c(-Inf, breaks, Inf),labels = labls)
  # 
  # if(grepl(pattern = "bund",size_value)){
  #   t_map = t_map %>% mutate(size_s = mean_abundance)}
  # if(grepl(pattern = "umber",size_value)){
  #   t_map = t_map %>% mutate(size_s = mean_n_surveys) 
  # }
  # if(grepl(pattern = "bserve",size_value)){
  #   t_map = t_map %>% mutate(size_s = mean_observed) 
  # }
  
  if(grepl(pattern = "pecies",size_value)){
    t_map = t_map %>% mutate(size_s = nspecies) 
  }
  if(grepl(pattern = "zero",size_value)){
    t_map = t_map %>% mutate(size_s = p_not_zero) 
  }
  
  
  fyr = unique(t_map$start_year)
  lyr = unique(t_map$end_year)
  
  # ptit = paste(sp,"trends",fyr,"-",lyr)
  # 
  
  bb = st_bbox(hex_map)
  xl = c(bb["xmin"],bb["xmax"])
  yl = c(bb["ymin"],bb["ymax"])
  
  tmap = ggplot()+
    
    geom_sf(data = strata_map,alpha = 0,colour = grey(0.8))+
    geom_sf(data = hex_map,alpha = 0,colour = grey(0.9))+
    geom_sf(data = t_map,aes(colour = Tplot,size = size_s),alpha = 1)+
    labs(title = paste0("Mean regional variation in trend ",tlab))+
    coord_sf(xlim = xl,ylim = yl)+
    theme_bw()+
    scale_colour_manual(values = map_palette, aesthetics = c("colour"),
                        guide = ggplot2::guide_legend(reverse=TRUE),
                        name = paste0(" Difference\n","Trend"))+
    scale_size(name = size_value,range = c(3,9))
  
  
  return(tmap)
  
}#end function








