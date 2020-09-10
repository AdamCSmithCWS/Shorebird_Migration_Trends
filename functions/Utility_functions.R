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
      cl2 = str_extract(dat[,cl],pattern = pat2)
      
      d = max(nchar(cl2))
      
      pat = paste0(pat,"[:digit:]{1,",d,"}[:punct:]")
    }
  }
  
  
  pat = paste0(pat,")[:digit:]+")
  dds = as.integer(str_extract(dat[,cl],pattern = pat))
  return(dds)
  
}



texp <- function(x,ny = 2019-1974){
  (x^(1/ny)-1)*100
}


extr_inds <- function(param = "n_s",
                      sumtable = sums,
                      regions = TRUE){
  
  pat <- paste0(param,"[")
  wn_s <- grep(pattern = pat,
               x = sumtable$Parameter,
               fixed = TRUE)
  inds <- sumtable[wn_s,]
  
  if(regions){
    inds[,"s"] <- jags_dim(dim = 1,
                         var = param,
                         dat = inds)
    inds[,"y"] <- jags_dim(dim = 2,
                           var = param,
                           dat = inds)
    inds <- left_join(inds,strats,by = "s")
  
    }else{
  
  inds[,"y"] <- jags_dim(dim = 1,
                         var = param,
                         dat = inds)
    }
  inds$year <- inds$y +1973

  return(inds)
}



ItoT <- function(inds = NSamples,
                 start = 1974,
                 end = 2019,
                 regions = FALSE,
                 qs = 95){
  
  lq = (1-(qs/100))/2
  uq = ((qs/100))+lq
  nyrs = end-start
  
  indt <- inds %>% filter(year %in% c(start,end)) %>% 
    ungroup %>% 
    select(-y) %>% 
    pivot_wider(names_from = year,
                values_from = .value)
  
  indt[,"start"] <- indt[,as.character(start)] 
  indt[,"end"] <- indt[,as.character(end)] 
  
  
  tt <- indt %>% group_by(.iteration) %>% 
    summarise(t = texp(end/start,ny = nyrs),.groups = "keep") %>%
    ungroup() %>% 
    summarise(trend = mean(t),
              lci = quantile(t,lq,names = FALSE),
              uqi = quantile(t,uq,names = FALSE))
  
  if(regions){
    
    indt <- inds %>% filter(year %in% c(start,end)) %>% 
      ungroup %>% 
      select(-y,-s) %>% 
      group_by(Region) %>% 
      pivot_wider(names_from = year,
                  values_from = .value)
    
    indt[,"start"] <- indt[,as.character(start)] 
    indt[,"end"] <- indt[,as.character(end)] 
    
    
    tt1 <- indt %>% group_by(.iteration,Region) %>% 
      summarise(t = texp(end/start,ny = nyrs),.groups = "keep") %>%
      ungroup() %>% 
      group_by(Region) %>% 
      summarise(trend = mean(t),
                lci = quantile(t,lq,names = FALSE),
                uqi = quantile(t,uq,names = FALSE),
                .groups = "keep")
    
    tt2 <- indt %>% group_by(.iteration) %>% 
      summarise(end = sum(end),
                start = sum(start),.groups = "keep") %>% 
      summarise(t = texp(end/start,ny = nyrs),.groups = "keep") %>%
      ungroup() %>% 
      summarise(trend = mean(t),
                lci = quantile(t,lq,names = FALSE),
                uqi = quantile(t,uq,names = FALSE))
    
    tt2$Region <- "Composite"     
    tt <- bind_rows(tt2,tt1)
    
  }
  
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
      select(Region,year) %>% 
      group_by(Region,year) %>% 
      slice_sample(prop = 0.2)
    # summarizing observed values ---------------------------------------------
    
    
    obs <- raw %>% 
      group_by(strat,yr) %>% 
      summarise(n_counts = n(),
                mean_counts = mean(count))
    obs$s = obs$strat
    
    obs <- left_join(obs,strats,by = "s")
    obs$year <- obs$yr+1973
    
  }else{
    annotobs = filter(obs,year == 1980)
    
    
    ss <- raw %>% 
      select(year) %>% 
      group_by(year) %>% 
      slice_sample(prop = 0.1)
    
    obs <- raw %>% 
      group_by(yr) %>% 
      summarise(n_counts = n(),
                mean_counts = mean(count))
    
    obs$year <- obs$yr+1973
    
  
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
      ggplot2::geom_line(data = inds, ggplot2::aes(x = year, y = median,group = Region),colour = cls[1]) +
      ggplot2::geom_ribbon(data = inds, ggplot2::aes(x = year, ymin = lci, ymax = uci,group = Region),alpha = 0.3,fill = cls[1])+
      
      ggplot2::scale_x_continuous(breaks = yys)+
      ggplot2::scale_y_continuous(limits = c(0,NA))
    
    if(!is.null(smooth_inds)){
      
      
      p <- p+ ggplot2::geom_line(data = smooth_inds, ggplot2::aes(x = year, y = median,group = Region),colour = cls[2]) +
        ggplot2::geom_ribbon(data = smooth_inds, ggplot2::aes(x = year, ymin = lci, ymax = uci,group = Region),alpha = 0.3,fill = cls[2])
      
      
    }
    if(add_samplesize){
      
      p <- p + ggplot2::geom_dotplot(data = ss,mapping = ggplot2::aes(x = year,group = Region),drop = TRUE,binaxis = "x", stackdir = "up",method = "histodot",binwidth = 1,width = 0.2,inherit.aes = FALSE,fill = grDevices::grey(0.6),colour = grDevices::grey(0.6),alpha = 0.2,dotsize = 0.3)
      
    }
    if(add_observed){
      p <- p + ggplot2::geom_point(data = obs,ggplot2::aes(x = year,y = mean_counts,group = Region),colour = grDevices::grey(0.6))
        # ggplot2::annotate(geom = "text",x = annotobs$year,y = annotobs$mean_counts,label = "Observed means",colour = grDevices::grey(0.6))
        # 
    }
    
    p <- p + facet_wrap(facets = ~ Region,nrow = 3,scales = "free")
    
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















