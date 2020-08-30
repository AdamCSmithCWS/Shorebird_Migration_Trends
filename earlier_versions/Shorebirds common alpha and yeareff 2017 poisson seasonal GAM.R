### shorebird fitted season model, using multiple counts at each site within seasons
### seasonal pattern not estimated, assume each count is a representation of the average count within each season
### 
setwd("c:/Shorebirds")

library(rjags)

library(RColorBrewer)
library(sp)
library(maptools)
library(rgdal)
gpclibPermit()
fyear <- 1974
# sp2 <- "LEYE"
#sp2 <- "LESA"
#sp2 <- "DUNL"
#sp2 <- "BBPL"

# df <- read.csv("input/Combined_data_2017_horiz.csv",
#                na.strings = "#DIV/0!")
# 
# storun <- names(df)[-c(1:7)]
# for(j in storun){
#   df[is.na(df[,j]),j] <- 0
# }
df <- read.csv("input/Combined_data_2017 vertical.csv")



sps <- read.csv("input/species list.csv")
sps[,"dirnames"] <- gsub(sps$column.names,
                         pattern = ".",
                         replacement = "",
                         fixed = T)
# sps$Species %in% unique(df$Species)
# unique(df$Species) %in% sps$Species


### match the sites to the regions
regions <- readOGR(dsn = "mapping",
                   layer = "region_boundaries_FeatureToP2")

#plot(regions)
sitessp <- unique(df[,c("Locality","Lat","Lon","Country","Prov_State")])  
sitessp$Locality <- as.character(sitessp$Locality)
sitessp[,"Locality2"] <- as.character(sitessp$Locality)
dupsites <- sitessp[duplicated(sitessp$Locality),"Locality"]
sitessp[which(sitessp$Locality %in% dupsites),]


for(s in dupsites){
  tmp <- which(sitessp$Locality %in% s)
  i = 1
  for(j in tmp){
    sitessp[j,"Locality2"] <- paste0(sitessp[j,"Locality"],i)
    i = i+1
  }
  
}

df <- merge(df,sitessp[,c("Locality","Locality2","Lat","Lon")],
            by = c("Locality","Lat","Lon"))

df[,"survey"] <- paste(df[,"Locality2"], df[,"ddmmyyyy"],sep = "---")



sitessp.o <- sitessp
sitessp[,"Region"] <- ""
coordinates(sitessp) <- c("Lon","Lat")
#plot(sitessp, pch = 1)


sitessp@data[,"col"] <- "black"

for (j in 1:length(regions)) {
  p <- regions@data[j,"Region"]
  poly1 <- regions@polygons[[j]]
  in.t <- which(point.in.polygon(point.x = coordinates(sitessp)[,1], point.y = coordinates(sitessp)[,2], pol.x = poly1@Polygons[[1]]@coords[,1], pol.y = poly1@Polygons[[1]]@coords[,2]) %in% c(1,2))
  sitessp@data[in.t,"Region"] <- as.character(p)
  sitessp@data[in.t,"col"] <- rainbow(6)[j]
  sitessp.o[in.t,"Region"] <- p
}


sitessp@data[which(sitessp@data[,"Region"] == "" & sitessp@data[,"Prov_State"] == "US-TX"),"col"] <- rainbow(6)[1]

sitessp@data[which(sitessp@data[,"Region"] == "" & sitessp@data[,"Prov_State"] == "US-TX"),"Region"] <- "Texas Coastal"

sitessp@data[which(sitessp@data[,"Region"] == "" & sitessp@data[,"Prov_State"] == "US-LA"),"col"] <- rainbow(6)[2]

sitessp@data[which(sitessp@data[,"Region"] == "" & sitessp@data[,"Prov_State"] == "US-LA"),"Region"] <- "Southeast Coastal"

# plot(regions,col = "grey")
# plot(sitessp,pch = 19, col = sitessp@data[,"col"], add = T)

sitessp@data[,"Region2"] <- as.character(sitessp@data[,"Region"])
sitessp@data[which(sitessp@data[,"Region"] == "East Inland" &
                     sitessp@data[,"Prov_State"] == "ON"),"Region2"] <- "Ontario"
sitessp@data[which(sitessp@data[,"Region"] == "Northeast Coastal" &
                     sitessp@data[,"Country"] == "Canada"),"Region2"] <- "Atlantic Canada"

sitessp@data[which(sitessp@data[,"Region"] == "Northeast Coastal" &
                     sitessp@data[,"Country"] == "US"),"Region2"] <- "Northeast US Coastal"

i <- 0
for(j in unique(sitessp@data[,"Region2"])){
  i = i+1
  sitessp@data[which(sitessp@data[,"Region2"] == j),"col2"] <- rainbow(9)[i]
}

# x11()
# plot(regions,col = "grey")
# plot(sitessp,pch = 19, col = sitessp@data[,"col2"], add = T)
### test of region assignment

sites <- sitessp@data
sites <- sites[which(sites$Region2 != ""),]

df2 <- merge(df,
             sites[,c("Locality2","Region","Region2")],
             by = c("Locality2"))
### above merge removes data from 20 sites that are outside of the regions
### mostly Alaska and Hawaii, but also 4 sites from Northern Labrador


# sp2 <- storun[26]
#ww <- 22


sspwt <- read.csv("archived species results/Long-term continental shorebird trends.csv")
spwt <- unique(sspwt$species)
spwt <- gsub(spwt,pattern = ".",replacement = "",fixed = T)
spwt <- gsub(spwt,pattern = " ",replacement = "",fixed = T)


storun <- as.character(sps[which(sps$dirnames %in% spwt),"Species"])


# 
# ww = 25
# storuntemp = c("Lesser Yellowlegs",
#                "Whimbrel","Willet",
#                "Black-necked Stilt",
#                "Spotted Sandpiper",
#                "Wilsons Phalarope")
# 



#### generate a table of species*strata combinations
##fill columns indicating the number of counts, and sites, dropped at each step of the data
## exclusion process.

for (sp2 in storun[(ww):(ww+5)]) {
  print(paste(sp2,ww,"starting",Sys.time()))
  
  sp2d <- sps[which(sps$Species == sp2),"dirnames"]
  #sp1 <- snrun[sp2]
  #dir.spsp <- paste0("output/",sp2d)
  dir.spsp <- paste0("output/",sp2d,"poisson seasonal GAM")
  dir.create(dir.spsp)
  dtf <- df2[which(df2$Species == sp2),]
  
  
  ## make list of sites with the species
  siteswsp <- unique(dtf$Locality2)
  ### and then create a complete list of surveys at those sites
  surveysatspsites <- unique(df2[which(df2$Locality2 %in% siteswsp),-which(names(df2) %in% c("Species","Count"))])
  ### then merge the counts and add in all the missing zeros
  fulldtf <- merge(surveysatspsites,
                   dtf[,c("survey","Species","Count")],
                   by = "survey",
                   all.x = T)
  ### replace NAs w zeros
  fulldtf[which(is.na(fulldtf[,"Count"])),"Count"] <- 0
  
  ### remove non-fall migration surveys and surveys before first year
  ### day 182 is July 1
  fulldtf <-  fulldtf[which(fulldtf$Year >= fyear),]
  
  ### replace all the non-fall surveys with zeros
  fulldtf[which(fulldtf$juldate < 182),"Count"] <- 0
  
  
  
  
  meansbyday <- tapply(fulldtf$Count,fulldtf[,c("juldate","Region2")],mean)
  pdf(paste0(dir.spsp,"//mean counts by day during fall ",sp2d,".pdf"))
  par(mfrow = c(4,2))
  for(reg in names(meansbyday[1,])){
    plot(meansbyday[,reg],
         main = reg)
  }
  dev.off()
  
  
  
  ### identifying each regions 95% migration window
  ## and removing counts from outside those windows
  csums <- list()
  length(csums) <- ncol(meansbyday)
  names(csums) <- names(meansbyday[1,])
  migwindows <- data.frame(strat = names(meansbyday[1,]),
                           firstday = 0,
                           lastday = 0)
  for(reg in names(meansbyday[1,])){
    meansbyday[which(is.na(meansbyday[,reg])),reg] <- 0
    csums[[reg]] <- cumsum(meansbyday[,reg])
    mx <- max(csums[[reg]])
    ll <- mx*0.025
    ul <- mx*0.975
    frstd <- min(which(csums[[reg]] >= ll))
    lstd <- max(which(csums[[reg]] <= ul))
    migwindows[which(migwindows$strat == reg),"firstday"] <- frstd
    migwindows[which(migwindows$strat == reg),"lastday"] <- lstd
    
    tmp <- which(fulldtf$Region2 == reg &
                   fulldtf$juldate >= frstd &
                   fulldtf$juldate <= lstd)
    if(reg == names(meansbyday[1,])[1]){
      tokeep <- tmp
    }else{
      tokeep <- c(tokeep,tmp)
    }
    
  }
  fulldtf <- fulldtf[tokeep,]
  
  write.csv(migwindows,
            paste0(dir.spsp,"/migration windows ",sp2d,".csv"),
            row.names = F)
  
  
  
  ### removing sites with less than X years of observations
  minn <- 2 #minimum number of years with nonzero observations  to retain a site - this could be modified to select only sites with some threshold number of years observed.
  
  
  dtf2 <- fulldtf[which(fulldtf$Count > 0),]
  yearbysit <- table(dtf2$Locality2,dtf2$Year)
  yearbysit[which(yearbysit > 0)] <- 1
  
  nyrbysit <- rowSums(yearbysit)
  sitswcrity <- names(nyrbysit)[which(nyrbysit > (minn-1))]  #removing sites with too few years
  sitswtoofewy <- names(nyrbysit)[which(nyrbysit < (minn))]
  
  dts <- fulldtf[which(fulldtf$Locality2 %in% sitswcrity),]
  
  
  
  
  ####### removing sites with too short a span of time covered
  minspan <- 5 # minimum number of years between first and last years surveyed
  rangebysit <- tapply(dts$Year,dts$Locality2,range)
  spanbysit <- integer(length = length(rangebysit))
  names(spanbysit) <- names(rangebysit)
  
  for(st in names(rangebysit)){
    spanbysit[st] <- diff(rangebysit[[st]])
  }
  
  
  
  sitswcritspan <- names(spanbysit)[which(spanbysit >= minspan)]
  
  
  dts <- dts[which(dts$Locality2 %in% sitswcritspan),]
  
  # meansbyday <- tapply(dts$Count,dts[,c("juldate","Region2")],mean)
  # x11()
  # par(mfrow = c(4,2))
  # for(reg in names(meansbyday[1,])){
  #   plot(meansbyday[,reg],
  #        main = reg)
  # }
  # ##
  
  
  

  
  ### idenfity which regions to keep for this species
  
  minsitsbyreg <- 10
  sitsinreg <- unique(dts[,c("Locality2","Region2")])
  nsitsbyreg <- table(sitsinreg$Region2)
  
  regs <- names(nsitsbyreg)[which(nsitsbyreg >= minsitsbyreg)]
  
  dts <- dts[which(dts$Region2 %in% regs),]
  
  
  
  #### regions with a minimum number of years with the species observed
  
  
  dts2 <- dts[which(dts$Count > 0),] 
  nonzerocountsbyregyear <- table(dts[,c("Region2","Year")])
  nonzerocountsbyregyear[which(nonzerocountsbyregyear > 0)] <- 1
  
  nyearswnonzerobyreg <- rowSums(nonzerocountsbyregyear)
  
  minyears <- 15 ### minimum nubmer of years with nonzero counts in strata
  regs <- regs[which(regs %in% names(nyearswnonzerobyreg)[which(nyearswnonzerobyreg >= minyears)])]
  dts <- dts[which(dts$Region2 %in% regs),]
  
  
  
  
  
  stratkeep <- regs
  # this line keeps strata where the species has been observed in more than X years of the time series
  
  
  
  
  if (length(stratkeep) > 0){
    
    
    names(dts)[which(names(dts) == "Locality2")] <- "site"
    names(dts)[which(names(dts) == "Region2")] <- "strat"
    
    
    
    # 
    # #meancountbystratyr <- meancountbystratyr[order(meancountbystratyr[,c("strat","year")]),]
    # stcols <- brewer.pal(8,"Dark2")
    # names(stcols) <- unique(meancountbystratyr$strat)
    # pdf(file = paste(sp.dir,"observed mean counts by strata year NEW2",sp1,".pdf"))
    # par(mar = c(3,4,1,8))
    # for(st in unique(meancountbystratyr$strat)) {
    # 
    # tmp <- meancountbystratyr[which(meancountbystratyr$strat == st),]
    # tmp <- tmp[order(tmp$year),]
    # if (st == unique(meancountbystratyr$strat)[1]) {
    #   if(st %in% stratkeep) {plot(tmp$year,log(tmp$meancount+1,10),type = "l",col = stcols[st],ylim = c(0,4.2),lwd = 3)}else{
    #     plot(tmp$year,log(tmp$meancount+1,10),type = "n",col = stcols[st],ylim = c(0,4.2),lwd = 3)
    #   }}else{  if (st %in% stratkeep){
    # lines(tmp$year,log(tmp$meancount+1,10),col = stcols[st],lwd = 3)}}
    # if (st %in% stratkeep) {
    #   mtext(text = st,side = 4,at = log(tmp$meancount+1,10)[length(tmp$meancount)],line = 0.5,col = stcols[st],las = 1)
    # rm(list = "tmp") }
    # }
    # legend("topleft",legend = paste("-",stratkeep))
    # dev.off()
    # 
    # #if (length(stratkeep > 0){
    #   write.csv(meancountbysityr[which(meancountbysityr$strat %in% stratkeep),],paste(sp.dir,"observed mean counts by site year NEW3",sp1,".csv"))#}else{
    #     #write.csv(meancountbysityr[which(meancountbysityr$strat %in% stratkeep)],paste(sp.dir,"observed mean counts by site year NEW2",sp1,".csv"))
    # #}
    # 
    # #}#temporary end of species loop
    # 
    # 
    # 
    # 
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    #pdf(file = paste(sp.dir,"observed max counts by strata year",sp1,".pdf"))
    #x11()
    #par(mar = c(3,4,1,8))
    #for(st in unique(maxcountbystratyr$strat)) {
    #tmp <- maxcountbystratyr[which(maxcountbystratyr$strat == st),]
    #tmp <- tmp[order(tmp$year),]
    #if (st == unique(maxcountbystratyr$strat)[1]) {plot(tmp$year,log(tmp$maxcount+1,10),type = "l",col = stcols[st],ylim = c(-1,8),lwd = 3)}else{
    #lines(tmp$year,log(tmp$maxcount+1,10),col = stcols[st],lwd = 3)}
    #mtext(text = st,side = 4,at = log(tmp[which.max(tmp$year),"maxcount"]+1,10),line = 0.5,col = stcols[st],las = 1)
    #points(tmp[which.max(tmp$year),"year"],log(tmp[which.max(tmp$year),"maxcount"]+1,10),col = stcols[st],pch = 20)
    ##rm(list = "tmp") 
    #}
    #dev.off()
    
    
    
    #}
    
    
    
    
    
    
    
    if (length(stratkeep) > 2){ ### only runs bugs model for species with 3 or more strata
      
      
      
      ### set up factors to track site and stratum identifications through winbugs
      
      dts[,"stratf"] <- factor(dts[,"strat"])
      dts[,"sitef"] <- factor(dts[,"site"])
      
      dts[,"stratx"] <- as.integer(unclass(dts[,"stratf"]))
      dts[,"sitex"] <- as.integer(unclass(dts[,"sitef"]))
      dts[,"yearx"] <- dts[,"Year"]-(min(dts[,"Year"])-1)
      # dts[,"dateS1"] <- dts[,"date"]-(min(dts[,"date"])-1)
      # dts[,"dateS"] <- (dts[,"dateS1"]-floor(diff(range(dts[,"dateS1"]))/2))/10
      
      dmin <- tapply(dts$juldate,dts$stratx,min)
      for(stt in unique(dts$stratx)){
        dts[which(dts$stratx == stt),"date"] <- (dts[which(dts$stratx == stt),"juldate"]-dmin[as.character(stt)])+1
      }
      
      dmax <- tapply(dts$date,dts$stratx,max)
      
      midyear <- floor(max(dts$yearx)/2)
      
      stratsitlout <- unique(dts[,c("strat","site","stratf","sitef","stratx","sitex")])
      for(sr in unique(stratsitlout$stratx)){
        tmp <- dts[which(dts$stratx == sr),]
        for(si in unique(tmp$sitex)){
          tmp2 <- tmp[which(tmp$sitex == si),]
          stratsitlout[which(stratsitlout$stratx == sr & stratsitlout$sitex == si),"nyears"] <- length(unique(tmp2$year))
          stratsitlout[which(stratsitlout$stratx == sr & stratsitlout$sitex == si),"ncounts"] <- nrow(tmp2)
          
        }
      }
      stratsitlout <- stratsitlout[order(stratsitlout$sitex),]
      write.csv(stratsitlout,paste(dir.spsp,"strata and site lists ",sp2d,".csv",sep = ""),row.names = F)
      ## above writes a text file that tracks site and stratum IDs through winbugs
      
      
      ncounts <- nrow(dts)
      nsites <- max(dts$sitex)
      nyears <- max(dts$yearx)
      nstrata <- max(dts$stratx)
      #nDays <- max(dts$dateS1)
      #days <- as.numeric(unique(dts[,"dateS"]))
      #midday <- floor(diff(range(dts[,"dateS1"]))/2)
      nsitbystr <- as.numeric(tapply(stratsitlout$sitex,stratsitlout$stratx,length))
      nknots = 7
      
      
      ############ GAM basis function
      X.basis = array(NA,dim = c(max(dmax),nknots,nstrata))
      for(stt in 1:nstrata){
        knotsX<- round(seq(1,dmax[stt],length=(nknots+2))[-c(1,nknots+2)])
        X_K<-(abs(outer(1:dmax[stt],knotsX,"-")))^3
        X_OMEGA_all<-(abs(outer(knotsX,knotsX,"-")))^3
        X_svd.OMEGA_all<-svd(X_OMEGA_all)
        X_sqrt.OMEGA_all<-t(X_svd.OMEGA_all$v  %*% (t(X_svd.OMEGA_all$u)*sqrt(X_svd.OMEGA_all$d)))
        X.basist<-t(solve(X_sqrt.OMEGA_all,t(X_K)))
        
        X.basis[1:dmax[stt],,stt] <- X.basist
      }#stt
      
      
      ### create list of data for JAGS
      #
      
      
      ## need to add nknots, x.basis, dmax
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      write.csv(dts,paste0(dir.spsp,"/original data table ",sp2d,".csv"),
                row.names = F)
      #data
      data.jags <- list(count = dts[,"Count"], 
                        year = dts[,"yearx"], 
                        site = dts[,"sitex"],
                        strat = dts[,"stratx"],
                        date = dts[,"date"],
                        ncounts = ncounts,
                        nsites = nsites,
                        nyears = nyears,
                        nstrata = nstrata,  
                        midyear = midyear,
                        nknots = nknots,
                        X.basis = X.basis,
                        dmax = dmax)
      
      
      #   
      # #inits
      # inits <- function(){
      #   
      #   sdyear = runif(1,1,2)
      #   sdbeta = runif(1,0.01,0.02)
      #   sdnoise = runif(1,1,2)
      #   sdsite = runif(1,1,2)
      #   taubetasit = runif(1,100,1000)
      #   alpha = runif(1,-0.5,-0.1)
      #   beta = runif(nstrata,-0.0001,0.0001)
      #   betaH = 0
      #   #bS = c(runif(floor(nknots/2),0.1,0.5),0,runif(floor(nknots/2),-0.5,-0.1))
      #   bet = c(runif(nsites,-0.005,0.005))
      #   #taudate = 1/(runif(1,0.5,1)^2)
      #   
      #   
      #   list(
      #     
      #     sdyear = sdyear,
      #     sdnoise = sdnoise,
      #     sdsite = sdsite,
      #     alpha = alpha,
      #     beta = beta,
      #     sdbetastr = sdbetastr,
      #     betaH = betaH,
      #     #bS1 = bS,
      #     #taudate = taudate,
      #     bet = bet,
      #     taubetasit = taubetasit)
      #   
      # }
      
      # params to monitor
      # params <- c("pchange1","pchange2","pchange3","pchange4","betaH",
      #             "beta","ste","size","ann","sdsite","sdyear","gof","gofnew","posdiff",
      #             "nzeronew","meannew","sizecont","pchange1c",
      #             "pchange2c","pchange3c","pchange4c",
      #             "sdbeta","alpha",
      #             "r") #"u95new",
      params <- c("betaH",
                  "beta","ste","ann","sdsite","sdyear","gof","gofnew","posdiff",
                  "nzeronew","meannew",
                  "sdbeta","sdbetagam","alpha",
                  "zi",
                  "sdnoise",
                  #"B.X",
                  "beta.X") #"u95new",
      mod <- paste(getwd(),"/Shorebirdmodel common alpha and yeareff Poisson seasonal GAM.txt",sep = "")
      
      
      adaptSteps = 500              # Number of steps to "tune" the samplers.
      burnInSteps = 70000            # Number of steps to "burn-in" the samplers.
      nChains = 3                 # Number of chains to run.
      numSavedSteps=2000           # Total number of steps to save.
      thinSteps=50                   # Number of steps to "thin" (1=keep every step).
      nIter = ceiling( ( numSavedSteps * thinSteps ) / nChains ) # Steps per chain.
      
      t1 = Sys.time()
      
      
      jagsMod = jags.model( mod, 
                            data= data.jags ,  
                            #inits= sp.inits,  
                            n.chains= nChains , 
                            n.adapt= adaptSteps )
      # adaptest <- F
      # while(adaptest == F){
      #   
      #   adaptest <- adapt(jagsMod,n.iter = adaptSteps)
      # }
      
      
      # Burn-in:
      cat( "Burning in the MCMC chain...\n" )
      update( jagsMod , n.iter=burnInSteps )
      
      initnew <- coef(jagsMod)
      save(list = "initnew",file = paste0(dir.spsp,"\\initnew.RData",sep = ""))
      
      # The saved MCMC chain:
      cd13 = coda.samples( jagsMod , variable.names=params ,
                           n.iter=nIter , thin=thinSteps )
      
      t2 = Sys.time()
      save(list = "cd13",file = paste(dir.spsp,"\\cd13.RData",sep = ""))
      
      t2-t1
      
      pdf(paste0(dir.spsp,"\\",sp2d,"mcmc diagnostic.pdf"))
      plot(cd13)
      dev.off()  
      
      newinits = coef(jagsMod)
      
      newparms = c("lambda")
      cd13lambda = coda.samples( jagsMod , variable.names= newparms,
                                 n.iter=500 , thin=50 )
      
      save(list = c("cd13lambda",
                    "newinits",
                    "dts",
                    "data.jags"), file = paste(dir.spsp,"\\cd13lambda.RData",sep = ""))
      
    }#end if else 
  }#end if else 
  
    print(paste(sp2,ww,"finished"))
  }#end sp species loop, if used (and uncommented)
  
  
  
################# 

### extract the posdif and other posterior predictive checks for the Neg binomial
### and the poisson w seasonal gam models

### run a negbinomial seasonal gam



  # 
  # 
  # par(mfrow = c(nstrata,1))
  # for(j in 1:nstrata){
  #   #rng = (1+sum(dmax[-c(j:nstrata)])):dmax[j]
  #   rng = paste0("zi[",1:dmax[j],",",j,"]")
  #   plot(invlogit(tmp2[rng,"Mean"]),type = "l")
  # }
  # 
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  