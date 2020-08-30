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

siteexc = expand.grid(region = c(unique(sites$Region2),"Continental"),
                      species = storun)


for (sp2 in storun) {
  print(paste(sp2,ww,"starting",Sys.time()))
  
  sp2d <- sps[which(sps$Species == sp2),"dirnames"]
  #sp1 <- snrun[sp2]
  #dir.spsp <- paste0("output/",sp2d)
  dir.spsp <- paste0("output/",sp2d,"poisson seasonal GAM")

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
  # pdf(paste0(dir.spsp,"//mean counts by day during fall ",sp2d,".pdf"))
  # par(mfrow = c(4,2))
  # for(reg in names(meansbyday[1,])){
  #   plot(meansbyday[,reg],
  #        main = reg)
  # }
  # dev.off()
  # 
  # 
  
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
  # 
  # write.csv(migwindows,
  #           paste0(dir.spsp,"/migration windows ",sp2d,".csv"),
  #           row.names = F)
  # 
  
  for(rg in unique(dtf$Region2)){
    ii = which(siteexc$region == rg & siteexc$species == sp2)
    tmpp = fulldtf[which(fulldtf$Region2 == rg),]
    hdr = "total"
    siteexc[ii,paste("n.counts",hdr)] = nrow(tmpp)
    siteexc[ii,paste("n.nonzero.counts",hdr)] = length(which(tmpp$Count > 0))
    siteexc[ii,paste("n.sites",hdr)] = length(unique(tmpp$Locality2))
    
    
  }
  rg = "Continental"
  ii = which(siteexc$region == rg & siteexc$species == sp2)
  tmpp = fulldtf
  hdr = "total"
  siteexc[ii,paste("n.counts",hdr)] = nrow(tmpp)
  siteexc[ii,paste("n.nonzero.counts",hdr)] = length(which(tmpp$Count > 0))
  siteexc[ii,paste("n.sites",hdr)] = length(unique(tmpp$Locality2))
  
  
  
  ### removing sites with less than X years of observations
  minn <- 2 #minimum number of years with nonzero observations  to retain a site - this could be modified to select only sites with some threshold number of years observed.
  
  
  dtf2 <- fulldtf[which(fulldtf$Count > 0),]
  yearbysit <- table(dtf2$Locality2,dtf2$Year)
  yearbysit[which(yearbysit > 0)] <- 1
  
  nyrbysit <- rowSums(yearbysit)
  sitswcrity <- names(nyrbysit)[which(nyrbysit > (minn-1))]  #removing sites with too few years
  sitswtoofewy <- names(nyrbysit)[which(nyrbysit < (minn))]
  
  dts <- fulldtf[which(fulldtf$Locality2 %in% sitswcrity),]
  
  for(rg in unique(dtf$Region2)){
    ii = which(siteexc$region == rg & siteexc$species == sp2)
    tmpp = dts[which(dts$Region2 == rg),]
    hdr = "removed after excl. sites w 1 observation"
    hdrl1 = "total"
    siteexc[ii,paste("n.counts",hdr)] = siteexc[ii,paste("n.counts",hdrl1)] - nrow(tmpp)
    siteexc[ii,paste("n.nonzero.counts",hdr)] = siteexc[ii,paste("n.nonzero.counts",hdrl1)] - length(which(tmpp$Count > 0))
    siteexc[ii,paste("n.sites",hdr)] = siteexc[ii,paste("n.sites",hdrl1)] - length(unique(tmpp$Locality2))
    
  }
  rg = "Continental"
  ii = which(siteexc$region == rg & siteexc$species == sp2)
  tmpp = dts
  hdr = "removed after excl. sites w 1 observation"
  hdrl1 = "total"
  siteexc[ii,paste("n.counts",hdr)] = siteexc[ii,paste("n.counts",hdrl1)] - nrow(tmpp)
  siteexc[ii,paste("n.nonzero.counts",hdr)] = siteexc[ii,paste("n.nonzero.counts",hdrl1)] - length(which(tmpp$Count > 0))
  siteexc[ii,paste("n.sites",hdr)] = siteexc[ii,paste("n.sites",hdrl1)] - length(unique(tmpp$Locality2))
  
  
  
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
  
  for(rg in unique(dtf$Region2)){
    ii = which(siteexc$region == rg & siteexc$species == sp2)
    tmpp = dts[which(dts$Region2 == rg),]
    hdrl2 = "removed after excl. sites w 1 observation"
    hdr = "removed after excl. sites w less than 5yr span"
    siteexc[ii,paste("n.counts",hdr)] = siteexc[ii,paste("n.counts",hdrl1)]-siteexc[ii,paste("n.counts",hdrl2)] - nrow(tmpp)
    siteexc[ii,paste("n.nonzero.counts",hdr)] = siteexc[ii,paste("n.nonzero.counts",hdrl1)] - siteexc[ii,paste("n.nonzero.counts",hdrl2)] - length(which(tmpp$Count > 0))
    siteexc[ii,paste("n.sites",hdr)] = siteexc[ii,paste("n.sites",hdrl1)] - siteexc[ii,paste("n.sites",hdrl2)] - length(unique(tmpp$Locality2))
    
    
  }
  rg = "Continental"
  ii = which(siteexc$region == rg & siteexc$species == sp2)
  tmpp = dts
  hdrl2 = "removed after excl. sites w 1 observation"
  hdr = "removed after excl. sites w less than 5yr span"
  siteexc[ii,paste("n.counts",hdr)] = siteexc[ii,paste("n.counts",hdrl1)]-siteexc[ii,paste("n.counts",hdrl2)] - nrow(tmpp)
  siteexc[ii,paste("n.nonzero.counts",hdr)] = siteexc[ii,paste("n.nonzero.counts",hdrl1)] - siteexc[ii,paste("n.nonzero.counts",hdrl2)] - length(which(tmpp$Count > 0))
  siteexc[ii,paste("n.sites",hdr)] = siteexc[ii,paste("n.sites",hdrl1)] - siteexc[ii,paste("n.sites",hdrl2)] - length(unique(tmpp$Locality2))
  
  

  
  for(rg in unique(dtf$Region2)){
    ii = which(siteexc$region == rg & siteexc$species == sp2)
    tmpp = dts[which(dts$Region2 == rg),]
    hdr = "final"
    siteexc[ii,paste("n.counts",hdr)] = nrow(tmpp)
    siteexc[ii,paste("n.nonzero.counts",hdr)] = length(which(tmpp$Count > 0))
    siteexc[ii,paste("n.sites",hdr)] = length(unique(tmpp$Locality2))
    
    
  }
  rg = "Continental"
  ii = which(siteexc$region == rg & siteexc$species == sp2)
  tmpp = dts
  hdr = "final"
  siteexc[ii,paste("n.counts",hdr)] = nrow(tmpp)
  siteexc[ii,paste("n.nonzero.counts",hdr)] = length(which(tmpp$Count > 0))
  siteexc[ii,paste("n.sites",hdr)] = length(unique(tmpp$Locality2))
  
  
  
  
  ### idenfity which regions to keep for this species
  
  minsitsbyreg <- 10
  sitsinreg <- unique(dts[,c("Locality2","Region2")])
  nsitsbyreg <- table(sitsinreg$Region2)
  
  regs <- names(nsitsbyreg)[which(nsitsbyreg >= minsitsbyreg)]
  
  dts <- dts[which(dts$Region2 %in% regs),]
  
  for(rg in unique(dtf$Region2)){
    ii = which(siteexc$region == rg & siteexc$species == sp2)
    tmpp = dts[which(dts$Region2 == rg),]
    hdr = "Region.excl.lt.10.sites"
    if(nrow(tmpp) > 0){
    siteexc[ii,paste(hdr)] = ""
    }else{
      siteexc[ii,paste(hdr)] = "Exluded"
}
    
  }
  
  
  #### regions with a minimum number of years with the species observed
  
  
  dts2 <- dts[which(dts$Count > 0),] 
  nonzerocountsbyregyear <- table(dts[,c("Region2","Year")])
  nonzerocountsbyregyear[which(nonzerocountsbyregyear > 0)] <- 1
  
  nyearswnonzerobyreg <- rowSums(nonzerocountsbyregyear)
  
  minyears <- 15 ### minimum nubmer of years with nonzero counts in strata
  regs <- regs[which(regs %in% names(nyearswnonzerobyreg)[which(nyearswnonzerobyreg >= minyears)])]
  dts <- dts[which(dts$Region2 %in% regs),]
  
  
  for(rg in unique(dtf$Region2)){
    ii = which(siteexc$region == rg & siteexc$species == sp2)
    tmpp = dts[which(dts$Region2 == rg),]
    hdr = "Region.excl.lt.15.years"
    if(nrow(tmpp) > 0){
      siteexc[ii,paste(hdr)] = ""
    }else{
      siteexc[ii,paste(hdr)] = "Exluded"
    }
    
  }
  
  
  
  stratkeep <- regs
  # this line keeps strata where the species has been observed in more than X years of the time series
  
  

  
  
  
  
  
  
  
  
  
}#end species loop
  
  
  write.csv(siteexc,
            "site and strata exclusions.csv")
  
  
  
  
  
  
  
  