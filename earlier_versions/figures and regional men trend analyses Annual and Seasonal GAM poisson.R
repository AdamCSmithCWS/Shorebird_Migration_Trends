


### add the observed mean counts in each year to the relevant trajectory plots




setwd("c:/Shorebirds")

library(rjags)

library(RColorBrewer)
library(sp)
library(maptools)
library(rgdal)
gpclibPermit()
library(RColorBrewer)

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


stratalistfull <- data.frame(stratname = unique(sites$Region2))
plotting <- T

nstrat <- length(stratalistfull$stratname)
col.st <- brewer.pal(nstrat, "Dark2")
names(col.st) <- stratalistfull$stratname
# 
############# arbitrary species used to retrieve a list of the strata names to generate a colour scheme that
############## is constant across species
# stlist <- read.csv(paste(paste(getwd(),"/","Wilson.s Phalarope","jan/",sep = ""),"strata and site lists ","Wilson.s Phalarope",".csv",sep = ""),stringsAsFactors = F)
# sts <- unique(stlist[,c("strat","stratx")])
# names(sts) <- c("stratname","strat")
# 

# rm(list = c("sts","stlist"))
# ############# end of arbitrary code to generate colour scheme

sspwt <- read.csv("archived species results/Long-term continental shorebird trends.csv")
spwt <- unique(sspwt$species)
spwt <- gsub(spwt,pattern = ".",replacement = "",fixed = T)
spwt <- gsub(spwt,pattern = " ",replacement = "",fixed = T)


storun <- as.character(sps[which(sps$dirnames %in% spwt),"Species"])
source("c:/Functions/transparency function.r")


trajcols <- c(grey(0.2),col.st)
trajcolst = trajcols
trajcolst[1] = transp.func(trajcols[i],0.3)

for(i in 2:length(trajcols)){
  trajcolst[i] = transp.func(trajcols[i],0.1)
}
names(trajcols) <- c("Continental",names(col.st))
names(trajcolst) <- c("Continental",names(col.st))




for (sp2 in storun) {

# sp2d <- snrun[sp2]


sp2d <- sps[which(sps$Species == sp2),"dirnames"]
#sp2d <- snrun[sp2]
#dir.spsp <- paste0("output/",sp2d)
sp.dir <- paste0("output/",sp2d,"poisson annual and seasonal GAM")
# dtf <- df2[which(df2$Species == sp2),]

#col.st <- brewer.pal(8, "Dark2")
#names(col.st) <- sts$stratname
#nstrat <- length(sts$stratname)

if (file.exists(paste(sp.dir,"/cd13.RData", sep = ""))) {
  
  ##original data
  dts <- read.csv(paste(sp.dir,"/original data table ",sp2d,".csv",sep = ""), stringsAsFactors = F)
  
  sts <- unique(dts[,c("strat","stratf","stratx")])
  names(sts) <- c("stratname","stratnamefact","strat")
  
  
#   ncz <- function(x){length(which(x == 0))/length(x)}  
# pzs = tapply(dts$Count,dts[,c("juldate","stratf")],ncz) 
#  
# x11()
# par(mfrow = c(3,3))
#   for(j in 1:ncol(pzs)){
#     plot(pzs[,j],
#          type = "l",
#          main = )}
# 
 
 
  
load(paste(sp.dir,"/cd13.RData", sep = ""))

############################################


#####################
# generate full posterior for each stratum-level index, 

# calculate the annual indices for each region, scaled to reflect the mean of the site-level counts
# use these indices to calculate a continental scale index to contrast with the "size" parameter




# calculate trends for each region and continental:
# the slope trends (betas)
# the end-pooint trends (change)
# just for the continental trend calculate the following...
# end-point trends after accounting for the mean abundance across sites in a region
# by adding the site-level abundances back into the annual indices
#
#####################

midyear <- floor(max(dts$yearx)/2)

parms <- names(cd13[[1]][1,])

walpha = "alpha"
alpha <- unlist(cd13[,walpha,drop = F])



stebystrat <- unique(dts[,c("stratx","sitex")])
ylist <- unique(dts[,c("yearx","Year")])
ylist <- ylist[order(ylist$yearx),]


dmins = tapply(dts$juldate,dts$stratx,min)
dmins = dmins-1
dmaxs = tapply(dts$date,dts$stratx,max)


### stratum level indices two ways
strindices = expand.grid(stratx = c(sts$strat,"continental"),
                      yearx = min(dts$yearx):max(dts$yearx),
                      species = sp2,
                      stringsAsFactors = F)

strsmooths = expand.grid(stratx = c(sts$strat),
                         date = min(dts$date):max(dts$date),
                         species = sp2,
                         stringsAsFactors = F)

postmatsm <- array(NA,
                   dim = c(length(sts$strat),
                           length(min(dts$date):max(dts$date)),
                           length(alpha)))

postmatii <- array(NA,
                   dim = c(length(sts$strat),
                           length(min(dts$yearx):max(dts$yearx)),
                           length(alpha)))
postmatsiz = postmatii
postmatslii = postmatii
postmatslsiz = postmatii

postmatsizcont = matrix(NA,
                        ncol = length(alpha),
                        nrow = length(min(dts$yearx):max(dts$yearx)))
postmatslsizcont = postmatsizcont



####### seasonal smooths
for(st in sts$strat){
  
for(d in 1:dmaxs[st]){
 D <- d+dmins[st]
  
  wz <- paste0("zi[",d,",",st,"]")
  zi <- unlist(cd13[,wz,drop = F])
  

    ST = sts[which(sts$strat == st),"stratname"] 
    wind = which(strsmooths$stratx == st & strsmooths$date == d)  
    sits = stebystrat[which(stebystrat$stratx == st),"sitex"]

 
    strsmooths[wind,"Julian date"] <- D
    strsmooths[wind,"Stratname"] <- ST
    
    
    strsmooths[wind,"smooth"] <- signif(median(zi),3)
    strsmooths[wind,"smooth.lci"] <- signif(quantile(zi,0.025,names = F),3)
    strsmooths[wind,"smooth.uci"] <- signif(quantile(zi,0.975,names = F),3)
    strsmooths[wind,"smooth.lci.90"] <- signif(quantile(zi,0.05,names = F),3)
    strsmooths[wind,"smooth.uci.90"] <- signif(quantile(zi,0.9,names = F),3)
    strsmooths[wind,"smooth.lci.50"] <- signif(quantile(zi,0.25,names = F),3)
    strsmooths[wind,"smooth.uci.50"] <- signif(quantile(zi,0.75,names = F),3)
    
    
    
    postmatsm[st,d,] <- zi
      rm(list = c("zi","wind"))
    
    
  }#d
  
}#st


strsmooths <- strsmooths[which(!is.na(strsmooths$smooth)),]
strsmooths <- strsmooths[order(strsmooths$stratx,strsmooths$date),]
wzimax = 0
zimax = ""
for(st in sts$strat){
  tmp = strsmooths[which(strsmooths$strat == st),]
  ww = which.max(tmp$smooth)
  wwd = tmp[ww,"date"]
  zimax[st] <- paste0("zi[",wwd,",",st,"]")
  wzimax[st] = wwd
}




################### plotting the smooths by region

pdf(height = 12,
    width = 8,
    file = paste0(sp.dir,"/seasonal smooths",sp2d,".pdf"))

par(mfcol = c(4,2))
ee = "smooth"
medc = ee
lcic = paste0(ee,".lci.50")
ucic = paste0(ee,".uci.50")
lcic9 = paste0(ee,".lci")
ucic9 = paste0(ee,".uci")
lcic6 = paste0(ee,".lci.90")
ucic6 = paste0(ee,".uci.90")
een = "Seasonal Smooth"


for(st in c(sts$stratname)){
  tmp = strsmooths[which(strsmooths$Stratname == st),]
  stt = sts[which(sts$stratname == st),"strat"]
  tmpobsmeancnts <- dts[which(dts$stratx == stt),]
  obsmeancnts = log(tapply(tmpobsmeancnts$Count,tmpobsmeancnts$juldate,mean,na.rm = T)+1)
  peakdays = c(((wzimax[stt]+dmins[stt])-7):((wzimax[stt]+dmins[stt])+7))
  obsmeancnts = obsmeancnts/quantile(obsmeancnts,0.95)#[which(names(obsmeancnts) %in% 
                                             #        peakdays)])
  obsmeancntstr = min(tmp[,lcic9])+(obsmeancnts*(diff(range(tmp$smooth))))
  yobsmn = as.integer(names(obsmeancntstr))
  tsmoo = loess(obsmeancntstr~yobsmn)
  sqrtncnts = sqrt(tapply(tmpobsmeancnts$Count,tmpobsmeancnts$juldate,length))
  sqrtncnts = sqrtncnts/mean(sqrtncnts)
  
yup = max(tmp[,ucic9])
ydown = min(tmp[,lcic9])
xlims = c(min(strsmooths[,"Julian date"]),max(strsmooths[,"Julian date"]))
xplot = tmp[,"Julian date"]
if(st == sts$stratname[1]){
  plot(1000,
     col = "white",
     xlim = xlims,
     ylim = c(ydown,yup),
     xlab = "",
     ylab = "",
     xaxs = "i",
     yaxs = "i",
     bty = "l",
     main = paste(st,"seasonal smooth on count")) 
}else{
  plot(1000,
       col = "white",
       xlim = xlims,
       ylim = c(ydown,yup),
       xlab = "",
       ylab = "",
       xaxs = "i",
       yaxs = "i",
       bty = "l",
       main = paste(st)) 
  
  
}
###### instead of the daily means, calculate weekly means, or plot a smoother...
  points(x = as.integer(names(obsmeancntstr)),
         y = obsmeancntstr,
         cex = sqrtncnts,
         pch = 19,
         col = grey(0.8))
  lines(y = predict(tsmoo),
                 x = as.integer(names(obsmeancntstr)),
                 cex = sqrtncnts,
         pch = 19,
         col = grey(0.8),
        lwd = 1.5)
  polygon(y = c(tmp[,lcic],rev(tmp[,ucic])),
          x = c(xplot,rev(xplot)),
          border = NA,
          col = transp.func(trajcols[st],0.1))
  polygon(y = c(tmp[,lcic9],rev(tmp[,ucic9])),
          x = c(xplot,rev(xplot)),
          border = NA,
          col = transp.func(trajcols[st],0.1))
  polygon(y = c(tmp[,lcic6],rev(tmp[,ucic6])),
          x = c(xplot,rev(xplot)),
          border = NA,
          col = transp.func(trajcols[st],0.1))
  lines(y = tmp[,medc],
        x = xplot,
        col = trajcols[st],
        lwd = 2)

  
}#st




dev.off()



























#### compiling the mean siteeffects

sits = stebystrat[,"sitex"]

siteffs = vector(length = length(sits))
names(siteffs) = as.character(sits)
for(si in sits){
  wsit = paste0("ste[",si,"]")
  sit = unlist(cd13[,wsit,drop = F])
  siteffs[as.character(si)] = mean(sit)
}#si







### indices
for(y in ylist$yearx){
  Y <- ylist[which(ylist$yearx == y),"Year"]
  

  for(st in sts$strat){

   ST = sts[which(sts$strat == st),"stratname"]
   
   wann <- paste0("ann[",y,",",st,"]")
   ann <- unlist(cd13[,wann,drop = F])
   
    wind = which(strindices$stratx == st & strindices$yearx == y) 
   
   sits = stebystrat[which(stebystrat$stratx == st),"sitex"]
    #wbetas = paste0("beta[",st,"]")
    #b = unlist(cd13[,wbetas,drop = F])
    zi = unlist(cd13[,zimax[st],drop = F])
    
    ii <- matrix(NA,
                 nrow = length(ann),
                 ncol = length(sits))
   
    
    j = 1
    for(si in sits){
      wsit = paste0("ste[",si,"]")
      sit = unlist(cd13[,wsit,drop = F])
      ii[,j] <- exp(alpha + ann + sit + zi)
      
        j = j+1
    }#si
    
    iiv <- rowMeans(ii)
    
    rm(list = c("ii"))
    siz = exp(alpha + ann + zi)
    
    strindices[wind,"Year"] <- Y
    strindices[wind,"Stratname"] <- ST
    
    
    strindices[wind,"size"] <- signif(median(siz),3)
    strindices[wind,"size.lci"] <- signif(quantile(siz,0.025,names = F),3)
    strindices[wind,"size.uci"] <- signif(quantile(siz,0.975,names = F),3)
    strindices[wind,"size.lci.90"] <- signif(quantile(siz,0.05,names = F),3)
    strindices[wind,"size.uci.90"] <- signif(quantile(siz,0.9,names = F),3)
    strindices[wind,"size.lci.50"] <- signif(quantile(siz,0.25,names = F),3)
    strindices[wind,"size.uci.50"] <- signif(quantile(siz,0.75,names = F),3)
 
    
    # strindices[wind,"slopeind"] <- signif(median(sislp),3)
    # strindices[wind,"slopeind.lci"] <- signif(quantile(sislp,0.025,names = F),3)
    # strindices[wind,"slopeind.uci"] <- signif(quantile(sislp,0.975,names = F),3)
    # strindices[wind,"slopeind.lci.90"] <- signif(quantile(sislp,0.05,names = F),3)
    # strindices[wind,"slopeind.uci.90"] <- signif(quantile(sislp,0.9,names = F),3)
    # strindices[wind,"slopeind.lci.50"] <- signif(quantile(sislp,0.25,names = F),3)
    # strindices[wind,"slopeind.uci.50"] <- signif(quantile(sislp,0.75,names = F),3)
    # 
    #   strindices[wind,"slopeind2"] <- signif(median(iislv),3)
    # strindices[wind,"slopeind2.lci"] <- signif(quantile(iislv,0.025,names = F),3)
    # strindices[wind,"slopeind2.uci"] <- signif(quantile(iislv,0.975,names = F),3)
    # strindices[wind,"slopeind2.lci.90"] <- signif(quantile(iislv,0.05,names = F),3)
    # strindices[wind,"slopeind2.uci.90"] <- signif(quantile(iislv,0.9,names = F),3)
    # strindices[wind,"slopeind2.lci.50"] <- signif(quantile(iislv,0.25,names = F),3)
    # strindices[wind,"slopeind2.uci.50"] <- signif(quantile(iislv,0.75,names = F),3)
    # 
    
    strindices[wind,"index"] <- signif(median(iiv),3)
    strindices[wind,"index.lci"] <- signif(quantile(iiv,0.025,names = F),3)
    strindices[wind,"index.uci"] <- signif(quantile(iiv,0.975,names = F),3)
    strindices[wind,"index.lci.90"] <- signif(quantile(iiv,0.05,names = F),3)
    strindices[wind,"index.uci.90"] <- signif(quantile(iiv,0.9,names = F),3)
    strindices[wind,"index.lci.50"] <- signif(quantile(iiv,0.25,names = F),3)
    strindices[wind,"index.uci.50"] <- signif(quantile(iiv,0.75,names = F),3)

    wdts = which(dts$Year == Y & dts$stratx == st)
    sitsthisyear = as.character(dts[wdts,"sitex"])
    strindices[wind,"meansiteeff"] <- mean(siteffs[sitsthisyear]) 
    strindices[wind,"obsmeancount"] <- mean(dts[wdts,"Count"]) 
      strindices[wind,"obsmediancount"] <- median(dts[wdts,"Count"])
      strindices[wind,"pzerocounts"] <- length(which(dts[wdts,"Count"] == 0))/length(wdts) 
      strindices[wind,"ncounts"] <- length(wdts) 
      strindices[wind,"nsites"] <- length(unique(sitsthisyear)) 
      
      
    postmatii[st,y,] <- iiv
    postmatsiz[st,y,] <- siz
    # postmatslii[st,y,] <- iislv
    # postmatslsiz[st,y,] <- sislp
    rm(list = c("iiv","siz","iislv","sislp","wind"))

    
  }#st
  
  
  # wbetas = paste0("betaH")
  # bh = unlist(cd13[,wbetas,drop = F])
  # 
  # ctslsiz <- exp(alpha + bh*(y-midyear))
  # ctsz <- exp(alpha + bh*(y-midyear) + ann)
  # 
  # postmatsizcont[y,] <- ctsz
  # postmatslsizcont[y,] <- ctslsiz
  # 
  
  ctii <- colMeans(postmatii[,y,])
  #ctslii <- colMeans(postmatslii[,y,])
  
  wind = which(strindices$stratx == "continental" & strindices$yearx == y)
  
  strindices[wind,"Year"] <- Y
  strindices[wind,"Stratname"] <- "Continental"
  
  
  # strindices[wind,"size"] <- signif(median(ctsz),3)
  # strindices[wind,"size.lci"] <- signif(quantile(ctsz,0.025,names = F),3)
  # strindices[wind,"size.uci"] <- signif(quantile(ctsz,0.975,names = F),3)
  # strindices[wind,"size.lci.90"] <- signif(quantile(ctsz,0.05,names = F),3)
  # strindices[wind,"size.uci.90"] <- signif(quantile(ctsz,0.9,names = F),3)
  # strindices[wind,"size.lci.50"] <- signif(quantile(ctsz,0.25,names = F),3)
  # strindices[wind,"size.uci.50"] <- signif(quantile(ctsz,0.75,names = F),3)
  # 
  
  # strindices[wind,"slopeind"] <- signif(median(ctslsiz),3)
  # strindices[wind,"slopeind.lci"] <- signif(quantile(ctslsiz,0.025,names = F),3)
  # strindices[wind,"slopeind.uci"] <- signif(quantile(ctslsiz,0.975,names = F),3)
  # strindices[wind,"slopeind.lci.90"] <- signif(quantile(ctslsiz,0.05,names = F),3)
  # strindices[wind,"slopeind.uci.90"] <- signif(quantile(ctslsiz,0.9,names = F),3)
  # strindices[wind,"slopeind.lci.50"] <- signif(quantile(ctslsiz,0.25,names = F),3)
  # strindices[wind,"slopeind.uci.50"] <- signif(quantile(ctslsiz,0.75,names = F),3)
  # 
  # strindices[wind,"slopeind2"] <- signif(median(ctslii),3)
  # strindices[wind,"slopeind2.lci"] <- signif(quantile(ctslii,0.025,names = F),3)
  # strindices[wind,"slopeind2.uci"] <- signif(quantile(ctslii,0.975,names = F),3)
  # strindices[wind,"slopeind2.lci.90"] <- signif(quantile(ctslii,0.05,names = F),3)
  # strindices[wind,"slopeind2.uci.90"] <- signif(quantile(ctslii,0.9,names = F),3)
  # strindices[wind,"slopeind2.lci.50"] <- signif(quantile(ctslii,0.25,names = F),3)
  # strindices[wind,"slopeind2.uci.50"] <- signif(quantile(ctslii,0.75,names = F),3)
  # 
  
  strindices[wind,"index"] <- signif(median(ctii),3)
  strindices[wind,"index.lci"] <- signif(quantile(ctii,0.025,names = F),3)
  strindices[wind,"index.uci"] <- signif(quantile(ctii,0.975,names = F),3)
  strindices[wind,"index.lci.90"] <- signif(quantile(ctii,0.05,names = F),3)
  strindices[wind,"index.uci.90"] <- signif(quantile(ctii,0.9,names = F),3)
  strindices[wind,"index.lci.50"] <- signif(quantile(ctii,0.25,names = F),3)
  strindices[wind,"index.uci.50"] <- signif(quantile(ctii,0.75,names = F),3)
  
  
  wdts = which(dts$Year == Y)
  sitsthisyear = as.character(dts[wdts,"sitex"])
  strindices[wind,"meansiteeff"] <- mean(siteffs[sitsthisyear]) 
  strindices[wind,"obsmeancount"] <- mean(dts[wdts,"Count"]) 
  strindices[wind,"obsmediancount"] <- median(dts[wdts,"Count"])
  strindices[wind,"pzerocounts"] <- length(which(dts[wdts,"Count"] == 0))/length(wdts) 
  strindices[wind,"ncounts"] <- length(wdts) 
  strindices[wind,"nsites"] <- length(unique(sitsthisyear)) 
  
  
  
  rm(list = c("wind"))
  
}#y

strindices = strindices[order(strindices$stratx,strindices$Year),]


### trend calculations:
syrs = c(1974,1974,1995,2006)
eyrs = c(2016,1995,2016,2016)

timspns = paste(syrs,eyrs,sep = "-")

trendsp <- expand.grid(strat = as.character(c(sts$strat,"continental")),
                       species = sp2,
                       timespan = timspns,
                       stringsAsFactors = F)
for(eyi in 1:4){
  ey1 = eyrs[eyi]
  sy1 = syrs[eyi]
  
  ey = ylist[which(ylist$Year == ey1),"yearx"]
  sy = ylist[which(ylist$Year == sy1),"yearx"]
  
  eyv <- matrix(NA,
                ncol = length(sts$strat),
                nrow = length(alpha))
  syv = eyv
  spn = timspns[eyi]
  
  yint = sy:ey
  
  
  for(st in sts$strat){
    ST = sts[which(sts$strat == st),"stratname"]
    wst = which(trendsp$strat == st &
                  trendsp$timespan == spn)
    
    trendsp[wst,"stratname"] <- ST
    
    trendsp[wst,"startyear"] <- sy1
    trendsp[wst,"endyear"] <- ey1
    
    
    sz1 <- postmatsiz[st,sy,]
    sz2 <- postmatsiz[st,ey,]
    csz <- (((sz2/sz1)^(1/(ey-sy)))-1)*100
    
    trendsp[wst,"trend.size"] <- signif(mean(csz),3)
    trendsp[wst,"trend.size.lci"] <- signif(quantile(csz,0.025,names = F),3)
    trendsp[wst,"trend.size.uci"] <- signif(quantile(csz,0.975,names = F),3)
    trendsp[wst,"trend.size.lci.90"] <- signif(quantile(csz,0.05,names = F),3)
    trendsp[wst,"trend.size.uci.90"] <- signif(quantile(csz,0.9,names = F),3)
    trendsp[wst,"trend.size.lci.50"] <- signif(quantile(csz,0.25,names = F),3)
    trendsp[wst,"trend.size.uci.50"] <- signif(quantile(csz,0.75,names = F),3)
    rm("csz")
  
    ###slope trend ignoring relative abundance
    slps = vector(length = length(alpha))
    for(yi in 1:length(alpha)){
    sz1 <- log(postmatsiz[st,yint,yi])
    m = lm(sz1~yint)
    slps[yi] = (exp(coef(m)["yint"])-1)*100
    }#yi
    trendsp[wst,"trend.size.slope"] <- signif(mean(slps),3)
    trendsp[wst,"trend.size.slope.lci"] <- signif(quantile(slps,0.025,names = F),3)
    trendsp[wst,"trend.size.slope.uci"] <- signif(quantile(slps,0.975,names = F),3)
    trendsp[wst,"trend.size.slope.lci.90"] <- signif(quantile(slps,0.05,names = F),3)
    trendsp[wst,"trend.size.slope.uci.90"] <- signif(quantile(slps,0.9,names = F),3)
    trendsp[wst,"trend.size.slope.lci.50"] <- signif(quantile(slps,0.25,names = F),3)
    trendsp[wst,"trend.size.slope.uci.50"] <- signif(quantile(slps,0.75,names = F),3)
    rm("slps")


      
    # sz1 <- postmatslsiz[st,sy,]
    # sz2 <- postmatslsiz[st,ey,]
    # csz <- (((sz2/sz1)^(1/(ey-sy)))-1)*100
    # 
    # trendsp[wst,"trend.beta"] <- signif(mean(csz),3)
    # trendsp[wst,"trend.beta.lci"] <- signif(quantile(csz,0.025,names = F),3)
    # trendsp[wst,"trend.beta.uci"] <- signif(quantile(csz,0.975,names = F),3)
    # trendsp[wst,"trend.beta.lci.90"] <- signif(quantile(csz,0.05,names = F),3)
    # trendsp[wst,"trend.beta.uci.90"] <- signif(quantile(csz,0.9,names = F),3)
    # trendsp[wst,"trend.beta.lci.50"] <- signif(quantile(csz,0.25,names = F),3)
    # trendsp[wst,"trend.beta.uci.50"] <- signif(quantile(csz,0.75,names = F),3)
    # rm("csz")
    
   
    # sz1 <- postmatslii[st,sy,]
    # sz2 <- postmatslii[st,ey,]
    # csz <- (((sz2/sz1)^(1/(ey-sy)))-1)*100
    # 
    # trendsp[wst,"trend.beta2"] <- signif(mean(csz),3)
    # trendsp[wst,"trend.beta2.lci"] <- signif(quantile(csz,0.025,names = F),3)
    # trendsp[wst,"trend.beta2.uci"] <- signif(quantile(csz,0.975,names = F),3)
    # trendsp[wst,"trend.beta2.lci.90"] <- signif(quantile(csz,0.05,names = F),3)
    # trendsp[wst,"trend.beta2.uci.90"] <- signif(quantile(csz,0.9,names = F),3)
    # trendsp[wst,"trend.beta2.lci.50"] <- signif(quantile(csz,0.25,names = F),3)
    # trendsp[wst,"trend.beta2.uci.50"] <- signif(quantile(csz,0.75,names = F),3)
    # rm("csz")
    
    
    
    
    sz1 <- postmatii[st,sy,]
    sz2 <- postmatii[st,ey,]
    csz <- (((sz2/sz1)^(1/(ey-sy)))-1)*100
    
    trendsp[wst,"trend.index"] <- signif(mean(csz),3)
    trendsp[wst,"trend.index.lci"] <- signif(quantile(csz,0.025,names = F),3)
    trendsp[wst,"trend.index.uci"] <- signif(quantile(csz,0.975,names = F),3)
    trendsp[wst,"trend.index.lci.90"] <- signif(quantile(csz,0.05,names = F),3)
    trendsp[wst,"trend.index.uci.90"] <- signif(quantile(csz,0.9,names = F),3)
    trendsp[wst,"trend.index.lci.50"] <- signif(quantile(csz,0.25,names = F),3)
    trendsp[wst,"trend.index.uci.50"] <- signif(quantile(csz,0.75,names = F),3)
    rm("csz")
    
    
    
    ###slope trend accounting for relative abundance
    slps = vector(length = length(alpha))
    for(yi in 1:length(alpha)){
      sz1 <- log(postmatii[st,yint,yi])
      m = lm(sz1~yint)
      slps[yi] = (exp(coef(m)["yint"])-1)*100
    }#yi
    trendsp[wst,"trend.index.slope"] <- signif(mean(slps),3)
    trendsp[wst,"trend.index.slope.lci"] <- signif(quantile(slps,0.025,names = F),3)
    trendsp[wst,"trend.index.slope.uci"] <- signif(quantile(slps,0.975,names = F),3)
    trendsp[wst,"trend.index.slope.lci.90"] <- signif(quantile(slps,0.05,names = F),3)
    trendsp[wst,"trend.index.slope.uci.90"] <- signif(quantile(slps,0.9,names = F),3)
    trendsp[wst,"trend.index.slope.lci.50"] <- signif(quantile(slps,0.25,names = F),3)
    trendsp[wst,"trend.index.slope.uci.50"] <- signif(quantile(slps,0.75,names = F),3)
    rm("slps")
    
    
    
    
  }
  
  #trends with continental summed indices accounting for mean counts by regions
  wst = which(trendsp$strat == "continental" &
                trendsp$timespan == spn)
  

  trendsp[wst,"startyear"] <- sy1
  trendsp[wst,"endyear"] <- ey1
  trendsp[wst,"stratname"] <- "Continental"
  
  
  # sz1 <- postmatsizcont[sy,]
  # sz2 <- postmatsizcont[ey,]
  # csz <- (((sz2/sz1)^(1/(ey-sy)))-1)*100
  # 
  # trendsp[wst,"trend.size"] <- signif(mean(csz),3)
  # trendsp[wst,"trend.size.lci"] <- signif(quantile(csz,0.025,names = F),3)
  # trendsp[wst,"trend.size.uci"] <- signif(quantile(csz,0.975,names = F),3)
  # trendsp[wst,"trend.size.lci.90"] <- signif(quantile(csz,0.05,names = F),3)
  # trendsp[wst,"trend.size.uci.90"] <- signif(quantile(csz,0.9,names = F),3)
  # trendsp[wst,"trend.size.lci.50"] <- signif(quantile(csz,0.25,names = F),3)
  # trendsp[wst,"trend.size.uci.50"] <- signif(quantile(csz,0.75,names = F),3)
  # rm("csz")
  # 
  # 
  # ###slope trend ignoring relative abundance
  # slps = vector(length = length(alpha))
  # for(yi in 1:length(alpha)){
  #   sz1 <- log(postmatsizcont[yint,yi])
  #   m = lm(sz1~yint)
  #   slps[yi] = (exp(coef(m)["yint"])-1)*100
  # }#yi
  # trendsp[wst,"trend.size.slope"] <- signif(mean(slps),3)
  # trendsp[wst,"trend.size.slope.lci"] <- signif(quantile(slps,0.025,names = F),3)
  # trendsp[wst,"trend.size.slope.uci"] <- signif(quantile(slps,0.975,names = F),3)
  # trendsp[wst,"trend.size.slope.lci.90"] <- signif(quantile(slps,0.05,names = F),3)
  # trendsp[wst,"trend.size.slope.uci.90"] <- signif(quantile(slps,0.9,names = F),3)
  # trendsp[wst,"trend.size.slope.lci.50"] <- signif(quantile(slps,0.25,names = F),3)
  # trendsp[wst,"trend.size.slope.uci.50"] <- signif(quantile(slps,0.75,names = F),3)
  # rm("slps")
  # 
  # 
  # 
  # sz1 <- postmatslsizcont[sy,]
  # sz2 <- postmatslsizcont[ey,]
  # csz <- (((sz2/sz1)^(1/(ey-sy)))-1)*100
  # 
  # trendsp[wst,"trend.beta"] <- signif(mean(csz),3)
  # trendsp[wst,"trend.beta.lci"] <- signif(quantile(csz,0.025,names = F),3)
  # trendsp[wst,"trend.beta.uci"] <- signif(quantile(csz,0.975,names = F),3)
  # trendsp[wst,"trend.beta.lci.90"] <- signif(quantile(csz,0.05,names = F),3)
  # trendsp[wst,"trend.beta.uci.90"] <- signif(quantile(csz,0.9,names = F),3)
  # trendsp[wst,"trend.beta.lci.50"] <- signif(quantile(csz,0.25,names = F),3)
  # trendsp[wst,"trend.beta.uci.50"] <- signif(quantile(csz,0.75,names = F),3)
  # rm("csz")
  # 
  
  # sz1 <- colSums(postmatslii[,sy,])
  # sz2 <- colSums(postmatslii[,ey,])
  # csz <- (((sz2/sz1)^(1/(ey-sy)))-1)*100
  # 
  # trendsp[wst,"trend.beta2"] <- signif(mean(csz),3)
  # trendsp[wst,"trend.beta2.lci"] <- signif(quantile(csz,0.025,names = F),3)
  # trendsp[wst,"trend.beta2.uci"] <- signif(quantile(csz,0.975,names = F),3)
  # trendsp[wst,"trend.beta2.lci.90"] <- signif(quantile(csz,0.05,names = F),3)
  # trendsp[wst,"trend.beta2.uci.90"] <- signif(quantile(csz,0.9,names = F),3)
  # trendsp[wst,"trend.beta2.lci.50"] <- signif(quantile(csz,0.25,names = F),3)
  # trendsp[wst,"trend.beta2.uci.50"] <- signif(quantile(csz,0.75,names = F),3)
  # rm("csz")
  
  
  
  sz1 <- colSums(postmatii[,sy,])
  sz2 <- colSums(postmatii[,ey,])
  csz <- (((sz2/sz1)^(1/(ey-sy)))-1)*100
  
  trendsp[wst,"trend.index"] <- signif(mean(csz),3)
  trendsp[wst,"trend.index.lci"] <- signif(quantile(csz,0.025,names = F),3)
  trendsp[wst,"trend.index.uci"] <- signif(quantile(csz,0.975,names = F),3)
  trendsp[wst,"trend.index.lci.90"] <- signif(quantile(csz,0.05,names = F),3)
  trendsp[wst,"trend.index.uci.90"] <- signif(quantile(csz,0.9,names = F),3)
  trendsp[wst,"trend.index.lci.50"] <- signif(quantile(csz,0.25,names = F),3)
  trendsp[wst,"trend.index.uci.50"] <- signif(quantile(csz,0.75,names = F),3)
  
  
  ###slope trend accounting for relative abundance
  slps = vector(length = length(alpha))
  for(yi in 1:length(alpha)){
    sz1 <- log(colSums(postmatii[,yint,yi]))
    m = lm(sz1~yint)
    slps[yi] = (exp(coef(m)["yint"])-1)*100
  }#yi
  trendsp[wst,"trend.index.slope"] <- signif(mean(slps),3)
  trendsp[wst,"trend.index.slope.lci"] <- signif(quantile(slps,0.025,names = F),3)
  trendsp[wst,"trend.index.slope.uci"] <- signif(quantile(slps,0.975,names = F),3)
  trendsp[wst,"trend.index.slope.lci.90"] <- signif(quantile(slps,0.05,names = F),3)
  trendsp[wst,"trend.index.slope.uci.90"] <- signif(quantile(slps,0.9,names = F),3)
  trendsp[wst,"trend.index.slope.lci.50"] <- signif(quantile(slps,0.25,names = F),3)
  trendsp[wst,"trend.index.slope.uci.50"] <- signif(quantile(slps,0.75,names = F),3)
  rm("slps")
  
  
  
}#eyi



# ### trend comparisons
# for(tsp in timspns){
# pdf(file = paste0(sp.dir,"/",sp2d,"trend comparisons ",tsp,".pdf"))
# par(mfrow = c(3,3))
# for(st in c("Continental",sts$stratname)){
# plot(999,
#      xlim = c(0.5,6.5),
#      ylim = c(-10,10),
#      main = st,
#      ylab = "trend",
#      xlab = "",
#      xaxt = "n",
#      bty = "l")
# wst = which(trendsp$stratname == st & trendsp$timespan == tsp)
# 
# abline(h = 0,lty = 2,col = grey(0.8))
# j = 1
# for(ttp in c("trend.size","trend.size.slope","trend.beta","trend.beta2","trend.index.slope","trend.index")){
# 
#   
#   arrows(y0 = trendsp[wst,paste0(ttp,".lci")],
#        y1 = trendsp[wst,paste0(ttp,".uci")],
#        x0 = j,
#        x1 = j,
#        length = 0,
#        col = transp.func(trajcols[st],0.3),
#        lwd = 2)
#   arrows(y0 = trendsp[wst,paste0(ttp,".lci.90")],
#          y1 = trendsp[wst,paste0(ttp,".uci.90")],
#          x0 = j,
#          x1 = j,
#          length = 0,
#          col = transp.func(trajcols[st],0.4),
#          lwd = 2)
#   
#   arrows(y0 = trendsp[wst,paste0(ttp,".lci.50")],
#          y1 = trendsp[wst,paste0(ttp,".uci.50")],
#          x0 = j,
#          x1 = j,
#          length = 0,
#          col = transp.func(trajcols[st],0.5),
#          lwd = 2)
#   
#   points(y = trendsp[wst,ttp],
#          x = j,
#          col = trajcols[st])
#   
#   j = j+1
# }
# 
# axis(side = 1,
#      at = c(1,2,3,4,5,6),
#      labels = c("size","lmsize","slope","slope2","lmindex","index"),
#      las = 3)
# 
# 
# 
# }
# dev.off()
# 
# }#tsp
# 



# ##### overplotting the indices
# 
# estimates = c("size", "index")
# names(estimates) <- c("size","index")
# 
# 
# pdf(paste0(sp.dir,"/index overplots",sp2d,".pdf"))
# 
# for(ee in estimates){
#   medc = ee
#   lcic = paste0(ee,".lci.50")
#   ucic = paste0(ee,".uci.50")
#   lcic9 = paste0(ee,".lci")
#   ucic9 = paste0(ee,".uci")
#   lcic6 = paste0(ee,".lci.90")
#   ucic6 = paste0(ee,".uci.90")
#  een = names(ee) 
# 
#  yup = min(c(max(strindices[,medc])*2,max(strindices[,ucic])))
#  xplot = 1974:2016
# if(ee %in% estimates[3:4]){
#   plot(1000,
#       col = "white",
#       xlim = c(1974,2025),
#       ylim = c(0.0001,yup),
#       xlab = "",
#       ylab = "",
#       xaxs = "i",
#       yaxs = "i",
#       bty = "l",
#       main = ee,
#       log = "y") 
#   for(st in c(sts$stratname,"Continental")){
#     tmp = strindices[which(strindices$Stratname == st),]
#     # polygon(y = c(tmp[,lcic],rev(tmp[,ucic])),
#     #         x = c(xplot,rev(xplot)),
#     #         border = NA,
#     #         col = trajcolst[st])
#     # polygon(y = c(tmp[,lcic9],rev(tmp[,ucic9])),
#     #         x = c(xplot,rev(xplot)),
#     #         border = NA,
#     #         col = transp.func(trajcols[st],0.1))
#     polygon(y = c(tmp[,lcic6],rev(tmp[,ucic6])),
#             x = c(xplot,rev(xplot)),
#             border = NA,
#             col = transp.func(trajcols[st],0.1))
#     lines(y = tmp[,medc],
#           x = xplot,
#           col = trajcols[st],
#           lwd = 2)
#     text(st,
#          x = 2016,
#          y = tmp[which(tmp$Year == 2016),medc],
#          pos = 4,
#          col = trajcols[st])
#     
#     
#   }#st
#   
#   
# }
#  
#     plot(1000,
#          col = "white",
#          xlim = c(1974,2025),
#          ylim = c(0,yup),
#          xlab = "",
#          ylab = "",
#          xaxs = "i",
#          yaxs = "i",
#          bty = "l",
#          main = ee) 
# 
#  for(st in c(sts$stratname,"Continental")){
#    tmp = strindices[which(strindices$Stratname == st),]
#    # polygon(y = c(tmp[,lcic],rev(tmp[,ucic])),
#    #         x = c(xplot,rev(xplot)),
#    #         border = NA,
#    #         col = trajcolst[st])
#    # polygon(y = c(tmp[,lcic9],rev(tmp[,ucic9])),
#    #         x = c(xplot,rev(xplot)),
#    #         border = NA,
#    #         col = transp.func(trajcols[st],0.1))
#    polygon(y = c(tmp[,lcic6],rev(tmp[,ucic6])),
#            x = c(xplot,rev(xplot)),
#            border = NA,
#            col = transp.func(trajcols[st],0.1))
#    lines(y = tmp[,medc],
#          x = xplot,
#          col = trajcols[st],
#          lwd = 2)
#    text(st,
#         x = 2016,
#         y = tmp[which(tmp$Year == 2016),medc],
#         pos = 4,
#         col = trajcols[st])
#    
#    
#  }#st
#  
#  
#    
#   
# }#ee
# 
# 
# dev.off()
# 




pdf(file = paste0(sp.dir,"/overall trajectory",sp2d,".pdf"))


ee = estimates[2]
medc = ee
lcic = paste0(ee,".lci.50")
ucic = paste0(ee,".uci.50")
lcic9 = paste0(ee,".lci")
ucic9 = paste0(ee,".uci")
lcic6 = paste0(ee,".lci.90")
ucic6 = paste0(ee,".uci.90")
een = names(ee) 

yup = min(c(max(strindices[which(strindices$Stratname == "Continental"),medc])*2,
            max(strindices[which(strindices$Stratname == "Continental"),ucic6])))
xplot = 1974:2016

  plot(0.001,
       col = "white",
       xlim = c(1974,2025),
       ylim = c(0,yup),
       xlab = "",
       ylab = "",
       xaxs = "i",
       yaxs = "i",
       bty = "l",
       main = ee) 
  

for(st in c(sts$stratname)){
  tmp = strindices[which(strindices$Stratname == st),]
  # polygon(y = c(tmp[,lcic],rev(tmp[,ucic])),
  #         x = c(xplot,rev(xplot)),
  #         border = NA,
  #         col = transp.func(trajcols[st],0.1))
  # polygon(y = c(tmp[,lcic9],rev(tmp[,ucic9])),
  #         x = c(xplot,rev(xplot)),
  #         border = NA,
  #         col = transp.func(trajcols[st],0.1))
  # polygon(y = c(tmp[,lcic6],rev(tmp[,ucic6])),
  #         x = c(xplot,rev(xplot)),
  #         border = NA,
  #         col = transp.func(trajcols[st],0.1))
  lines(y = tmp[,medc],
        x = xplot,
        col = trajcols[st],
        lwd = 2)
  text(st,
       x = 2016,
       y = tmp[which(tmp$Year == 2016),medc],
       pos = 4,
       col = transp.func(trajcols[st],0.3))
  
  
}#st

  ee = estimates[1]
  
st = "Continental"
medc = ee
lcic = paste0(ee,".lci.50")
ucic = paste0(ee,".uci.50")
lcic9 = paste0(ee,".lci")
ucic9 = paste0(ee,".uci")
lcic6 = paste0(ee,".lci.90")
ucic6 = paste0(ee,".uci.90")
een = names(ee) 

tmp = strindices[which(strindices$Stratname == st),]
polygon(y = c(tmp[,lcic9],rev(tmp[,ucic9])),
        x = c(xplot,rev(xplot)),
        border = NA,
        col = trajcolst[st])
# polygon(y = c(tmp[,lcic9],rev(tmp[,ucic9])),
#         x = c(xplot,rev(xplot)),
#         border = NA,
#         col = transp.func(trajcols[st],0.1))
# polygon(y = c(tmp[,lcic6],rev(tmp[,ucic6])),
#         x = c(xplot,rev(xplot)),
#         border = NA,
#         col = transp.func(trajcols[st],0.1))
lines(y = tmp[,medc],
      x = xplot,
      col = trajcols[st],
      lwd = 2)
text(st,
     x = 2016,
     y = tmp[which(tmp$Year == 2016),medc],
     pos = 4,
     col = trajcols[st])




dev.off()



# 
# 
# 
# 
# pdf(file = paste0(sp.dir,"/alternate overall trajectory",sp2d,".pdf"))
# 
# 
# ee = estimates[3]
# medc = ee
# lcic = paste0(ee,".lci.50")
# ucic = paste0(ee,".uci.50")
# lcic9 = paste0(ee,".lci")
# ucic9 = paste0(ee,".uci")
# lcic6 = paste0(ee,".lci.90")
# ucic6 = paste0(ee,".uci.90")
# een = names(ee) 
# 
# yup = min(c(max(strindices[which(strindices$Stratname == "Continental"),medc])*2,
#             max(strindices[which(strindices$Stratname == "Continental"),ucic6])))
# xplot = 1974:2016
# 
# plot(0.001,
#      col = "white",
#      xlim = c(1974,2025),
#      ylim = c(0,yup),
#      xlab = "",
#      ylab = "",
#      xaxs = "i",
#      yaxs = "i",
#      bty = "l",
#      main = paste(sp2d,"alternate continental trajectory")) 
# 
# 
# for(st in c(sts$stratname)){
#   tmp = strindices[which(strindices$Stratname == st),]
#   # polygon(y = c(tmp[,lcic],rev(tmp[,ucic])),
#   #         x = c(xplot,rev(xplot)),
#   #         border = NA,
#   #         col = transp.func(trajcols[st],0.1))
#   # polygon(y = c(tmp[,lcic9],rev(tmp[,ucic9])),
#   #         x = c(xplot,rev(xplot)),
#   #         border = NA,
#   #         col = transp.func(trajcols[st],0.1))
#   # polygon(y = c(tmp[,lcic6],rev(tmp[,ucic6])),
#   #         x = c(xplot,rev(xplot)),
#   #         border = NA,
#   #         col = transp.func(trajcols[st],0.1))
#   lines(y = tmp[,medc],
#         x = xplot,
#         col = trajcols[st],
#         lwd = 2)
#   text(st,
#        x = 2016,
#        y = tmp[which(tmp$Year == 2016),medc],
#        pos = 4,
#        col = transp.func(trajcols[st],0.3))
#   
#   
# }#st
# 
# ee = estimates[4]
# 
# st = "Continental"
# medc = ee
# lcic = paste0(ee,".lci.50")
# ucic = paste0(ee,".uci.50")
# lcic9 = paste0(ee,".lci")
# ucic9 = paste0(ee,".uci")
# lcic6 = paste0(ee,".lci.90")
# ucic6 = paste0(ee,".uci.90")
# een = names(ee) 
# 
# tmp = strindices[which(strindices$Stratname == st),]
# polygon(y = c(tmp[,lcic9],rev(tmp[,ucic9])),
#         x = c(xplot,rev(xplot)),
#         border = NA,
#         col = trajcolst[st])
# # polygon(y = c(tmp[,lcic9],rev(tmp[,ucic9])),
# #         x = c(xplot,rev(xplot)),
# #         border = NA,
# #         col = transp.func(trajcols[st],0.1))
# # polygon(y = c(tmp[,lcic6],rev(tmp[,ucic6])),
# #         x = c(xplot,rev(xplot)),
# #         border = NA,
# #         col = transp.func(trajcols[st],0.1))
# lines(y = tmp[,medc],
#       x = xplot,
#       col = trajcols[st],
#       lwd = 2)
# text(st,
#      x = 2016,
#      y = tmp[which(tmp$Year == 2016),medc],
#      pos = 4,
#      col = trajcols[st])
# 
# 
# 
# 
# dev.off()
# 





strindices.o = strindices

#strindices = strindices[-which(is.na(strindices$index)),]




pdf(file = paste0(sp.dir,"/panel plot trajectories ",sp2d,".pdf"),
    width = 9,
    height = 9)

par(mfrow = c(3,3))
ee = estimates[1]
medc = ee
lcic = paste0(ee,".lci.50")
ucic = paste0(ee,".uci.50")
lcic9 = paste0(ee,".lci")
ucic9 = paste0(ee,".uci")
lcic6 = paste0(ee,".lci.90")
ucic6 = paste0(ee,".uci.90")
een = names(ee) 

ee2 = estimates[2]
lcic92 = paste0(ee2,".lci")
ucic92 = paste0(ee2,".uci")


# ee2 = estimates[3]
# lcic92 = paste0(ee2,".lci")
# ucic92 = paste0(ee2,".uci")

st = "Continental"


yup = min(c(max(strindices[which(strindices$Stratname == st),ee2])*3,
            max(strindices[which(strindices$Stratname == st),ucic92])))
xplot = 1974:2016

plot(0.001,
     col = "white",
     xlim = c(1974,2017),
     ylim = c(0,yup),
     xlab = "",
     ylab = "",
     xaxs = "i",
     yaxs = "i",
     bty = "l",
     main = paste(st,sp2)) 

ymid = yup/2
yrng1 = diff(c(ymid*0.8,ymid*1.2))

tmp = strindices[which(strindices$Stratname == st),]
nct = tmp$ncounts/max(tmp$ncounts)
### ncounts plotting
points(x = tmp$Year,
       y = (nct*yup),
       col = transp.func("darkorange",0.3),
       pch = 3)
axis(side = 4,
     at = seq(0,yup,length = 4),
     labels = round(seq(0,max(tmp$ncounts),length = 4)),
     col = "darkorange",
     col.ticks = "darkorange",
     col.axis = "darkorange")
mtext(side = 4,
      at = yup*0.5,
      "Ncounts",
      col = "darkorange",
      cex = 0.5)


polygon(y = c(tmp[,lcic],rev(tmp[,ucic])),
        x = c(xplot,rev(xplot)),
        border = NA,
        col = transp.func(trajcols[st],0.1))
polygon(y = c(tmp[,lcic9],rev(tmp[,ucic9])),
        x = c(xplot,rev(xplot)),
        border = NA,
        col = transp.func(trajcols[st],0.1))
polygon(y = c(tmp[,lcic6],rev(tmp[,ucic6])),
        x = c(xplot,rev(xplot)),
        border = NA,
        col = transp.func(trajcols[st],0.1))
lines(y = tmp[,medc],
      x = xplot,
      col = trajcols[st],
      lwd = 1)

polygon(y = c(tmp[,lcic92],rev(tmp[,ucic92])),
        x = c(xplot,rev(xplot)),
        border = NA,
        col = transp.func(trajcols[st],0.1))

lines(y = tmp[,ee2],
      x = xplot,
      col = trajcols[st],
      lwd = 1)

legend("topright",
       pch = c(3,1,18),
       col = c(transp.func("darkorange",0.3),
               transp.func("darkgreen",0.3),
               transp.func("red",0.5)),
       legend = c("ncounts",
                  "mean site effects",
                  "mean obs count"),
       cex = 0.5,
       bty = "n")

# tplot = trendsp[which(trendsp$stratname == st & trendsp$timespan == "1974-2016"),"trend.beta"]
# tplotl = trendsp[which(trendsp$stratname == st & trendsp$timespan == "1974-2016"),"trend.beta.lci"]
# tplotu = trendsp[which(trendsp$stratname == st & trendsp$timespan == "1974-2016"),"trend.beta.uci"]
# 
# mtext(side = 3,
#       at = 1995,
#       line = 0,
#       cex = 0.6,
#       paste("trend =",
#             round(tplot,1),
#             ":",
#             round(tplotl,1),
#             "-",
#             round(tplotu,1)))
# 

ee = estimates[2]
medc = ee
lcic = paste0(ee,".lci.50")
ucic = paste0(ee,".uci.50")
lcic9 = paste0(ee,".lci")
ucic9 = paste0(ee,".uci")
lcic6 = paste0(ee,".lci.90")
ucic6 = paste0(ee,".uci.90")
een = names(ee) 

ee2 = estimates[1]
lcic92 = paste0(ee2,".lci")
ucic92 = paste0(ee2,".uci")


# ee2 = estimates[3]
# lcic92 = paste0(ee2,".lci")
# ucic92 = paste0(ee2,".uci")

for(st in c("Continental",sts$stratname)){
  
yup = min(c(max(strindices[which(strindices$Stratname == st),medc])*3,
            max(strindices[which(strindices$Stratname == st),ucic6])))
xplot = 1974:2016

if(st == "Continental"){
plot(0.001,
     col = "white",
     xlim = c(1974,2017),
     ylim = c(0,yup),
     xlab = "",
     ylab = "",
     xaxs = "i",
     yaxs = "i",
     bty = "l",
     main = paste(st,"weight reg. abundance")) 
}else{
  plot(0.001,
       col = "white",
       xlim = c(1974,2017),
       ylim = c(0,yup),
       xlab = "",
       ylab = "",
       xaxs = "i",
       yaxs = "i",
       bty = "l",
       main = paste(st)) 
  
}
ymid = yup/2
yrng1 = diff(c(ymid*0.8,ymid*1.2))

  tmp = strindices[which(strindices$Stratname == st),]
 
  
  nct = tmp$ncounts/max(tmp$ncounts)
  ### ncounts plotting
  points(x = tmp$Year,
         y = (nct*yup),
         col = transp.func("darkorange",0.3),
         pch = 3)
  axis(side = 4,
       at = seq(0,yup,length = 4),
       labels = round(seq(0,max(tmp$ncounts),length = 4)),
       col = "darkorange",
       col.ticks = "darkorange",
       col.axis = "darkorange")
  mtext(side = 4,
        at = yup*0.5,
        "Ncounts",
        col = "darkorange",
        cex = 0.5)
  
  ### siteeffect plotting
  points(x = tmp$Year,
         y = (tmp$meansiteeff*yrng1)+ymid,
         col = transp.func("darkgreen",0.3),
         cex = 0.5)
  # abline(h = c(ymid,ymid*0.8),
  #        col = transp.func("darkgreen",0.3),
  #        lty = c(1,3))
  # 
  ### mean count plotting
  points(x = tmp$Year,
         y = tmp$obsmeancount,
         col = transp.func("red",0.6),
         pch = 18)
  if(length(which(tmp$obsmeancount > yup)) > 0){
    points(tmp[which(tmp$obsmeancount > yup),"Year"],
         y = rep(yup*0.99,length(which(tmp$obsmeancount > yup))),
         col = transp.func("red",0.5),
         pch = 17) 
  }
  
  polygon(y = c(tmp[,lcic],rev(tmp[,ucic])),
          x = c(xplot,rev(xplot)),
          border = NA,
          col = transp.func(trajcols[st],0.1))
  polygon(y = c(tmp[,lcic9],rev(tmp[,ucic9])),
          x = c(xplot,rev(xplot)),
          border = NA,
          col = transp.func(trajcols[st],0.1))
  polygon(y = c(tmp[,lcic6],rev(tmp[,ucic6])),
          x = c(xplot,rev(xplot)),
          border = NA,
          col = transp.func(trajcols[st],0.1))
  lines(y = tmp[,medc],
        x = xplot,
        col = trajcols[st],
        lwd = 1)

  polygon(y = c(tmp[,lcic92],rev(tmp[,ucic92])),
          x = c(xplot,rev(xplot)),
          border = NA,
          col = transp.func(trajcols[st],0.1))
  
  lines(y = tmp[,ee2],
        x = xplot,
        col = trajcols[st],
        lwd = 1)
  
  
  # tplot = trendsp[which(trendsp$stratname == st & trendsp$timespan == "1974-2016"),"trend.beta2"]
  # tplotl = trendsp[which(trendsp$stratname == st & trendsp$timespan == "1974-2016"),"trend.beta2.lci"]
  # tplotu = trendsp[which(trendsp$stratname == st & trendsp$timespan == "1974-2016"),"trend.beta2.uci"]
  # 
  # mtext(side = 3,
  #       at = 1995,
  #       line = 0,
  #       cex = 0.6,
  #       paste("trend =",
  #             round(tplot,1),
  #             ":",
  #             round(tplotl,1),
  #             "-",
  #             round(tplotu,1)))
  # 
  # 
}#st


dev.off()




















if(sp2 == storun[1]){
  indicesout <- strindices
  
  trendsoutall <- trendsp
  
}else{
  indicesout <- rbind(indicesout,strindices)
  trendsoutall <- rbind(trendsoutall,trendsp)
}



write.csv(indicesout,"shorebird new indices poisson seasonal GAM.csv")
write.csv(trendsoutall,"shorebird new trends poisson seasonal GAM.csv")

}#temp end if loop

}#temp end species loop





















#####################################

# calcualte the pooling factors so that there is some discussion about how much shrinkage there is in the annual indices
### consider fitting the year effects separately among regions so that there is more shrinkage
### consider how much weight to place in the year-effects if the shrinkage is very small
### if the shrinkage is small also determine if it is consistent across regions or dominated by one region (which should
### indicate whether the year effects should be include in the estimate of long-term trends: if the year effects are large
### and dominated by process in one or a few regions, then the slope term is likely a better indication of long-term change)
### also, consider calculating a change in the trend estimate from the sizecont[,] parameters. should provide some information on 
### how likely it is that the declines are ameliorating...




##### calculate rolling ten-year trends and look for patterns
### look at the slope of a regression through the year-effects.
##############################
setwd("C:/Shorebirds")
library(rjags)
library(RColorBrewer)
library(R2WinBUGS)
library(R.utils)#for GString function below
#### plot of paired (by species) slope and interval trend estimates

trendsoutall = read.csv("shorebird new trends poisson annual and seasonal GAM.csv", stringsAsFactors = F)
indicesout = read.csv("shorebird new indices poisson annual and seasonal GAM.csv", stringsAsFactors = F)

#################### on Friday morning: recreate the following basic plots (possibly not the pooling factor calcluations)
### also import the 2013 trends and compare the trend values



contt <- trendsoutall[which(trendsoutall$stratname == "Continental" &
                              trendsoutall$timespan == "1974-2016"),]


contt <- contt[order(contt[,"trend.index.slope"]),]
splot <- unique(contt$species)
contt$species = factor(contt$species,
                       ordered = T,
                       levels = splot)
splot <- levels(contt$species)

#ylims <- c(min(contt[,"trend lci"]),max(contt[,"trend uci"]))
ylims <- c(-15,15)
offs <- 0.2

pdf(width = 7,height = 6,file = "figure 1 annual GAM.pdf")
par(mar = c(11,4,1,1))
plot(x = c(1:length(splot)),
     y = rep(0,length(splot)),
     ylim = ylims,
     ylab = "Continental trend (%/year)",
     xlab = "",
     xaxt = "n",
     xlim = c(0.5,length(splot)+0.5),
     type = "n",
     bty = "l")
abline(h = 0,lty = 3)
axis(side = 1,at = 1:length(splot),labels = splot,las = 3)
points(x = 1:length(splot)-offs,y = contt[,"trend.index.slope"],col = "black",pch = 20,cex = 0.7)
points(x = 1:length(splot)+offs,y = contt[,"trend.index"],col = grey(0.5),pch = 20,cex = 0.7)

arrows(x0 = 1:length(splot)-offs,
       x1 = 1:length(splot)-offs,
       y0 = contt[,"trend.index.slope.lci"],
       y1 = contt[,"trend.index.slope.uci"],
       length = 0,
       col = "black")
arrows(x0 = 1:length(splot)+offs,
       x1 = 1:length(splot)+offs,
       y0 = contt[,"trend.index.lci"],
       y1 = contt[,"trend.index.uci"],
       length = 0,
       col = grey(0.5))
legend("topleft",
       legend = c("Mean survey-wide long-term slope (1974-2016)",
                  "End-point trend by size"),
       lwd = c(2,2),
       col = c("black",
               grey(0.5)),
       pch = 20,
       bty = "n",
       cex = 0.8)
dev.off()





#### figure 1 using the lm slope trends
# 
# 
# #ylims <- c(min(contt[,"trend lci"]),max(contt[,"trend uci"]))
# ylims <- c(-15,15)
# offs <- 0.2
# 
# pdf(width = 7,height = 6,file = "figure 1 test.pdf")
# par(mar = c(11,4,1,1))
# plot(x = c(1:length(splot)),
#      y = rep(0,length(splot)),
#      ylim = ylims,
#      ylab = "Continental trend (%/year)",
#      xlab = "",
#      xaxt = "n",
#      xlim = c(0.5,length(splot)+0.5),
#      type = "n",
#      bty = "l")
# abline(h = 0,lty = 3)
# axis(side = 1,at = 1:length(splot),labels = splot,las = 3)
# points(x = 1:length(splot)-offs,y = contt[,"trend.index.slope"],col = "black",pch = 20,cex = 0.7)
# points(x = 1:length(splot)+offs,y = contt[,"trend.size.slope"],col = grey(0.5),pch = 20,cex = 0.7)
# 
# arrows(x0 = 1:length(splot)-offs,
#        x1 = 1:length(splot)-offs,
#        y0 = contt[,"trend.index.slope.lci"],
#        y1 = contt[,"trend.index.slope.uci"],
#        length = 0,
#        col = "black")
# arrows(x0 = 1:length(splot)+offs,
#        x1 = 1:length(splot)+offs,
#        y0 = contt[,"trend.size.slope.lci"],
#        y1 = contt[,"trend.size.slope.uci"],
#        length = 0,
#        col = grey(0.5))
# legend("topleft",
#        legend = c("Mean survey-wide long-term slope (1974-2016)",
#                   "lm slope trend on annual indices"),
#        lwd = c(2,2),
#        col = c("black",
#                grey(0.5)),
#        pch = 20,
#        bty = "n",
#        cex = 0.8)
# dev.off()
# 






#### plotting the slope trends against the slope trends weighted by abundance
# 
# pdf(width = 7,height = 5,file = "figure 1 comp slope vs slope abund.pdf")
# par(mar = c(11,4,1,1))
# plot(x = c(1:length(splot)),
#      y = rep(0,length(splot)),
#      ylim = ylims,
#      ylab = "Continental trend (%/year)",
#      xlab = "",
#      xaxt = "n",
#      xlim = c(0.5,length(splot)+0.5),
#      type = "n",
#      bty = "l")
# abline(h = 0,lty = 3)
# axis(side = 1,at = 1:length(splot),labels = splot,las = 3)
# points(x = 1:length(splot)-offs,y = contt[,"trend.index.slope"],col = "black",pch = 20,cex = 0.7)
# points(x = 1:length(splot)+offs,y = contt[,"trend.index.slope2"],col = grey(0.5),pch = 20,cex = 0.7)
# 
# arrows(x0 = 1:length(splot)-offs,
#        x1 = 1:length(splot)-offs,
#        y0 = contt[,"trend.index.slope.lci"],
#        y1 = contt[,"trend.index.slope.uci"],
#        length = 0,
#        col = "black")
# arrows(x0 = 1:length(splot)+offs,
#        x1 = 1:length(splot)+offs,
#        y0 = contt[,"trend.index.slope2.lci"],
#        y1 = contt[,"trend.index.slope2.uci"],
#        length = 0,
#        col = grey(0.5))
# legend("topleft",
#        legend = c("Mean survey-wide long-term slope (1974-2016)",
#                   "Mean slope weighted by stratum abundance"),
#        lwd = c(2,2),
#        col = c("black",
#                grey(0.5)),
#        pch = 20,
#        bty = "n",
#        cex = 0.8)
# dev.off()
# 
# 









contt2 <- trendsoutall[which(trendsoutall$stratname == "Continental" &
                               trendsoutall$timespan == "2006-2016"),]

contt2$species = factor(contt2$species,
                        ordered = T,
                        levels = splot)

contt2 <- contt2[order(contt2$species),]

#### plotting the endpoint long-term trends against the endpoint short-term trends
# 
# pdf(width = 7,height = 5,file = "figure 1 comp endpoint long vs short.pdf")
# par(mar = c(11,4,1,1))
# plot(x = c(1:length(splot)),
#      y = rep(0,length(splot)),
#      ylim = ylims,
#      ylab = "Continental trend (%/year)",
#      xlab = "",
#      xaxt = "n",
#      xlim = c(0.5,length(splot)+0.5),
#      type = "n",
#      bty = "l")
# abline(h = 0,lty = 3)
# axis(side = 1,at = 1:length(splot),labels = splot,las = 3)
# points(x = 1:length(splot)-offs,y = contt[,"trend.size"],col = "black",pch = 20,cex = 0.7)
# points(x = 1:length(splot)+offs,y = contt2[,"trend.size"],col = grey(0.5),pch = 20,cex = 0.7)
# 
# arrows(x0 = 1:length(splot)-offs,
#        x1 = 1:length(splot)-offs,
#        y0 = contt[,"trend.size.lci"],
#        y1 = contt[,"trend.size.uci"],
#        length = 0,
#        col = "black")
# arrows(x0 = 1:length(splot)+offs,
#        x1 = 1:length(splot)+offs,
#        y0 = contt2[,"trend.size.lci"],
#        y1 = contt2[,"trend.size.uci"],
#        length = 0,
#        col = grey(0.5))
# legend("topleft",
#        legend = c("Survey-wide long-term endpoint (1974-2016)",
#                   "Survey-wide short-term endpoint (2006-2016)"),
#        lwd = c(2,2),
#        col = c("black",
#                grey(0.5)),
#        pch = 20,
#        bty = "n",
#        cex = 0.8)
# dev.off()
# 





##### long vs short comparison using lm slope trends

pdf(width = 7,height = 5,file = "figure 1 comp lm slope long vs short annual GAM.pdf")
par(mar = c(11,4,1,1))
plot(x = c(1:length(splot)),
     y = rep(0,length(splot)),
     ylim = ylims,
     ylab = "Continental trend (%/year)",
     xlab = "",
     xaxt = "n",
     xlim = c(0.5,length(splot)+0.5),
     type = "n",
     bty = "l")
abline(h = 0,lty = 3)
axis(side = 1,at = 1:length(splot),labels = splot,las = 3)
points(x = 1:length(splot)-offs,y = contt[,"trend.index.slope"],col = "black",pch = 20,cex = 0.7)
points(x = 1:length(splot)+offs,y = contt2[,"trend.index.slope"],col = grey(0.5),pch = 20,cex = 0.7)

arrows(x0 = 1:length(splot)-offs,
       x1 = 1:length(splot)-offs,
       y0 = contt[,"trend.index.slope.lci"],
       y1 = contt[,"trend.index.slope.uci"],
       length = 0,
       col = "black")
arrows(x0 = 1:length(splot)+offs,
       x1 = 1:length(splot)+offs,
       y0 = contt2[,"trend.index.slope.lci"],
       y1 = contt2[,"trend.index.slope.uci"],
       length = 0,
       col = grey(0.5))
legend("topleft",
       legend = c("Survey-wide long-term lm slope (1974-2016)",
                  "Survey-wide short-term lm slope (2006-2016)"),
       lwd = c(2,2),
       col = c("black",
               grey(0.5)),
       pch = 20,
       bty = "n",
       cex = 0.8)
dev.off()







contt3 <- trendsoutall[which(trendsoutall$stratname == "Continental" &
                               trendsoutall$timespan == "1995-2016"),]

contt3$species = factor(contt3$species,
                        ordered = T,
                        levels = splot)

contt3 <- contt3[order(contt3$species),]



contt4 <- trendsoutall[which(trendsoutall$stratname == "Continental" &
                               trendsoutall$timespan == "1974-1995"),]

contt4$species = factor(contt4$species,
                        ordered = T,
                        levels = splot)

contt4 <- contt4[order(contt4$species),]

#### plotting the endpoint long-term trends against the endpoint short-term trends
# 
# pdf(width = 7,height = 5,file = "figure 1 comp endpoint first20 vs last 20.pdf")
# par(mar = c(11,4,1,1))
# plot(x = c(1:length(splot)),
#      y = rep(0,length(splot)),
#      ylim = ylims,
#      ylab = "Continental trend (%/year)",
#      xlab = "",
#      xaxt = "n",
#      xlim = c(0.5,length(splot)+0.5),
#      type = "n",
#      bty = "l")
# abline(h = 0,lty = 3)
# axis(side = 1,at = 1:length(splot),labels = splot,las = 3)
# points(x = 1:length(splot)-offs,y = contt4[,"trend.size"],col = "black",pch = 20,cex = 0.7)
# points(x = 1:length(splot)+offs,y = contt3[,"trend.size"],col = grey(0.5),pch = 20,cex = 0.7)
# 
# arrows(x0 = 1:length(splot)-offs,
#        x1 = 1:length(splot)-offs,
#        y0 = contt4[,"trend.size.lci"],
#        y1 = contt4[,"trend.size.uci"],
#        length = 0,
#        col = "black")
# arrows(x0 = 1:length(splot)+offs,
#        x1 = 1:length(splot)+offs,
#        y0 = contt3[,"trend.size.lci"],
#        y1 = contt3[,"trend.size.uci"],
#        length = 0,
#        col = grey(0.5))
# legend("topleft",
#        legend = c("Survey-wide first phase endpoint (1974-1995)",
#                   "Survey-wide last phase endpoint (1995-2016)"),
#        lwd = c(2,2),
#        col = c("black",
#                grey(0.5)),
#        pch = 20,
#        bty = "n",
#        cex = 0.8)
# dev.off()
# 
# 
# 



### same as above using the lm slope estimates



pdf(width = 7,height = 5,file = "figure 1 comp lm slope first20 vs last 20 annual GAM.pdf")
par(mar = c(11,4,1,1))
plot(x = c(1:length(splot)),
     y = rep(0,length(splot)),
     ylim = ylims,
     ylab = "Continental trend (%/year)",
     xlab = "",
     xaxt = "n",
     xlim = c(0.5,length(splot)+0.5),
     type = "n",
     bty = "l")
abline(h = 0,lty = 3)
axis(side = 1,at = 1:length(splot),labels = splot,las = 3)
points(x = 1:length(splot)-offs,y = contt4[,"trend.index.slope"],col = "black",pch = 20,cex = 0.7)
points(x = 1:length(splot)+offs,y = contt3[,"trend.index.slope"],col = grey(0.5),pch = 20,cex = 0.7)

arrows(x0 = 1:length(splot)-offs,
       x1 = 1:length(splot)-offs,
       y0 = contt4[,"trend.index.slope.lci"],
       y1 = contt4[,"trend.index.slope.uci"],
       length = 0,
       col = "black")
arrows(x0 = 1:length(splot)+offs,
       x1 = 1:length(splot)+offs,
       y0 = contt3[,"trend.index.slope.lci"],
       y1 = contt3[,"trend.index.slope.uci"],
       length = 0,
       col = grey(0.5))
legend("topleft",
       legend = c("Survey-wide first phase lm slope (1974-1995)",
                  "Survey-wide last phase lm slope (1995-2016)"),
       lwd = c(2,2),
       col = c("black",
               grey(0.5)),
       pch = 20,
       bty = "n",
       cex = 0.8)
dev.off()





# 
# 
# 
regtta <- trendsoutall[-which(trendsoutall$stratname == "Continental"),]
# 
# ####### calculate pooling factors on year effects for each species and add them to contt file
# for (sp2 in storun) {
#   sp2d <- snrun[sp2]
#   sp <- gsub(sp2d,pattern = " ",replace = "_")
#  sp.dir <- paste(getwd(),"/",sp2d,"5/",sep = "")  
#   
#   if (file.exists(paste(sp.dir,"cd13.RData", sep = ""))) {
#   load(paste(sp.dir,"cd13.RData", sep = ""))  
#   wann <- grep(names(cd13[[1]][1,]),pattern = "ann[",fixed = T)
#   ann <- matrix(NA,nrow = length(cd13[[1]][,1])*3,ncol = length(wann))
#   e.ann <- ann
#   sl.ann <- ann[,1]
#   for(j in 1:length(wann)){
#     ji <- wann[j]
#   ann[,j] <- c(cd13[[1]][,ji],cd13[[2]][,ji],cd13[[3]][,ji])    #rows = mcmc iterations, cols = years
#   } 
#   for(j in 1:length(ann[,1])){  #for each mcmc iteration
#    
#     e.ann[j,] <- ann[j,]-mean(ann[j,])  # for each j-iteration, difference between each year's value and the mean across years (which should be close to zero, but is not exactly zero).
#  a <- ann[j,]
#  y = 1:length(a)
#  m.ann <- lm(a~y)
#  sl.ann[j] <- as.numeric(coef(m.ann)["y"])
#  
#  
#   } 
#   
#   
#   #aa <- apply (ann, 2, mean)
#   
#   lambda.a <- 1 - var (apply (e.ann, 2, mean))/ mean (apply (e.ann, 1, var))#pooling factor - see 21.6_Summarizing the amount of partial pooling.R in the folder M:\My Documents\R\Functions\pooling factor calculation from Gelman and Hill  
#   lambda.a2 <- 1 - var (apply (ann, 2, mean))/ mean (apply (ann, 1, var))
#     # var (apply (ann, 2, mean)) = variance of posterior means of the annual effects
#   #mean (apply (ann, 1, var)) = posterior mean of the variance of year effects
# #  rsqa <- 1 - mean (apply (e.ann, 1, var))/ mean (apply (ann, 1, var))#pooling factor - see 21.6_Summarizing the amount of partial pooling.R in the folder M:\My Documents\R\Functions\pooling factor calculation from Gelman and Hill  
#   contt[which(contt$sp == sp2),"year.pooling"] <- lambda.a
#   contt[which(contt$sp == sp2),"year.pooling2"] <- lambda.a2
#   contt[which(contt$sp == sp2),"year.slope"] <- mean(sl.ann)
#   contt[which(contt$sp == sp2),"year.slope.se"] <- round(sd(sl.ann),4)
#   contt[which(contt$sp == sp2),"year.slope.lci"] <- round(quantile(sl.ann,probs = 0.025,names = F),4)
#   contt[which(contt$sp == sp2),"year.slope.uci"] <- round(quantile(sl.ann,probs = 0.975,names = F),4)
#  
#  m.a <- apply (ann, 2, mean)
#  uci.a <-  apply (ann, 2, quantile,probs = 0.975)
#  lci.a <-  apply (ann, 2, quantile,probs = 0.025)
#   
#  png(file = paste("year effects in time",sp2,".png"),height = 8*72,width = 8*72)
# plot(y = m.a,x = 1:length(m.a),xaxt = "n",ylim = c(min(lci.a),max(uci.a)),main = "year effects over time")
# 
# #slt <- sample(1:length(ann[,1]),size = 100,replace = F)
# #for(jj in 1:100){
# #
# #abline(a = mean(ann[jj,]),b = sl.ann[slt[jj]],col = grey(0.9)) 
# #
# #}#end jj random samples of the slopes estimated for year effects
# ##abline(a = m.a[floor(mean(1:length(m.a)))],b = mean(sl.ann)) 
# #abline(a = mean(m.a),b = mean(sl.ann),col = "blue") 
# abline(h = 0)
# arrows(y0 = lci.a,y1 = uci.a,x0 = 1:length(m.a),x1 = 1:length(m.a),length = 0,col = grey(0.5))
# points(y = m.a,x = 1:length(m.a))
# axis(side = 1,at = c(1,8,18,28,38),labels = c("1974",as.character(seq(from = 1980,to = 2010,by = 10))))  
# 
# dev.off()

for (sp2 in storun) {
  
sp2d <- sps[which(sps$Species == sp2),"dirnames"]
#sp2d <- snrun[sp2]
#dir.spsp <- paste0("output/",sp2d)
sp.dir <- paste0("output/",sp2d,"poisson seasonal GAM")
# dtf <- df2[which(df2$Species == sp2),]

#col.st <- brewer.pal(8, "Dark2")
#names(col.st) <- sts$stratname
#nstrat <- length(sts$stratname)

if (file.exists(paste(sp.dir,"/cd13.RData", sep = ""))) {
  
  
 #  sp2d <- snrun[sp2]
 #  sp <- gsub(sp2d,pattern = " ",replace = "_")
 # sp.dir <- paste(getwd(),"/",sp2d,"5/",sep = "")
 # 
 #  if (file.exists(paste(sp.dir,"cd13.RData", sep = ""))) {
  load(paste(sp.dir,"/cd13.RData", sep = ""))

######### calculate the posterior distribution of the difference between the regional trend and the global trend
  betas <- grep(names(cd13[[1]][1,]),pattern = "beta[",fixed = T)
  betasM <- matrix(NA,nrow = length(cd13[[1]][,1])*3,ncol = length(betas))
  betaH <- grep(names(cd13[[1]][1,]),pattern = "betaH",fixed = T)
  betaHM <- matrix(c(cd13[[1]][,betaH],cd13[[2]][,betaH],cd13[[3]][,betaH]),nrow = length(cd13[[1]][,1])*3,ncol = 1)
 
  for(j in 1:length(betas)){
    ji <- betas[j]
    betasM[,j] <- c(cd13[[1]][,ji],cd13[[2]][,ji],cd13[[3]][,ji])-betaHM
    
    regtta[which(regtta$species == sp2 & regtta$strat == j),"trenddifH"] <- mean(betasM[,j])
    regtta[which(regtta$species == sp2 & regtta$strat == j),"trenddifHlci"] <- quantile(betasM[,j],0.025)
    regtta[which(regtta$species == sp2 & regtta$strat == j),"trenddifHuci"] <- quantile(betasM[,j],0.975)
    regtta[which(regtta$species == sp2 & regtta$strat == j),"trenddifprec"] <- 1/var(betasM[,j])
    
  } 
  
  ####### pooling factor on the betas

  
  
  
  cl4 <- "betaH"
  ri4 <- c(as.numeric(cd13[[1]][,cl4]),as.numeric(cd13[[2]][,cl4]),as.numeric(cd13[[3]][,cl4]))
  Bet <- round(mean(ri4,na.rm = T),4)
  
  for (st in 1:length(betas)) {
    #j <- which(pool.p$strat == st)
    wregtt = which(regtta$strat == st & regtta$species == sp2)
    
    regtmp <- unique(regtta[wregtt,c("stratname")])
    
    cl <- paste("beta[",st,"]",sep = "") 
    ri <- c(as.numeric(cd13[[1]][,cl]),as.numeric(cd13[[2]][,cl]),as.numeric(cd13[[3]][,cl]))
    
    cl2 <- paste("sdbeta",sep = "") 
    ri2 <- c(as.numeric(cd13[[1]][,cl2]),as.numeric(cd13[[2]][,cl2]),as.numeric(cd13[[3]][,cl2]))
    
    
    # cl3 <- paste("strata[",st,"]",sep = "") 
    # ri3 <- c(as.numeric(cd13[[1]][,cl3]),as.numeric(cd13[[2]][,cl3]),as.numeric(cd13[[3]][,cl3]))
    # #       
    # pooling factor (pg 478, Gelman and Hill) for each stratum and year = [in R code] (sd(yearcar[y,j])/sdyear[y])^2
    regtta[wregtt,"beta"] <- mean(ri,na.rm = T)
    regtta[wregtt,"beta.se"] <- sd(ri,na.rm = T)
    regtta[wregtt,"beta.lci"] <- quantile(ri,probs = 0.025,na.rm = T)
    regtta[wregtt,"beta.uci"] <- quantile(ri,probs = 0.975,na.rm = T)
    
    regtta[wregtt,"pool"] <- min((sd(ri-ri4,na.rm = T)/mean(ri2,na.rm = T))^2,1)

  }#st
  
  
  
  
  
  
    }#if cd13 exists
    
}#sp2



write.csv(regtta,"regional trends w pooling factor beta.csv")





regtta[,"stratnamef"] <- factor(regtta[,"stratname"],
                              levels = c("Atlantic Canada",
                                         "Northeast US Coastal",
                                         "Ontario",
                                         "East Inland",
                                         "Midcontinental",
                                         "Pacific and Intermountain",
                                         "Southeast Coastal",
                                         "Texas Coastal"),
                              ordered = T)
############## plot of regional slope trends as differences from the global trends



mregt <- tapply(regtta$trenddifH,regtta$stratname,mean)
sdregt <- tapply(regtta$trenddifH,regtta$stratname,sd)

sdregt
mregt
boxplot(trenddifH~stratname,data = regtta)

######### full trend comparison model (incorporating the distribution of hyperparameters among regions)
Tsum <- data.frame(stratnamef = levels(regtta$stratnamef),mean = NA,sd = NA, lci = NA, lqu = NA,med = NA,uqu = NA,uci = NA)

nspecies <- NA 
betahat1 <- tapply(regtta$trenddifH,regtta[,c("stratnamef","species")],mean)
prec.betahat1 <-tapply(regtta$trenddifprec,regtta[,c("stratnamef","species")],mean)
betahat <- matrix(NA,nrow = length(levels(regtta$stratnamef)),ncol = length(unique(regtta$species)))
prec.betahat <- betahat
i = 0
for(rr in levels(regtta$stratnamef)[-8]){
  i = i+1
  betah <- regtta[which(regtta$stratnamef == rr &
                          regtta$timespan == "1974-2016"),"trenddifH"]
  betahat[i,1:length(betah)] <- betah
  prec.betah <- regtta[which(regtta$stratnamef == rr &
                               regtta$timespan == "1974-2016"),"trenddifprec"]
  prec.betahat[i,1:length(prec.betah)] <- prec.betah  
}

for(i in 1:nrow(betahat)){
  nspecies[i]<- length(which(!is.na(betahat[i,]))) 
}

nspecies# 
nspecies = nspecies[-length(nspecies)]
betahat = betahat[-nrow(betahat),]
prec.betahat = prec.betahat[-nrow(prec.betahat),]

#    species1<- which(!is.na(betahat[1,]))
# species2<- which(!is.na(betahat[2,]))
# species3<- which(!is.na(betahat[3,]))
# species4<- which(!is.na(betahat[4,]))
# species5<- which(!is.na(betahat[5,]))
# species6<- which(!is.na(betahat[6,]))
# species7<- which(!is.na(betahat[7,]))
# species8<- which(!is.na(betahat[8,]))
# nspecies1 <- length(species1)
# nspecies2 <- length(species2)
# nspecies3 <- length(species3)
# nspecies4 <- length(species4)
# nspecies5 <- length(species5)
# nspecies6 <- length(species6)
# nspecies7 <- length(species7)
# nspecies8 <- length(species8)
# 
# mspecies1<- which(is.na(betahat[1,]))
# mspecies2<- which(is.na(betahat[2,]))
# mspecies3<- which(is.na(betahat[3,]))
# mspecies4<- which(is.na(betahat[4,]))
# mspecies5<- which(is.na(betahat[5,]))
# mspecies6<- which(is.na(betahat[6,]))
# mspecies7<- which(is.na(betahat[7,]))
# mspecies8<- which(is.na(betahat[8,]))
# nmspecies1 <- length(mspecies1)
# nmspecies2 <- length(mspecies2)
# nmspecies3 <- length(mspecies3)
# nmspecies4 <- length(mspecies4)
# nmspecies5 <- length(mspecies5)
# nmspecies6 <- length(mspecies6)
# nmspecies7 <- length(mspecies7)
# nmspecies8 <- length(mspecies8)


  d.y <- list(nregions = nrow(betahat),
#               nspecies1 = nspecies1,
#               nspecies2 = nspecies2,
#               nspecies3 = nspecies3,
#               nspecies4 = nspecies4,
#               nspecies5 = nspecies5,
#               nspecies6 = nspecies6,
#               nspecies7 = nspecies7,
#               nspecies8 = nspecies8,
#               species1 = species1,
#               species2 = species2,
#               species3 = species3,
#               species4 = species4,
#               species5 = species5,
#               species6 = species6,
#               species7 = species7,
#               species8 = species8,
#               mspecies1 = mspecies1,
#               #mspecies2 = mspecies2,
#               mspecies3 = mspecies3,
#               mspecies4 = mspecies4,
#               mspecies5 = mspecies5,
#               mspecies6 = mspecies6,
#               mspecies7 = mspecies7,
#               mspecies8 = mspecies8,
#               nmspecies1 = nmspecies1,
#               #nmspecies2 = nmspecies2, # there are no missing species in region 2
#               #nmspecies3 = nmspecies3,
#               nmspecies4 = nmspecies4,
#               nmspecies5 = nmspecies5,
#               #nmspecies6 = nmspecies6,
#               nmspecies7 = nmspecies7,
#               nmspecies8 = nmspecies8,
              nspecies = nspecies,
              betahat=betahat, 
              prec.betahat = prec.betahat)
  param <- c("Tmu","T","sd.beta","sd.mu.beta","mu.beta","dift","pt")
  model.f <- "model_mubeta2.txt" #laptop and desktop
  initial <- list(list(prec.mu.beta = 1, prec.beta = rep(1,nrow(betahat)),mu.beta = runif(nrow(betahat),-0.05,0.05)))
  
  sim2 <- bugs(d.y,initial, param, model.f, n.chains = 1, n.iter = 50000,n.burnin = 20000,n.thin = 1, 
               bugs.directory="C:\\WinBUGS14",program="WinBUGS", debug = F)    #my desktop
  

save(file = paste0("sim2 all regions.RData"),list = c("sim2","d.y","betahat1","prec.betahat1"))
pt <- sim2$mean$pt  


regtta[,"trenddifCI2"] <- (((exp(regtta$trenddifHuci)-1)*100-(exp(regtta$trenddifHlci)-1)*100)/2)
range(regtta$trenddifCI2)
regtta[,"trendCI2"] <- (regtta[,"trend.beta.uci"]-regtta[,"trend.beta.lci"])/2
range(regtta$trendCI2)

i = 0  
for(rr in levels(regtta$stratnamef)[-8]){
 i = i+1 
 cl <- paste0("Tmu[",i,"]")
     Tsum[which(Tsum$stratnamef == rr),2:8] <- sim2$summary[cl,]

  Tsum[which(Tsum$stratnamef == rr),"rawmean"] <- mean(betahat[i,],na.rm = T)
  Tsum[which(Tsum$stratnamef == rr),"rawsd"] <- sd(betahat[i,],na.rm = T)
  Tsum[which(Tsum$stratnamef == rr),"nspecies"] <- nspecies[i]
  png(height = 7*72,width = 7*72, file = paste0("histogram raw species trends ",rr,".png"))
  hist((exp(betahat[i,])-1)*100,freq = T,breaks = c(seq(from = -100,to = 100,by = 1)),xlim = c(-10,10))
  
  dev.off()
  tshrnk <- sim2$summary[grep(names(sim2$summary[,1]),pattern = paste0("T[",i),fixed = T),"mean"]
  tshrnkCIw <- (sim2$summary[grep(names(sim2$summary[,1]),pattern = paste0("T[",i),fixed = T),"97.5%"]-sim2$summary[grep(names(sim2$summary[,1]),pattern = paste0("T[",i),fixed = T),"2.5%"])/2
  
  png(height = 7*72,width = 7*72, file = paste0("histogram shrunk species trends ",rr,".png"))
  hist(sim2$summary[grep(names(sim2$summary[,1]),pattern = paste0("T[",i),fixed = T),"mean"],freq = T,breaks = c(seq(from = -100,to = 100,by = 1)),xlim = c(-10,10))
  
  dev.off()
 
 
 
 
 betahat.y <- regtta[which(regtta$stratnamef == rr &
                             regtta$timespan == "1974-2016"),]  
  betahat.y[,"trenddifCIw"] <- (((exp(betahat.y$trenddifHuci)-1)*100-(exp(betahat.y$trenddifHlci)-1)*100)/2)/8
  betahat.y[,"trendCIw"] <- (betahat.y[,"trend.beta.uci"]-betahat.y[,"trend.beta.lci"])/2
# betahat.y[,"trenddifCI2"] <- (((exp(betahat.y$trenddifHuci)-1)*100-(exp(betahat.y$trenddifHlci)-1)*100)/2)
 
 betahat.y[which(betahat.y[,"trenddifCIw"] > 0.7 & betahat.y[,"trenddifCIw"] < 1.5),"trenddifCIw"] <- 0.7
  betahat.y[which(betahat.y[,"trenddifCIw"] >= 1.5),"trenddifCIw"] <- 0.9
  betahat.y[which(betahat.y[,"trenddifCIw"] < 0.4),"trenddifCIw"] <- 0.2
 betahat.y[,"trendCI3"] <- betahat.y[,"trendCI2"]/8
 betahat.y[which(betahat.y[,"trendCI3"] > 1),"trendCI3"] <- 1
 
 
  
  pdf(height = 7,width = 10, file = paste0("shrinkage plot species trends 2 ",rr,".pdf"))
  plot(y = rep(1,nrow(betahat.y)),x = (exp(betahat.y$trenddifH)-1)*100,ylim = c(-0.5,2),
       xlab = "Difference between regional and continental trend (%/year)",
       ylab = "",yaxt = "n",pch = 20,cex = 0.8,
       col = grey(betahat.y$trenddifCIw),
       xlim = range((exp(betahat.y$trenddifH)-1)*100),
       bty = "n")
  abline(v = 0,lty = 3,col = grey(0.4))
  arrows(x0 = (exp(betahat.y$trenddifH)-1)*100,x1 = (exp(betahat.y$trenddifH)-1)*100,
         y0 = 1,#1-(betahat.y[,"trendCIw"]/max(betahat.y[,"trendCIw"])),
         #y1 = 1.5,
          y1 = 1+(betahat.y[,"trendCI3"]),
         length = 0,
         col = grey(betahat.y$trenddifCIw))
  points(y = rep(0.5,nrow(betahat.y)),x = tshrnk,pch = 20,cex = 0.8,
         col = grey(betahat.y$trenddifCIw))
  arrows(x0 = tshrnk,x1 = tshrnk,
         y0 = 0.5,
         y1 = 0.5-(betahat.y[,"trendCI3"]),
         # y1 = 1+(betahat.y[,"trendCI3"]),
         length = 0,
         col = grey(betahat.y$trenddifCIw))
  
  arrows(x0 = tshrnk,x1 = (exp(betahat.y$trenddifH)-1)*100,
         y0 = 0.5,
          y1 = 1,
         length = 0,
         col = grey(betahat.y$trenddifCIw))
  axis(side = 2,at = c(-0.5,0,0.5),labels = c("8+","4","0"),tcl = -0.25,las = 1,lwd = 2)
 axis(side = 2,at = c(1,1.5,2),labels = c("0","4","8+"),tcl = -0.25,las = 1,lwd = 2)
 axis(side = 2,at = seq(from = -0.5,to = 0.5,length = 9),labels = F,tcl = -0.1,lwd = 2)
 axis(side = 2,at = seq(from = 1,to = 2,length = 9),labels = F,tcl = -0.1,lwd = 2)
 mtext(side = 2, at = c(0,1.5),text = rep("Half-width of trend CI (%/year)",2),line = 2)
 
 axis(side = 1,at = Tsum[which(Tsum$stratnamef == rr),2],labels = F,tcl = 0.5,lwd = 3)
 axis(side = 1,at = Tsum[which(Tsum$stratnamef == rr),2],labels = F,tcl = -0.5,lwd = 3)
 axis(side = 1,at = c(-20,20),labels = F,tcl = 0.5)
 
  dev.off() 
  
}#rr




# save.image("figures and regional mean trend analyses April.Rdata")
# 
# write.csv(contt,"continental trends w pooling factor year.csv")
# 
# load("figures and regional mean trend analyses April.Rdata")

# oropt <- options()
# options(scipen=10)
# 
# jjc = 0
# for(jj in c(0,9,18,27)){
#   jjc = jjc+1
#   
#   precsp <- contt[1:9+jj,"sp"]
#   
#   
#   
#   
#   inds3 <- indout[which(indout$stratname == "Continental" & indout$sp %in% precsp),]
#   
# mag = 8  
# #  pdf(width = 7,height = 7,file = paste0("Annual indices panel ",jjc," ",GString$getBuiltinTime(format="%H-%M-%S"),".pdf"))
#   png(width = 7*72*mag,height = 7*72*mag,file = paste0("Annual indices panel ",jjc," ",GString$getBuiltinTime(format="%H-%M-%S"),".png"))
#   
#   par(mfrow = c(3,3),mar = c(1,1,1,1)*mag,oma = c(4,mag+3,1,1))
#   
#   for(ssp in precsp){
#     ssp2 <- contt[jj,"species"]
#   inds2 <- inds3[which(inds3$sp == ssp),]  
#   #for(ss in sts$stratname) {
#   #  lines(indst[which(indst$stratname == ss),"year"],indst[which(indst$stratname == ss),"indsmooth"],col = col.st[ss],lwd = 1)
#   #  mtext(ss,side = 4,at = indst[which(indst$stratname == ss),"indsmooth"][40],line = 0.5,col = col.st[ss],adj = 0,las = 1)
#   #if(ssp == precsp[4]){ 
#                                
#     plot(inds2$YEAR,inds2[,"50%"],ylim = c(0,max(inds2[,"50%"])*2),type = "l",
#                                     lwd = 3,xlab = "",xaxt = "n", 
#                                     ylab = "",cex.axis = mag)  
#      mtext(ssp,side = 3,at = 1993,line = -1.3*mag,cex = mag)                        #}  
#   #}#ss
#   axis(side = 1,at = seq(1973,2013,by = 10),labels = seq(1973,2013,by = 10),cex = mag,cex.axis = mag,line = mag-2,tick = F)
#     axis(side = 1,at = seq(1973,2013,by = 10),labels = F,tick = T,tcl = -0.5)
#   lines(inds2$YEAR,inds2[,"X2.5."],lty = 1,lwd = 1,col = grey(0.5),cex = mag)
#   lines(inds2$YEAR,inds2[,"X97.5."],lty = 1,lwd = 1,col = grey(0.5),cex = mag)
#   lines(inds2$YEAR,inds2[,"50%"],lwd = 1.5,cex = mag)
#   points(inds2$YEAR,inds2[,"50%"],cex = mag/2)
#   }
#   
#   mtext(side = 2,"Continental Annual Index of Abundance",line = 2,outer = T,cex = mag,cex.axis = mag)
#   #legend("topright",legend = names(col.st)[1:nstrat],col = col.st[1:nstrat],lty = 1,lwd = 1,bty = "n",text.col = col.st[1:nstrat])
#   dev.off()
#   
# }#jjc
# 
# 
# options(oropt)

# ######### simple trend comparison model (ignoring the distribution of hyperparameters among regions)
# Tsum <- data.frame(stratname = levels(regtta$stratname),mean = NA,sd = NA, lci = NA, lqu = NA,med = NA,uqu = NA,uci = NA)
# for (rr in levels(regtta$stratname)) {
#   
#   
#   betahat.y <- regtta[which(regtta$stratname == rr),]
#   d.y <- list(nspecies = nrow(betahat.y), betahat=round(betahat.y$trenddifH,4), prec.betahat = round(betahat.y$trenddifprec,4))
#   param <- c("Tmu","T","sd.beta")
#   model.f <- "model_mubeta.txt" #laptop and desktop
#   initial <- list(list(prec.beta = 1, mu.beta = 0))
#   
#   sim2 <- bugs(d.y,initial, param, model.f, n.chains = 1, n.iter = 50000,n.burnin = 20000,n.thin = 1, bugs.directory="C:\\WinBUGS14",program="WinBUGS", debug = F)    #my desktop
# save(file = paste0("sim2 ",rr,".RData"),list = "sim2")
# 
#   
#   Tsum[which(Tsum$stratname == rr),2:8] <- sim2$summary["Tmu",]
# 
#   Tsum[which(Tsum$stratname == rr),"rawmean"] <- mean(betahat.y$trenddifH)
#   Tsum[which(Tsum$stratname == rr),"rawsd"] <- sd(betahat.y$trenddifH)
#   Tsum[which(Tsum$stratname == rr),"nspecies"] <- nrow(betahat.y)
#   png(height = 7*72,width = 7*72, file = paste0("histogram raw species trends ",rr,".png"))
#   hist((exp(betahat.y$trenddifH)-1)*100,freq = T,breaks = c(seq(from = -100,to = 100,by = 1)),xlim = c(-10,10))
#   
#   dev.off()
#   tshrnk <- sim2$summary[grep(names(sim2$summary[,1]),pattern = "T[",fixed = T),"mean"]
#   tshrnkCIw <- (sim2$summary[grep(names(sim2$summary[,1]),pattern = "T[",fixed = T),"97.5%"]-sim2$summary[grep(names(sim2$summary[,1]),pattern = "T[",fixed = T),"2.5%"])/2
#   
#   png(height = 7*72,width = 7*72, file = paste0("histogram shrunk species trends ",rr,".png"))
#   hist(sim2$summary[grep(names(sim2$summary[,1]),pattern = "T[",fixed = T),"mean"],freq = T,breaks = c(seq(from = -100,to = 100,by = 1)),xlim = c(-10,10))
#   
#   dev.off()
# 
#   betahat.y[,"trenddifCIw"] <- (((exp(betahat.y$trenddifHuci)-1)*100-(exp(betahat.y$trenddifHlci)-1)*100)/2)/8
#   betahat.y[,"trendCIw"] <- (betahat.y[,"trend uci"]-betahat.y[,"trend lci"])/2
#   betahat.y[which(betahat.y[,"trenddifCIw"] > 0.7 & betahat.y[,"trenddifCIw"] < 1.5),"trenddifCIw"] <- 0.7
#   betahat.y[which(betahat.y[,"trenddifCIw"] >= 1.5),"trenddifCIw"] <- 0.9
#   betahat.y[which(betahat.y[,"trenddifCIw"] < 0.4),"trenddifCIw"] <- 0.2
#   
#   
#   pdf(height = 7,width = 10, file = paste0("shrinkage plot species trends ",rr,".pdf"))
#   plot(y = rep(1,nrow(betahat.y)),x = (exp(betahat.y$trenddifH)-1)*100,ylim = c(-0.25,1.75),
#        xlab = "Difference between regional and continental trend (%/year)",
#        ylab = "",yaxt = "n",pch = 20,cex = 0.5,
#        col = grey(betahat.y$trenddifCIw),
#        xlim = range((exp(betahat.y$trenddifH)-1)*100),
#        bty = "l")
#   abline(v = Tsum[which(Tsum$stratname == rr),2],lty = 3)
#   arrows(x0 = (exp(betahat.y$trenddifH)-1)*100,x1 = (exp(betahat.y$trenddifH)-1)*100,
#          y0 = 1,#1-(betahat.y[,"trendCIw"]/max(betahat.y[,"trendCIw"])),
#          y1 = 1.5,
#          # y1 = 1+(betahat.y[,"trendCIw"]/max(betahat.y[,"trendCIw"])),
#          length = 0,
#          col = grey(betahat.y$trenddifCIw))
#   points(y = rep(0.5,nrow(betahat.y)),x = tshrnk,pch = 20,cex = 0.5,
#          col = grey(betahat.y$trenddifCIw))
#   arrows(x0 = tshrnk,x1 = tshrnk,
#          y0 = 0.5,
#          y1 = 0,
#          # y1 = 0.5-(tshrnkCIw/max(betahat.y[,"trendCIw"])),
#          length = 0,
#          col = grey(betahat.y$trenddifCIw))
# 
#   arrows(x0 = tshrnk,x1 = (exp(betahat.y$trenddifH)-1)*100,
#          y0 = 0.5,
#          y1 = 1,
#          length = 0,
#          col = grey(betahat.y$trenddifCIw))
# #axis(side = 2,at = c(0.5,1),labels = c("Shrunk","Raw"))
#   
#   dev.off() 
#   
# }#rr
# 
# 
# 


Tsum = Tsum[-8,]

pdf(width = 7,height = 7,file = paste0("mean dif btw regional and continental trends 2 ",".pdf"))
par(mar = c(13,4,1,1))
plot(1:nrow(Tsum),Tsum$mean,pch = 20,ylab = "Mean difference from continental trend (%/year)",
     xlab = "",xaxt = "n",bty = "l",ylim = c(min(Tsum$lci),max(Tsum$uci)))
abline(h = 0,col = grey(0.7))
arrows(x0 = 1:nrow(Tsum),x1 = 1:nrow(Tsum),y0 = Tsum$lci,y1 = Tsum$uci,length = 0,col = grey(0.4))
points(1:nrow(Tsum),Tsum$mean,pch = 20)
#points(1:nrow(Tsum),(exp(Tsum$rawmean)-1)*100)
axis(side = 1,at = 1:nrow(Tsum),labels = paste0(levels(regtta$stratnamef)[-8]," (",Tsum$nspecies,")"),las = 3,cex = 1.2)
dev.off()







##################### compare trend estimates among different models

#2013 vs 2017 Poisson model

t13 = read.csv("archived species results/all trendout.csv", stringsAsFactors = F)
i13 = read.csv("archived species results/all indices.csv", stringsAsFactors = F)

t17 = read.csv("shorebird new trends poisson seasonal GAM.csv", stringsAsFactors = F)
i17 = read.csv("shorebird new indices poisson seasonal GAM.csv", stringsAsFactors = F)

t13$splink = gsub(gsub(t13$species,
                  pattern = "[[:punct:]]",
                  replacement = ""),
                  pattern = "[[:blank:]]",
                  replacement = "")

i13$splink = gsub(gsub(t13$species,
                       pattern = "[[:punct:]]",
                       replacement = ""),
                  pattern = "[[:blank:]]",
                  replacement = "")


t17$splink = gsub(gsub(t17$species,
                       pattern = "[[:punct:]]",
                       replacement = ""),
                  pattern = "[[:blank:]]",
                  replacement = "")

i17$splink = gsub(gsub(i17$species,
                       pattern = "[[:punct:]]",
                       replacement = ""),
                  pattern = "[[:blank:]]",
                  replacement = "")



sp13 = unique(t13$splink)
sp17 = unique(t17$splink)

st13 = unique(t13$stratname)
st17 = unique(t17$stratname)

t13[which(t13$stratname == "Northeast US Coast"),"stratname"] <- "Northeast US Coastal"
t13[which(t13$stratname == "Southeast US Coast"),"stratname"] <- "Southeast Coastal"


i13[which(i13$stratname == "Northeast US Coast"),"stratname"] <- "Northeast US Coastal"
i13[which(i13$stratname == "Southeast US Coast"),"stratname"] <- "Southeast Coastal"

st13 = unique(t13$stratname)
st17 = unique(t17$stratname)


t17 = t17[which(t17$timespan == "1974-2016"),]


tmerg = merge(t17[,c("stratname",
                      "splink",
                      "species",
                      "trend.beta",
                      "trend.beta.lci",
                      "trend.beta.uci")],
              t13[,c("stratname",
                     "sp",
                      "splink",
                      "species",
                      "trend.lci",
                      "trend.",
                      "trend.uci")],
              by = c("stratname",
                     "splink"))

pdf(file = "trend comparisons 2013 vs 2017.pdf")
trange = range(c(tmerg$trend.beta,tmerg$trend.))

for(st in c("Continental",
            unique(tmerg$stratname)[-which(unique(tmerg$stratname) == "Continental")])){

  tmp = tmerg[which(tmerg$stratname == st),]
  if(st == "Continental"){
    plot(x = tmp$trend.,
         y = tmp$trend.beta,
         pch = 19,
         ylim = c(-10,5),
         xlim = c(-10,5),
         main = st,
         xlab = "2013 slope trend",
         ylab = "2017 slope trend",
         col = grey(0.4))
  }else{
  plot(x = tmp$trend.,
       y = tmp$trend.beta,
       pch = 19,
       ylim = trange,
       xlim = trange,
       main = st,
       xlab = "2013 slope trend",
       ylab = "2017 slope trend",
       col = grey(0.4))}
  abline(0,1,
         col = grey(0.5),
         lwd = 0.5)
  abline(h = 0,
         col = grey(0.5),
         lwd = 0.5)
  abline(v = 0,
         col = grey(0.5),
         lwd = 0.5)
  
  arrows(x0 = tmp$trend.lci,
         x1 = tmp$trend.uci,
         y0 = tmp$trend.beta,
         y1 = tmp$trend.beta,
         length = 0,
         col = transp.func("black",0.1))
  arrows(x0 = tmp$trend.,
         x1 = tmp$trend.,
         y0 = tmp$trend.beta.lci,
         y1 = tmp$trend.beta.uci,
         length = 0,
         col = transp.func("black",0.1))
  text(tmp$sp,
       y = tmp$trend.beta,
       x = tmp$trend.,
       cex = 0.5,
       pos = 1)
  
  
}#st
dev.off()












tmerg = merge(t17[,c("stratname",
                     "splink",
                     "species",
                     "trend.beta2",
                     "trend.beta2.lci",
                     "trend.beta2.uci")],
              t13[,c("stratname",
                     "sp",
                     "splink",
                     "species",
                     "trend.lci",
                     "trend.",
                     "trend.uci")],
              by = c("stratname",
                     "splink"))

pdf(file = "trend comparisons 2013 vs 2017 slope weighted by abundance.pdf")
trange = range(c(tmerg$trend.beta,tmerg$trend.))

for(st in c("Continental",
            unique(tmerg$stratname)[-which(unique(tmerg$stratname) == "Continental")])){
  
  tmp = tmerg[which(tmerg$stratname == st),]
  if(st == "Continental"){
    plot(x = tmp$trend.,
         y = tmp$trend.beta2,
         pch = 19,
         ylim = c(-10,5),
         xlim = c(-10,5),
         main = st,
         xlab = "2013 slope trend",
         ylab = "2017 slope trend",
         col = grey(0.4))
  }else{
    plot(x = tmp$trend.,
         y = tmp$trend.beta2,
         pch = 19,
         ylim = trange,
         xlim = trange,
         main = st,
         xlab = "2013 slope trend",
         ylab = "2017 slope trend",
         col = grey(0.4))}
  abline(0,1,
         col = grey(0.5),
         lwd = 0.5)
  abline(h = 0,
         col = grey(0.5),
         lwd = 0.5)
  abline(v = 0,
         col = grey(0.5),
         lwd = 0.5)
  
  arrows(x0 = tmp$trend.lci,
         x1 = tmp$trend.uci,
         y0 = tmp$trend.beta2,
         y1 = tmp$trend.beta2,
         length = 0,
         col = transp.func("black",0.1))
  arrows(x0 = tmp$trend.,
         x1 = tmp$trend.,
         y0 = tmp$trend.beta2.lci,
         y1 = tmp$trend.beta2.uci,
         length = 0,
         col = transp.func("black",0.1))
  text(tmp$sp,
       y = tmp$trend.beta2,
       x = tmp$trend.,
       cex = 0.5,
       pos = 1)
  
  
}#st
dev.off()






############## panel plot of the slope trends by regions (comparison of slopes among regions for a selected group of species with relatively precise trend estimates)
# 
# 
# jjc = 0
# for(jj in c(0,9,18,27)){
#   jjc = jjc+1
# 
# precsp <- contt[1:9+jj,"sp"]
# 
# regtt <- regtta[which(regtta$sp %in% precsp),]
# #regplot <- unique(regtt[,"stratname"])
# regtt[,"stratname"] <- factor(regtt[,"stratname"],
#                   levels = c("Atlantic Canada","Northeast US Coast","Southeast US Coast",
#                              "Texas Coastal","Ontario","East Inland","Midcontinental",
#                              "Pacific and Intermountain"),ordered = T)
# 
# 
# 
# pdf(width = 7,height = 7,file = paste0("regional trends panel ",jjc," ",GString$getBuiltinTime(format="%H-%M-%S"),".pdf"))
# par(mfrow = c(3,3),mar = c(1,1,1,1),oma = c(12,4,1,1))
# 
#   for(j in precsp){
# 
#   tmp <- regtt[which(regtt$sp == j),]
#   plot(y = c(min(regtt[,"trend lci"]),max(regtt[,"trend uci"])),
#        x = c(1,length(levels(tmp$stratname))),type = "n",
#        main = unique(regtt[which(regtt$sp == j),"species"]),xaxt = "n")
#   ### add in the continental trend estimate for the species
#   abline(h = contt[which(contt$sp == j),"trend "],col = gray(0.3),lty = 3)
#   if(j %in% precsp[c((length(precsp)-2):(length(precsp)))]) {
#     axis(side = 1,at = c(1:8),labels = levels(regtt$stratname),las = 3,cex = 1.2)
#   }
#   cj = 0
#   for(i in levels(tmp$stratname)){
#     cj <- cj+1
#   if(length(which(tmp$stratname == i))>0){
#   tmp2 <- tmp[which(tmp$stratname == i),]
#   points(y = tmp2[,"trend "],x = cj,pch = 20)
#   arrows(x0 = cj,x1 = cj,y0 = tmp2[,"trend lci"],y1 = tmp2[,"trend uci"],length = 0,col = gray(0.7))
#   }
#   }#i
# }#j
# dev.off()
#   
# }#jjc
# 
#     
#     
#     
#     
#     
#     
#   
#   
#   
#   
#   
#   
#   
#   
#   
#   
#   
#   
#   
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
    
    #################### SOCB style composite indicators for each region













# 
# 
# 
# indcout <- read.csv("all continental indices.csv",stringsAsFactors = F)
# spls <- read.csv("Copy of ShorebirdSotb2014Analysis_PAS.csv",stringsAsFactors = F)
# 
# 
# data <- indcout
# data[,"year"] <- 1973+data[,"year"]
# 
# 
# names(data)[c(6,8)] <- c("lci","uci")
# 
# 
# base.yr <- 1974
# data[,"index.s"] <- NA
# data[,"cv"] <- NA
# data[,"cvar.s"] <- NA
# data[,"se.s"] <- NA
# data[,"prec.logthetahat"] <- NA
# data[,"se.logthetahat"] <- NA
# data[,"logindex.plusse"] <- NA
# data[,"logthetahat"] <- NA
# 
# base.i <- rep(NA, length = length(unique(data$species)))
# names(base.i) <- unique(data$species)
# base.se.sp <- rep(NA, length = length(unique(data$species)))
# names(base.se.sp) <- unique(data$species)
# 
# mn.i <- vector(length = length(unique(data$species)))
# names(mn.i) <- unique(data$species)
# 
# 
# for (s in names(base.i)) {
# 
# ind <- "median"
# se <- "sd"
# r <- which(data[,"species"] == s)
# 
# base.s <- data[which(data$species == s & data$year == base.yr),ind] #stores the base index value - 1990
# base.se <- data[which(data$species == s & data$year == base.yr),se] #stores the base se value - 1990
# 
# mn.i[s]    <- mean(data[which(data$species == s),ind],na.rm = T)
#    for (y in r) {
# 
# data[y,"index.s"] <- (data[y,ind])/base.s
# data[y,"cv"] <- data[y,se]/data[y,ind]
# data[y,"cvar.s"] <- (((data[y,se]^2)/(data[y,ind]^2))+((base.se^2)/(base.s^2))) #this line sums the ratios of the variance over the index for year y and for the base year -
# ### the line above assumes that index estimates are independent among years, which is clearly false, however it is reasonable given that:
# ## there are no estimates of the covariance among years for indices from any of the data sources, and
# ## the scaling of variances in the final analysis are all relative not absolute and in the absence of more information, assuming that the covariance of annual indices is the same among different species seems reasonable.
# data[y,"se.s"] <- (sqrt(data[y,"cvar.s"])*data[y,"index.s"])
# data[y,"logthetahat"] <- log(data[y,"index.s"])
# 
# #### to convert the variance from the linear scale into the log scale, I've used the CV (as a porportion) to re-scale the SE estimates - column data$se.s
# #### then I have calculated the log of (the standardized index (shore$index.s) plus one standard error (shore$se.s))
# #### for many index values it might be better to calculate 1/4 of the distance between the logs of the upper and lower CI
# #### however for many of the se and index combinations this is not possible because the lower CI is negative
# data[y,"logindex.plusse"] <- log(data[y,"index.s"]+data[y,"se.s"])
# #### then calculated the difference between the logindex shore[y,"logthetahat"] and the log of the index plus one SE shore[y,"logthetahat.se"]
# data[y,"se.logthetahat"] <- data[y,"logindex.plusse"]-data[y,"logthetahat"]
# data[y,"prec.logthetahat"] <- 1/(data[y,"se.logthetahat"]^2)
# #data[y,"source"] <- d
# ############ taken from aerial insectivore changepoint analysis code ("R code to fit trend model with trend only multispecies full survey seprate FC and SSN.r")
# #y.lvar <- log(1+((((y.s$uci-y.s$lci)/4)^2)/((y.s$median)^2)))         # definition of the variance of a lognormal distribution, by the mean and variance of the non-logged distribution
# ##y.lse <- sqrt(y.lvar)
# #y.prec[,j] <- 1/y.lvar
#  data[,"test.prec"] <- 1/log(1+((data[,"se.s"]^2)/(data[,"index.s"]^2)))
# #
# #
# #
#   }
#   print(s)
# }
# 
# 
# ##################### read in list of long-distance migrants and arctic breeders, from Paul Smith, email dated June 2, 2014.
# 
# ldmig <- read.csv("LDmigrants.csv",stringsAsFactors = F)
# 
#  data.o <- data
#  spsob <- spls[which(spls$Species %in% ldmig[,1]),"code"]
#  spsob <- spsob[-which(spsob == "WESA")] # species selection for SOCB indicator calculation below
# 
#  data <- data[which(data$sp %in% spsob),]
# 
# 
#  sim.store <- list()
# 
# sp.list <- unique(data$sp)
# 
# 
# region = "continental"
# subgroup = "shorebirds"
# start.year = 1974
# end.year = 2013
# 
# 
# 
# 
# d.l <- list()
# length(d.l) <- (end.year-start.year)+1
# names(d.l) <- as.character(start.year:end.year)
#  indices <- data.frame(year = c(start.year:end.year),median = NA,sd = NA,x2.5 = NA,x25 = NA,x50 = NA,x75 = NA,x97.5 = NA)
#  indices.tr <- data.frame(year = c(start.year:end.year),median = NA,sd = NA,x2.5 = NA,x25 = NA,x50 = NA,x75 = NA,x97.5 = NA)
#  sim.out <- list()
# length(sim.out) <- (end.year-start.year)+1
# names(sim.out) <- as.character(start.year:end.year)
# sp.out <- list()
# length(sim.out) <- (end.year-start.year)+1
# names(sim.out) <- as.character(start.year:end.year)
# 
# 
# 
# 
# 
#  data.b <- data
#  library(R2WinBUGS)
# # library(BRugs)
#  for (y in start.year:end.year) {
# 
# 
#  lthat.y <- data.b[which(data.b$year == y & !is.na(data.b$logthetahat)),c("species","logthetahat","prec.logthetahat")]
# #d.y <- list(nspecies = nrow(lthat.y), logthetahat=round(lthat.y$logthetahat,4), prec.logthetahat = round(lthat.y$prec.logthetahat,4), fid = -0.693)
# d.y <- list(nspecies = nrow(lthat.y), logthetahat=round(lthat.y$logthetahat,4), prec.logthetahat = round(lthat.y$prec.logthetahat,4))
#   param <- c("expmu","logtheta")#,"gt.fifty","lt.fifty","pos")#,,"mu.logtheta","tau","theta","meantr","numpos","ranking","stable")
# #model.f <- "C:/Adam/SOTB/model_logtheta.txt" #Jon's computer
# model.f <- "model_logtheta.txt" #laptop and desktop
#   initial <- list(list(prec.logtheta = 1, mu.logtheta = 0, logtheta = rep(0,nrow(lthat.y))))
# 
#   sim2 <- bugs(d.y,initial, param, model.f, n.chains = 1, n.iter = 50000,n.burnin = 20000,n.thin = 1, bugs.directory="C:\\WinBUGS14",program="WinBUGS", debug = F)    #my desktop
# # sim2 <- bugs(d.y,initial, param, model.f, n.chains = 1, n.iter = 10000,n.burnin = 2000,n.thin = 1,bugs.directory="C:\\Program Files\\OpenBUGS\\OpenBUGS312",program="OpenBUGS", debug = T)#,working.directory=paste("EastBoreal\\Bugs\\",y,"\\", sep = ""),)
#  indices[indices$year == y,2:8] <- sim2$summary["expmu",]
#  indices[indices$year == y,"median"] <- sim2$median$expmu
# 
# # indices.tr[indices.tr$year == y,2:8] <- sim2$summary["meantr",]
# #  indices.tr[indices.tr$year == y,"median"] <- sim2$median$meantr
# # indices.th[indices.th$year == y,2:8] <- sim2$summary["theta",]
# yr <- as.character(y)
#  d.l[[yr]] <- d.y
#  sim.out[[yr]] <- sim2$sims.list
#  sp.out[[yr]] <- lthat.y$species
# 
# 
# img <- paste("bugsout/",y,region,".RData",sep = "")
# 
# #
# if (y == end.year) {save.image(img) }
# #     save(list = ls(envir = environment(), all.names = TRUE),
# #     file = img, envir = environment())
# rmve <- ls()
# rmve <- rmve[-which(rmve %in% c("indices","sim.out","sp.out","indices.tr","d.l","data","y","data.b","start.year","end.year","region","subgroup","bugs.run","sp.list"))]
# #rmve <- rmve[-which(rmve %in% c("indices","indices.tr","d.l","data","y","data.b","start.year","end.year","region","subgroup","bugs.run","sp.list"))]
# rm(list = rmve)
#  print(y)
# }
# 
# 
# yup <- max(indices[,c("median","x2.5","x97.5")])
# ydn <- min(indices[,c("median","x2.5","x97.5")])
# 
# plot(indices$year,indices$median,ylim = c(ydn,yup),ylab = "index wrt 1974 (1974 = 1.0)",type = "l",lwd = 3,xaxt = "n",bty = "l")
# lines(indices$year,indices[,"x2.5"],col = grey(0.5),lty = 3,lwd = 2)
# lines(indices$year,indices[,"x97.5"],col = grey(0.5),lty = 3,lwd = 2)
# abline(h = 1,col = grey(0.5))
# lines(indices$year,indices[,"median"],lwd = 3)
# axis(side = 1,at = c(1974,seq(1980,2010,by = 5),2013),labels = c(1974,seq(1980,2010,by = 5),2013))
# 
# write.csv(indices,"State of Birds indicator LDmigrant shorebirds.csv")
# 
# 
# #source("index plotting function.r")
# #index.plot(index = national.nat.all.compindex, region = "national",subgroup = "nat.all", error.plot = T)
# ##index.plot(index = national.nat.all.compindex, region = "national",subgroup = "nat.all")
# #save.image("national\\national.nat.all2.postbugs.RData")
# #


