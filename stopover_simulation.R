 nbirds = 1000 #nbirds
 ncounts = 10 #ncounts/year
 ndaysSepCounts = 7 #5 days between subsequent counts
 countdays = seq(5,5+(ncounts-1)*ndaysSepCounts,by = ndaysSepCounts)
 nyears = 40 #
 deltastopy = 0.98 #yearl proportional change
 deltastop = deltastopy^40 #total change in stopover duration
 stop1 = 10 #starting stopover duration
 stopy = c(stop1,stop1*deltastopy^(1:(nyears-1)))
 sdstopy = stopy*0.2
 
 ndseason = 80
 meanArrival = 30
 arrivalDays = rnorm(nbirds,meanArrival,5)
 if(any(arrivalDays < 1)){
   arrivalDays[which(arrivalDays < 1)] <- 1
 }
 beta = 1 #true multiplicative annual change
 nbirdsy = c(nbirds,nbirds*beta^(1:(nyears-1)))
 
 
 birdsPresent = matrix(0,nrow = nyears,ncol = ndseason)
 counts = matrix(0,nrow = nyears,ncol = length(countdays))
 for(y in 1:nyears){
   birdocc = matrix(0,nrow = nbirdsy[y],ncol = ndseason)
   d1s = sample(arrivalDays,nbirdsy[y])
   for(b in 1:nbirdsy[y]){
     birdocc[b,c(ceiling(d1s[b]):ceiling(d1s[b]+rnorm(1,stopy[y],sdstopy[y])))] <- 1
     
   }
   
   birdsPresent[y,] <- colSums(birdocc)
   
   counts[y,] <- birdsPresent[y,countdays]
 }
 
 
 
 