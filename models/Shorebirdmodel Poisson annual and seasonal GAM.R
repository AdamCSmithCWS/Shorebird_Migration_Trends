## WinBUGS model used to estimate annual indices of abundance and population trends from shorebird migration data
## alternative model from Smith and Smith designed for use in State of Canada's Birds 2019
#	
model
{
for(i in 1:ncounts){

	loglambda[i] <- alpha + ste[site[i]] + ann[year[i],strat[i]] + zi[date[i],strat[i]] + noise[i]  ### common intercept, varying slopes, so that the site effect accounts for all of the variation in abundance.
	
	noise[i] ~ dt(0,taunoise,nu)
	log(lambda[i]) <- loglambda[i]
	
	count[i] ~ dpois(lambda[i])


}#i



	sdnoise <- 1/pow(taunoise,0.5)
	taunoise ~ dscaled.gamma(0.5,50)
  nu ~ dgamma(2,0.2)

alpha~dnorm(0,1)







###########COMPUTING GAMs for seasonal patterns (fixed GAMs by strata)
# Addapted from Crainiceanu, C. M., Ruppert, D. & Wand, M. P. (2005). Bayesian Analysis for Penalized Spline Regression Using WinBUGS. Journal of Statistical Softare, 14 (14), 1-24.
# X.basis is data computed in R

#tauX~dgamma(1.0E-3,1.0E-7) #alternate prior, original from Cainiceanu et al.
#tauX <- 1/pow(sdX,2) # prior on precision of gam hyperparameters
#sdX <- 1/(pow(tauX,0.5)) # ~ dunif(0,5)
taubetagam <- 1/pow(sdbetagam,2) # prior on precision of gam coefficients
sdbetagam ~ dunif(0,5)

   for(j in 1:nknots){ # Computation of GAM components
         #B.X[j] ~ dnorm(0,tauX)

                
	for(k in 1:nstrata){
		beta.X[k,j] ~ dnorm(0,taubetagam )
		

         for ( i in 1: dmax[k] )
         {
             X.part[i,j,k] <- beta.X[k,j]*(X.basis[i,j,k])
           
         }#i

}#k
    }#j

	for(k in 1:nstrata){
    for (i in 1: dmax[k])
    {

        zi[i,k] <- sum(X.part[i,1:nknots,k])
    }#k
    }#i


       




###########COMPUTING GAMs for year effects (hierarchical GAMs by strata)
# Addapted from Crainiceanu, C. M., Ruppert, D. & Wand, M. P. (2005). Bayesian Analysis for Penalized Spline Regression Using WinBUGS. Journal of Statistical Softare, 14 (14), 1-24.
# Y.basis is data computed in R


#tauY~dgamma(1.0E-2,1.0E-4) 
### alternate prior, original from Cainiceanu et al. 
### second gamma parameter == 0.0001 << (abs(mean(B.Y[]))^2)/2, mean(B.Y[]) ~ 0.2

tauY <- 1/pow(sdY,2) # prior on precision of gam hyperparameters
sdY ~ dunif(0,5) #<- 1/(pow(tauY,0.5)) # #alternate prior
taubetagamY <- 1/pow(sdbetagamY,2) # prior on precision of gam coefficients
sdbetagamY ~ dunif(0,5)

   for(j in 1:nknotsy){ # Computation of GAM components
         B.Y[j] ~ dnorm(0,tauY)

                
	for(k in 1:nstrata){
		beta.Y[k,j] ~ dnorm(B.Y[j],taubetagamY )
		

         for ( i in 1: nyears)
         {
             Y.part[i,j,k] <- beta.Y[k,j]*(Y.basis[i,j])
           
         }#i

}#k
    }#j

	for(k in 1:nstrata){
    for (i in 1: nyears)
    {

        ann[i,k] <- sum(Y.part[i,1:nknotsy,k])
    }#k
    }#i


       







######################
## site effects

for (j in 1:nsites){
	ste[j]~dnorm(0, tausite)
	}#s
sdsite ~ dunif(0,20)
tausite <- 1/pow(sdsite,2)


######################


}#end model
