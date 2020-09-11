## JAGS model used to estimate annual indices of abundance and population trends from shorebird migration data
## alternative model from Smith and Smith designed for use in State of Canada's Birds 2019
#	
model
{
for(i in 1:ncounts){

	loglambda[i] <- alpha[strat[i]] + ste[site[i],strat[i]] + beta[strat[i]]*(yr[i]-midyear) + year_effect[yr[i],strat[i]] + gam.sm_season[i] + noise[i]  ### common intercept, varying slopes, so that the site effect accounts for all of the variation in abundance.
	
	noise[i] ~ dt(0,taunoise,nu)
	#noise[i] ~ dnorm(0,taunoise)
	log(lambda[i]) <- loglambda[i]
	
	count[i] ~ dpois(lambda[i])


}#i



	# sdnoise <- 1/pow(taunoise,0.5)
	# taunoise ~ dscaled.gamma(0.5,50)
  taunoise <- 1/pow(sdnoise,2)
  sdnoise ~ dt(0, 1, 4)T(0,) # half-t prior on sd (chung et al. 2013) DOI: 10.1007/S11336-013-9328-2 places ~95% of the prior < 3.0
	
  nu ~ dgamma(2,0.2)







###########COMPUTING GAMs for seasonal patterns hierarchical smooth across strata
# Addapted from Crainiceanu, C. M., Ruppert, D. & Wand, M. P. (2005). Bayesian Analysis for Penalized Spline Regression Using WinBUGS. Journal of Statistical Softare, 14 (14), 1-24.

taugam_season <- 1/pow(sdgam_season,2)
sdgam_season ~ dt(0, 1, 10)T(0,) # half-t prior on sd (chung et al. 2013) DOI: 10.1007/S11336-013-9328-2 places ~95% of the prior < 1.0


for(k in 1:nknots_season){
  beta_season[k] ~ dnorm(0,taugam_season)
}

gam.sm_season <-	season_basis %*% beta_season


##### derived parameters to visualize the smooth

vis.sm_season <-	season_basispred %*% beta_season

mn_sm_season <- mean(vis.sm_season)



###########COMPUTING slopes


tau_beta <- 1/pow(sd_beta,2)
sd_beta ~ dt(0, 1, 4)T(0,) # half-t prior on sd (chung et al. 2013) DOI: 10.1007/S11336-013-9328-2 places ~95% of the prior < 3.0

B ~ dnorm(0,100) #regularizing prior putting 95% of the prior with absolute value of trends < 20%/year
  for(s in 1:nstrata){
    
    b[s] ~ dnorm(0,tau_beta)
    beta[s] <- B + b[s]
    }


##############Year Effects
tau_year <- 1/pow(sd_year,2) #variance in continental year effects, controls shrinkage towards the line
sd_year ~ dt(0, 1, 10)T(0,) # half-t prior on sd (chung et al. 2013) DOI: 10.1007/S11336-013-9328-2 places ~95% of the prior < 2

for(s in 1:nstrata){
  tau_ye[s] <- 1/pow(sd_ye[s],2) #variance in strata-level departures from the continental year-effects - further shrinkage
  sd_ye[s] ~ dt(0, 1, 10)T(0,) # half-t prior on sd (chung et al. 2013) DOI: 10.1007/S11336-013-9328-2 places ~95% of the prior < 2
  
}
for(y in 1:nyears){
  YE[y] ~ dnorm(0,tau_year)
  for(s in 1:nstrata){
  ye[y,s] ~ dnorm(0,tau_ye[s])
    year_effect[y,s] <- YE[y]+ye[y,s]
      }
}




######################
## site effects
for(s in 1:nstrata){
  alpha[s] ~ dnorm(0,1)
for (j in 1:nsites[s]){
	ste[j,s]~dnorm(0,tausite[s])
	}#j

# sdsite[s] <- 1/pow(tausite[s],0.5)
# tausite[s] ~ dscaled.gamma(2,10) #complexity penalty

  tausite[s] <- 1/pow(sdsite[s],2)
sdsite[s] ~ dt(0, 1, 4)T(0,) # half-t prior on sd (chung et al. 2013) DOI: 10.1007/S11336-013-9328-2 places ~95% of the prior < 3.0


retrans_js[s] <- 0.5*(1/tausite[s])
}
retrans_j <- mean(retrans_js[1:nstrata]) 

######################
##  derived parameters
#retrans <- 0.5*(1/taunoise)
retrans <- 0.5*(1/taunoise)/nu_ret
nu_ret <- (1.422*nu^0.906)/(1+(1.422*nu^0.906)) #approximate retransformation to equate a t-distribution to a normal distribution - see appendix of Link et al. 2020 BBS model selection paper



for(y in 1:nyears){
  for(s in 1:nstrata){
    for(j in 1:nsites[s]){
      n_sj[j,s,y] <- exp(alpha[s] + ste[j,s] + beta[s]*(y-midyear) + year_effect[y,s] + mn_sm_season + retrans) #site-level predictions including strata-level yearly smooths
    }#j
    n_s[s,y] <- mean(n_sj[1:nsites[s],s,y])#stratum predictions including strata-level yearly smooths and scaled to mean across stratum sites
    n_s_scaled[s,y] <- exp(alpha[s] + beta[s]*(y-midyear) + year_effect[y,s] + mn_sm_season + retrans + retrans_j) #stratum predictions including strata-level yearly smooths and on a common scale (visualisation only)
  }#s
  N[y] <- exp(B*(y-midyear) + YE[y] + mn_sm_season + retrans + retrans_j) #continental predictions including only the hyperparameter smooth
}

}#end model
