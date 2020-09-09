## JAGS model used to estimate annual indices of abundance and population trends from shorebird migration data
## alternative model from Smith and Smith designed for use in State of Canada's Birds 2019
#	
model
{
for(i in 1:ncounts){

	loglambda[i] <- alpha[strat[i]] + ste[site[i],strat[i]] + sm_year[yr[i],strat[i]] + gam.sm_season[i] + noise[i]  ### common intercept, varying slopes, so that the site effect accounts for all of the variation in abundance.
	
	#noise[i] ~ dt(0,taunoise,nu)
	noise[i] ~ dnorm(0,taunoise)
	log(lambda[i]) <- loglambda[i]
	
	count[i] ~ dpois(lambda[i])


}#i



	# sdnoise <- 1/pow(taunoise,0.5)
	# taunoise ~ dscaled.gamma(0.5,50)
  taunoise <- 1/pow(sdnoise,2)
  sdnoise ~ dt(0, 1, 4)T(0,) # half-t prior on sd (chung et al. 2013) DOI: 10.1007/S11336-013-9328-2 places ~95% of the prior < 3.0
	
  #nu ~ dgamma(2,0.2)







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



###########COMPUTING GAMs for year effects - Hyperparameters


taugam_year <- 1/pow(sdgam_year,2)
sdgam_year ~ dt(0, 1, 4)T(0,) # half-t prior on sd (chung et al. 2013) DOI: 10.1007/S11336-013-9328-2 places ~95% of the prior < 3.0

# sdgam_year <- 1/pow(taugam_year,0.5)
# taugam_year ~ dscaled.gamma(0.5,50) #complexity penalty
### random effects
# sdgam_year_b <- 1/pow(taugam_year_b,0.5)
# taugam_year_b ~ dscaled.gamma(0.1,50) #shrinkage to hyperparameter

# taugam_year_b <- 1/pow(sdgam_year_b,2)
# sdgam_year_b ~ dt(0, 0.5, 10)T(0,) # half-t prior on sd (chung et al. 2013) DOI: 10.1007/S11336-013-9328-2 places ~95% of the prior < 1.0


for(k in 1:nknots_year){
  taugam_year_b[k] <- 1/pow(sdgam_year_b[k],2)
  sdgam_year_b[k] ~ dt(0, 0.5, 10)T(0,) # half-t prior on sd (chung et al. 2013) DOI: 10.1007/S11336-013-9328-2 places ~95% of the prior < 1.0
  
  B_year[k] ~ dnorm(0,taugam_year)
  for(s in 1:nstrata){
    
    b_year[k,s] ~ dnorm(0,taugam_year_b[k])
    beta_year[k,s] <- B_year[k]
    #beta_year[k,s] <- b_year[k,s]+B_year[k]
    }
}

#strata smooths
for(s in 1:nstrata){
  
  sm_year[1:nyears,s] <-	year_basispred %*% beta_year[1:nknots_year,s]
  
}#s

##### hyperparameter smooth

 
sm_year_B[1:nyears] <-	year_basispred %*% B_year[1:nknots_year]




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
retrans <- 0.5*(1/taunoise)
# retrans <- 0.5*(1/taunoise)/nu_ret 
# nu_ret <- (1.422*nu^0.906)/(1+(1.422*nu^0.906)) #approximate retransformation to equate a t-distribution to a normal distribution - see appendix of Link et al. 2020 BBS model selection paper



for(y in 1:nyears){
  for(s in 1:nstrata){
    for(j in 1:nsites[s]){
      n_sj[j,s,y] <- exp(alpha[s] + ste[j,s] + sm_year[y,s] + vis.sm_season[76] + retrans) #site-level predictions including strata-level yearly smooths
    }#j
    n_s[s,y] <- mean(n_sj[1:nsites[s],s,y])#stratum predictions including strata-level yearly smooths and scaled to mean across stratum sites
    n_s_scaled[s,y] <- exp(alpha[s] + sm_year[y,s] + vis.sm_season[76] + retrans + retrans_j) #stratum predictions including strata-level yearly smooths and on a common scale (visualisation only)
  }#s
  N[y] <- exp( sm_year_B[y] + vis.sm_season[76] + retrans + retrans_j) #continental predictions including only the hyperparameter smooth
}

}#end model
