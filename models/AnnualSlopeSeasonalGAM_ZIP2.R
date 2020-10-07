## JAGS model used to estimate annual indices of abundance and population trends from shorebird migration data
## alternative model from Smith and Smith designed for use in State of Canada's Birds 2019
#	
model
{
for(i in 1:ncounts){

	elambda1[i] <- alpha[strat[i]] + ste[site[i],strat[i]] + beta[strat[i]]*(yr[i]-midyear) + year_effect[yr[i],strat[i]] + gam.sm_season[i] + noise[i]  ### common intercept, varying slopes, so that the site effect accounts for all of the variation in abundance.
	
	noise[i] ~ dt(0,taunoise[strat[i]],nu[strat[i]])
	#noise[i] ~ dnorm(0,taunoise)
	log(lambda1[i]) <- elambda1[i]
	
	
	lambda[i] <- lambda1[i]*z[i] + 0.00001 ## hack required for JAGS -- otherwise 'incompatible'-error
	z[i] ~ dbern(psi[site[i],strat[i]])#psi[year[i]]) #psi = proportion of non-zeros for each year
	### how tomodel the zip component...at what scale do the zeros vary? sites, strata, years? probably sites given the crazy count distributions at some sites, but maybe just overall, given the fact that many sites have extreme variations
	
	
	count[i] ~ dpois(lambda[i])


}#i

 
  
  for(s in 1:nstrata){
   
      # sdnoise <- 1/pow(taunoise,0.5)
      # taunoise ~ dscaled.gamma(0.5,50)
      taunoise[s] <- 1/pow(sdnoise[s],2)
      sdnoise[s] ~ dt(0, 0.5, 4)T(0,) # half-t prior on sd (chung et al. 2013) DOI: 10.1007/S11336-013-9328-2 places ~99% of the prior < 2.25
      # this is a relatively informative prior to avoid large estimates of sdnoise that have a strong influence on teh scaling of the stratum level estimates
      nu[s] ~ dgamma(2,0.5) #puts the mean of the prior at ~ 4 and ~99% of the prior < 13
      
    }
    
    
    




  
  ###########COMPUTING GAMs for seasonal patterns hierarchical smooth across strata
  # Addapted from Crainiceanu, C. M., Ruppert, D. & Wand, M. P. (2005). Bayesian Analysis for Penalized Spline Regression Using WinBUGS. Journal of Statistical Softare, 14 (14), 1-24.
  
  # taugam_season <- 1/pow(sdgam_season,2)
  # sdgam_season ~ dt(0, 1, 10)T(0,) # half-t prior on sd (chung et al. 2013) DOI: 10.1007/S11336-013-9328-2 places ~95% of the prior < 1.0
  # 
  # 
  # for(k in 1:nknots_season){
  #   B_season[k] ~ dnorm(0,taugam_season)
  #   
  #   taugam_season_s[k] <- 1/pow(sdgam_season_s[k],2)
  #   sdgam_season_s[k] ~ dt(0, 1, 10)T(0,) # half-t prior on sd (chung et al. 2013) DOI: 10.1007/S11336-013-9328-2 places ~95% of the prior < 1.0
  #   
  #   for(s in 1:nstrata){
  #     
  #     b_season[k,s] ~ dnorm(0,taugam_season_s[k])
  #     
  #     beta_season[k,s] <- b_season[k,s]+B_season[k]
  #   }
  # }
  # 
  # 
  # 
  # #strata smooths
  # for(s in 1:nstrata){
  #   
  #   gam.sm_season[1:ncounts,s] <-	season_basis %*% beta_season[1:nknots_season,s]
  #   ##### derived parameters to visualize the smooth
  #   
  #   vis.sm_season[1:ndays,s] <-	season_basispred %*% beta_season[1:nknots_season,s]
  #   mn_sm_season[s] <- mean(vis.sm_season[,s])
  #   
  # }#s
  
#   
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



PSI ~ dbeta(1,1) #mean site-level zip parameter
ZIP <- logit(PSI)
tau_zip_js <- 1/pow(sd_zip_js,2)
sd_zip_js ~ dt(0, 0.2, 4)T(0,) #half-t prior on sd (chung et al. 2013) DOI: 10.1007/S11336-013-9328-2 places ~95% of the prior < 1

  for(s in 1:nstrata){
    alpha[s] ~ dnorm(0,1)

  for (j in 1:nsites[s]){
    ste[j,s]~dnorm(0,tausite)
    
    zip_js[j,s] ~ dnorm(0,tau_zip_js)
    zip[j,s] <- zip_js[j,s] + ZIP
    logit(psi[j,s]) <- zip[j,s]
  }#j
}#s

tausite <- 1/pow(sdsite,2)
sdsite ~ dt(0, 1, 4)T(0,) # half-t prior on sd (chung et al. 2013) DOI: 10.1007/S11336-013-9328-2 places ~95% of the prior < 3
#nu_site ~ dgamma(2,0.2) #

retrans_j <- (0.5*(1/tausite))

#nu_ret_j <- (1.422*nu_site^0.906)/(1+(1.422*nu_site^0.906)) #approximate retransformation to equate a t-distribution to a normal distribution - see appendix of Link et al. 2020 BBS model selection paper

# 
# 
# ######################   Consider changing this to a single parameter across all strata, and using the informative priors on sdnoise
# ## site effects
# for(s in 1:nstrata){
#   alpha[s] ~ dnorm(0,1)
#   for (j in 1:nsites[s]){
#     ste[j,s]~dt(0,tausite[s],nu_site[s])
#   }#j
#   
#   # sdsite[s] <- 1/pow(tausite[s],0.5)
#   # tausite[s] ~ dscaled.gamma(2,10) #complexity penalty
#   
#   tausite[s] <- 1/pow(sdsite[s],2)
#   sdsite[s] ~ dt(0, 0.3, 4)T(0,) # half-t prior on sd (chung et al. 2013) DOI: 10.1007/S11336-013-9328-2 places ~99% of the prior < 1.5
#   # this is a relatively informative prior to avoid large estimates of sdsite that have a strong influence on teh scaling of the stratum level estimates
#   nu_site[s] ~ dgamma(2,0.5) #puts the mean of the prior at ~ 4 and ~99% of the prior < 13
#   
#   retrans_js[s] <- (0.5*(1/tausite[s]))/nu_ret_j[s]
#   
#   nu_ret_j[s] <- (1.422*nu_site[s]^0.906)/(1+(1.422*nu_site[s]^0.906)) #approximate retransformation to equate a t-distribution to a normal distribution - see appendix of Link et al. 2020 BBS model selection paper
#   
# }
# retrans_j <- mean(retrans_js[1:nstrata]) 

######################
##  derived parameters
for(s in 1:nstrata){
  retrans_a1[s] <- 0.5*(1/taunoise[s])
  retrans[s] <- 0.5*(1/taunoise[s])/nu_ret[s]
  nu_ret[s] <- (1.422*nu[s]^0.906)/(1+(1.422*nu[s]^0.906)) #approximate retransformation to equate a t-distribution to a normal distribution - see appendix of Link et al. 2020 BBS model selection paper
}
retrans_m <- mean(retrans[1:nstrata])
alpha_m <- mean(alpha[1:nstrata])


for(y in 1:nyears){
  for(s in 1:nstrata){
    for(j in 1:nsites[s]){
      n_sj[j,s,y] <- exp(alpha[s] + ste[j,s] + beta[s]*(y-midyear) + year_effect[y,s] + mn_sm_season + retrans[s])* psi[j,s] #site-level predictions including strata-level yearly smooths
      # n_sj_a1[j,s,y] <- exp(alpha[s] + ste[j,s] + beta[s]*(y-midyear) + year_effect[y,s] + mn_sm_season[s] + retrans_a1[s])* psi #site-level predictions including only the normal noise component (ignoring the t-distributed error)
      # n_sj_a2[j,s,y] <- exp(alpha[s] + ste[j,s] + beta[s]*(y-midyear) + year_effect[y,s] + mn_sm_season[s])* psi #site-level predictions excluding the log-normal retransformation completely
    }#j
    n_s[s,y] <- mean(n_sj[1:nsites[s],s,y])#stratum predictions including strata-level yearly smooths and scaled to mean across stratum sites
    # n_s_a1[s,y] <- mean(n_sj_a1[1:nsites[s],s,y])#stratum predictions including strata-level yearly smooths and scaled to mean across stratum sites
    # n_s_a2[s,y] <- mean(n_sj_a2[1:nsites[s],s,y])#stratum predictions including strata-level yearly smooths and scaled to mean across stratum sites
    # 
     # n_s_scaled[s,y] <- exp(alpha[s] + beta[s]*(y-midyear) + year_effect[y,s] + mn_sm_season[s] + retrans[s] + retrans_j)* psi #stratum predictions including strata-level yearly smooths and on a common scale (visualisation only)
     # n_s_scaled2[s,y] <- exp(alpha[s] + beta[s]*(y-midyear) + year_effect[y,s] + mn_sm_season[s] + retrans[s] + retrans_j)* psi #stratum predictions including strata-level yearly smooths and on a common scale (visualisation only)
     # 
      }#s
  N[y] <- exp(alpha_m + B*(y-midyear) + YE[y] + mn_sm_season + retrans_m + retrans_j)*PSI #continental predictions including only the hyperparameter smooth
  N_comp[y] <- mean(n_s[1:nstrata,y]) #continental predictions 
  # N_sc[y] <- mean(n_s_scaled[1:nstrata,y]) #coh
  # N_sc2[y] <- mean(n_s_scaled2[1:nstrata,y]) #c
  # 
  }

}#end model
