### a function to calculate the basis function(s) for a GAM component of a JAGS or Stan model
### orig.preds = predictor variable to be modeled as a smooth, column from the dataset being modeled, must match with the indexing used in teh main likelihood definition (i.e., one value for each statistical unit)
### nknots = number of knots in the GAM smooth, (careful! you need to know what you're doing :-) 
### random = logical; TRUE if the GAM parameters are to be estimated as random effects (e.g., separate smooths within strata, centered on a mean smooth across all strata)
### npredpoints = number of points at which to teh user wants to visualize predictions for the smooth
### standardize = text; "z" if the predictor should be standardized (centered and re-scaled) to a mean = 0 and SD = 1, "range" if the user prefers a standard value wit hmean = mid-point of the range, and SD ~ 1 depending on the shape of the distribution (this approach can be helpful if predictor values are strongly skewed)
### even_gaps = T; logical, if TRUE the teh knots are spaced evenly along the range fot he predictor, if FALSE the data quantiles are used
### sm_name = text, text string to uniquely identify the predictor variable and its related parameters in the smooth


gam.basis.func <- function(orig.preds = dts[,"yr"],
                           nknots = 6,
                           standardize = "z",
                           random = F,
                           npredpoints = 100,
                           even_gaps = FALSE,
                           sm_name = ""){
  
  
  
  
  if(any(is.na(orig.preds) == T)){
    stop("This GAM formulation cannot handle missing values in the predictor")
  }
  

  vmin = min(orig.preds) 
  vmax = max(orig.preds)
  
  
  n = length(orig.preds)
  
  predpoints = seq(vmin,vmax,length.out = npredpoints)
  
  if(standardize == "range"){ 
    recenter = floor(diff(c(vmin,vmax))/2) #centering to the middle of the range
    #ignores the distribution of the data points, which may be dangerous if the data are not well distributed
    rescale = ((vmax-vmin)+1)/4 # rescales to a unit value = 1/4 of the range, similar to SD = 1 for BBS data
    vscale = (orig.preds-recenter)/rescale
  }
  if(standardize == "z"){
    vscale = (orig.preds-mean(orig.preds))/sd(orig.preds)
  }
  
  
  
  
  yminsc = min(vscale,na.rm = T)
  ymaxsc = max(vscale,na.rm = T)
  
  ############ GAM basis function
  
  if(even_gaps){
    knotsgamx<- seq(yminsc,ymaxsc,length=(nknots+2))[-c(1,nknots+2)]
  }else{
    qqs = seq(0,1,length = (nknots+2))[-c(1,nknots+2)]
    knotsgamx<- quantile(vscale,qqs,names = F)
  }
  
  gamx_OMEGA_all<-(abs(outer(knotsgamx,knotsgamx,"-")))^3
  gamx_svd.OMEGA_all<-svd(gamx_OMEGA_all)
  gamx_sqrt.OMEGA_all<-t(gamx_svd.OMEGA_all$v  %*% (t(gamx_svd.OMEGA_all$u)*sqrt(gamx_svd.OMEGA_all$d)))
  
  
  gamx_K<-(abs(outer(seq(yminsc,ymaxsc,length = npredpoints),knotsgamx,"-")))^3
  gamx.basispred<-t(solve(gamx_sqrt.OMEGA_all,t(gamx_K)))
  
  
  gamx_K1<-(abs(outer(vscale,knotsgamx,"-")))^3
  gamx.basis<-t(solve(gamx_sqrt.OMEGA_all,t(gamx_K1)))
  
  if(random){
    # mod.code <- paste0("
    # ###########insert the following text into an existing JAGS model##############
    # ###########making sure to link the relevant parameters into section##############
    # ###########of the model that estimates the full likelihood ##############
    # ###################################################################
    # # the GAM smooth is predicted as intercepts for each whole value of the
    # # predictor (e.g., each year in a temporal smooth)
    # # so the likelihood line will need to include a component such as the following for each data point-i
    # # y[i] <- ..... + gam.sm[i] +....
    # 
    # for ( i in 1:n ){
    # for(k in 1:nknots){
    # gamx.part[i,k] <- beta.gamx[group[i],k]*(gamx.basis[i,k])
    # }
    # gam.sm[i] <- sum(gamx.part[i,1:nknots])
    # }#i
    # 
    # 
    # ###################################################################
    # ##### this text assumes the GAM parameters (betas) are estimated as random effects########
    # ##### your will need to provide the variable ngroups (an integer indicating the number of groups in the random effect)####
    # ##### and ensure that your R-script retains the correct order of groups in relation to the predicted values indexed by the 1:ngroups section below ######
    # 
    # # Following Crainiceanu, C. M., Ruppert, D. & Wand, M. P. (2005). Bayesian Analysis for Penalized Spline Regression Using WinBUGS. Journal of Statistical Softare, 14 (14), 1-24.
    # 
    # taugamx~dgamma(1.0E-2,1.0E-4) #alternate prior, original from Cainiceanu et al. second gamma parameter == 0.0001 << (abs(mean(B.gamx[]))^2)/2, mean(B.gamx[]) ~ 0.2
    # sdgamx <- 1/(pow(taugamx,0.5)) # 
    # taubeta <- 1/pow(sdbeta,2) # prior on precision of gam coefficients(
    # sdbeta ~ dunif(0,5)
    # 
    # for(k in 1:nknots)
    # { # Computation of GAM components
    # B.gamx[k] ~ dnorm(0,taugamx)
    # 
    # for(j in 1:ngroups)
    # {
    # beta.gamx[j,k] ~ dnorm(B.gamx[k],taubeta)
    # 
    # 
    # for ( i in 1:npredpoints ){ #this allows the user to model just the GAM component at particular values of the predictor
    # gamx.partpred[i,k,j] <- beta.gamx[j,k]*(gamx.basispred[i,k])
    # 
    # }#i
    # 
    # }#k
    # }#j
    # 
    # for (i in 1:npredpoints){
    # for(j in 1:ngroups)
    # {
    # gam.smooth[i,j] <- sum(gamx.part[i,1:nknots,j])
    # }#k
    # }#i
    # ")
  }else{
    mod.code = paste0("
                      ###########insert the following text into an existing JAGS model##############
                      ###########making sure to link the relevant parameters into section##############
                      ###########of the model that estimates the full likelihood ##############
                      ###################################################################
                      ###################################################################
                      # the GAM smooth is predicted as intercepts for each whole value of the
                      # predictor (e.g., each year in a temporal smooth)
                      # so the likelihood line will need to include a component such as the following for each data point-i
                      # for (i in 1:n){
                      # y[i] <- ..... + gam.sm_",sm_name,"[i] +....
                      #}
                      
                      ######
                      #the effect of this prior is worth exploring
                      #alternatives such as uniform priors on the SD or half-Cauchy may perform better
                      taugam_",sm_name,"~dgamma(1.0E-2,1.0E-4) #original from Cainiceanu et al. second gamma parameter == 0.0001 << (abs(mean(B.gamx[]))^2)/2, mean(B.gamx[]) ~ 0.2
                      sdgam_",sm_name," <- 1/(pow(taugam_",sm_name,",0.5)) # 
                      #or
                      sdgam_",sm_name," <- 1/pow(taunoise,0.5)
                      taugam_",sm_name," ~ dscaled.gamma(0.5,50)
                      # Computation of GAM components
                      for(k in 1:nknots_",sm_name,"){
                      beta_",sm_name,"[k] ~ dnorm(0,taugam_",sm_name,")
                      }
                      
                      gam.sm_",sm_name," <-	",sm_name,"_basis %*% beta_",sm_name,"
                      
                      
                      ##### derived parameters to visualize the smooth
                      
                      vis.sm_",sm_name," <-	",sm_name,"_basispred %*% beta_",sm_name,"
                      
                      ##### suggested components to add to annual index calculation
                      #### value of vis.sm_",sm_name," below is just the middle of the values used to visualize the smooth, i.e., predpoints[round(0.5*npredpoints)]
                      #### if there is a better value then insert the index for the relevant value of predpoints
                      # for(y in 1:nyears){
                      # for(s in 1:nstrata){
                      # n[i,t] <- exp(strata[i] + beta[i]*(t-fixedyear) + yeareffect[i,t] + vis.sm_",sm_name,"[",round(0.5*npredpoints),"] ) 
                      # }
                      # }
                      # 	
                      #----------------------------------#
                      ") 
  }
  
  writeLines(mod.code,con = "functions/text_to_add_to_jags_model.txt")
  outlist <- list(gamx.basis = gamx.basis,
                  gamx.basispred = gamx.basispred,
                  vscale = vscale,
                  orig.preds = orig.preds,
                  predpoints = predpoints,
                  nknots = nknots,
                  npredpoints = npredpoints,
                  mod.code = mod.code,
                  knotsgamx = knotsgamx)
  names(outlist) <- c(paste0(sm_name,"_basis"),
                      paste0(sm_name,"_basispred"),
                      paste0(sm_name,"_scaled_values"),
                      "original_predictor_values",
                      paste0(sm_name,"_visualized_predictor_values"),
                      paste0("nknots_",sm_name),
                      paste0("npredpoints_",sm_name),
                      paste0(sm_name,"_model_code"),
                      paste0(sm_name,"_knots"))
  
  
  return(outlist)
  
  
  
}