# ================================================================================================
# Name: A3
# Purpose:    Propensity score-based power priors for the simulation  
# Programmer: Mark Warden
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Input:
#     trial - dataframe that contains the simulated trial data
#     HxN_size - a way to manually set the minimum number of matched controls. In this method, a 
#     critical matching step is required. In the paper, we allowed the matching pool to be large
#     but the final number of Hx controls to be the same as the other compared methods
#     Hx_rel=1 - a second way to define the minimum number of matched controls relative to the size of the
#                trial's experimental treatment arm (0-1 as a percent)
# Output:
#     jags.coda - an MCMC list that contains the posterior samples for all relevant parameters
# ================================================================================================


# A3 function programming
library(rjags)
library(R2jags)
library(optmatch)
# remotes::install_url("https://cran.r-project.org/src/contrib/Archive/optmatch/optmatch_0.9-17.tar.gz")
library(RItools)
# library(MatchIt)
library(tidyverse)
library(HDInterval)
library(magrittr)

A3<-function(trial,HxN_size,Hx_rel=1){
  ### Corresponds to Matching Scheme 1 in Lin
  
  # Matching on Link-level fitted values using Mahalanobis distance
  # Calculates the distance matrix
  current_distance<-match_on(E ~ linkfits,
                              controls=2,
                              data=trial[which(trial$Hx==1|trial$E==1),],
                              method="mahalanobis")
  
  # Creates matched 1:1 pairs
  pair_match_result<-pairmatch(E ~ linkfits,
                               controls=2,
                               data=trial[which(trial$Hx==1|trial$E==1),],
                               distance=current_distance,
                               remove.unmatchables = FALSE)
  
  # Adds the matches to the dataset
  trial$pair<-NA
  trial[which(trial$Hx==1|trial$E==1),"pair"]<-pair_match_result
  
  # Subsets to trial with only Active trial participants and matched pairs
  # trial<-rbind(trial[which(trial$E==1),],
  #              trial[sample(which(trial$Hx!=1 & trial$E==0),size=ceiling(nrow(trial[which(trial$E==1),])/2)),],
  #              trial[sample(trial[which(trial$E==1),"pair"],size=ceiling(nrow(trial[which(trial$E==1),])/2)),]
  # )
  
  trial<-rbind(trial[which(trial$E==1),],
               trial[which(trial$Hx!=1 & trial$E==0),],
               trial[sample(which(trial$Hx==1 & !is.na(trial$pair)),size=min(HxN_size,Hx_rel*sum(trial$E))),]
  )
  trial %<>% filter(!is.na(Hx))

  # ## Model with uniform prior for d
  # modelstring <- "
  # model{
  # 
  # for (i in 1:(N_E+N_C)) {
  # y1[i] ~ dbern(p[i])
  # logit(p[i]) <- alpha + d*z1[i]
  # }
  # 
  # alpha ~ dnorm(prior.alpha.mean, prior.alpha.prec)
  # 
  # prior.alpha.mean ~ dnorm(0,0.001)
  # prior.alpha.prec ~ dgamma(0.001, 0.001)
  # 
  # k<-10000
  # for( i in 1:(Ntot-N)) {
  # l[i]<- a02[i]*(y2[i]*log(ilogit(alpha))+ (1-y2[i])*log((1-ilogit(alpha))))
  # phi[i]<-  -l[i]+k
  # zeros[i]~dpois(phi[i])
  # }
  # 
  # d~dunif(-100,100)
  # 
  # alpha_ilogit <- ilogit(alpha)
  # diff_ilogit <- ilogit(alpha+d) - alpha_ilogit
  # ratio_ilogit <- ilogit(alpha+d) / alpha_ilogit
  # HR <- (-log(1-ilogit(alpha+d))) / (-log(1-alpha_ilogit))
  # be<-log(HR)
  # 
  # }
  # 
  # 
  # "
  # model_file <- "JAGS_model_A3.bug"
  # writeLines(modelstring,con=model_file)
  
  ### INPUTS to JAGS
  #parameters    <- c("diff_ilogit","ratio_ilogit","d")
  parameters    <- c("be")
  data_list     <- list( N_E=sum(trial$E==1), N_C=sum(trial$E==0 & trial$Hx==0), Ntot=nrow(trial), N=sum(trial$Hx==0),
                         y1=trial[which(trial$Hx==0),"y"], 
                         y2=trial[which(trial$Hx==1),"y"],
                         z1=trial[which(trial$Hx==0),"E"], 
                         a02=trial[which(trial$Hx==1),"p"],
                         zeros=rep(0,sum(trial$Hx==1)))
  
  inits_list    <- list(alpha=0.5,
                        prior.alpha.mean=0,
                        d=0)
  
  jags.fit      <- jags.model(file='JAGS_model_A3.bug',data=data_list,inits=inits_list,n.chains=1,n.adapt=20000)
  jags.fit.out1 <- coda.samples(jags.fit, parameters, n.iter=120000)
  jags.coda     <- as.mcmc.list(jags.fit.out1) 

  
  
  # par(mfrow=c(3,1))
  # traceplot(jags.coda)
  # gelman.plot(jags.coda) # converge at around 20000 iteration
  # densplot(jags.coda)
  
  return(S=jags.coda)
}