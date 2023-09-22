# ================================================================================================
# Name: A2
# Purpose:    Simulation implementation of Commensurate priors   
# Programmer: Mark Warden
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Input:
#     trial - dataframe that contains a trial simulated using the D1 function
#
# Output:
#     jags.coda - an MCMC list that contains the posterior samples for all relevant parameters
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# A2 function programming
library(rjags)
library(R2jags)
library(RItools)
library(tidyverse)
library(HDInterval)

A2<-function(trial){

# ### Commensurate Prior Model
#   modelstring <- "
#   model{
#   
#   for (i in 1:N) {
#   y[i] ~ dpois(lambda[i])
#   log(lambda[i]) <- mu[i]
#   mu[i] <- b0[study[i]] + b1[study[i]]*x1[i] + b2[study[i]]*x2[i] + be*E[i] + log(time[i])
#   }
#   
#   # Note for the study variable 1=historical 2=active
#   # Therefore beta[1] are historical parameters and betas[2] are active trial parameters
# 
#   # Spike and slab hyper priors for tau
#   pi[1] ~ dbern(0.5)
#   pi[2] ~ dbern(0.5)
#   pi[3] ~ dbern(0.5)
#   
#   slab[1] ~ dunif(0.005,30)
#   slab[2] ~ dunif(0.005,100)
#   slab[3] ~ dunif(0.005,100)
#   
#   # Tau now is distributed as a mixture of the spike (at 200) and slab
#   # Even though tau is determinisitic (<-) it has the desired hyper prior distribution
#   tau[1] <- slab[1]*pi[1] + 30*(1-pi[1]) 
#   tau[2] <- slab[2]*pi[2] + 300*(1-pi[2]) 
#   tau[3] <- slab[3]*pi[3] + 300*(1-pi[3])
#   
#   # Vague normal prior's for historical parameters and treatment parameter (Note: [1] is Hx and [2] is active trial)
#   b0[1] ~ dnorm(0,1/10000)  
#   b1[1] ~ dnorm(0,1/10000)  # Note that 1/10000 is precision or the inverse of the variance
#   b2[1] ~ dnorm(0,1/10000)  # Meaning the variance is 10,000 and standard deviation is 100
#   be ~ dnorm(0,1/10000)
#   
#   # Normal distribution around historical beta with tau precision (Note: [1] is Hx and [2] is active trial)
#   b0[2] ~ dnorm(b0[1],tau[1])
#   b1[2] ~ dnorm(b1[1],tau[2])
#   b2[2] ~ dnorm(b2[1],tau[3])
#   
#   HR <- exp(be)
#   }
#   "
#   model_file <- "JAGS_model_A2.bug"
#   writeLines(modelstring,con=model_file)
  
  ### INPUTS to JAGS
  parameters    <- c("be")
  data_list     <- list(N=nrow(trial),
                        y=trial$y, 
                        study=ifelse(trial$Hx==1,1,2),
                        x1=trial$X1,
                        x2=trial$X2,
                        E=trial$E,
                        time=trial$t)
  
  inits_list    <- list(pi=c(1,1,1),
                        slab=c(1,1,1),
                        b0=c(0,0),
                        b1=c(0,0),
                        b2=c(0,0),
                        be=0)
  
  jags.fit      <- jags.model(file='JAGS_model_A2.bug',data=data_list,inits=inits_list,n.chains=1,n.adapt=20000)
  jags.fit.out1 <- coda.samples(jags.fit, parameters, n.iter=120000)
  jags.coda     <- as.mcmc.list(jags.fit.out1) 
  
  # par(mfrow=c(3,1))
  # traceplot(jags.coda)
  # gelman.plot(jags.coda) # converge at around 20000 iteration
  # densplot(jags.coda)  
  
}