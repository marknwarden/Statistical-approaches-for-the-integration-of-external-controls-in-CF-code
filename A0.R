# ================================================================================================
# Name: A0
# Purpose:  Propensity Score Estimation      
# Programmer: Mark Warden
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 
# Input:
#       data - Dataframe from which to estimate a propensity score
#       p_model - the propensity score model to be fit (logistic for trial membership)
#       diagnostics - Will print a plot of the fit propensity scores for visual inspection
#       keep_link_fits - controls if propensity scores on the link scale should also be saved to the output
#
# Output:
#       trial - the input dataframe with the propensity scores and inverse probability weights appended
#                as column variables "fit" and "ipw_austin" (referencing the Austin et. al. paper that describes
#                how to fit the weights). Additionally, propensity scores on the link scale (log odds) can be saved
#                as well under the variable name "linkfits"
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


# 0. Loading Packages -----------------------------------------------------
library(ggplot2)


# 1. Propensity Score Estimation ------------------------------------------

A0<-function(data,p_model=(1-Hx) ~ X1 + X2,diagnostics=F,keep_link_fits=T){
  # (1-Hx) predicts membership in the active trial
  propensity<-glm(formula=p_model,data=data,family="binomial")
  fits<-predict(propensity,type="response")
  lfits<-predict(propensity,type="link")
  
  # Austin IPW weight formula
  ipw_austin<-function(z,p){
    # Z is the treatment assignment and p is the propensity score (i.e. the probability of Z=1|X)
    weights<-ifelse(z==1,1/p,1/(1-p))
    return(weights)
  }
  ipw_austin<-Vectorize(ipw_austin)
  
  propensity<-data.frame(p=fits,w=ipw_austin(z=1-data[,"Hx"],p=fits))  #Again, 1-trial[,"Hx"] is to define z as the active trial
  trial<-data.frame(data,propensity)
  
  if(keep_link_fits==T){trial<-data.frame(data,propensity,linkfits=lfits)}
  
  if(diagnostics==T){
    print("Includes a 0.01 propensity score jitter because of the binary variables.")
    print(ggplot(data = trial, aes(x = factor(Hx,levels=c(1,0),labels=c("Historical Trial","Current Trial")), y = p)) +
      geom_jitter(width = 0.05, height = 0.01) +
      labs(x = "Study Source", y = "Propensity Score") +
      coord_flip() +
      theme_bw())
    
    print("Weights Summary:")
    print(tapply(trial[,"w"],factor(trial[,"Hx"],levels=c(1,0),labels=c("Historical Trial","Current Trial")),summary))
    
    print("Propensity score Summary:")
    print(tapply(trial[,"p"],factor(trial[,"Hx"],levels=c(1,0),labels=c("Historical Trial","Current Trial")),summary))
  }
  
  return(trial)
}
