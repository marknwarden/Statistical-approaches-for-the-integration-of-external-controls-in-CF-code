# ================================================================================================
# Name: A1
# Purpose:  Poisson Regression with Inverse probability weighting      
# Programmer: Mark Warden
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Input:
#       Data - dataframe that includes trial data, propensity scores, and inverse probability weights
#       O_model - the outcome model written in the form of a Poisson regression
#       full_info - whether all covariate model information shoudl be outputted or only treatment effect
# Output:
#       output - Calculates model terms, standard errors, p-values, robust standard errors and p-values,
#                and confidence intervals from Poisson Regression with Inverse probability weighting
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


# 0. Loading Needed Packages ----------------------------------------------
library(magrittr)
library(broom)
library(MASS)
library(sandwich)

# 1. Running the IPW Poisson Regression -----------------------------------
A1<-function(data,O_model=y ~ E + X1 + X2 + offset(log(t)),full_info=F){

outcome_model<-glm(formula=O_model,data=data,weights=w,family = "poisson")
output<-setNames((outcome_model %>% coef),c("A1_Int",paste0("A1_",names(coef(outcome_model))[-1])))
output %<>% c(setNames(coef(summary(outcome_model))[,"Std. Error"],c("A1_Int_naive_SE",paste0(paste0("A1_",names(coef(outcome_model))[-1]),"_naive_SE"))))
outcome_profile<-profile(outcome_model)
output %<>% c(setNames(confint(outcome_profile)[,1],c("A1_Int_LCI",paste0(paste0("A1_",names(coef(outcome_model))[-1]),"_LCI"))))
output %<>% c(setNames(confint(outcome_profile)[,2],c("A1_Int_UCI",paste0(paste0("A1_",names(coef(outcome_model))[-1]),"_UCI"))))
output %<>% c(setNames(coef(summary(outcome_model))[,"Pr(>|z|)"],c("A1_Int_profile_pv",paste0(paste0("A1_",names(coef(outcome_model))[-1]),"_profile_pv"))))

# Calculating Huber-White Sandwich estimates for the SE
std.err<-sqrt(diag(vcovHC(outcome_model,type="HC0")))
robust<-cbind(Estimate= coef(outcome_model), "Robust_SE" = std.err,
      "Pr(>|z|)" = 2 * pnorm(abs(coef(outcome_model)/std.err), lower.tail=FALSE),
      LL = coef(outcome_model) + qnorm(0.025) * std.err,
      UL = coef(outcome_model) + qnorm(0.975) * std.err)


output %<>% c(setNames(robust[,"Robust_SE"],c("A1_Int_robust_SE",paste0(paste0("A1_",names(coef(outcome_model))[-1]),"_robust_SE"))))
output %<>% c(setNames(robust[,"LL"],c("A1_Int_robust_LCI",paste0(paste0("A1_",names(coef(outcome_model))[-1]),"_robust_LCI"))))
output %<>% c(setNames(robust[,"UL"],c("A1_Int_robust_UCI",paste0(paste0("A1_",names(coef(outcome_model))[-1]),"_robust_UCI"))))
output %<>% c(setNames(robust[,"Pr(>|z|)"],c("A1_Int_robust_pv",paste0(paste0("A1_",names(coef(outcome_model))[-1]),"_robust_pv"))))

if(full_info==T){return(output)}

output<-output[grep("A1_E",names(output))]
return(output)
}
