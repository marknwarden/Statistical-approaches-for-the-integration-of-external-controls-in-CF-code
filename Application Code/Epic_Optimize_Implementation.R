# ================================================================================================
# Name: Epic_Optimize_implementations
# Purpose:    To implement IPW, Propensity score-based power priors, and commensurate priors in Epic/Optimize   
# Programmer: Mark Warden
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Changes and Edits:
# 
# Input:
#
# Output:
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# 0. Loading needed packages ----------------------------------------------

# Propensity Score Base Power Priors
library(rjags)
library(R2jags)
library(optmatch)
# remotes::install_url("https://cran.r-project.org/src/contrib/Archive/optmatch/optmatch_0.9-17.tar.gz")
library(RItools)
# library(MatchIt)
library(tidyverse)
library(HDInterval)
library(magrittr)

# Commensurate priors
library(rjags)
library(R2jags)
library(RItools)
library(tidyverse)
library(HDInterval)
library(runjags)
library(bootstrap)

# 1. Making small needed modifications to the propensity score and --------
#    Analytic datasets

# Saving new analytic dataset and calculating needed variables with imputation
trial<-Epic_Opt_clean

# Generating a propensity score value for those with missing variables (at most 10 participants)
links<-predict(multivariable)
trial$links<-links
rm(links)

# Simplified propensity score for those with missing values
multivariable_imp<-lrm(formula=source2 ~ agecat+sex_txt+NonWhite+genocat+dalfa+hs+asthma+PANUM,
                   data=propensity)

links<-predict(multivariable_imp)
fits<-predict(multivariable_imp,type="fitted")

# trial[which(is.na(trial$fits)),] %>% dplyr::select(source,subject,agecat,sex_txt,NonWhite,genocat,fevppcat2_txt,haz,dalfa,hs,asthma,pabase,PANUM)
# trial[which(is.na(trial$fits)),"weights"]<-ipw_austin(trial[which(is.na(trial$fits)),"az2",drop=T],fits[which(is.na(trial$fits))])
trial[which(is.na(trial$fits)),"fits"]<-fits[which(is.na(trial$fits))]
trial[which(is.na(trial$links)),"links"]<-links[which(is.na(trial$links))]

rm(fits)
rm(links)

# summary(coxph(formula=Surv(time,pe1_event) ~ az2+agecat_txt,
#       data=trial,
#       weights=weights,
#       robust=F))

# Modifying the bcanon function to also calculate a p-value
bcanon_p<-function (x, nboot, theta, ..., alpha = c(0.025, 0.05, 0.1, 0.16, 
                                          0.84, 0.9, 0.95, 0.975)) {
  if (!all(alpha < 1) || !all(alpha > 0)) 
    stop("All elements of alpha must be in (0,1)")
  alpha_sorted <- sort(alpha)
  if (nboot <= 1/min(alpha_sorted[1], 1 - alpha_sorted[length(alpha_sorted)])) 
    warning("nboot is not large enough to estimate your chosen alpha.")
  call <- match.call()
  n <- length(x)
  thetahat <- theta(x, ...)
  bootsam <- matrix(sample(x, size = n * nboot, replace = TRUE), 
                    nrow = nboot)
  thetastar <- apply(bootsam, 1, theta, ...)
  z0 <- qnorm(sum(thetastar < thetahat)/nboot)
  u <- rep(0, n)
  for (i in 1:n) {
    u[i] <- theta(x[-i], ...)
  }
  uu <- mean(u) - u
  acc <- sum(uu * uu * uu)/(6 * (sum(uu * uu))^1.5)
  zalpha <- qnorm(alpha)
  tt <- pnorm(z0 + (z0 + zalpha)/(1 - acc * (z0 + zalpha)))
  confpoints <- quantile(x = thetastar, probs = tt, type = 1)
  names(confpoints) <- NULL
  confpoints <- cbind(alpha, confpoints,exp(confpoints))
  dimnames(confpoints)[[2]] <- c("alpha", "bca point","Confint")
  pvalue<-2*((1+sum(thetastar>=0))/(nboot+1))
  percentile_0<-ecdf(thetastar)(0)
  pvalueBCa<-2*(1-uniroot(f=function(x){pnorm(z0 + (z0 + qnorm(x))/(1 - acc * (z0 + qnorm(x))))-percentile_0},c(0.001,0.999))$root)
  return(list(confpoints = confpoints, z0 = z0, acc = acc, 
              p = pvalue,pvalue_BCa=pvalueBCa,estimate=c(link=thetahat,HR=exp(thetahat)),
              bootstrap_sample<-thetastar,call = call))
}

# 2. [Augmenting, Bootstrap] IPW -----------------------------------------------------
# Saving new analytic dataset and calculating needed variables with imputation
trial<-propensity

# Generating a propensity score value for those with missing variables
links<-predict(multivariable)
trial$links<-links
rm(links)

# Simplified propensity score for those with missing values
multivariable_imp<-lrm(formula=source2 ~ agecat+sex_txt+NonWhite+genocat+dalfa+hs+asthma+PANUM,
                       data=propensity)

links<-predict(multivariable_imp)
fits<-predict(multivariable_imp,type="fitted")

# trial[which(is.na(trial$fits)),] %>% dplyr::select(source,subject,agecat,sex_txt,NonWhite,genocat,fevppcat2_txt,haz,dalfa,hs,asthma,pabase,PANUM)
# trial[which(is.na(trial$fits)),"weights"]<-ipw_austin(trial[which(is.na(trial$fits)),"az2",drop=T],fits[which(is.na(trial$fits))])
trial[which(is.na(trial$fits)),"fits"]<-fits[which(is.na(trial$fits))]
trial[which(is.na(trial$links)),"links"]<-links[which(is.na(trial$links))]

rm(fits)
rm(links)

# Creating the 100 samples with N=18 Optimize Placebo
set.seed(1121992)
{aug_ipw<-list((trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
              (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
              (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
              (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
              (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
              (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
              (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
              (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
              (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
              (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
              (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
              (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
              (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
              (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
              (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
              (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
              (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
              (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
              (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
              (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
              (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
              (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
              (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
              (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
              (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
              (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
              (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
              (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
              (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
              (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
              (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
              (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
              (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
              (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
              (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
              (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
              (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
              (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
              (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
              (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
              (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
              (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
              (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
              (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
              (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
              (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
              (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
              (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
              (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
              (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
              (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
              (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
              (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
              (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
              (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
              (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
              (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
              (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
              (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
              (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
              (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
              (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
              (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
              (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
              (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
              (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
              (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
              (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
              (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
              (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
              (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
              (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
              (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
              (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
              (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
              (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
              (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
              (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
              (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
              (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
              (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
              (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
              (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
              (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
              (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
              (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
              (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
              (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
              (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
              (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
              (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
              (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
              (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
              (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
              (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
              (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
              (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
              (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
              (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
              (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)))}
aug_ipw %<>% lapply(function(z){rbind(z,(trial %>% filter(TRT!=6)))})

# Bootstrapping - A single random draw of half the controls as a check
treat<-function(x,data){coef(glm(formula=pe1_event ~ az2+agecat_txt+offset(log(time)),
                                 data=data[x,],
                                 weights=weights,
                                 family = "poisson"))["az2TRUE"]}
# bcanon_p(1:nrow(aug_ipw[[1]]),10000,treat,data=aug_ipw[[1]],alpha=c(0.025,0.975))[1:6]

# Doing one massive bootstrap
aug_ipw %<>% lapply(function(z){bcanon_p(1:nrow(z),5000,treat,data=z,alpha=c(0.025,0.975))[1:7]})

apply(Reduce(rbind,lapply(aug_ipw,function(x){x$estimate})),2,mean)
apply(Reduce(rbind,lapply(aug_ipw,function(x){unlist(x[1])})),2,summary)
apply(Reduce(rbind,lapply(aug_ipw,function(x){unlist(x[4:5])})),2,mean)

thetahat<-apply(Reduce(rbind,lapply(aug_ipw,function(x){x$estimate[1]})),2,mean)
z0<-apply(Reduce(rbind,lapply(aug_ipw,function(x){x$z0})),2,mean)
acc<-apply(Reduce(rbind,lapply(aug_ipw,function(x){x$acc})),2,mean)
thetastar<-unlist(lapply(aug_ipw,function(x){x[[7]]}))
zalpha <- qnorm(c(0.025,0.975))
tt <- pnorm(z0 + (z0 + zalpha)/(1 - acc * (z0 + zalpha)))
confpoints <- quantile(x = thetastar, probs = tt, type = 1)
names(confpoints) <- NULL
confpoints <- cbind(c(0.025,0.975), confpoints,exp(confpoints))
dimnames(confpoints)[[2]] <- c("alpha", "bca point","Confint")
pvalue<-2*((1+sum(thetastar>=0))/(length(thetastar)+1))
percentile_0<-ecdf(thetastar)(0)
pvalueBCa<-2*(1-uniroot(f=function(x){pnorm(z0 + (z0 + qnorm(x))/(1 - acc * (z0 + qnorm(x))))-percentile_0},c(0.001,0.999))$root)
list(confpoints = confpoints, z0 = z0, acc = acc, 
     p = pvalue,pvalue_BCa=pvalueBCa,estimate=setNames(c(thetahat,exp(thetahat)),c("link","HR")))
rm(list = c("thetahat","z0","acc","thetastar","zalpha","tt","confpoints","pvalue","percentile_0","pvalueBCa"))

# Creating the 100 samples with N=35 Optimize Placebo
set.seed(1121993)
{aug_ipw<-list((trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)))}
aug_ipw %<>% lapply(function(z){rbind(z,(trial %>% filter(TRT!=6)))})

# Doing one massive bootstrap
aug_ipw %<>% lapply(function(z){bcanon_p(1:nrow(z),5000,treat,data=z,alpha=c(0.025,0.975))[1:7]})

apply(Reduce(rbind,lapply(aug_ipw,function(x){x$estimate})),2,mean)
apply(Reduce(rbind,lapply(aug_ipw,function(x){unlist(x[1])})),2,summary)
apply(Reduce(rbind,lapply(aug_ipw,function(x){unlist(x[4:5])})),2,mean)

thetahat<-apply(Reduce(rbind,lapply(aug_ipw,function(x){x$estimate[1]})),2,mean)
z0<-apply(Reduce(rbind,lapply(aug_ipw,function(x){x$z0})),2,mean)
acc<-apply(Reduce(rbind,lapply(aug_ipw,function(x){x$acc})),2,mean)
thetastar<-unlist(lapply(aug_ipw,function(x){x[[7]]}))
zalpha <- qnorm(c(0.025,0.975))
tt <- pnorm(z0 + (z0 + zalpha)/(1 - acc * (z0 + zalpha)))
confpoints <- quantile(x = thetastar, probs = tt, type = 1)
names(confpoints) <- NULL
confpoints <- cbind(c(0.025,0.975), confpoints,exp(confpoints))
dimnames(confpoints)[[2]] <- c("alpha", "bca point","Confint")
pvalue<-2*((1+sum(thetastar>=0))/(length(thetastar)+1))
percentile_0<-ecdf(thetastar)(0)
pvalueBCa<-2*(1-uniroot(f=function(x){pnorm(z0 + (z0 + qnorm(x))/(1 - acc * (z0 + qnorm(x))))-percentile_0},c(0.001,0.999))$root)
list(confpoints = confpoints, z0 = z0, acc = acc, 
     p = pvalue,pvalue_BCa=pvalueBCa,estimate=setNames(c(thetahat,exp(thetahat)),c("link","HR")))
rm(list = c("thetahat","z0","acc","thetastar","zalpha","tt","confpoints","pvalue","percentile_0","pvalueBCa"))

# Creating the 100 samples with N=43 Optimize Placebo
set.seed(1121994)
{aug_ipw<-list((trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)))}
aug_ipw %<>% lapply(function(z){rbind(z,(trial %>% filter(TRT!=6)))})

# Doing one massive bootstrap
aug_ipw %<>% lapply(function(z){bcanon_p(1:nrow(z),5000,treat,data=z,alpha=c(0.025,0.975))[1:7]})

apply(Reduce(rbind,lapply(aug_ipw,function(x){x$estimate})),2,mean)
apply(Reduce(rbind,lapply(aug_ipw,function(x){unlist(x[1])})),2,summary)
apply(Reduce(rbind,lapply(aug_ipw,function(x){unlist(x[4:5])})),2,mean)

thetahat<-apply(Reduce(rbind,lapply(aug_ipw,function(x){x$estimate[1]})),2,mean)
z0<-apply(Reduce(rbind,lapply(aug_ipw,function(x){x$z0})),2,mean)
acc<-apply(Reduce(rbind,lapply(aug_ipw,function(x){x$acc})),2,mean)
thetastar<-unlist(lapply(aug_ipw,function(x){x[[7]]}))
zalpha <- qnorm(c(0.025,0.975))
tt <- pnorm(z0 + (z0 + zalpha)/(1 - acc * (z0 + zalpha)))
confpoints <- quantile(x = thetastar, probs = tt, type = 1)
names(confpoints) <- NULL
confpoints <- cbind(c(0.025,0.975), confpoints,exp(confpoints))
dimnames(confpoints)[[2]] <- c("alpha", "bca point","Confint")
pvalue<-2*((1+sum(thetastar>=0))/(length(thetastar)+1))
percentile_0<-ecdf(thetastar)(0)
pvalueBCa<-2*(1-uniroot(f=function(x){pnorm(z0 + (z0 + qnorm(x))/(1 - acc * (z0 + qnorm(x))))-percentile_0},c(0.001,0.999))$root)
list(confpoints = confpoints, z0 = z0, acc = acc, 
     p = pvalue,pvalue_BCa=pvalueBCa,estimate=setNames(c(thetahat,exp(thetahat)),c("link","HR")))
rm(list = c("thetahat","z0","acc","thetastar","zalpha","tt","confpoints","pvalue","percentile_0","pvalueBCa"))


# Creating the 100 samples with N=52 Optimize Placebo
set.seed(1121995)
{aug_ipw<-list((trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)))}
aug_ipw %<>% lapply(function(z){rbind(z,(trial %>% filter(TRT!=6)))})

# Doing one massive bootstrap
aug_ipw %<>% lapply(function(z){bcanon_p(1:nrow(z),5000,treat,data=z,alpha=c(0.025,0.975))[1:7]})

apply(Reduce(rbind,lapply(aug_ipw,function(x){x$estimate})),2,mean)
apply(Reduce(rbind,lapply(aug_ipw,function(x){unlist(x[1])})),2,summary)
apply(Reduce(rbind,lapply(aug_ipw,function(x){unlist(x[4:5])})),2,mean)

thetahat<-apply(Reduce(rbind,lapply(aug_ipw,function(x){x$estimate[1]})),2,mean)
z0<-apply(Reduce(rbind,lapply(aug_ipw,function(x){x$z0})),2,mean)
acc<-apply(Reduce(rbind,lapply(aug_ipw,function(x){x$acc})),2,mean)
thetastar<-unlist(lapply(aug_ipw,function(x){x[[7]]}))
zalpha <- qnorm(c(0.025,0.975))
tt <- pnorm(z0 + (z0 + zalpha)/(1 - acc * (z0 + zalpha)))
confpoints <- quantile(x = thetastar, probs = tt, type = 1)
names(confpoints) <- NULL
confpoints <- cbind(c(0.025,0.975), confpoints,exp(confpoints))
dimnames(confpoints)[[2]] <- c("alpha", "bca point","Confint")
pvalue<-2*((1+sum(thetastar>=0))/(length(thetastar)+1))
percentile_0<-ecdf(thetastar)(0)
pvalueBCa<-2*(1-uniroot(f=function(x){pnorm(z0 + (z0 + qnorm(x))/(1 - acc * (z0 + qnorm(x))))-percentile_0},c(0.001,0.999))$root)
list(confpoints = confpoints, z0 = z0, acc = acc, 
     p = pvalue,pvalue_BCa=pvalueBCa,estimate=setNames(c(thetahat,exp(thetahat)),c("link","HR")))
rm(list = c("thetahat","z0","acc","thetastar","zalpha","tt","confpoints","pvalue","percentile_0","pvalueBCa"))

# Creating the 100 samples with N=69 Optimize Placebo
set.seed(1121996)
{aug_ipw<-list((trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
               (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)))}
aug_ipw %<>% lapply(function(z){rbind(z,(trial %>% filter(TRT!=6)))})

# Doing one massive bootstrap
aug_ipw %<>% lapply(function(z){bcanon_p(1:nrow(z),5000,treat,data=z,alpha=c(0.025,0.975))[1:7]})

apply(Reduce(rbind,lapply(aug_ipw,function(x){x$estimate})),2,mean)
apply(Reduce(rbind,lapply(aug_ipw,function(x){unlist(x[1])})),2,summary)
apply(Reduce(rbind,lapply(aug_ipw,function(x){unlist(x[4:5])})),2,mean)

thetahat<-apply(Reduce(rbind,lapply(aug_ipw,function(x){x$estimate[1]})),2,mean)
z0<-apply(Reduce(rbind,lapply(aug_ipw,function(x){x$z0})),2,mean)
acc<-apply(Reduce(rbind,lapply(aug_ipw,function(x){x$acc})),2,mean)
thetastar<-unlist(lapply(aug_ipw,function(x){x[[7]]}))
zalpha <- qnorm(c(0.025,0.975))
tt <- pnorm(z0 + (z0 + zalpha)/(1 - acc * (z0 + zalpha)))
confpoints <- quantile(x = thetastar, probs = tt, type = 1)
names(confpoints) <- NULL
confpoints <- cbind(c(0.025,0.975), confpoints,exp(confpoints))
dimnames(confpoints)[[2]] <- c("alpha", "bca point","Confint")
pvalue<-2*((1+sum(thetastar>=0))/(length(thetastar)+1))
percentile_0<-ecdf(thetastar)(0)
pvalueBCa<-2*(1-uniroot(f=function(x){pnorm(z0 + (z0 + qnorm(x))/(1 - acc * (z0 + qnorm(x))))-percentile_0},c(0.001,0.999))$root)
list(confpoints = confpoints, z0 = z0, acc = acc, 
     p = pvalue,pvalue_BCa=pvalueBCa,estimate=setNames(c(thetahat,exp(thetahat)),c("link","HR")))
rm(list = c("thetahat","z0","acc","thetastar","zalpha","tt","confpoints","pvalue","percentile_0","pvalueBCa"))

# 2.a [Pooling] IPW -------------------------------------------------------
# Pooling
bcanon_p(1:nrow(trial),1000,treat,data=trial,alpha=c(0.025,0.975))[1:6]
bcanon_p(1:nrow(trial),5000,treat,data=trial,alpha=c(0.025,0.975))[1:6]
bcanon_p(1:nrow(trial),10000,treat,data=trial,alpha=c(0.025,0.975))[1:6]

# 3. [Augmenting, Multi-chain] PSBPP ---------------------------------------------------
# Saving new analytic dataset and calculating needed variables with imputation
trial<-propensity

# Generating a propensity score value for those with missing variables
links<-predict(multivariable)
trial$links<-links
rm(links)

# Simplified propensity score for those with missing values
multivariable_imp<-lrm(formula=source2 ~ agecat+sex_txt+NonWhite+genocat+dalfa+hs+asthma+PANUM,
                       data=propensity)

links<-predict(multivariable_imp)
fits<-predict(multivariable_imp,type="fitted")

# trial[which(is.na(trial$fits)),] %>% dplyr::select(source,subject,agecat,sex_txt,NonWhite,genocat,fevppcat2_txt,haz,dalfa,hs,asthma,pabase,PANUM)
# trial[which(is.na(trial$fits)),"weights"]<-ipw_austin(trial[which(is.na(trial$fits)),"az2",drop=T],fits[which(is.na(trial$fits))])
trial[which(is.na(trial$fits)),"fits"]<-fits[which(is.na(trial$fits))]
trial[which(is.na(trial$links)),"links"]<-links[which(is.na(trial$links))]

rm(fits)
rm(links)

# Creating the 100 samples with N=18 Optimize Placebo
set.seed(8231989)
{aug_pspp<-list((trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)))}
aug_pspp %<>% lapply(function(z){rbind(z,(trial %>% filter(TRT!=6)))})

# Creating the Analytic Samples
aug_pspp %<>% lapply(function(z){
  
  # Matching on Link-level fitted values using Mahalanobis distance
  # Calculates the distance matrix
  current_distance<-match_on(source2 ~ links,
                             controls=3,
                             data=z[which(z$source_std==1|z$az2==TRUE),],
                             method="mahalanobis")
  
  # Creates matched 1:1 pairs
  pair_match_result<-pairmatch(source2 ~ links,
                               controls=3,
                               data=z[which(z$source_std==1|z$az2==TRUE),],
                               distance=current_distance,
                               remove.unmatchables = FALSE)
  
  # Adds the matches to the dataset
  z$pair<-NA
  z[which(z$source_std==1|z$az2==TRUE),"pair"]<-pair_match_result
  
  z<-rbind(z[which(z$az2==1),],
           z[which(z$source_std!=1 & z$az2==0),],
           z[sample(which(z$source_std==1 & !is.na(z$pair)),size=3*sum(z$az2)),]
  )
  z %<>% filter(!is.na(source_std))
  return(z)
})
# Running the Bayesian Analysis
aug_pspp %<>% lapply(function(z){
  parameters    <- c("be")
  data_list     <- list( N_E=sum(z$az2==1), N_C=sum(z$az2==0 & z$source_std==2), Ntot=nrow(z), N=sum(z$source_std==2),
                         y1=z[which(z$source_std==2),"pe1_event",drop=TRUE], 
                         y2=z[which(z$source_std==1),"pe1_event",drop=TRUE],
                         z1=as.integer(z[which(z$source_std==2),"az2",drop=TRUE]), 
                         a02=z[which(z$source_std==1),"fits",drop=TRUE],
                         zeros=rep(0,sum(z$source_std==1)))
  
  inits_list    <- list(alpha=0.5,
                        prior.alpha.mean=0,
                        d=0)
  
  jags.fit      <- jags.model(file='JAGS_model_A3.bug',data=data_list,inits=inits_list,n.chains=4,n.adapt=20000)
  jags.fit.out1 <- coda.samples(jags.fit, parameters, n.iter=120000, thin = 2)
  jags.coda     <- as.mcmc.list(jags.fit.out1)
  
  # Generating Results
  return(jags.coda)
})

# Combine 100 4-chain MCMCs into a single object
aug_pspp %<>% combine.mcmc(thin=2)
aug_pspp_chains_combined <- combine.mcmc(aug_pspp)
aug_pspp_table_less0<-table(aug_pspp_chains_combined<0)

# Final summaries
saveRDS(c(Mode=density(aug_pspp_chains_combined)$x[which.max(density(aug_pspp_chains_combined)$y)],
  Two_sided_P_value=2*min(ifelse(is.na(prop.table(aug_pspp_table_less0)["FALSE"]),0,prop.table(aug_pspp_table_less0)["FALSE"]),
                          ifelse(is.na(prop.table(aug_pspp_table_less0)["TRUE"]),0,prop.table(aug_pspp_table_less0)["TRUE"])),
  One_sided_P_Value=min(ifelse(is.na(prop.table(aug_pspp_table_less0)["FALSE"]),0,prop.table(aug_pspp_table_less0)["FALSE"]),
                        ifelse(is.na(prop.table(aug_pspp_table_less0)["TRUE"]),0,prop.table(aug_pspp_table_less0)["TRUE"])),
  Lower_CI=hdi(aug_pspp)[1],Upper_Ci=hdi(aug_pspp)[2]),"N18_psbpp.rds")

# Creating the 100 samples with N=35 Optimize Placebo
set.seed(8231990)
{aug_pspp<-list((trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)))}
aug_pspp %<>% lapply(function(z){rbind(z,(trial %>% filter(TRT!=6)))})

# Creating the Analytic Samples
aug_pspp %<>% lapply(function(z){
  
  # Matching on Link-level fitted values using Mahalanobis distance
  # Calculates the distance matrix
  current_distance<-match_on(source2 ~ links,
                             controls=3,
                             data=z[which(z$source_std==1|z$az2==TRUE),],
                             method="mahalanobis")
  
  # Creates matched 1:1 pairs
  pair_match_result<-pairmatch(source2 ~ links,
                               controls=3,
                               data=z[which(z$source_std==1|z$az2==TRUE),],
                               distance=current_distance,
                               remove.unmatchables = FALSE)
  
  # Adds the matches to the dataset
  z$pair<-NA
  z[which(z$source_std==1|z$az2==TRUE),"pair"]<-pair_match_result
  
  z<-rbind(z[which(z$az2==1),],
           z[which(z$source_std!=1 & z$az2==0),],
           z[sample(which(z$source_std==1 & !is.na(z$pair)),size=3*sum(z$az2)),]
  )
  z %<>% filter(!is.na(source_std))
  return(z)
})
# Running the Bayesian Analysis
aug_pspp %<>% lapply(function(z){
  parameters    <- c("be")
  data_list     <- list( N_E=sum(z$az2==1), N_C=sum(z$az2==0 & z$source_std==2), Ntot=nrow(z), N=sum(z$source_std==2),
                         y1=z[which(z$source_std==2),"pe1_event",drop=TRUE], 
                         y2=z[which(z$source_std==1),"pe1_event",drop=TRUE],
                         z1=as.integer(z[which(z$source_std==2),"az2",drop=TRUE]), 
                         a02=z[which(z$source_std==1),"fits",drop=TRUE],
                         zeros=rep(0,sum(z$source_std==1)))
  
  inits_list    <- list(alpha=0.5,
                        prior.alpha.mean=0,
                        d=0)
  
  jags.fit      <- jags.model(file='JAGS_model_A3.bug',data=data_list,inits=inits_list,n.chains=4,n.adapt=20000)
  jags.fit.out1 <- coda.samples(jags.fit, parameters, n.iter=120000, thin = 2)
  jags.coda     <- as.mcmc.list(jags.fit.out1)
  
  # Generating Results
  return(jags.coda)
})

# Combine 100 4-chain MCMCs into a single object
aug_pspp %<>% combine.mcmc(thin=2)
aug_pspp_chains_combined <- combine.mcmc(aug_pspp)
aug_pspp_table_less0<-table(aug_pspp_chains_combined<0)

# Final summaries
saveRDS(c(Mode=density(aug_pspp_chains_combined)$x[which.max(density(aug_pspp_chains_combined)$y)],
  Two_sided_P_value=2*min(ifelse(is.na(prop.table(aug_pspp_table_less0)["FALSE"]),0,prop.table(aug_pspp_table_less0)["FALSE"]),
                          ifelse(is.na(prop.table(aug_pspp_table_less0)["TRUE"]),0,prop.table(aug_pspp_table_less0)["TRUE"])),
  One_sided_P_Value=min(ifelse(is.na(prop.table(aug_pspp_table_less0)["FALSE"]),0,prop.table(aug_pspp_table_less0)["FALSE"]),
                        ifelse(is.na(prop.table(aug_pspp_table_less0)["TRUE"]),0,prop.table(aug_pspp_table_less0)["TRUE"])),
  Lower_CI=hdi(aug_pspp)[1],Upper_Ci=hdi(aug_pspp)[2]),"N35_psbpp.rds")

# Creating the 100 samples with N=43 Optimize Placebo
set.seed(8231991)
{aug_pspp<-list((trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)))}
aug_pspp %<>% lapply(function(z){rbind(z,(trial %>% filter(TRT!=6)))})

# Creating the Analytic Samples
aug_pspp %<>% lapply(function(z){
  
  # Matching on Link-level fitted values using Mahalanobis distance
  # Calculates the distance matrix
  current_distance<-match_on(source2 ~ links,
                             controls=3,
                             data=z[which(z$source_std==1|z$az2==TRUE),],
                             method="mahalanobis")
  
  # Creates matched 1:1 pairs
  pair_match_result<-pairmatch(source2 ~ links,
                               controls=3,
                               data=z[which(z$source_std==1|z$az2==TRUE),],
                               distance=current_distance,
                               remove.unmatchables = FALSE)
  
  # Adds the matches to the dataset
  z$pair<-NA
  z[which(z$source_std==1|z$az2==TRUE),"pair"]<-pair_match_result
  
  z<-rbind(z[which(z$az2==1),],
           z[which(z$source_std!=1 & z$az2==0),],
           z[sample(which(z$source_std==1 & !is.na(z$pair)),size=3*sum(z$az2)),]
  )
  z %<>% filter(!is.na(source_std))
  return(z)
})
# Running the Bayesian Analysis
aug_pspp %<>% lapply(function(z){
  parameters    <- c("be")
  data_list     <- list( N_E=sum(z$az2==1), N_C=sum(z$az2==0 & z$source_std==2), Ntot=nrow(z), N=sum(z$source_std==2),
                         y1=z[which(z$source_std==2),"pe1_event",drop=TRUE], 
                         y2=z[which(z$source_std==1),"pe1_event",drop=TRUE],
                         z1=as.integer(z[which(z$source_std==2),"az2",drop=TRUE]), 
                         a02=z[which(z$source_std==1),"fits",drop=TRUE],
                         zeros=rep(0,sum(z$source_std==1)))
  
  inits_list    <- list(alpha=0.5,
                        prior.alpha.mean=0,
                        d=0)
  
  jags.fit      <- jags.model(file='JAGS_model_A3.bug',data=data_list,inits=inits_list,n.chains=4,n.adapt=20000)
  jags.fit.out1 <- coda.samples(jags.fit, parameters, n.iter=120000, thin = 2)
  jags.coda     <- as.mcmc.list(jags.fit.out1)
  
  # Generating Results
  return(jags.coda)
})

# Combine 100 4-chain MCMCs into a single object
aug_pspp %<>% combine.mcmc(thin=2)
aug_pspp_chains_combined <- combine.mcmc(aug_pspp)
aug_pspp_table_less0<-table(aug_pspp_chains_combined<0)

# Final summaries
saveRDS(c(Mode=density(aug_pspp_chains_combined)$x[which.max(density(aug_pspp_chains_combined)$y)],
  Two_sided_P_value=2*min(ifelse(is.na(prop.table(aug_pspp_table_less0)["FALSE"]),0,prop.table(aug_pspp_table_less0)["FALSE"]),
                          ifelse(is.na(prop.table(aug_pspp_table_less0)["TRUE"]),0,prop.table(aug_pspp_table_less0)["TRUE"])),
  One_sided_P_Value=min(ifelse(is.na(prop.table(aug_pspp_table_less0)["FALSE"]),0,prop.table(aug_pspp_table_less0)["FALSE"]),
                        ifelse(is.na(prop.table(aug_pspp_table_less0)["TRUE"]),0,prop.table(aug_pspp_table_less0)["TRUE"])),
  Lower_CI=hdi(aug_pspp)[1],Upper_Ci=hdi(aug_pspp)[2]),"N43_psbpp.rds")

# Creating the 100 samples with N=52 Optimize Placebo
set.seed(8231992)
{aug_pspp<-list((trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)))}
aug_pspp %<>% lapply(function(z){rbind(z,(trial %>% filter(TRT!=6)))})

# Creating the Analytic Samples
aug_pspp %<>% lapply(function(z){
  
  # Matching on Link-level fitted values using Mahalanobis distance
  # Calculates the distance matrix
  current_distance<-match_on(source2 ~ links,
                             controls=3,
                             data=z[which(z$source_std==1|z$az2==TRUE),],
                             method="mahalanobis")
  
  # Creates matched 1:1 pairs
  pair_match_result<-pairmatch(source2 ~ links,
                               controls=3,
                               data=z[which(z$source_std==1|z$az2==TRUE),],
                               distance=current_distance,
                               remove.unmatchables = FALSE)
  
  # Adds the matches to the dataset
  z$pair<-NA
  z[which(z$source_std==1|z$az2==TRUE),"pair"]<-pair_match_result
  
  z<-rbind(z[which(z$az2==1),],
           z[which(z$source_std!=1 & z$az2==0),],
           z[sample(which(z$source_std==1 & !is.na(z$pair)),size=3*sum(z$az2)),]
  )
  z %<>% filter(!is.na(source_std))
  return(z)
})
# Running the Bayesian Analysis
aug_pspp %<>% lapply(function(z){
  parameters    <- c("be")
  data_list     <- list( N_E=sum(z$az2==1), N_C=sum(z$az2==0 & z$source_std==2), Ntot=nrow(z), N=sum(z$source_std==2),
                         y1=z[which(z$source_std==2),"pe1_event",drop=TRUE], 
                         y2=z[which(z$source_std==1),"pe1_event",drop=TRUE],
                         z1=as.integer(z[which(z$source_std==2),"az2",drop=TRUE]), 
                         a02=z[which(z$source_std==1),"fits",drop=TRUE],
                         zeros=rep(0,sum(z$source_std==1)))
  
  inits_list    <- list(alpha=0.5,
                        prior.alpha.mean=0,
                        d=0)
  
  jags.fit      <- jags.model(file='JAGS_model_A3.bug',data=data_list,inits=inits_list,n.chains=4,n.adapt=20000)
  jags.fit.out1 <- coda.samples(jags.fit, parameters, n.iter=120000, thin = 2)
  jags.coda     <- as.mcmc.list(jags.fit.out1)
  
  # Generating Results
  return(jags.coda)
})

# Combine 100 4-chain MCMCs into a single object
aug_pspp %<>% combine.mcmc(thin=2)
aug_pspp_chains_combined <- combine.mcmc(aug_pspp)
aug_pspp_table_less0<-table(aug_pspp_chains_combined<0)

# Final summaries
saveRDS(c(Mode=density(aug_pspp_chains_combined)$x[which.max(density(aug_pspp_chains_combined)$y)],
  Two_sided_P_value=2*min(ifelse(is.na(prop.table(aug_pspp_table_less0)["FALSE"]),0,prop.table(aug_pspp_table_less0)["FALSE"]),
                          ifelse(is.na(prop.table(aug_pspp_table_less0)["TRUE"]),0,prop.table(aug_pspp_table_less0)["TRUE"])),
  One_sided_P_Value=min(ifelse(is.na(prop.table(aug_pspp_table_less0)["FALSE"]),0,prop.table(aug_pspp_table_less0)["FALSE"]),
                        ifelse(is.na(prop.table(aug_pspp_table_less0)["TRUE"]),0,prop.table(aug_pspp_table_less0)["TRUE"])),
  Lower_CI=hdi(aug_pspp)[1],Upper_Ci=hdi(aug_pspp)[2]),"N52_psbpp.rds")

# Creating the 100 samples with N=69 Optimize Placebo
set.seed(8231993)
{aug_pspp<-list((trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)))}
aug_pspp %<>% lapply(function(z){rbind(z,(trial %>% filter(TRT!=6)))})
# Creating the Analytic Samples
aug_pspp %<>% lapply(function(z){
  
  # Matching on Link-level fitted values using Mahalanobis distance
  # Calculates the distance matrix
  current_distance<-match_on(source2 ~ links,
                             controls=3,
                             data=z[which(z$source_std==1|z$az2==TRUE),],
                             method="mahalanobis")
  
  # Creates matched 1:1 pairs
  pair_match_result<-pairmatch(source2 ~ links,
                               controls=3,
                               data=z[which(z$source_std==1|z$az2==TRUE),],
                               distance=current_distance,
                               remove.unmatchables = FALSE)
  
  # Adds the matches to the dataset
  z$pair<-NA
  z[which(z$source_std==1|z$az2==TRUE),"pair"]<-pair_match_result
  
  z<-rbind(z[which(z$az2==1),],
           z[which(z$source_std!=1 & z$az2==0),],
           z[sample(which(z$source_std==1 & !is.na(z$pair)),size=3*sum(z$az2)),]
  )
  z %<>% filter(!is.na(source_std))
  return(z)
})
# Running the Bayesian Analysis
aug_pspp %<>% lapply(function(z){
  parameters    <- c("be")
  data_list     <- list( N_E=sum(z$az2==1), N_C=sum(z$az2==0 & z$source_std==2), Ntot=nrow(z), N=sum(z$source_std==2),
                         y1=z[which(z$source_std==2),"pe1_event",drop=TRUE], 
                         y2=z[which(z$source_std==1),"pe1_event",drop=TRUE],
                         z1=as.integer(z[which(z$source_std==2),"az2",drop=TRUE]), 
                         a02=z[which(z$source_std==1),"fits",drop=TRUE],
                         zeros=rep(0,sum(z$source_std==1)))
  
  inits_list    <- list(alpha=0.5,
                        prior.alpha.mean=0,
                        d=0)
  
  jags.fit      <- jags.model(file='JAGS_model_A3.bug',data=data_list,inits=inits_list,n.chains=4,n.adapt=20000)
  jags.fit.out1 <- coda.samples(jags.fit, parameters, n.iter=120000, thin = 2)
  jags.coda     <- as.mcmc.list(jags.fit.out1)
  
  # Generating Results
  return(jags.coda)
})

# Combine 100 4-chain MCMCs into a single object
aug_pspp %<>% combine.mcmc(thin=2)
aug_pspp_chains_combined <- combine.mcmc(aug_pspp)
aug_pspp_table_less0<-table(aug_pspp_chains_combined<0)

# Final summaries
saveRDS(c(Mode=density(aug_pspp_chains_combined)$x[which.max(density(aug_pspp_chains_combined)$y)],
  Two_sided_P_value=2*min(ifelse(is.na(prop.table(aug_pspp_table_less0)["FALSE"]),0,prop.table(aug_pspp_table_less0)["FALSE"]),
                          ifelse(is.na(prop.table(aug_pspp_table_less0)["TRUE"]),0,prop.table(aug_pspp_table_less0)["TRUE"])),
  One_sided_P_Value=min(ifelse(is.na(prop.table(aug_pspp_table_less0)["FALSE"]),0,prop.table(aug_pspp_table_less0)["FALSE"]),
                        ifelse(is.na(prop.table(aug_pspp_table_less0)["TRUE"]),0,prop.table(aug_pspp_table_less0)["TRUE"])),
  Lower_CI=hdi(aug_pspp)[1],Upper_Ci=hdi(aug_pspp)[2]),"N69_psbpp.rds")

# 3.a [Pooling] PSBPP -----------------------------------------------------------

# Matching on Link-level fitted values using Mahalanobis distance
# Calculates the distance matrix
current_distance<-match_on(source2 ~ links,
                           controls=3,
                           data=trial[which(trial$source_std==1|trial$az2==TRUE),],
                           method="mahalanobis")

# Creates matched 1:1 pairs
pair_match_result<-pairmatch(source2 ~ links,
                             controls=3,
                             data=trial[which(trial$source_std==1|trial$az2==TRUE),],
                             distance=current_distance,
                             remove.unmatchables = FALSE)

# Adds the matches to the dataset
trial$pair<-NA
trial[which(trial$source_std==1|trial$az2==TRUE),"pair"]<-pair_match_result

# Subsets to trial with only Active trial participants and matched pairs
# trial<-rbind(trial[which(trial$E==1),],
#              trial[sample(which(trial$Hx!=1 & trial$E==0),size=ceiling(nrow(trial[which(trial$E==1),])/2)),],
#              trial[sample(trial[which(trial$E==1),"pair"],size=ceiling(nrow(trial[which(trial$E==1),])/2)),]
# )

trial<-rbind(trial[which(trial$az2==1),],
             trial[which(trial$source_std!=1 & trial$az2==0),],
             trial[sample(which(trial$source_std==1 & !is.na(trial$pair)),size=3*sum(trial$az2)),]
)
trial %<>% filter(!is.na(source_std))


# Performing the Bayesian analysis

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
# zeros[i]~dpois(phi[i])
# phi[i]<-  -l[i]+k
# l[i]<- a02[i]*(y2[i]*log(ilogit(alpha))+ (1-y2[i])*log((1-ilogit(alpha))))
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
data_list     <- list( N_E=sum(trial$az2==1), N_C=sum(trial$az2==0 & trial$source_std==2), Ntot=nrow(trial), N=sum(trial$source_std==2),
                       y1=trial[which(trial$source_std==2),"pe1_event",drop=TRUE], 
                       y2=trial[which(trial$source_std==1),"pe1_event",drop=TRUE],
                       z1=as.integer(trial[which(trial$source_std==2),"az2",drop=TRUE]), 
                       a02=trial[which(trial$source_std==1),"fits",drop=TRUE],
                       zeros=rep(0,sum(trial$source_std==1)))

inits_list    <- list(alpha=0.5,
                      prior.alpha.mean=0,
                      d=0)

jags.fit      <- jags.model(file='JAGS_model_A3.bug',data=data_list,inits=inits_list,n.chains=4,n.adapt=20000)
jags.fit.out1 <- coda.samples(jags.fit, parameters, n.iter=120000)
jags.coda     <- as.mcmc.list(jags.fit.out1)


# Generating Results
pool_pspp<-c(summary(jags.coda)$statistics,summary(jags.coda)$quantiles,
             Mode=density(jags.coda[[1]])$x[which.max(density(jags.coda[[1]])$y)],
             Two_sided_P_value=2*min(ifelse(is.na(prop.table(table(jags.coda[[1]]<0))["FALSE"]),0,prop.table(table(jags.coda[[1]]<0))["FALSE"]),
                                     ifelse(is.na(prop.table(table(jags.coda[[1]]<0))["TRUE"]),0,prop.table(table(jags.coda[[1]]<0))["TRUE"])),
             One_sided_P_Value=min(ifelse(is.na(prop.table(table(jags.coda[[1]]<0))["FALSE"]),0,prop.table(table(jags.coda[[1]]<0))["FALSE"]),
                                   ifelse(is.na(prop.table(table(jags.coda[[1]]<0))["TRUE"]),0,prop.table(table(jags.coda[[1]]<0))["TRUE"])),
             Lower_CI=hdi(jags.coda)[1],Upper_Ci=hdi(jags.coda)[2])

pool_pspp
exp(pool_pspp[c("Mean","Lower_CI","Upper_Ci")])

rm(list=c("current_distance","pair_match_result","parameters","data_list","inits_list","jags.fit","jags.fit.out1","jags.coda","pool_pspp"))

# 4. [Augmenting, Multi-chain] CM ------------------------------------------------------

### Reseting the Analytic dataset
# Saving new analytic dataset and calculating needed variables with imputation
trial<-propensity
# Generating a propensity score value for those with missing variables
links<-predict(multivariable)
trial$links<-links
rm(links)
# Simplified propensity score for those with missing values
multivariable_imp<-lrm(formula=source2 ~ agecat+sex_txt+NonWhite+genocat+dalfa+hs+asthma+PANUM,
                       data=propensity)
links<-predict(multivariable_imp)
fits<-predict(multivariable_imp,type="fitted")

trial[which(is.na(trial$fits)),"fits"]<-fits[which(is.na(trial$fits))]
trial[which(is.na(trial$links)),"links"]<-links[which(is.na(trial$links))]   
rm(fits)
rm(links)

trial %<>% filter(!is.na(fevppcat2_txt))
trial %<>% filter(!is.na(haz))
trial %<>% filter(!is.na(pabase))

# Creating the 100 subsets for N=18 OPTIMIZE placebo Arm
set.seed(2121956)
{aug_cm<-list((trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=18)))}
aug_cm %<>% lapply(function(z){rbind(z,(trial %>% filter(TRT!=6)))})

# Doing the Bayesian Analysis
aug_cm %<>% lapply(function(z){
  ### INPUTS to JAGS
  parameters    <- c("be")
  data_list     <- list(N=nrow(z),
                        y=z$pe1_event, 
                        study=ifelse(z$source_std==1,1,2),
                        x1=ifelse(z$agecat=="0-3",1,0),
                        x2=ifelse(z$agecat==">6 to 12",1,0),
                        x3=ifelse(z$sex_txt=="Male",1,0),
                        x4=ifelse(z$NonWhite==1,1,0),
                        x5=ifelse(z$genocat=="Delta F508 Heterozygous",1,0),
                        x6=ifelse(z$genocat=="Other",1,0),
                        x7=ifelse(z$genocat=="Unknown",1,0),
                        x8=ifelse(z$fevppcat2_txt=="<90% predicted",1,0),
                        x9=ifelse(z$fevppcat2_txt==">=90% predicted",1,0),
                        x10=z$haz,
                        x11=z$dalfa,
                        x12=z$hs,
                        x13=z$asthma,
                        x14=z$pabase,
                        x15=ifelse(z$PANUM==2,1,0),
                        E=as.integer(z$az2),
                        time=z$time)
  
  inits_list    <- list(pi=rep(1,16),
                        slab=rep(1,16),
                        b0=c(0,0),
                        b1=c(0,0),
                        b2=c(0,0),
                        b3=c(0,0),
                        b4=c(0,0),
                        b5=c(0,0),
                        b6=c(0,0),
                        b7=c(0,0),
                        b8=c(0,0),
                        b9=c(0,0),
                        b10=c(0,0),
                        b11=c(0,0),
                        b12=c(0,0),
                        b13=c(0,0),
                        b14=c(0,0),
                        b15=c(0,0),
                        be=0)
  
  jags.fit      <- jags.model(file='JAGS_model_Epic_Optimize_CP.bug',data=data_list,inits=inits_list,n.chains=3,n.adapt=20000)
  
  return(jags.fit)
  
})
parameters    <- c("be")
aug_cm_results<-lapply(aug_cm,function(z){
  jags.fit.out1 <- coda.samples(z, parameters, n.iter=120000,thin=2)
  jags.coda     <- as.mcmc.list(jags.fit.out1)
  return(jags.coda)
})

# Combine 100 4-chain MCMCs into a single object
aug_cm_results %<>% combine.mcmc(thin=2)
aug_cm_chains_combined <- combine.mcmc(aug_cm_results)
aug_cm_table_less0<-table(aug_cm_chains_combined<0)

# Final summaries
saveRDS(c(Mode=density(aug_cm_chains_combined)$x[which.max(density(aug_cm_chains_combined)$y)],
  Two_sided_P_value=2*min(ifelse(is.na(prop.table(aug_cm_table_less0)["FALSE"]),0,prop.table(aug_cm_table_less0)["FALSE"]),
                          ifelse(is.na(prop.table(aug_cm_table_less0)["TRUE"]),0,prop.table(aug_cm_table_less0)["TRUE"])),
  One_sided_P_Value=min(ifelse(is.na(prop.table(aug_cm_table_less0)["FALSE"]),0,prop.table(aug_cm_table_less0)["FALSE"]),
                        ifelse(is.na(prop.table(aug_cm_table_less0)["TRUE"]),0,prop.table(aug_cm_table_less0)["TRUE"])),
  Lower_CI=hdi(aug_cm_chains_combined)[1],Upper_Ci=hdi(aug_cm_chains_combined)[2]),"N18_cm.rds")


# Creating the 100 subsets for N=35 OPTIMIZE placebo Arm
set.seed(2121957)
{aug_cm<-list((trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=35)))}
aug_cm %<>% lapply(function(z){rbind(z,(trial %>% filter(TRT!=6)))})

# Doing the Bayesian Analysis
aug_cm %<>% lapply(function(z){
  ### INPUTS to JAGS
  parameters    <- c("be")
  data_list     <- list(N=nrow(z),
                        y=z$pe1_event, 
                        study=ifelse(z$source_std==1,1,2),
                        x1=ifelse(z$agecat=="0-3",1,0),
                        x2=ifelse(z$agecat==">6 to 12",1,0),
                        x3=ifelse(z$sex_txt=="Male",1,0),
                        x4=ifelse(z$NonWhite==1,1,0),
                        x5=ifelse(z$genocat=="Delta F508 Heterozygous",1,0),
                        x6=ifelse(z$genocat=="Other",1,0),
                        x7=ifelse(z$genocat=="Unknown",1,0),
                        x8=ifelse(z$fevppcat2_txt=="<90% predicted",1,0),
                        x9=ifelse(z$fevppcat2_txt==">=90% predicted",1,0),
                        x10=z$haz,
                        x11=z$dalfa,
                        x12=z$hs,
                        x13=z$asthma,
                        x14=z$pabase,
                        x15=ifelse(z$PANUM==2,1,0),
                        E=as.integer(z$az2),
                        time=z$time)
  
  inits_list    <- list(pi=rep(1,16),
                        slab=rep(1,16),
                        b0=c(0,0),
                        b1=c(0,0),
                        b2=c(0,0),
                        b3=c(0,0),
                        b4=c(0,0),
                        b5=c(0,0),
                        b6=c(0,0),
                        b7=c(0,0),
                        b8=c(0,0),
                        b9=c(0,0),
                        b10=c(0,0),
                        b11=c(0,0),
                        b12=c(0,0),
                        b13=c(0,0),
                        b14=c(0,0),
                        b15=c(0,0),
                        be=0)
  
  jags.fit      <- jags.model(file='JAGS_model_Epic_Optimize_CP.bug',data=data_list,inits=inits_list,n.chains=3,n.adapt=20000)
  
  return(jags.fit)
  
})
parameters    <- c("be")
aug_cm_results<-lapply(aug_cm,function(z){
  jags.fit.out1 <- coda.samples(z, parameters, n.iter=120000,thin=2)
  jags.coda     <- as.mcmc.list(jags.fit.out1)
  return(jags.coda)
})

# Combine 100 4-chain MCMCs into a single object
aug_cm_results %<>% combine.mcmc(thin=2)
aug_cm_chains_combined <- combine.mcmc(aug_cm_results)
aug_cm_table_less0<-table(aug_cm_chains_combined<0)

# Final summaries
saveRDS(c(Mode=density(aug_cm_chains_combined)$x[which.max(density(aug_cm_chains_combined)$y)],
          Two_sided_P_value=2*min(ifelse(is.na(prop.table(aug_cm_table_less0)["FALSE"]),0,prop.table(aug_cm_table_less0)["FALSE"]),
                                  ifelse(is.na(prop.table(aug_cm_table_less0)["TRUE"]),0,prop.table(aug_cm_table_less0)["TRUE"])),
          One_sided_P_Value=min(ifelse(is.na(prop.table(aug_cm_table_less0)["FALSE"]),0,prop.table(aug_cm_table_less0)["FALSE"]),
                                ifelse(is.na(prop.table(aug_cm_table_less0)["TRUE"]),0,prop.table(aug_cm_table_less0)["TRUE"])),
          Lower_CI=hdi(aug_cm_chains_combined)[1],Upper_Ci=hdi(aug_cm_chains_combined)[2]),"N35_cm.rds")

# Creating the 100 subsets for N=43 OPTIMIZE placebo Arm
set.seed(2121958)
{aug_cm<-list((trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=43)))}
aug_cm %<>% lapply(function(z){rbind(z,(trial %>% filter(TRT!=6)))})

# Doing the Bayesian Analysis
aug_cm %<>% lapply(function(z){
  ### INPUTS to JAGS
  parameters    <- c("be")
  data_list     <- list(N=nrow(z),
                        y=z$pe1_event, 
                        study=ifelse(z$source_std==1,1,2),
                        x1=ifelse(z$agecat=="0-3",1,0),
                        x2=ifelse(z$agecat==">6 to 12",1,0),
                        x3=ifelse(z$sex_txt=="Male",1,0),
                        x4=ifelse(z$NonWhite==1,1,0),
                        x5=ifelse(z$genocat=="Delta F508 Heterozygous",1,0),
                        x6=ifelse(z$genocat=="Other",1,0),
                        x7=ifelse(z$genocat=="Unknown",1,0),
                        x8=ifelse(z$fevppcat2_txt=="<90% predicted",1,0),
                        x9=ifelse(z$fevppcat2_txt==">=90% predicted",1,0),
                        x10=z$haz,
                        x11=z$dalfa,
                        x12=z$hs,
                        x13=z$asthma,
                        x14=z$pabase,
                        x15=ifelse(z$PANUM==2,1,0),
                        E=as.integer(z$az2),
                        time=z$time)
  
  inits_list    <- list(pi=rep(1,16),
                        slab=rep(1,16),
                        b0=c(0,0),
                        b1=c(0,0),
                        b2=c(0,0),
                        b3=c(0,0),
                        b4=c(0,0),
                        b5=c(0,0),
                        b6=c(0,0),
                        b7=c(0,0),
                        b8=c(0,0),
                        b9=c(0,0),
                        b10=c(0,0),
                        b11=c(0,0),
                        b12=c(0,0),
                        b13=c(0,0),
                        b14=c(0,0),
                        b15=c(0,0),
                        be=0)
  
  jags.fit      <- jags.model(file='JAGS_model_Epic_Optimize_CP.bug',data=data_list,inits=inits_list,n.chains=3,n.adapt=20000)
  
  return(jags.fit)
  
})
parameters    <- c("be")
aug_cm_results<-lapply(aug_cm,function(z){
  jags.fit.out1 <- coda.samples(z, parameters, n.iter=120000,thin=2)
  jags.coda     <- as.mcmc.list(jags.fit.out1)
  return(jags.coda)
})

# Combine 100 4-chain MCMCs into a single object
aug_cm_results %<>% combine.mcmc(thin=2)
aug_cm_chains_combined <- combine.mcmc(aug_cm_results)
aug_cm_table_less0<-table(aug_cm_chains_combined<0)

# Final summaries
saveRDS(c(Mode=density(aug_cm_chains_combined)$x[which.max(density(aug_cm_chains_combined)$y)],
          Two_sided_P_value=2*min(ifelse(is.na(prop.table(aug_cm_table_less0)["FALSE"]),0,prop.table(aug_cm_table_less0)["FALSE"]),
                                  ifelse(is.na(prop.table(aug_cm_table_less0)["TRUE"]),0,prop.table(aug_cm_table_less0)["TRUE"])),
          One_sided_P_Value=min(ifelse(is.na(prop.table(aug_cm_table_less0)["FALSE"]),0,prop.table(aug_cm_table_less0)["FALSE"]),
                                ifelse(is.na(prop.table(aug_cm_table_less0)["TRUE"]),0,prop.table(aug_cm_table_less0)["TRUE"])),
          Lower_CI=hdi(aug_cm_chains_combined)[1],Upper_Ci=hdi(aug_cm_chains_combined)[2]),"N43_cm.rds")

# Creating the 100 subsets for N=52 OPTIMIZE placebo Arm
set.seed(2121959)
{aug_cm<-list((trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=52)))}
aug_cm %<>% lapply(function(z){rbind(z,(trial %>% filter(TRT!=6)))})

# Doing the Bayesian Analysis
aug_cm %<>% lapply(function(z){
  ### INPUTS to JAGS
  parameters    <- c("be")
  data_list     <- list(N=nrow(z),
                        y=z$pe1_event, 
                        study=ifelse(z$source_std==1,1,2),
                        x1=ifelse(z$agecat=="0-3",1,0),
                        x2=ifelse(z$agecat==">6 to 12",1,0),
                        x3=ifelse(z$sex_txt=="Male",1,0),
                        x4=ifelse(z$NonWhite==1,1,0),
                        x5=ifelse(z$genocat=="Delta F508 Heterozygous",1,0),
                        x6=ifelse(z$genocat=="Other",1,0),
                        x7=ifelse(z$genocat=="Unknown",1,0),
                        x8=ifelse(z$fevppcat2_txt=="<90% predicted",1,0),
                        x9=ifelse(z$fevppcat2_txt==">=90% predicted",1,0),
                        x10=z$haz,
                        x11=z$dalfa,
                        x12=z$hs,
                        x13=z$asthma,
                        x14=z$pabase,
                        x15=ifelse(z$PANUM==2,1,0),
                        E=as.integer(z$az2),
                        time=z$time)
  
  inits_list    <- list(pi=rep(1,16),
                        slab=rep(1,16),
                        b0=c(0,0),
                        b1=c(0,0),
                        b2=c(0,0),
                        b3=c(0,0),
                        b4=c(0,0),
                        b5=c(0,0),
                        b6=c(0,0),
                        b7=c(0,0),
                        b8=c(0,0),
                        b9=c(0,0),
                        b10=c(0,0),
                        b11=c(0,0),
                        b12=c(0,0),
                        b13=c(0,0),
                        b14=c(0,0),
                        b15=c(0,0),
                        be=0)
  
  jags.fit      <- jags.model(file='JAGS_model_Epic_Optimize_CP.bug',data=data_list,inits=inits_list,n.chains=3,n.adapt=20000)
  
  return(jags.fit)
  
})
parameters    <- c("be")
aug_cm_results<-lapply(aug_cm,function(z){
  jags.fit.out1 <- coda.samples(z, parameters, n.iter=120000,thin=2)
  jags.coda     <- as.mcmc.list(jags.fit.out1)
  return(jags.coda)
})

# Combine 100 4-chain MCMCs into a single object
aug_cm_results %<>% combine.mcmc(thin=2)
aug_cm_chains_combined <- combine.mcmc(aug_cm_results)
aug_cm_table_less0<-table(aug_cm_chains_combined<0)

# Final summaries
saveRDS(c(Mode=density(aug_cm_chains_combined)$x[which.max(density(aug_cm_chains_combined)$y)],
          Two_sided_P_value=2*min(ifelse(is.na(prop.table(aug_cm_table_less0)["FALSE"]),0,prop.table(aug_cm_table_less0)["FALSE"]),
                                  ifelse(is.na(prop.table(aug_cm_table_less0)["TRUE"]),0,prop.table(aug_cm_table_less0)["TRUE"])),
          One_sided_P_Value=min(ifelse(is.na(prop.table(aug_cm_table_less0)["FALSE"]),0,prop.table(aug_cm_table_less0)["FALSE"]),
                                ifelse(is.na(prop.table(aug_cm_table_less0)["TRUE"]),0,prop.table(aug_cm_table_less0)["TRUE"])),
          Lower_CI=hdi(aug_cm_chains_combined)[1],Upper_Ci=hdi(aug_cm_chains_combined)[2]),"N52_cm.rds")

# Creating the 100 subsets for N=69 OPTIMIZE placebo Arm
set.seed(2121960)
{aug_cm<-list((trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)),
                (trial %>% filter(source_std==2) %>% filter(TRT==6) %>% slice_sample(n=69)))}
aug_cm %<>% lapply(function(z){rbind(z,(trial %>% filter(TRT!=6)))})

# Doing the Bayesian Analysis
aug_cm %<>% lapply(function(z){
  ### INPUTS to JAGS
  parameters    <- c("be")
  data_list     <- list(N=nrow(z),
                        y=z$pe1_event, 
                        study=ifelse(z$source_std==1,1,2),
                        x1=ifelse(z$agecat=="0-3",1,0),
                        x2=ifelse(z$agecat==">6 to 12",1,0),
                        x3=ifelse(z$sex_txt=="Male",1,0),
                        x4=ifelse(z$NonWhite==1,1,0),
                        x5=ifelse(z$genocat=="Delta F508 Heterozygous",1,0),
                        x6=ifelse(z$genocat=="Other",1,0),
                        x7=ifelse(z$genocat=="Unknown",1,0),
                        x8=ifelse(z$fevppcat2_txt=="<90% predicted",1,0),
                        x9=ifelse(z$fevppcat2_txt==">=90% predicted",1,0),
                        x10=z$haz,
                        x11=z$dalfa,
                        x12=z$hs,
                        x13=z$asthma,
                        x14=z$pabase,
                        x15=ifelse(z$PANUM==2,1,0),
                        E=as.integer(z$az2),
                        time=z$time)
  
  inits_list    <- list(pi=rep(1,16),
                        slab=rep(1,16),
                        b0=c(0,0),
                        b1=c(0,0),
                        b2=c(0,0),
                        b3=c(0,0),
                        b4=c(0,0),
                        b5=c(0,0),
                        b6=c(0,0),
                        b7=c(0,0),
                        b8=c(0,0),
                        b9=c(0,0),
                        b10=c(0,0),
                        b11=c(0,0),
                        b12=c(0,0),
                        b13=c(0,0),
                        b14=c(0,0),
                        b15=c(0,0),
                        be=0)
  
  jags.fit      <- jags.model(file='JAGS_model_Epic_Optimize_CP.bug',data=data_list,inits=inits_list,n.chains=3,n.adapt=20000)
  
  return(jags.fit)
  
})
parameters    <- c("be")
aug_cm_results<-lapply(aug_cm,function(z){
  jags.fit.out1 <- coda.samples(z, parameters, n.iter=120000,thin=2)
  jags.coda     <- as.mcmc.list(jags.fit.out1)
  return(jags.coda)
})

# Combine 100 4-chain MCMCs into a single object
aug_cm_results %<>% combine.mcmc(thin=2)
aug_cm_chains_combined <- combine.mcmc(aug_cm_results)
aug_cm_table_less0<-table(aug_cm_chains_combined<0)

# Final summaries
saveRDS(c(Mode=density(aug_cm_chains_combined)$x[which.max(density(aug_cm_chains_combined)$y)],
          Two_sided_P_value=2*min(ifelse(is.na(prop.table(aug_cm_table_less0)["FALSE"]),0,prop.table(aug_cm_table_less0)["FALSE"]),
                                  ifelse(is.na(prop.table(aug_cm_table_less0)["TRUE"]),0,prop.table(aug_cm_table_less0)["TRUE"])),
          One_sided_P_Value=min(ifelse(is.na(prop.table(aug_cm_table_less0)["FALSE"]),0,prop.table(aug_cm_table_less0)["FALSE"]),
                                ifelse(is.na(prop.table(aug_cm_table_less0)["TRUE"]),0,prop.table(aug_cm_table_less0)["TRUE"])),
          Lower_CI=hdi(aug_cm_chains_combined)[1],Upper_Ci=hdi(aug_cm_chains_combined)[2]),"N69_cm.rds")


# 4.a [Pooling] CM --------------------------------------------------------

### Reseting the Analytic dataset

# Saving new analytic dataset and calculating needed variables with imputation
trial<-propensity
# Generating a propensity score value for those with missing variables
links<-predict(multivariable)
trial$links<-links
rm(links)
# Simplified propensity score for those with missing values
multivariable_imp<-lrm(formula=source2 ~ agecat+sex_txt+NonWhite+genocat+dalfa+hs+asthma+PANUM,
                       data=propensity)
links<-predict(multivariable_imp)
fits<-predict(multivariable_imp,type="fitted")

trial[which(is.na(trial$fits)),"fits"]<-fits[which(is.na(trial$fits))]
trial[which(is.na(trial$links)),"links"]<-links[which(is.na(trial$links))]   
rm(fits)
rm(links)

trial %<>% filter(!is.na(fevppcat2_txt))
trial %<>% filter(!is.na(haz))
trial %<>% filter(!is.na(pabase))

# # ### Commensurate Prior Model
  modelstring <- "
  model{

  for (i in 1:N) {
  y[i] ~ dpois(lambda[i])
  log(lambda[i]) <- mu[i]
  mu[i] <- b0[study[i]] + be*E[i] + b1[study[i]]*x1[i] + b2[study[i]]*x2[i] + b3[study[i]]*x3[i] + b4[study[i]]*x4[i] + b5[study[i]]*x5[i] + b6[study[i]]*x6[i] + b7[study[i]]*x7[i] + b8[study[i]]*x8[i] + b9[study[i]]*x9[i] + b10[study[i]]*x10[i] + b11[study[i]]*x11[i] + b12[study[i]]*x12[i] + b13[study[i]]*x13[i] + b14[study[i]]*x14[i] + b15[study[i]]*x15[i] + log(time[i])
  }

  # Note for the study variable 1=historical 2=active
  # Therefore beta[1] are historical parameters and betas[2] are active trial parameters

  # Spike and slab hyper priors for tau
  pi[1] ~ dbern(0.5)
  pi[2] ~ dbern(0.5)
  pi[3] ~ dbern(0.5)
  pi[4] ~ dbern(0.5)
  pi[5] ~ dbern(0.5)
  pi[6] ~ dbern(0.5)
  pi[7] ~ dbern(0.5)
  pi[8] ~ dbern(0.5)
  pi[9] ~ dbern(0.5)
  pi[10] ~ dbern(0.5)
  pi[11] ~ dbern(0.5)
  pi[12] ~ dbern(0.5)
  pi[13] ~ dbern(0.5)
  pi[14] ~ dbern(0.5)
  pi[15] ~ dbern(0.5)
  pi[16] ~ dbern(0.5)


  slab[1] ~ dunif(0.005,7)
  slab[2] ~ dunif(0.005,9)
  slab[3] ~ dunif(0.005,4)
  slab[4] ~ dunif(0.005,19)
  slab[5] ~ dunif(0.005,4)
  slab[6] ~ dunif(0.005,17)
  slab[7] ~ dunif(0.005,7)
  slab[8] ~ dunif(0.005,2)
  slab[9] ~ dunif(0.005,3)
  slab[10] ~ dunif(0.005,3)
  slab[11] ~ dunif(0.005,64)
  slab[12] ~ dunif(0.005,18)
  slab[13] ~ dunif(0.005,3)
  slab[14] ~ dunif(0.005,13)
  slab[15] ~ dunif(0.005,19)
  slab[16] ~ dunif(0.005,13)


  # Tau now is distributed as a mixture of the spike (at 200) and slab
  # Even though tau is determinisitic (<-) it has the desired hyper prior distribution
  tau[1] <- slab[1]*pi[1] + 14*(1-pi[1])
  tau[2] <- slab[2]*pi[2] + 17*(1-pi[2])
  tau[3] <- slab[3]*pi[3] + 7*(1-pi[3])
  tau[4] <- slab[4]*pi[4] + 38*(1-pi[4])
  tau[5] <- slab[5]*pi[5] + 7*(1-pi[5])
  tau[6] <- slab[6]*pi[6] + 33*(1-pi[6])
  tau[7] <- slab[7]*pi[7] + 13*(1-pi[7])
  tau[8] <- slab[8]*pi[8] + 3*(1-pi[8])
  tau[9] <- slab[9]*pi[9] + 5*(1-pi[9])
  tau[10] <- slab[10]*pi[10] + 5*(1-pi[10])
  tau[11] <- slab[11]*pi[11] + 127*(1-pi[11])
  tau[12] <- slab[12]*pi[12] + 36*(1-pi[12])
  tau[13] <- slab[13]*pi[13] + 5*(1-pi[13])
  tau[14] <- slab[14]*pi[14] + 26*(1-pi[14])
  tau[15] <- slab[15]*pi[15] + 38*(1-pi[15])
  tau[16] <- slab[16]*pi[16] + 26*(1-pi[16])


  # Vague normal prior's for historical parameters and treatment parameter (Note: [1] is Hx and [2] is active trial)
  b0[1] ~ dnorm(0,1)
  b1[1] ~ dnorm(0,1)  
  b2[1] ~ dnorm(0,1)  
  b3[1] ~ dnorm(0,1)  # 5-31-23 Changing to weakly informative priors (dnorm(0,1)) per Gelman prior blog advice
  b4[1] ~ dnorm(0,1)
  b5[1] ~ dnorm(0,1)
  b6[1] ~ dnorm(0,1)
  b7[1] ~ dnorm(0,1)
  b8[1] ~ dnorm(0,1)
  b9[1] ~ dnorm(0,1)
  b10[1] ~ dnorm(0,1)
  b11[1] ~ dnorm(0,1)
  b12[1] ~ dnorm(0,1)
  b13[1] ~ dnorm(0,1)
  b14[1] ~ dnorm(0,1)
  b15[1] ~ dnorm(0,1)
  be ~ dnorm(0,1)

  # Normal distribution around historical beta with tau precision (Note: [1] is Hx and [2] is active trial)
  b0[2] ~ dnorm(b0[1],tau[1])
  b1[2] ~ dnorm(b1[1],tau[2])
  b2[2] ~ dnorm(b2[1],tau[3])
  b3[2] ~ dnorm(b3[1],tau[4])
  b4[2] ~ dnorm(b4[1],tau[5])
  b5[2] ~ dnorm(b5[1],tau[6])
  b6[2] ~ dnorm(b6[1],tau[7])
  b7[2] ~ dnorm(b7[1],tau[8])
  b8[2] ~ dnorm(b8[1],tau[9])
  b9[2] ~ dnorm(b9[1],tau[10])
  b10[2] ~ dnorm(b10[1],tau[11])
  b11[2] ~ dnorm(b11[1],tau[12])
  b12[2] ~ dnorm(b12[1],tau[13])
  b13[2] ~ dnorm(b13[1],tau[14])
  b14[2] ~ dnorm(b14[1],tau[15])
  b15[2] ~ dnorm(b15[1],tau[16])

  HR <- exp(be)
  }
#   "
#   model_file <- "JAGS_model_Epic_Optimize_CP.bug"
#   writeLines(modelstring,con=model_file)

### INPUTS to JAGS
parameters    <- c("be")
data_list     <- list(N=nrow(trial),
                      y=trial$pe1_event, 
                      study=ifelse(trial$source_std==1,1,2),
                      x1=ifelse(trial$agecat=="0-3",1,0),
                      x2=ifelse(trial$agecat==">6 to 12",1,0),
                      x3=ifelse(trial$sex_txt=="Male",1,0),
                      x4=ifelse(trial$NonWhite==1,1,0),
                      x5=ifelse(trial$genocat=="Delta F508 Heterozygous",1,0),
                      x6=ifelse(trial$genocat=="Other",1,0),
                      x7=ifelse(trial$genocat=="Unknown",1,0),
                      x8=ifelse(trial$fevppcat2_txt=="<90% predicted",1,0),
                      x9=ifelse(trial$fevppcat2_txt==">=90% predicted",1,0),
                      x10=trial$haz,
                      x11=trial$dalfa,
                      x12=trial$hs,
                      x13=trial$asthma,
                      x14=trial$pabase,
                      x15=ifelse(trial$PANUM==2,1,0),
                      E=as.integer(trial$az2),
                      time=trial$time)

inits_list    <- list(pi=rep(1,16),
                      slab=rep(1,16),
                      b0=c(0,0),
                      b1=c(0,0),
                      b2=c(0,0),
                      b3=c(0,0),
                      b4=c(0,0),
                      b5=c(0,0),
                      b6=c(0,0),
                      b7=c(0,0),
                      b8=c(0,0),
                      b9=c(0,0),
                      b10=c(0,0),
                      b11=c(0,0),
                      b12=c(0,0),
                      b13=c(0,0),
                      b14=c(0,0),
                      b15=c(0,0),
                      be=0)

jags.fit      <- jags.model(file='JAGS_model_Epic_Optimize_CP.bug',data=data_list,inits=inits_list,n.chains=4,n.adapt=20000)
jags.fit.out1 <- coda.samples(jags.fit, parameters, n.iter=120000)
jags.coda     <- as.mcmc.list(jags.fit.out1)

# par(mfrow=c(3,1))
# traceplot(jags.coda)
# gelman.plot(jags.coda) # converge at around 20000 iteration
# densplot(jags.coda)

c(summary(jags.coda)$statistics,summary(jags.coda)$quantiles,
  Mode=density(jags.coda[[1]])$x[which.max(density(jags.coda[[1]])$y)],
  Two_sided_P_value=2*min(ifelse(is.na(prop.table(table(jags.coda[[1]]<0))["FALSE"]),0,prop.table(table(jags.coda[[1]]<0))["FALSE"]),
                          ifelse(is.na(prop.table(table(jags.coda[[1]]<0))["TRUE"]),0,prop.table(table(jags.coda[[1]]<0))["TRUE"])),
  One_sided_P_Value=min(ifelse(is.na(prop.table(table(jags.coda[[1]]<0))["FALSE"]),0,prop.table(table(jags.coda[[1]]<0))["FALSE"]),
                        ifelse(is.na(prop.table(table(jags.coda[[1]]<0))["TRUE"]),0,prop.table(table(jags.coda[[1]]<0))["TRUE"])),
  Lower_CI=hdi(jags.coda)[1],Upper_Ci=hdi(jags.coda)[2])

# 5. [Optimize Only] IPW --------------------------------------------------
# No Weighting
hold<-glm(formula=pe1_event ~ AZ2+agecat_txt+offset(log(time)),
    data=left_join(outcome,dm,by=c("subject","TRT","source","source_std")) %>% filter(source_std==2) %>%
      mutate(AZ2=ifelse(TRT==5,1,0)),family = "poisson")

summary(hold)
hold %>% coef %>% exp %>% round(2)
hold %>% confint %>% exp %>% round(2)

bcanon_p(1:nrow(trial),10000,function(x,xdata){coef(glm(formula=pe1_event ~ AZ2+agecat_txt+offset(log(time)),
                                                       data=xdata[x,],
                                                       family = "poisson"))["AZ2"]},
         xdata=left_join(outcome,dm,by=c("subject","TRT","source","source_std")) %>% filter(source_std==2) %>%
           mutate(AZ2=ifelse(TRT==5,1,0)),alpha=c(0.025,0.975))[1:6]

# 6. [Optimize Only] PSBPP ------------------------------------------------
# Saving new analytic dataset and calculating needed variables with imputation
trial<-propensity

# Generating a propensity score value for those with missing variables
links<-predict(multivariable)
trial$links<-links
rm(links)

# Simplified propensity score for those with missing values
multivariable_imp<-lrm(formula=source2 ~ agecat+sex_txt+NonWhite+genocat+dalfa+hs+asthma+PANUM,
                       data=propensity)

links<-predict(multivariable_imp)
fits<-predict(multivariable_imp,type="fitted")

# trial[which(is.na(trial$fits)),] %>% dplyr::select(source,subject,agecat,sex_txt,NonWhite,genocat,fevppcat2_txt,haz,dalfa,hs,asthma,pabase,PANUM)
# trial[which(is.na(trial$fits)),"weights"]<-ipw_austin(trial[which(is.na(trial$fits)),"az2",drop=T],fits[which(is.na(trial$fits))])
trial[which(is.na(trial$fits)),"fits"]<-fits[which(is.na(trial$fits))]
trial[which(is.na(trial$links)),"links"]<-links[which(is.na(trial$links))]

rm(fits)
rm(links)

trial<-(trial %>% filter(source_std==2) %>%
          mutate(AZ2=ifelse(TRT==5,1,0)))

# Performing the Bayesian analysis

# ## Model with uniform prior for d
# modelstring <- "
# model{
# 
# for (i in 1:N) {
# y1[i] ~ dbern(p[i])
# logit(p[i]) <- alpha + d*z1[i]
# }
# 
# alpha ~ dnorm(prior.alpha.mean, prior.alpha.prec)
# 
# prior.alpha.mean ~ dnorm(0,0.001)
# prior.alpha.prec ~ dgamma(0.001, 0.001)
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
# model_file <- "JAGS_model_A3_opt_only.bug"
# writeLines(modelstring,con=model_file)

parameters    <- c("be")
data_list     <- list(N=nrow(trial),
                       y1=trial[which(trial$source_std==2),"pe1_event",drop=TRUE], 
                       z1=as.integer(trial[which(trial$source_std==2),"az2",drop=TRUE]))

inits_list    <- list(alpha=0.5,
                      prior.alpha.mean=0,
                      d=0)

jags.fit      <- jags.model(file='JAGS_model_A3_opt_only.bug',data=data_list,inits=inits_list,n.chains=4,n.adapt=20000)
jags.fit.out1 <- coda.samples(jags.fit, parameters, n.iter=120000)
jags.coda     <- as.mcmc.list(jags.fit.out1)

opt_pspp<-c(summary(jags.coda)$statistics,summary(jags.coda)$quantiles,
             Mode=density(jags.coda[[1]])$x[which.max(density(jags.coda[[1]])$y)],
             Two_sided_P_value=2*min(ifelse(is.na(prop.table(table(jags.coda[[1]]<0))["FALSE"]),0,prop.table(table(jags.coda[[1]]<0))["FALSE"]),
                                     ifelse(is.na(prop.table(table(jags.coda[[1]]<0))["TRUE"]),0,prop.table(table(jags.coda[[1]]<0))["TRUE"])),
             One_sided_P_Value=min(ifelse(is.na(prop.table(table(jags.coda[[1]]<0))["FALSE"]),0,prop.table(table(jags.coda[[1]]<0))["FALSE"]),
                                   ifelse(is.na(prop.table(table(jags.coda[[1]]<0))["TRUE"]),0,prop.table(table(jags.coda[[1]]<0))["TRUE"])),
             Lower_CI=hdi(jags.coda)[1],Upper_Ci=hdi(jags.coda)[2])

opt_pspp
exp(opt_pspp[c("Mean","Lower_CI","Upper_Ci")])

# 7. [Optimize Only] CM ---------------------------------------------------

### Reseting the Analytic dataset

# Saving new analytic dataset and calculating needed variables with imputation
trial<-propensity
# Generating a propensity score value for those with missing variables
links<-predict(multivariable)
trial$links<-links
rm(links)
# Simplified propensity score for those with missing values
multivariable_imp<-lrm(formula=source2 ~ agecat+sex_txt+NonWhite+genocat+dalfa+hs+asthma+PANUM,
                       data=propensity)
links<-predict(multivariable_imp)
fits<-predict(multivariable_imp,type="fitted")

trial[which(is.na(trial$fits)),"fits"]<-fits[which(is.na(trial$fits))]
trial[which(is.na(trial$links)),"links"]<-links[which(is.na(trial$links))]   
rm(fits)
rm(links)

trial %<>% filter(!is.na(fevppcat2_txt))
trial %<>% filter(!is.na(haz))
trial %<>% filter(!is.na(pabase))

trial<-(trial %>% filter(source_std==2) %>%
          mutate(AZ2=ifelse(TRT==5,1,0)))

# # ### Commensurate Prior Model
#   modelstring <- "
# model{
#   
#   for (i in 1:N) {
#     y[i] ~ dpois(lambda[i])
#     log(lambda[i]) <- mu[i]
#     mu[i] <- b0 + be*E[i] + b1*x1[i] + b2*x2[i] + b3*x3[i] + b4*x4[i] + b5*x5[i] + b6*x6[i] + b7*x7[i] + b8*x8[i] + b9*x9[i] + b10*x10[i] + b11*x11[i] + b12*x12[i] + b13*x13[i] + b14*x14[i] + b15*x15[i] + log(time[i])
#   }
#   
#   # Note for the study variable 1=historical 2=active
#   # Therefore beta[1] are historical parameters and betas[2] are active trial parameters
#   
#   # Weakly informatie normal prior's for historical parameters and treatment parameter (Note: [1] is Hx and [2] is active trial)
#   b0 ~ dnorm(0,1)
#   b1 ~ dnorm(0,1)  # Note that 1 is precision or the inverse of the variance
#   b2 ~ dnorm(0,1)  # Meaning the variance is 1 and standard deviation is 1
#   b3 ~ dnorm(0,1)
#   b4 ~ dnorm(0,1)
#   b5 ~ dnorm(0,1)
#   b6 ~ dnorm(0,1)
#   b7 ~ dnorm(0,1)
#   b8 ~ dnorm(0,1)
#   b9 ~ dnorm(0,1)
#   b10 ~ dnorm(0,1)
#   b11 ~ dnorm(0,1)
#   b12 ~ dnorm(0,1)
#   b13 ~ dnorm(0,1)
#   b14 ~ dnorm(0,1)
#   b15 ~ dnorm(0,1)
#   be ~ dnorm(0,1)
#   
#   HR <- exp(be)
# }
#   "
#   model_file <- "JAGS_model_Epic_Optimize_CP_opt_only.bug"
#   writeLines(modelstring,con=model_file)

### INPUTS to JAGS
parameters    <- c("be")
data_list     <- list(N=nrow(trial),
                      y=trial$pe1_event, 
                      x1=ifelse(trial$agecat=="0-3",1,0),
                      x2=ifelse(trial$agecat==">6 to 12",1,0),
                      x3=ifelse(trial$sex_txt=="Male",1,0),
                      x4=ifelse(trial$NonWhite==1,1,0),
                      x5=ifelse(trial$genocat=="Delta F508 Heterozygous",1,0),
                      x6=ifelse(trial$genocat=="Other",1,0),
                      x7=ifelse(trial$genocat=="Unknown",1,0),
                      x8=ifelse(trial$fevppcat2_txt=="<90% predicted",1,0),
                      x9=ifelse(trial$fevppcat2_txt==">=90% predicted",1,0),
                      x10=trial$haz,
                      x11=trial$dalfa,
                      x12=trial$hs,
                      x13=trial$asthma,
                      x14=trial$pabase,
                      x15=ifelse(trial$PANUM==2,1,0),
                      E=as.integer(trial$az2),
                      time=trial$time)

inits_list    <- list(b0=0,
                      b1=0,
                      b2=0,
                      b3=0,
                      b4=0,
                      b5=0,
                      b6=0,
                      b7=0,
                      b8=0,
                      b9=0,
                      b10=0,
                      b11=0,
                      b12=0,
                      b13=0,
                      b14=0,
                      b15=0,
                      be=0)

jags.fit      <- jags.model(file='JAGS_model_Epic_Optimize_CP_opt_only.bug',data=data_list,inits=inits_list,n.chains=4,n.adapt=20000)
jags.fit.out1 <- coda.samples(jags.fit, parameters, n.iter=120000)
jags.coda     <- as.mcmc.list(jags.fit.out1)

# par(mfrow=c(3,1))
# traceplot(jags.coda)
# gelman.plot(jags.coda) # converge at around 20000 iteration
# densplot(jags.coda)

c(summary(jags.coda)$statistics,summary(jags.coda)$quantiles,
  Mode=density(jags.coda[[1]])$x[which.max(density(jags.coda[[1]])$y)],
  Two_sided_P_value=2*min(ifelse(is.na(prop.table(table(jags.coda[[1]]<0))["FALSE"]),0,prop.table(table(jags.coda[[1]]<0))["FALSE"]),
                          ifelse(is.na(prop.table(table(jags.coda[[1]]<0))["TRUE"]),0,prop.table(table(jags.coda[[1]]<0))["TRUE"])),
  One_sided_P_Value=min(ifelse(is.na(prop.table(table(jags.coda[[1]]<0))["FALSE"]),0,prop.table(table(jags.coda[[1]]<0))["FALSE"]),
                        ifelse(is.na(prop.table(table(jags.coda[[1]]<0))["TRUE"]),0,prop.table(table(jags.coda[[1]]<0))["TRUE"])),
  Lower_CI=hdi(jags.coda)[1],Upper_Ci=hdi(jags.coda)[2])


    
