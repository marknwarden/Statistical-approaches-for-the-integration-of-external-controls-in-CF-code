# ================================================================================================
# Name: D1
# Purpose:  Data generation function for simulation     
# Programmer: Mark Warden
# Date started: 2021_6_23
# Date completed: 2022_4_27    
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 
# Input: - No external files required
#        Arguments:
#           N - The number of active trial participants
#           F1 - The cumulative event rate in 1 year
#           HxN - The number of Historical Controls
#           AR - Allocation Ratio (currently the code only supports a ratio of 0.5)
#           thetaHR - the treatment effect on the Hazard Ratio Scale
#           delta - The difference in prevalence between the active trial and hx trial.
#                       I.e., delta=c(-0.1,0.2) means that if the active trial has prevalence
#                       equal to 75% and 25% for X1 and X2 respectively, then the historical trial
#                       has (0.75,0.25) + (-0.1,0.2) = (0.65,0.45), or 
#                       65% prevalence for X1 and 45% prevalence for X2.
#           mu - The prevalence of the two covariates int he ACTIVE trial population
#                 I.e., mu=c(0.75,0.25) would mean that the binary X1 has a prevalence of 75%
#                       in the active trial and the binary X2 has a prevalence of 25% in the
#                       active trial.
#           sigma - The covariance matrix for the latent variables that generate X1 and X2 (also the
#                   correlation matrix because the variances are set at 1).
#                   I.e., sigma=matrix(c(1,0.2,0.2,1),2,2) would correspond to the covariance/correlation
#                   matrix    { 1   0.2} which has variance of 1 for X1,X2 and covariance of 0.2. Because
#                             {0.2   1 } both varainces are 1, this is also the correlation matrix.
#                   Note: The covariance only roughly corresponds to scales that measure binary correlation.
#           betaHR - The covariate effects on the Hazard ratio scale (i.e. betaHR=c(1.1,0.8) would mean
#                     X1=1 is associated with a 1.1 HR and X2=1 is associated with a HR of 0.8)
#           drop_c - A percentage of the active controls to drop ranging from 0 to 1 (such as 0.75)
#           dispersion - A multiplier for the variance in the outcome generation (1=Poisson distribution)
#           
# Output:
# Dataset that includes the following variables
#               i - the recruitment order for the active trial participants (NA if hx participant)
#               E - Whether the participant recieved the treatment
#               Hx - Whether the participant is from the historical trial (1) or the active trial (0)
#               X1 - the value of the x1 covariate (1=yes, 0 = no)
#               X2 - the value of the x2 covariate (1=yes, 0 = no)
#               t - time to the event or the end of the study (max is t=1)
#               y - Whether the outcome occured (1) or the participant reached the end of the study (0)
#     
#
# ================================================================================================


# 0. D1 ----------------------------------------------
library(rstream)
library(extraDistr)
library(mvtnorm)

# Setting Working Directory
# setwd("S:\\BCDM\\Statistical Research\\Projects\\Historical Controls\\Code\\Simulation Code\\Output")

N_power_calc<-function(HR,F1,D=1,alpha=0.05,power=0.80,drop=0,phi=1){
  # HR is hazard ratio of the treatment
  # F1 is cumulative event rate in 1 year
  # D is study duration in years
  # Alpha is type 1 error
  # Power is type 2 error.
  # Assumes a constant hazard function and by extension a linear cumulative hazard
  # See Lachin's formula detailed in Fundamentals of Clinical Trials by
  # Friedman, 4th edition, "Sample Size Calculations for 'Time to Failure'" pg. 152
  lambda_C=(-log(1-F1))  # Control lambda with no covariate modification
  # cov_dist=t(t(c(1-mu[2],mu[2])))  %*%  c(1-mu[1],mu[1])
  # Correcting the control rate average using the knowledge of the covariate distribution
  # lambda_C=(cov_dist[1,1])*lambda_C*(1*1)+(cov_dist[1,2])*lambda_C*(betaHR[1])+(cov_dist[2,1])*lambda_C*(betaHR[2])+(cov_dist[2,2])*lambda_C*(betaHR[1]*betaHR[2])
  lambda_I=lambda_C*HR   # Intervention lambda with no covariate modification
  # Correcting the control rate average using the knowledge of the covariate distribution
  # lambda_I=(cov_dist[1,1])*lambda_I*(1*1)+(cov_dist[1,2])*lambda_I*(betaHR[1])+(cov_dist[2,1])*lambda_I*(betaHR[2])+(cov_dist[2,2])*lambda_I*(betaHR[1]*betaHR[2])
  phi_C = lambda_C^2 / (1-exp(-lambda_C*D))
  phi_I = lambda_I^2 / (1-exp(-lambda_I*D))
  Z_alpha = qnorm(1-alpha/2)
  Z_beta = qnorm(power)
  N=ceiling(2*((Z_alpha+Z_beta)^2)*(phi*phi_C+phi*phi_I)/((lambda_C-lambda_I)^2))
  N_wdrop=ceiling(N/(1-drop))
  Nevents<-ceiling((1-exp(-lambda_C*D))*(N_wdrop/2)*(1-drop)+(1-exp(-lambda_I*D))*(N_wdrop/2)*(1-drop))
  return(c(N=N,N_wdrop=N_wdrop,Nevents=Nevents))
}

D1<-function(N=NA,F1=0.5,HxN,AR=0.5,thetaHR,delta,mu,sigma,betaHR,drop_c=0,dispersion=1, ...){
# 1. Generate i, Hx, and E ---------------------------------------
  # Generates an order of recruitment (i), treatment indicator (E), and
  # Historical participant indicator variable (Hx) when given the number of
  # total data rows (N), the number of Historical data rows (HxN), the
  # allocation ratio (AR, defaults to 0.5) for the active data rows, and the
  # random number seed (Rng).

# Calculates N based on sample size/power calculation function
if(is.na(N) & AR==0.5){N=N_power_calc(HR=thetaHR,F1,phi=dispersion,...)[2]}
if(is.na(N) & AR!=0.5){stop("The N_power_calc function is currently not able to calculate N if AR not equal to 0.5")}

# Checks if HxN is 0 and if it is handles the situation.
if (HxN==0){NoHx=1;
  HxN=1} else {NoHx=0}

trial<-data.frame(i=c(rep(NA,HxN),1:N),  # Note that order of recruitment is missing (NA) for historical controls
                    E=c(rep(0,HxN),rbern(N,prob=AR)), # Note that Exposure is bernouli.
                    Hx=c(rep(1,HxN),rep(0,N)))

# 2. Generate covariates --------------------------------------------------
  # Generates covariates (X_1,X_2) for each data row in (i,E,Hx). The simulation
  # characteristics and the random number seed (Rng). Uses a multivariate
  # normal distribution and probit model as the data generating mechanisms.


# Multivaraible normal with covariance set by sigma
covariates<-rbind(rmvnorm(HxN,mean=qnorm(mu+delta),sigma=sigma),
                    rmvnorm(N,mean=qnorm(mu),sigma=sigma))

# Converted to binary variables with p=mu.
covariates<-data.frame((covariates>0)*1)
names(covariates)<-c("X1","X2")

# Appending to trial data
trial<-data.frame(trial,covariates)
rm(covariates)

# 3. Generate outcomes ----------------------------------------------------
  # Generates outcome (t_i,y_i) for each data row. Requires (i,E,Hx,X_1,X_2).
  # It also requires specification of the over dispersion parameter,
  # study duration (D), but these are constant throughout the entire
  # simulation.

# Specification of b0 (Beta 0), D (study duration), and phi (over dispersion parameter)
lambda <-(-log(1-F1))   # This corresponds the hazard rate with a proportion experiencing the event (default to 0.5)
D<-1               # Study duration is set at 1 year for simplicity.

b0 <- rep(log(lambda),N+HxN)

# Phi is set to 1 by default
phi<-dispersion

r<-as.matrix(cbind(base=b0,trial[,c(4:5,2)])) %*% matrix(c(1,log(betaHR),log(thetaHR)),4,1)
r<-as.vector(r) # r is the modeled rate

# We originally wanted to consider dispersed data simulated using a random gamma distribution.
# However, we decided to instead set phi=1 such that the random gamma distribution simplifies to an
# exponential distribution. We left the random gamma distribution in in case we ever wanted to simulate
# overdispersed data.
t<-sapply(r,function(x){rgamma(1,shape=1/phi,scale=phi/exp(x))})
y<-(t<D)*1
t<-ifelse(t<=D,t,D)

trial<-data.frame(trial,t=t,y=y)

if (NoHx==1) {trial<-trial[which(trial$Hx==0),]}

if (drop_c!=0) {drop_controls<-sample(which(trial$Hx==0 & trial$E==0),size=(ceiling(drop_c*length(which(trial$Hx==0 & trial$E==0)))))
  trial<-trial[-drop_controls,]
  trial[which(trial$Hx==0),"i"]<-1:(length(which(trial$Hx==0)))
  }

return(trial)
}

