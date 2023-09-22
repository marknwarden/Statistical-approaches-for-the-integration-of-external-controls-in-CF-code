# ================================================================================================
# Name: L1
# Purpose:  Loop within 1 parameter set      
# Programmer: Mark Warden
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Input:
#       param_vector - a vector of simulation parameters that contain the following:
#           N - The number of active trial participants
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
#           Hx_relative_size - Controls the size of the historical control arm relative to the experimental
#                              arm (see A3)
#           file_set - an index for which parameter set is being considered
#       nsim - number of simulations to run
#       time_required_1 - the time in seconds for 1 simulation to run (used to print out estimates of progress)
# Output:
#       Output - writes to a csv file in the working directory that contains the following:
#   The Poisson Regression with inverse probability weighting model information:
#   "A1_Int","A1_E","A1_X1","A1_X2","A1_Int_naive_SE","A1_E_naive_SE","A1_X1_naive_SE","A1_X2_naive_SE","A1_Int_LCI","A1_E_LCI","A1_X1_LCI","A1_X2_LCI",
#   "A1_Int_UCI","A1_E_UCI","A1_X1_UCI","A1_X2_UCI","A1_Int_naive_pv","A1_E_naive_pv","A1_X1_naive_pv","A1_X2_naive_pv",
#   "A1_Int_robust_SE","A1_E_robust_SE","A1_X1_robust_SE","A1_X2_robust_SE","A1_Int_robust_LCI","A1_E_robust_LCI","A1_X1_robust_LCI","A1_X2_robust_LCI",
#   "A1_Int_robust_UCI","A1_E_robust_UCI","A1_X1_robust_UCI","A1_X2_robust_UCI","A1_Int_robust_pv","A1_E_robust_pv","A1_X1_robust_pv","A1_X2_robust_pv"
#
#   The commensurate prior model information:
#   "A2_E_mean","A2_E_SD","A2_E_NaiveSE","A2_Time_series_E_SE","A2_Posterior_2.5","A2_Posterior_25","A2_Posterior_50",
#   "A2_Posterior_75","A2_Posterior_97.5","A2_E_mode","A2_E_pv","A2_HDI_lower","A2_HDI_upper",
#
#   The propensity score based power prior model information
#   "A3_Posterior_E_mean","A3_Posterior_E_SD","A3_Posterior_E_NaiveSE","A3_Time_series_E_SE","A3_Posterior_2.5","A3_Posterior_25","A3_Posterior_50",
#   "A3_Posterior_75","A3_Posterior_97.5","A3_Posterior_E_mode","A3_E_pv","A3_HDI_lower","A3_HDI_upper",
#   
#   Simulation specific information for tracking and reproducibility (including the random seed)
#   "E","C","Hx","N",
#   "rseed1","rseed2","rseed3","rseed4","rseed5","rseed6","rseed7","node","procedure_time","completion_time"
# ================================================================================================
# Notes:
# 
# 
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# L1 ----------------------------------------------------------------------
library(lubridate)
library(HDInterval)
RNGkind("L'Ecuyer-CMRG")

L1<-function(param_vector,nsim=500,time_required_1=NA,...){
  N<-param_vector[1]
  HxN<-param_vector[2]
  AR<-param_vector[3]
  thetaHR<-param_vector[4]
  delta<-param_vector[c(5,6)]
  mu<-param_vector[c(7,8)]
  sigma<-matrix(c(1,param_vector[9],param_vector[9],1),2,2)
  betaHR<-param_vector[c(10,11)]
  Hx_relative_size<-param_vector[12]
  file_set<-param_vector[13]

  # Example of Result
# L1(as.numeric(structure(list(N = 572, HxN = 286, AR = 0.5, thetaHR = 1, delta_1 = 0,
#                   delta_2 = 0, mu_1 = 0.75, mu_2 = 0.2, sigma = 0.2, betaHR_1 = 1.2,
#                   betaHR_2 = 0.65, Hx_size = 1, file_set = 1L), row.names = 1L, class = "data.frame")),
#    nsim = 5,time_required_1 = 27.21)
    
# Step up the loop for step 2.
  
start_time <- Sys.time()

# Pre-allocating output matrix for all 5,000 runs
Output<-matrix(0,nrow=nsim,ncol=c(20+16+13+13+4+10))
colnames(Output)<-c("A1_Int","A1_E","A1_X1","A1_X2","A1_Int_naive_SE","A1_E_naive_SE","A1_X1_naive_SE","A1_X2_naive_SE","A1_Int_LCI","A1_E_LCI","A1_X1_LCI","A1_X2_LCI",
                    "A1_Int_UCI","A1_E_UCI","A1_X1_UCI","A1_X2_UCI","A1_Int_naive_pv","A1_E_naive_pv","A1_X1_naive_pv","A1_X2_naive_pv",
                    "A1_Int_robust_SE","A1_E_robust_SE","A1_X1_robust_SE","A1_X2_robust_SE","A1_Int_robust_LCI","A1_E_robust_LCI","A1_X1_robust_LCI","A1_X2_robust_LCI",
                    "A1_Int_robust_UCI","A1_E_robust_UCI","A1_X1_robust_UCI","A1_X2_robust_UCI","A1_Int_robust_pv","A1_E_robust_pv","A1_X1_robust_pv","A1_X2_robust_pv",
                    "A2_E_mean","A2_E_SD","A2_E_NaiveSE","A2_Time_series_E_SE","A2_Posterior_2.5","A2_Posterior_25","A2_Posterior_50",
                    "A2_Posterior_75","A2_Posterior_97.5","A2_E_mode","A2_E_pv","A2_HDI_lower","A2_HDI_upper",
                    "A3_Posterior_E_mean","A3_Posterior_E_SD","A3_Posterior_E_NaiveSE","A3_Time_series_E_SE","A3_Posterior_2.5","A3_Posterior_25","A3_Posterior_50",
                    "A3_Posterior_75","A3_Posterior_97.5","A3_Posterior_E_mode","A3_E_pv","A3_HDI_lower","A3_HDI_upper","E","C","Hx","N",
                    "rseed1","rseed2","rseed3","rseed4","rseed5","rseed6","rseed7","node","procedure_time","completion_time")

for (i in 1:nsim){
  # Step 1 - Within 1 dataset
  # Generates 1 dataset using D1. Conducts aim 1 analyses (A1-3).
  
  # Generates the data
  if (is.na(N)){
  trial<-D1(HxN=HxN,thetaHR=thetaHR,delta=delta,mu=mu,sigma=sigma,betaHR=betaHR,drop_c=0.5,dispersion=1)
  } else {
    trial<-D1(N=N,HxN=HxN,thetaHR=thetaHR,delta=delta,mu=mu,sigma=sigma,betaHR=betaHR,drop_c=0.5,dispersion=1)  
  }
  # Estimates the propensity score and calculates the weights.
  trial<-A0(trial,keep_link_fits=T)
  
  # Conducts Inverse probability weighting and saves estimates (20 outputs + 16 robust outputs)
  A1_out<-A1(trial,full_info = T)
  
  # Conducts Commensurate priors and saves estimates (13 outputs)
  A2_mcmcs<-A2(trial)
  A2_out<-c(summary(A2_mcmcs[[1]])$statistics,summary(A2_mcmcs[[1]])$quantiles,
            A2_E_mode=density(A2_mcmcs[[1]])$x[which.max(density(A2_mcmcs[[1]])$y)],
                  2*min(ifelse(is.na(prop.table(table(A2_mcmcs[[1]]<0))["FALSE"]),0,prop.table(table(A2_mcmcs[[1]]<0))["FALSE"]),
                        ifelse(is.na(prop.table(table(A2_mcmcs[[1]]<0))["TRUE"]),0,prop.table(table(A2_mcmcs[[1]]<0))["TRUE"])))
  names(A2_out)[11]<-"A2_E_pv"
  A2_out<-c(A2_out,hdi(A2_mcmcs[[1]])[1],hdi(A2_mcmcs[[1]])[2])
  names(A2_out)[12:13]<-c("A2_HDI_lower","A2_HDI_upper")
  
  # Conducts Propensity score-based power priors (13 outputs for one set, 26 if both Lin method's are run)
  
  #Expands HxN pool for A3
  pool<-D1(N=N,HxN=3*HxN,thetaHR=thetaHR,delta=delta,mu=mu,sigma=sigma,betaHR=betaHR)
  pool<-A0(pool,keep_link_fits=T)
  pool<-pool[which(pool$Hx==1),]
  
  A3_mcmcs<-A3(rbind(trial,pool),HxN_size=HxN,Hx_rel=Hx_relative_size)
  A3_out<-c(summary(A3_mcmcs)$statistics,summary(A3_mcmcs)$quantiles,
                  Mode=density(A3_mcmcs[[1]])$x[which.max(density(A3_mcmcs[[1]])$y)],
                  2*min(ifelse(is.na(prop.table(table(A3_mcmcs[[1]]<0))["FALSE"]),0,prop.table(table(A3_mcmcs[[1]]<0))["FALSE"]),
                        ifelse(is.na(prop.table(table(A3_mcmcs[[1]]<0))["TRUE"]),0,prop.table(table(A3_mcmcs[[1]]<0))["TRUE"])),
                  hdi(A3_mcmcs)[1],hdi(A3_mcmcs)[2])
  names(A3_out)[11]<-"pv"
  names(A3_out)[12:13]<-c("HDI_lower","HDI_upper")
  
  # Calculates process times
  
  end_time <- Sys.time()
  procedure_time<-as.numeric(as.duration(interval(start_time,end_time)),"seconds") 
  
# Step 2  - For 1 parameter set
# Runs step 1 nsim number of times and saves output into a data frame
  # Saving output (4 additional outputs)
  addition<-c(A1_out,A2_out,A3_out,sum(trial$E,na.rm=T),length(which(trial$E==0 & trial$Hx==0)),sum(trial$Hx,na.rm=T),nrow(trial))
  
  # Random seed, node, procedure time, and run time (10 additional outputs)
  addition<-c(addition,.Random.seed,node,as.numeric(procedure_time),as.numeric(Sys.time()))
  Output[i,]<-addition
  
  # End of For loop
  print("Simluations Completed:")
  print(paste0(i," of ",nsim," for parameter set ",file_set))
if (!is.na(time_required_1)) {print(paste0("Estimated time remaining is ",duration(time_required_1*(nsim-i)),"."))}
  }
rm(addition)
Output<-cbind(matrix(rep(param_vector[c(3:13)],nsim), ncol = 11,byrow=T),Output)
write.table(Output,file=paste0("parameter_set_",file_set,"_node_",node,".csv"),row.names=F,append=T,quote=F,col.names = F,sep = ",")
return("Done")
}