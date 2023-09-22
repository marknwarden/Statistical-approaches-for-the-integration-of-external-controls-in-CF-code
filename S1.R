# ================================================================================================
# Name: S1
# Purpose:  Summarize output from one Loop within 1 parameter set      
# Programmer: Mark Warden
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Changes and Edits:
# 
# Input: The output of L2
#
# Output: The summarized results in terms of performance measures
#
# ================================================================================================
# Notes:
#   
# 
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


S1<-function(output,...){

if (!all(output$thetaHR==output$thetaHR[1])) {stop("A single thetaHR should be present in one parameter set. Detecting more than 1.")}

thetaHR<-output$thetaHR[1]

# 1. Calculating Bias (and percent bias) ----------------------------------

# A1
A1_E_bar<-mean(output[,"A1_E"],na.rm=T)
A1_Bias<-(A1_E_bar - log(thetaHR))
A1_Perc_Bias<-ifelse(!is.infinite((A1_E_bar - log(thetaHR))/abs(log(thetaHR)))*100,
                     ((A1_E_bar - log(thetaHR))/abs(log(thetaHR)))*100,((exp(A1_E_bar) - thetaHR)/abs(thetaHR))*100)

# A2
A2_E_bar<-mean(output[,"A2_E_mode"],na.rm=T)
A2_Bias<-(A2_E_bar - log(thetaHR))
A2_Perc_Bias<-ifelse(!is.infinite((A2_E_bar - log(thetaHR))/abs(log(thetaHR))*100),
                     (A2_E_bar - log(thetaHR))/abs(log(thetaHR))*100,(exp(A2_E_bar) - thetaHR)/abs(thetaHR)*100)

# A3
A3_E_bar<-mean(output[,"A3_Posterior_E_mode"],na.rm=T)
A3_Bias<-(A3_E_bar - log(thetaHR))
A3_Perc_Bias<-ifelse(!is.infinite(((A3_E_bar - log(thetaHR))/abs(log(thetaHR)))*100),
                     ((A3_E_bar - log(thetaHR))/abs(log(thetaHR)))*100,((exp(A3_E_bar) - thetaHR)/abs(thetaHR))*100)

# 2. Calculating Empirical Standard Error ---------------------------------

# A1
A1_E_EmpSE<-sqrt((1/(sum(!is.na(output[,"A1_E"]))-1))*sum((output[,"A1_E"]-A1_E_bar)^2,na.rm=T))

# A2
A2_E_EmpSE<-sqrt((1/(sum(!is.na(output[,"A2_E_mode"]))-1))*sum((output[,"A2_E_mode"]-A2_E_bar)^2,na.rm=T))

# A3
A3_E_EmpSE<-sqrt((1/(sum(!is.na(output[,"A3_Posterior_E_mode"]))-1))*sum((output[,"A3_Posterior_E_mode"]-A3_E_bar)^2,na.rm=T))

# 3. Calculating Accuracy (MSE) -------------------------------------------

# A1
A1_MSE<-(1/sum(!is.na(output[,"A1_E"])))*sum((output[,"A1_E"]-log(thetaHR))^2,na.rm=T)

# A2
A2_MSE<-(1/sum(!is.na(output[,"A2_E_mode"])))*sum((output[,"A2_E_mode"]-log(thetaHR))^2,na.rm=T)

# A3
A3_MSE<-(1/sum(!is.na(output[,"A3_Posterior_E_mode"])))*sum((output[,"A3_Posterior_E_mode"]-log(thetaHR))^2,na.rm=T)

# 4. Calculating Coverage and average CI length -------------------------------------------------

# A1
A1_coverage<-(1/sum(!is.na(output[,"A1_E_robust_LCI"])))*sum(((output[,"A1_E_robust_LCI"]<log(thetaHR))&(output[,"A1_E_robust_UCI"]>log(thetaHR))),na.rm=T)
A1_Avg_CI_length<-(1/sum(!is.na(output[,"A1_E_robust_LCI"])))*sum(output[,"A1_E_robust_UCI"]-output[,"A1_E_robust_LCI"],na.rm=T)

# A2
A2_coverage<-(1/sum(!is.na(output[,"A2_HDI_lower"])))*sum(((output[,"A2_HDI_lower"]<log(thetaHR))&(output[,"A2_HDI_upper"]>log(thetaHR))),na.rm=T)
A2_Avg_CI_length<-(1/sum(!is.na(output[,"A2_HDI_lower"])))*sum(output[,"A2_HDI_upper"]-output[,"A2_HDI_lower"],na.rm=T)

# A3
A3_coverage<-(1/sum(!is.na(output[,"A3_HDI_lower"])))*sum(((output[,"A3_HDI_lower"]<log(thetaHR))&(output[,"A3_HDI_upper"]>log(thetaHR))),na.rm=T)
A3_Avg_CI_length<-(1/sum(!is.na(output[,"A3_HDI_lower"])))*sum(output[,"A3_HDI_upper"]-output[,"A3_HDI_lower"],na.rm=T)


# 5. Rejection Percentage --------------------------------------------------

# A1
A1_rejection<-(1/sum(!is.na(output[,"A1_E_robust_pv"])))*sum((output[,"A1_E_robust_pv"]<0.05),na.rm=T)

# A2
A2_rejection<-(1/sum(!is.na(output[,"A2_E_pv"])))*sum((output[,"A2_E_pv"]<0.05),na.rm=T)

# A3
A3_rejection<-(1/sum(!is.na(output[,"A3_E_pv"])))*sum((output[,"A3_E_pv"]<0.05),na.rm=T)

# 6. Monte Carlo bias SE of simulation estimates -------------------------

# A1
A1_MC_SE<-sqrt((1/(sum(!is.na(output[,"A1_E"]))*(sum(!is.na(output[,"A1_E"]))-1)))*sum((output[,"A1_E"]-A1_E_bar)^2,na.rm=T))

# A2
A2_MC_SE<-sqrt((1/(sum(!is.na(output[,"A2_E_mode"]))*(sum(!is.na(output[,"A2_E_mode"]))-1)))*sum((output[,"A2_E_mode"]-A2_E_bar)^2,na.rm=T))

# A3
A3_MC_SE<-sqrt((1/(sum(!is.na(output[,"A3_Posterior_E_mode"]))*(sum(!is.na(output[,"A3_Posterior_E_mode"]))-1)))*sum((output[,"A3_Posterior_E_mode"]-A3_E_bar)^2,na.rm=T))


# 7. Monte Carlo Coverage SE ----------------------------------------------
MC_Coverage_SE<-1.96*sqrt((0.8*0.2)/min(c(sum(!is.na(output[,"A1_E"])),
          sum(!is.na(output[,"A2_E_mode"])),
          sum(!is.na(output[,"A3_Posterior_E_mode"]))
      ),na.rm=T))


# 8. Monte Carlo Average CI SE --------------------------------------------
# A1
A1_AvgCI_MC_SE<-sqrt((1/(sum(!is.na(output[,"A1_E_robust_LCI"]))*(sum(!is.na(output[,"A1_E_robust_LCI"]))-1)))*sum(((output[,"A1_E_robust_UCI"]-output[,"A1_E_robust_LCI"])-A1_Avg_CI_length)^2,na.rm=T))

# A2
A2_AvgCI_MC_SE<-sqrt((1/(sum(!is.na(output[,"A2_HDI_lower"]))*(sum(!is.na(output[,"A2_HDI_lower"]))-1)))*sum(((output[,"A2_HDI_upper"]-output[,"A2_HDI_lower"])-A2_Avg_CI_length)^2,na.rm=T))


# A3
A3_AvgCI_MC_SE<-sqrt((1/(sum(!is.na(output[,"A3_HDI_lower"]))*(sum(!is.na(output[,"A3_HDI_lower"]))-1)))*sum(((output[,"A3_HDI_upper"]-output[,"A3_HDI_lower"])-A3_Avg_CI_length)^2,na.rm=T))


# 9. Monte Carlo EmpSE SE ------------------------------------------------
# A1
A1_EmpSE_MC_SE<-A1_E_EmpSE/sqrt(2*(sum(!is.na(output[,"A1_E"]))-1))

# A2
A2_EmpSE_MC_SE<-A2_E_EmpSE/sqrt(2*(sum(!is.na(output[,"A2_E_mode"]))-1))

# A3
A3_EmpSE_MC_SE<-A3_E_EmpSE/sqrt(2*(sum(!is.na(output[,"A3_Posterior_E_mode"]))-1))

# 10. Monte Carlo MSE SE --------------------------------------------------
# A1
A1_MSE_MC_SE<-sqrt((1/((sum(!is.na(output[,"A1_E"]))-1)*sum(!is.na(output[,"A1_E"]))))*sum(((output[,"A1_E"]-log(thetaHR))^2 - A1_MSE)^2,na.rm=T))

# A2
A2_MSE_MC_SE<-sqrt((1/((sum(!is.na(output[,"A2_E_mode"]))-1)*sum(!is.na(output[,"A2_E_mode"]))))*sum(((output[,"A2_E_mode"]-log(thetaHR))^2  - A2_MSE)^2,na.rm=T))

# A3
A3_MSE_MC_SE<-sqrt((1/((sum(!is.na(output[,"A3_Posterior_E_mode"]))-1)*sum(!is.na(output[,"A3_Posterior_E_mode"]))))*sum(((output[,"A3_Posterior_E_mode"]-log(thetaHR))^2  - A3_MSE)^2,na.rm=T))

# 11. Producing the summary output -----------------------------------------

result<-data.frame(cbind(cbind(cbind(c(A1_E_bar,A1_Bias,A1_Perc_Bias,A1_E_EmpSE,A1_MSE,A1_coverage,A1_Avg_CI_length,A1_rejection,A1_MC_SE,MC_Coverage_SE,A1_AvgCI_MC_SE,A1_EmpSE_MC_SE,A1_MSE_MC_SE),
              c(A2_E_bar,A2_Bias,A2_Perc_Bias,A2_E_EmpSE,A2_MSE,A2_coverage,A2_Avg_CI_length,A2_rejection,A2_MC_SE,MC_Coverage_SE,A2_AvgCI_MC_SE,A2_EmpSE_MC_SE,A2_MSE_MC_SE)),
              c(A3_E_bar,A3_Bias,A3_Perc_Bias,A3_E_EmpSE,A3_MSE,A3_coverage,A3_Avg_CI_length,A3_rejection,A3_MC_SE,MC_Coverage_SE,A3_AvgCI_MC_SE,A3_EmpSE_MC_SE,A3_MSE_MC_SE))),
              row.names = c("E_bar","Bias","Percent_Bias","Emperical_SE","Mean_Squared_Error","Coverage",
                            "Average_CI_length","Rejection_percentage","Monte_Carlo_bias_SE","Monte_Carlo_Coverage_SE","MC_SE_AvgCI","MC_SE_EmpSE","MC_SE_MSE"))
names(result)<-c("A1","A2","A3")

result<-rbind(setNames(data.frame(t(output[1:3,c("file","thetaHR","delta_1","delta_2","mu_1","mu_2","sigma","betaHR_1","betaHR_2","Hx_relative_size")])),
               c("A1","A2","A3")),
              set_rownames(setNames(data.frame(rbind(rep(ceiling(mean(output$N)),3),
                                        rep(ceiling(mean(output$E)),3),
                                        rep(ceiling(mean(output$C)),3),
                                        rep(ceiling(mean(output$Hx)),3))),
                       c("A1","A2","A3")),c("N_tot","E","C","Hx")),
      result
      )

return(result)

}
