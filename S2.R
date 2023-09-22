# ================================================================================================
# Name: S2
# Purpose:  Summarize output from all parameter sets    
# Programmer: Mark Warden
# 
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
library(stringr)

# Read all parameter sets into the workspace and name them correctly
parameter_sets<-mapply(function(x,name){setNames(read.csv(x,header = F),
                                        c("AR","thetaHR","delta_1","delta_2","mu_1","mu_2","sigma","betaHR_1","betaHR_2","Hx_relative_size","file",
                                          "A1_Int","A1_E","A1_X1","A1_X2","A1_Int_naive_SE","A1_E_naive_SE","A1_X1_naive_SE","A1_X2_naive_SE","A1_Int_LCI","A1_E_LCI","A1_X1_LCI","A1_X2_LCI",
                                          "A1_Int_UCI","A1_E_UCI","A1_X1_UCI","A1_X2_UCI","A1_Int_naive_pv","A1_E_naive_pv","A1_X1_naive_pv","A1_X2_naive_pv",
                                          "A1_Int_robust_SE","A1_E_robust_SE","A1_X1_robust_SE","A1_X2_robust_SE","A1_Int_robust_LCI","A1_E_robust_LCI","A1_X1_robust_LCI","A1_X2_robust_LCI",
                                          "A1_Int_robust_UCI","A1_E_robust_UCI","A1_X1_robust_UCI","A1_X2_robust_UCI","A1_Int_robust_pv","A1_E_robust_pv","A1_X1_robust_pv","A1_X2_robust_pv",
                                          "A2_E_mean","A2_E_SD","A2_E_NaiveSE","A2_Time_series_E_SE","A2_Posterior_2.5","A2_Posterior_25","A2_Posterior_50",
                                          "A2_Posterior_75","A2_Posterior_97.5","A2_E_mode","A2_E_pv","A2_HDI_lower","A2_HDI_upper",
                                          "A3_Posterior_E_mean","A3_Posterior_E_SD","A3_Posterior_E_NaiveSE","A3_Time_series_E_SE","A3_Posterior_2.5","A3_Posterior_25","A3_Posterior_50",
                                          "A3_Posterior_75","A3_Posterior_97.5","A3_Posterior_E_mode","A3_E_pv","A3_HDI_lower","A3_HDI_upper","E","C","Hx","N",
                                          "rseed1","rseed2","rseed3","rseed4","rseed5","rseed6","rseed7","node","procedure_time","completion_time"))
                                        },
       x=list.files(path="Output/Parameter set results",full.names=T),
       name=str_remove_all(list.files(path="Output/Parameter set results"),pattern=".csv"),
       SIMPLIFY = FALSE)

# Renaming parameter set list
names(parameter_sets)<-str_remove_all(str_remove_all(names(parameter_sets),pattern="Output/Parameter set results/"),
                                      pattern='.csv')
# Reordering parameter set list
parameter_sets<-parameter_sets[order(as.integer(str_remove_all(names(parameter_sets),pattern="parameter_set_")))]

parameter_sets<-lapply(parameter_sets,function(x){x<-x[,c("file","thetaHR","delta_1","delta_2","mu_1","mu_2","sigma","betaHR_1","betaHR_2","Hx_relative_size","N","E","C","Hx","AR",
                                                   "A1_Int","A1_E","A1_X1","A1_X2","A1_Int_naive_SE","A1_E_naive_SE","A1_X1_naive_SE","A1_X2_naive_SE","A1_Int_LCI","A1_E_LCI","A1_X1_LCI","A1_X2_LCI",
                                                   "A1_Int_UCI","A1_E_UCI","A1_X1_UCI","A1_X2_UCI","A1_Int_naive_pv","A1_E_naive_pv","A1_X1_naive_pv","A1_X2_naive_pv",
                                                   "A1_Int_robust_SE","A1_E_robust_SE","A1_X1_robust_SE","A1_X2_robust_SE","A1_Int_robust_LCI","A1_E_robust_LCI","A1_X1_robust_LCI","A1_X2_robust_LCI",
                                                   "A1_Int_robust_UCI","A1_E_robust_UCI","A1_X1_robust_UCI","A1_X2_robust_UCI","A1_Int_robust_pv","A1_E_robust_pv","A1_X1_robust_pv","A1_X2_robust_pv",
                                                   "A2_E_mean","A2_E_SD","A2_E_NaiveSE","A2_Time_series_E_SE","A2_Posterior_2.5","A2_Posterior_25","A2_Posterior_50",
                                                   "A2_Posterior_75","A2_Posterior_97.5","A2_E_mode","A2_E_pv","A2_HDI_lower","A2_HDI_upper",
                                                   "A3_Posterior_E_mean","A3_Posterior_E_SD","A3_Posterior_E_NaiveSE","A3_Time_series_E_SE","A3_Posterior_2.5","A3_Posterior_25","A3_Posterior_50",
                                                   "A3_Posterior_75","A3_Posterior_97.5","A3_Posterior_E_mode","A3_E_pv","A3_HDI_lower","A3_HDI_upper",
                                                   "rseed1","rseed2","rseed3","rseed4","rseed5","rseed6","rseed7","node","procedure_time","completion_time")]})

# final<-Reduce(rbind,parameter_sets)
parameter_set_summaries<-lapply(parameter_sets,S1)

# parameter_set_summaries<-lapply(parameter_sets,\(x){tryCatch(S1(x),error = function(e) NULL)})
# which(sapply(parameter_set_summaries,is.null)==TRUE)

# Extracting and saving A1 results
write.csv(A1_summary_results<-data.frame(Reduce(rbind,lapply(parameter_set_summaries,function(x){setNames(x[,"A1"],rownames(x))})),
                                         row.names = NULL),
          "A1_summary_results.csv",
          row.names = F)

write.csv(A2_summary_results<-data.frame(Reduce(rbind,lapply(parameter_set_summaries,function(x){setNames(x[,"A2"],rownames(x))})),row.names = NULL),
          "A2_summary_results.csv",
          row.names = F)

write.csv(A3_summary_results<-data.frame(Reduce(rbind,lapply(parameter_set_summaries,function(x){setNames(x[,"A3"],rownames(x))})),row.names = NULL),
          "A3_summary_results.csv",
          row.names = F)

