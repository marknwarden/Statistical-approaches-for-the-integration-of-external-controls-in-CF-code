# ================================================================================================
# Name: L2
# Purpose:  Loop over all parameter sets (114)     
# Programmer: Mark Warden
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Runs the whole simulation on the server using parallel programming.
#
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# L2 ----------------------------------------------------------------------
library(parallel)

# Setting work environment
# setwd("S:\\BCDM\\Statistical Research\\Projects\\Historical Controls\\Code\\Simulation Code\\Output")


# 1. Setting up the cluster -----------------------------------------------

# Making the cluster and loading all needed packages for each one
# 72 cores were available on the used server
cl <- makeCluster(72)

# Code sometimes needed to make sure JAGS worked correctly
# clusterCall(cl,function(){
#   Sys.setenv(JAGS_HOME="C:/Users/mwarden/AppData/Local/Programs/JAGS/JAGS-4.3.1")})


# Loads the required packages onto each of the socket clusters
clusterCall(cl, function(){
    needed_cluster_packages<-c("rstream","extraDistr","mvtnorm","ggplot2","magrittr","broom","MASS",
                             "RItools","tidyverse","HDInterval","optmatch","lubridate","sandwich")
    for (i in needed_cluster_packages){
    
    setwd("/home/rstudio/Dropbox/Dissertation/VM 1")
    library(i,character.only = TRUE)
    }
    return("Completed")})

clusterCall(cl,function(){
  library(rjags)
  library(R2jags)
  return("Completed")})

# Getting simulation written functions loaded to cluster
clusterExport(cl, varlist = c("D1","A0","A1","A2","A3","L1","N_power_calc"))

# Setting the RngStream to ensure simulation cluster independence
RNGkind("L'Ecuyer-CMRG")
clusterEvalQ(cl, { library(MASS); RNGkind("L'Ecuyer-CMRG") })
clusterEvalQ(cl, { library(MASS); RNGkind("L'Ecuyer-CMRG") })

# 2. Creating the parameter set data frame --------------------------------
# Input data frame

# Note: 0.7001 is placeholder's for thetaHR=1 for N calculation to be replaced later
param_set<-expand.grid(c(0.7001,0.5),c("A","C"),c("D","E"),c(0.2,0.6),c("F","G"),c(1,2))
names(param_set)<-c("thetaHR","delta","mu","sigma","betaHR","Hx_size")

# Fixing thetaHR after N calculation
param_set$N<-sapply(param_set$thetaHR,N_power_calc,F1=0.5,D=1,alpha=0.05,power=0.80,drop=0,phi=1)[1,]
param_set$thetaHR<-ifelse(param_set$thetaHR==0.7001,1,param_set$thetaHR)

# Setting HxN
param_set$HxN<-(ceiling((param_set$N/2))*param_set$Hx_size)

# Allocation Ratio
param_set$AR<-0.5

# Fixing delta
recoder <- Vectorize(vectorize.args = "a",
                 FUN = function(a) {
                   switch(as.character(a),
                          "A" = 0,
                          "B" = -0.2,
                          "C" = -0.4)})
param_set$delta_1<-recoder(as.character(param_set$delta))
recoder <- Vectorize(vectorize.args = "a",
                     FUN = function(a) {
                       switch(as.character(a),
                              "A" = 0,
                              "B" = 0.1,
                              "C" = 0.2)})
param_set$delta_2<-recoder(as.character(param_set$delta))

# Fixing mu
recoder <- Vectorize(vectorize.args = "a",
                     FUN = function(a) {
                       switch(as.character(a),
                              "D" = 0.75,
                              "E" = 0.85)})
param_set$mu_1<-recoder(as.character(param_set$mu))
recoder <- Vectorize(vectorize.args = "a",
                     FUN = function(a) {
                       switch(as.character(a),
                              "D" = 0.2,
                              "E" = 0.1)})
param_set$mu_2<-recoder(as.character(param_set$mu))

# Fixing betaHR
recoder <- Vectorize(vectorize.args = "a",
                     FUN = function(a) {
                       switch(as.character(a),
                              "F" = 1.2,
                              "G" = 1.1)})
param_set$betaHR_1<-recoder(as.character(param_set$betaHR))
recoder <- Vectorize(vectorize.args = "a",
                     FUN = function(a) {
                       switch(as.character(a),
                              "F" = 0.65,
                              "G" = 1.3)})
param_set$betaHR_2<-recoder(as.character(param_set$betaHR))

# Reordering data.frame
param_set<-param_set[,c("N","HxN","AR","thetaHR","delta_1","delta_2","mu_1","mu_2","sigma","betaHR_1","betaHR_2","Hx_size")]
param_set$file_set<-1:64

param_set<-data.frame(t(param_set))
param_set<-Reduce(cbind,list(param_set,param_set,param_set,param_set))

# Assignments were used to optimize computing time accross all cores (so all cores worked in tandem for a
# as long as possible in the computing process)
assignments<-readRDS("assignments_for_cluster_optimization.rds")

param_set<-param_set[,assignments]

# write.csv(param_set,file = "param_set.csv",row.names = F)

# 3. Setting up the random number generation ------------------------------

# This is based on the clusterSetRNGStream function from
# the parallel package, copyrighted by The R Core Team
getseeds <- function(ntasks, iseed) {
  RNGkind("L'Ecuyer-CMRG")
  set.seed(iseed)
  seeds <- vector("list", ntasks)
  seeds[[1]] <- .Random.seed
  for (i in seq_len(ntasks - 1)) {
    seeds[[i + 1]] <- nextRNGSubStream(seeds[[i]])
  }
  seeds
}

# VM_1
seeds<-getseeds(72,2201959)

# VM_2
# seeds<-getseeds(72,2121956)

# Assigning the seeds to the cluster

clusterApply(cl,seeds[1:72],function(x){assign(".Random.seed", x, envir=.GlobalEnv)})

# Creating node identifiers
clusterApply(cl,1:72, function(x){assign("node", x, envir=.GlobalEnv)})

# 4. Running the simulation and stopping the cluster ----------------------

# Running the cluster for the simulation

parLapply(cl, param_set, L1,nsim=450)

# Stopping the cluster
stopCluster(cl)
rm(cl)
Sys.sleep(5)
gc()
Sys.sleep(5)

system('sudo shutdown -h now', wait = FALSE)


