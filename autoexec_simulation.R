# Auto exec file for the AWS server

list.of.packages <- c("rstream","extraDistr","mvtnorm","ggplot2","magrittr","broom","MASS",
                      "RItools","tidyverse","HDInterval","lubridate","sandwich","rjags","R2jags","remotes")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library("rstream")
library("extraDistr")
library("mvtnorm")
library("ggplot2")
library("magrittr")
library("broom")
library("MASS")
library("RItools")
library("tidyverse")
library("HDInterval")
library("lubridate")
library("sandwich")
library("rjags")
library("R2jags")
library("remotes")

# An older version of the optmatch package is needed for compatibility with the AWS servers
remotes::install_url("https://cran.r-project.org/src/contrib/Archive/optmatch/optmatch_0.9-17.tar.gz")
library("optmatch")

rm(new.packages)
rm(list.of.packages)

