########################################################################
## Project: Alpha Part Variance
## Script purpose: Simulate a dairy cattle breeding programme
## Date: 2021-11-29
## Author: Thiago de Paula Oliveira
########################################################################
# SET OPTIONS ---------------------------------------
cat("SETTING OPTIONS... \n\n", sep = "")
#=======================================================================
# Packages
#=======================================================================
for(replicate in seq_len(30)){
  cat("Simulation:", replicate, "\n")
  library(tidyverse)
  library(AlphaSimR)
  #=======================================================================
  # Global Parameter
  #=======================================================================
  source("./simulation25Males/globalParameters.R")
  
  #=======================================================================
  # Creating founders
  #=======================================================================
  source("./simulation25Males/CreateFounders.R") # loading .RData
  #=======================================================================
  # Burnin
  #=======================================================================
  source("./simulation25Males/burnin.R")
  save.image(paste0("./Results/burnin25Males", replicate,".RData"))
  
  #=======================================================================
  # Scenario 1 - IS and Pheno Selection
  #=======================================================================
  source("./simulation25Males/Scenario1.R")
  #-----------------------------------------------------------------------
  # Save results
  #-----------------------------------------------------------------------
  saveRDS(RecSys, file=paste0("./Results/Scenario1_replicate25Males",replicate, ".rds"))
  saveRDS(RecVar, file=paste0("./Results/Scenario1Var_replicate25Males",replicate, ".rds"))
  
  #=======================================================================
  # Scenario 2 - IS and TBV Selection
  #=======================================================================
  source("./simulation25Males/Scenario2.R")
  #-----------------------------------------------------------------------
  # Save results
  #-----------------------------------------------------------------------
  saveRDS(RecSys, file=paste0("./Results/Scenario2_replicate25Males",replicate, ".rds"))
  saveRDS(RecVar, file=paste0("./Results/Scenario2Var_replicate25Males",replicate, ".rds"))
  rm(list = ls())
}

