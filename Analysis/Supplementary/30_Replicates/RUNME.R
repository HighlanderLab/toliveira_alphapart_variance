########################################################################
## Project: toliveira_alphapart_variance
## Script purpose: Analysis considering 30 replicates of the simulation
## Date: 2022-03-18
## Author: Thiago de Paula Oliveira
########################################################################
for(i in seq_len(30)){
  r <- i
  #=======================================================================
  # Pheno Scenario
  #=======================================================================
  source("./Analysis/Supplementary/30_Replicates/gibbs1f90Pheno.R")
  source("./Analysis/Supplementary/30_Replicates/AlphaPart_Gibbs_Pheno_Validation.R")
  setwd("./Analysis/Supplementary/30_Replicates/ModelOutput/")
  cat("Removing old files", "\n")
  file.remove(list.files())
  setwd("./../../../../")
  rm(list = ls()[!(ls() %in% c('i','r'))])
  gc()
  #=======================================================================
  # TBV scenario
  #=======================================================================
  source("./Analysis/Supplementary/30_Replicates/gibbs1f90TBV.R")
  source("./Analysis/Supplementary/30_Replicates/AlphaPart_Gibbs_TBV_Validation.R")
  setwd("./Analysis/Supplementary/30_Replicates/ModelOutputTBV/")
  cat("Removing old files", "\n")
  file.remove(list.files())
  setwd("./../../../../")
  rm(list = ls()[!(ls() %in% c('i','r'))])
  gc()
}
