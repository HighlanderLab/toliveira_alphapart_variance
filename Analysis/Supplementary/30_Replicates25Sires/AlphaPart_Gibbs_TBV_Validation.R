########################################################################
## Project: toliveira_alphapart_variance
## Script purpose: analyse using AlphaPart
## Date: 2022-01-18
## Author: Thiago de Paula Oliveira
########################################################################
#=======================================================================
# Loading packages and functions
#=======================================================================
source("./Analysis/CommonFunctions.R")

#=======================================================================
# Organizing gibbs1f90 database
#=======================================================================
ebvPath <- "./Analysis/Supplementary/30_Replicates25Sires/ModelOutputTBV/all_solutions"
pedPath <- "./Analysis/Supplementary/30_Replicates25Sires/ModelOutputTBV/renadd01.ped"
dataPath <- paste0("./Results/Scenario2_replicate25Males", r, ".rds")
dataEBV <- wrapperData(dataPath, pedPath, ebvPath)

#=======================================================================
# AlphaPart
#=======================================================================
m1 <- AlphaPart(dataEBV$data, colId = "ind", colFid = "father", 
                colMid = "mother", colPath = "path",
                colBV = dataEBV$colNames, 
                scaleEBV = list(center = TRUE, scale = TRUE))

m2 <- AlphaPart(dataEBV$data, colId = "ind", colFid = "father", 
                colMid = "mother", colPath = "path",
                colBV = "tbv", 
                scaleEBV = list(center = TRUE, scale = TRUE))

#=======================================================================
# Summary - Genetic Mean
#=======================================================================
SMean <- wrapperSummary(m1 = m1, m2 = m2, ncols = 1500,
                        levels=c("F", "M(N)", "M(S)", "Sum"))
save(SMean, 
     file = paste0("./Analysis/Supplementary/30_Replicates25Sires/Results/SMean_TBV_ValidationTBV", r,
                   ".rds"))

#=======================================================================
# Summary - Genetic Variance
#=======================================================================
SVar <- wrapperSummary(m1 = m1, m2 = m2, ncols = 1500,
                       FUN = var, cov = TRUE,
                       levels=c("F", "F:M(N)", "F:M(S)", "M(N)", 
                                "M(N):M(S)", "M(S)", "Sum"))
save(SVar, 
     file = paste0("./Analysis/Supplementary/30_Replicates25Sires/Results/SVar_TBV_ValidationTBV", r, 
                   ".rds"))