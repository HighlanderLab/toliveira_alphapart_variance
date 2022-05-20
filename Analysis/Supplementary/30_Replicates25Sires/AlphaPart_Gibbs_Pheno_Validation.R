########################################################################
## Project: toliveira_alphapart_variance
## Script purpose: Phenotypic selection - accounting for inbreeding
## Date: 2022-01-18
## Author: Thiago de Paula Oliveira
########################################################################
#=======================================================================
# Loading packages and functions
#=======================================================================
source("./Analysis/Supplementary/30_Replicates25Sires/CommonFunctions.R")

#=======================================================================
# Organizing gibbs1f90 database
#=======================================================================
cat("Reading data", "\n")
ebvPath <- "./Analysis/Supplementary/30_Replicates25Sires/ModelOutput/all_solutions"
pedPath <- "./Analysis/Supplementary/30_Replicates25Sires/ModelOutput/renadd01.ped"
dataPath <- paste0("./Results/Scenario1_replicate25Males", r, ".rds")
dataEBV <- wrapperData(dataPath, pedPath, ebvPath)

#=======================================================================
# AlphaPart
#=======================================================================
cat("Doing partition", "\n")
m1 <- AlphaPart(dataEBV$data, colId = "ind", colFid = "father", 
                colMid = "mother", colPath = "path",
                colBV = dataEBV$colNames, 
                scaleEBV = list(center = TRUE, scale = TRUE))

m2 <- AlphaPart(dataEBV$data, colId = "ind", colFid = "father", 
                colMid = "mother", colPath = "path",
                colBV = "tbv", 
                scaleEBV = list(center = TRUE, scale = TRUE))
rm(dataEBV)
gc()

#=======================================================================
# Summary - Genetic Mean
#=======================================================================
cat("Summarising", "\n")
SMean <- wrapperSummary(m1 = m1, m2 = m2, ncols = 1500,
                        levels=c("F", "M(N)", "M(S)", "Sum"))
save(SMean, 
     file = paste0("./Analysis/Supplementary/30_Replicates25Sires/Results/SMean_TBV_ValidationPheno", r,
                   ".rds"))

#=======================================================================
# Summary - Genetic Variance
#=======================================================================
SVar <- wrapperSummary(m1 = m1, m2 = m2, ncols = 1500,
                       FUN = var, cov = TRUE,
                       levels=c("F", "F:M(N)", "F:M(S)", "M(N)", 
                                "M(N):M(S)", "M(S)", "Sum"))
save(SVar, 
     file = paste0("./Analysis/Supplementary/30_Replicates25Sires/Results/SVar_TBV_ValidationPheno", r, 
                   ".rds"))
