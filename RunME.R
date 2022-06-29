########################################################################
## Project: AlphaPart Variance
## Script purpose: execute simulation, analysis and create ficgures
## Date: 2022-06-29
## Author: Thiago de Paula Oliveira
########################################################################

#=======================================================================
# Running simulation
#=======================================================================
source("./simulation5Males/breedingProgrammeScheme.R")

#=======================================================================
# Running Analysis -  1 replicate
#=======================================================================
# Descriptive Analysis
source("./Analysis/Descriptive.R")

# Analysis of true breeding values
source("./Analysis/AlphaPart_TruePartition.R")

# Running gibbs1f90
# Here I am assuming the blupf90 family of programmes is installed at 
# $HOME/bin/ directory
source("./Analysis/gibbs1f90.R")

# Running AlphaPart Analysis - Medium Accuracy Scenario
source("./Analysis/AlphaPart_Gibbs_Pheno_Validation.R")

# Running AlphaPart Analysis - High Accuracy Scenario
source("./Analysis/AlphaPart_Gibbs_TBV_Validation.R")


#=======================================================================
# Running Analysis -  20 replicate
#=======================================================================
# Analysis of true breeding values
source("./Analysis/AlphaPart_TruePartition30reps.R")

# Running gibbs1f90 and AlphaPart Analysis
# Here I am assuming the blupf90 family of programmes is installed at 
# $HOME/bin/ directory
source("./Analysis/Supplementary/30_Replicates/RUNME.R")



