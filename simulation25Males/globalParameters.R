########################################################################
## Project: Alpha Part Variance
## Script purpose: Global Parameters
## Date: 2021-11-29
## Author: Thiago de Paula Oliveira
########################################################################
nburnin <- 20
nfuture <- 20
#-----------------------------------------------------------------------
# Definitions
#-----------------------------------------------------------------------
nInd = 2000  # number of individuals
nQtl = 1000  # number of QTL's / Chromosome
nSnp = 1666  # number of SNPs / Chromosome
nChr = 30    # number of Chromosomes

# Founders parameters
meanG <- 0 #8000 # mean genetic value for the trait 
sigma2g <- 0.3 #1800^2 # genetic variance

# Herdability
h2 <- 0.3

# Error
(sigma2 <- sigma2g*(1/h2 - 1))

# Number of animals to be selected
nDams <- nInd/2 # no selection
nSires <- 25  

# Ne
Ne <- (4*nDams*nSires)/(nDams + nSires)
cat("Ne: ", Ne, "\n")
