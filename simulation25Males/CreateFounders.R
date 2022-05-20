########################################################################
## Project:
## Script purpose:
## Date: 2021-11-30
## Author: Thiago de Paula Oliveira
########################################################################
#-----------------------------------------------------------------------
# Haplotypes
#-----------------------------------------------------------------------
if(file.exists("./Results/haplotypes.RData")){
  load("./Results/haplotypes.RData")  
} else{
  founderPop <- runMacs(nInd = nInd,  nChr = nChr, segSites = nQtl + nSnp,
                        species = "CATTLE", ploidy = 2L)
  save(founderPop, file = "./Results/haplotypes.RData")  
}

#-----------------------------------------------------------------------
# Parameters
#-----------------------------------------------------------------------
SP <- SimParam$new(founderPop)

#-----------------------------------------------------------------------
# Quantitative trait (A = only have additive effects)
#-----------------------------------------------------------------------
SP$addTraitA(nQtlPerChr = nQtl,  mean = meanG,  var = sigma2g)

#-----------------------------------------------------------------------
# How sexes are determined in the simulation
## Sistematic assignment to ensure same number of males and females
#-----------------------------------------------------------------------
SP$setSexes("yes_sys")
#-----------------------------------------------------------------------
# Add SNPs Chip
#-----------------------------------------------------------------------
SP$addSnpChip(nSnp)

#-----------------------------------------------------------------------
# Attributing phenotypes - All individuals
#-----------------------------------------------------------------------
SP$setVarE(h2 = h2)

# New pop
founders <- newPop(founderPop)

