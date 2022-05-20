########################################################################
## Project: toliveira_alphapart_variance
## Script purpose: run gibbs1f90
## Date: 2022-01-17
## Author: Thiago de Paula Oliveira
########################################################################
library(tidyverse)
dataS2 <-  readRDS(paste0("./Results/Scenario2_replicate",r,".rds")) %>%
  dplyr::mutate(across(generation:mother, as.numeric)) %>%
  dplyr::mutate(across(c("sex", "type"), as.factor)) %>%
  arrange(generation, ind) %>%
  select(ind, father, mother, sex, type, generation, tbv, pheno)

dataS2$mean <- 1

setwd("./Analysis/Supplementary/30_Replicates/ModelOutputTBV/")

print("Saving phenotype file")
write.table(dataS2[, c("ind", "mean", "pheno")],
            "gibbsf90.dat", quote = FALSE, row.names = FALSE,
            col.names = FALSE, sep = " ", na = "0")

## Saving pedigree file
print("Saving pedigree file")
ped <- dataS2 %>%
  select(ind, father, mother)

write.table(ped, "gibbsf90.ped",
            quote = FALSE, row.names = FALSE, col.names = FALSE,
            sep = " ", na = "0")

## Creating parameter file
print("Saving parameter file")
sink("renumf90.par", type = "output")
writeLines("#renumf90 parameter file
DATAFILE
gibbsf90.dat
TRAITS # pheno
3
WEIGHT(S) \n
RESIDUAL_VARIANCE
0.7
EFFECT # ind as random
1 cross alpha
RANDOM
animal
FILE
gibbsf90.ped
INBREEDING
pedigree
(CO)VARIANCES
0.3 \n
EFFECT # overall mean
2 cross alpha
OPTION solution all 1
")
sink()
system(
  command = "echo renumf90.par | $HOME/bin/renumf90 | tee renum.log"
)


system(
  command = "$HOME/bin/gibbs1f90 renf90.par --rounds 80000 --burnin 20000 --thin 40 --thinprint 40"
)

setwd("./../../../../")

# postgibbsf90 renf90.par
