rm(list=ls())
load("./Results/burnin.RData")
for (generation in seq(nburnin + 1, nburnin + nfuture, 1)) {
  Program = "Scenario 2"
  cat("Evaluating", Program,": ", generation,"\n")
  
  year <- year + 1
  
  # Mate initial Sires and Dams
  pop_start <- randCross2(SelFemales, SelMales, 
                          nCrosses=SelFemales@nInd, nProgeny=1)
  # record pop_start
  pop_males <- pop_start[pop_start@sex=="M",]
  pop_females <- pop_start[pop_start@sex=="F",]
  
  RecSys <- RecSysMaleFemale(RecSys, pop_start)
  RecVar <- PullSumm(RecVar, pop_start)
  
  # Select next generation based on genetic value
  # Males
  SelMales <- 
    selectWithinFam(pop_start, nInd = 2, famType = "M", sex = "M") %>% 
    selectInd(nInd = nSires, use = "bv")
  RecSys$type[RecSys$ind %in% SelMales@id] <- "Selected"
  RecVar <- PullSumm(RecVar, SelMales, type ="Selected")
  
  # Females
  #SelFemales <- setPheno(SelFemales, h2 = h2)
  SelF <- selectInd(SelFemales, nDams/2, use = "bv")
  SelFemales = c(SelF, pop_females)
  
  RecVar <- PullSumm(RecVar, SelFemales, type ="Selected")
}

#-----------------------------------------------------------------------
# Save results
#-----------------------------------------------------------------------
saveRDS(RecSys, file="./Results/Scenario4.rds")
saveRDS(RecVar, file="./Results/Scenario4Var.rds")