rm(list=ls()[!ls() %in% c("replicate")])
load(paste0("./Results/burnin", replicate,".RData"))
for (generation in seq(nburnin + 1, nburnin + nfuture, 1)) {
  Program = "Scenario 1"
  cat("Evaluating", Program,": ", generation,"\n")
  
  # Mate initial Sires and Dams
  pop_start <- randCross2(SelFemales, SelMales, 
                          nCrosses=SelFemales@nInd, nProgeny=1)
  # record pop_start
  pop_males <- pop_start[pop_start@sex=="M",]
  pop_females <- pop_start[pop_start@sex=="F",]

  RecSys <- RecSysMaleFemale(RecSys, pop_start, gen = generation)
  RecVar <- PullSumm(RecVar, pop_start)

  
  # Select next generation based on genetic value
  # Males
  SelMales <- selectInd(pop_males, nSires, use = "gv")
  RecSys$type[RecSys$ind %in% SelMales@id] <- "Selected"
  RecVar <- PullSumm(RecVar, SelMales, type ="Selected")
  
  # Females
  SelF <- lastGenF 
  lastGenF <- pop_females
  SelFemales = c(SelF, pop_females)
  RecVar <- PullSumm(RecVar, SelFemales, type ="Selected")
  }