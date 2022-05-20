source("./simulation5Males/methods.R")
RecSys <- RecVar <- NULL
generation <- 0

RecSys <- RecSysMaleFemale(RecSys, founders, gen = generation)
RecVar <- PullSumm(RecVar, founders)

# No selection for females and top 10% males
SelMales <- selectInd(founders, nInd = nSires, 
                      use = "rand", sex = "M")
SelFemales <- selectInd(founders, nInd = nDams, 
                        use = "rand", sex = "F")
RecSys$type[RecSys$ind %in% SelMales@id] <- "Selected"
RecVar <- PullSumm(RecVar, SelMales, type = "Selected")
RecVar <- PullSumm(RecVar, SelFemales, type = "Selected")

lastGenF <- selectInd(SelFemales, nInd = nDams/2)
  
for (generation in seq_len(nburnin)) {
  Program = "Burnin"
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
  SelMales <- selectInd(pop_males, nSires, use = "pheno")
  RecSys$type[RecSys$ind %in% SelMales@id] <- "Selected"
  RecVar <- PullSumm(RecVar, SelMales, type ="Selected")
  
  # Females
  #SelFemales <- setPheno(SelFemales, h2 = h2)
  SelF <- lastGenF 
  lastGenF <- pop_females
  SelFemales = c(SelF, pop_females)
  RecVar <- PullSumm(RecVar, SelFemales, type ="Selected")
}