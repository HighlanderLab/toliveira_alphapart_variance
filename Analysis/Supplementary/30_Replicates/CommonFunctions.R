########################################################################
## Project: toliveira_alphapart_variance
## Script purpose: analyse using AlphaPart
## Date: 2022-01-18
## Author: Thiago de Paula Oliveira
########################################################################
#=======================================================================
# Library
#=======================================================================
library(tidyverse)
library(AlphaPart)
library(patchwork)

#=======================================================================
# Inbreeding analysis
#=======================================================================
wrapperDataInb <- function(dataPath, pedPath, ebvPath){
  data <- read.table(ebvPath)
  ped1 <- read_table(pedPath, 
                     col_names=FALSE, skip = 0,
                     col_types = cols(.default = col_character(),
                                      X1 = col_double()))
  data <- data %>%
    select(V2, V3)
  # Keep renumbered and original ID
  ped1 = ped1 %>%
    select("X1", "X10") %>%
    arrange(X1)
  colnames(ped1) = c("Level", "id")
  
  # Connect EBV with original animal ID (if RENUMF90 was used)
  data = inner_join(ped1, data, by=c("Level"="V3"))
  data$id <- as.numeric(data$id)
  data <- data[order(data$id), ]
  head(data)
  
  dataS2 <-  readRDS(dataPath) %>%
    dplyr::mutate(across(generation:mother, as.numeric)) %>%
    dplyr::mutate(across(c("sex", "type"), as.factor)) %>%
    arrange(generation, ind) %>%
    select(ind, father, mother, sex, type, generation, tbv, pheno)
  dataS2$path <- droplevels(dataS2$sex:dataS2$type)
  levels(dataS2$path)[1] <- "F"
  
  dataS2 <- dataS2[order(dataS2$ind), ]
  dataS2 <- cbind(dataS2, data)
  
  datPheno <- dataS2 %>%
    group_by(generation) %>%
    dplyr::mutate(generation = generation - 20) %>%
    summarise(
      Inbreeding = mean(V2),
      InbreedingSD = sd(V2)
    )
  
  datPheno$Lower <- datPheno$Inbreeding - 1.96*datPheno$InbreedingSD
  datPheno$Lower[datPheno$Lower<0]=0  
  datPheno$Upper <- datPheno$Inbreeding + 1.96*datPheno$InbreedingSD  
  return(datPheno)
}

#=======================================================================
# Preparing EBV data for AlphaPart
#=======================================================================
wrapperData <- function(dataPath, pedPath, ebvPath) {
  data <- read.table(ebvPath)
  data <- data %>%
    filter(V1 != -1) %>%
    select(V3, V4)
  
  data$Sample <- gl(n = nrow(data)/max(data$V3), k = max(data$V3))
  
  data <- data %>%
    pivot_wider(names_from = Sample, values_from = V4)
  
  colnames(data) <- c("Level", paste0("s", 1:(ncol(data)-1)))
  resp_names <- paste0("s", 1:(ncol(data)-1))
  
  #=======================================================================
  # Linking original data with gibbs
  #=======================================================================
  ped1 <- read_table(pedPath, 
                     col_names=FALSE, skip = 0,
                     col_types = cols(.default = col_character(),
                                      X1 = col_double()))
  # Keep renumbered and original ID
  ped2 = ped1 %>%
    select("X1", "X10") %>%
    arrange(X1)
  colnames(ped2) = c("Level", "id")
  
  # Connect EBV with original animal ID (if RENUMF90 was used)
  data = inner_join(ped2, data, by=c("Level"="Level"))
  data$id <- as.numeric(data$id)
  data <- data[order(data$id), ]
  head(data)
  
  #=======================================================================
  # Merging original and gibbs databases
  #=======================================================================
  dataS2 <-  readRDS(dataPath) %>%
    dplyr::mutate(across(generation:mother, as.numeric)) %>%
    dplyr::mutate(across(c("sex", "type"), as.factor)) %>%
    arrange(generation, ind) %>%
    select(ind, father, mother, sex, type, generation, tbv, pheno)
  dataS2$path <- droplevels(dataS2$sex:dataS2$type)
  levels(dataS2$path)[1] <- "F"
  
  dataS2 <- dataS2[order(dataS2$ind), ]
  dataS2 <- cbind(dataS2, data)
  
  # Test
  if(sum(dataS2$ind==dataS2$id)!=nrow(dataS2)){
    stop("Check line 114 of function wrapperData") # 42000 expected
  } 
  return(list(data = dataS2, colNames = resp_names))
}

#=======================================================================
# Summary AlphaPart
#=======================================================================
wrapperSummary <- function(m1, m2, ncols, FUN = mean, 
                           cov = FALSE, levels = NULL, nBurnin = 20){
  cat("Summarising partition MCMC", "\n")
  S2 <- summary(m1, by = "generation", FUN = FUN, cov = cov)
  cat("Summarising the true partition", "\n")
  S2True <- summary(m2, by = "generation", FUN = FUN, cov = cov)

  cat("Organizing data structure and computing 95% CI")
  getPatchCI <- NULL
  rangeCol = 3:ncol(S2$s1)
  for (i in rangeCol){
    mSamples <- matrix(NA, ncol = ncols, nrow = nrow(S2$s1))
    for (j in seq_len(ncols)){
      mSamples[,j] <- S2[[j]][,i]
      getName <- colnames(S2[[j]])[i]
    }
    mSamplesCI <- t(apply(mSamples, 1, quantile, probs = c(0.025, 0.5, 0.975)))
    colnames(mSamplesCI) <- c("Lower", "Median", "Upper")
    mSamplesCI <- as.data.frame(mSamplesCI)
    mSamplesCI$path <- getName
    getPatchCI <- rbind(getPatchCI, mSamplesCI)
  }
  
  getPatchCI <- data.frame(do.call("rbind", replicate(
    length(unique(getPatchCI$path)), 
    S2$s1[,1:2], simplify = FALSE)), 
    getPatchCI)
  getPatchCI$generation <- getPatchCI$generation-nBurnin
  getPatchCI$path <- as.factor(getPatchCI$path)
  if(!is.null(levels)){
    if(length(levels(getPatchCI$path)) != length(levels)){
      stop("Number of levels differ from Path variable")
    }
    levels(getPatchCI$path) <- levels
  }
  
  S2True2 <- S2True$tbv %>%
    pivot_longer(cols = all_of(rangeCol),
                 names_to = "path",
                 values_to = "Median") %>%
    dplyr::mutate(across(path, factor))
  
  if(!is.null(levels)){
    if(length(levels(getPatchCI$path)) != length(levels)){
      stop("Number of levels differ from Path variable")
    }
    levels(S2True2$path) <-  levels
  }
  
  S2True2$generation <- S2True2$generation-nBurnin
  
  if(identical(deparse(FUN),deparse(mean))){
    cat("Summarising the true Mendelian Sampling term", "\n")
    mstTrue <- m2$tbv %>%
      dplyr::mutate(generation = generation - nBurnin) %>%
      group_by(generation, path) %>%
      summarise(
        mstTrue = mean(tbv_w)
      )
    mstTrueSum <- m2$tbv %>%
      dplyr::mutate(generation = generation - nBurnin) %>%
      group_by(generation) %>%
      summarise(
        mstTrue = mean(tbv_w)
      ) %>%
      dplyr::mutate(path = as.factor("Sum")) %>%
      relocate(path, .before = mstTrue)
    mstTrue <- rbind(mstTrue, mstTrueSum)
    
    cat("Summarising the estimated Mendelian Sampling term", "\n")
    mSamples <- matrix(NA, ncol = ncols, nrow = nrow(m1$s1))    
    for (i in seq_len(ncols)){
      pos <- grep("_w", colnames(m1$s1))
      mSamples[,i] <- m1[[i]][,pos]
    }
    mSamples <- as.data.frame(mSamples)
    mSamples$path <- m1$s1$path
    mSamples$generation <- m1$s1$generation
    
    mSamplesSum <- mSamples %>%
      dplyr::mutate(generation = generation - nBurnin) %>%
      group_by(generation) %>%
      summarise_all(mean) %>%
      dplyr::mutate(path = as.factor("Sum")) %>%
      relocate(path, .after = generation)
    
    mSamples <- mSamples %>%
      dplyr::mutate(generation = generation - nBurnin) %>%
      group_by(generation, path) %>%
      summarise_all(mean)
    
    mSamples <- rbind(mSamples, mSamplesSum)
    
    mstEBV <- t(apply(mSamples[,-c(1:2)], 1, quantile, probs = c(0.025, 0.5, 0.975)))
    colnames(mstEBV) <- c("Lower", "Median", "Upper")
    mstEBV <- as.data.frame(mstEBV)
    mstEBV$generation <- mSamples$generation
    mstEBV$path <- mSamples$path
    if(!is.null(levels)){
      if(length(levels(mstEBV$path)) != length(levels)){
        stop("Number of levels differ from Path variable")
      }
      levels(mstTrue$path) <- levels
      levels(mstEBV$path) <- levels
    }
  }else{
    mstTrue <- mstEBV <- NULL
  }
  return(list(SummaryEBV = getPatchCI, 
              SummaryTrue = S2True2,
              mstTrue = mstTrue,
              mstEBV = mstEBV))
}
