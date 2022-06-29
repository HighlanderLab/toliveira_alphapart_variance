########################################################################
## Project: jobsteter_alphapart_variance
## Script purpose: Dairy Cattle BP Analysis
## Date: 2021-11-30
## Author: Thiago de Paula Oliveira
########################################################################
# Set working directory
#=======================================================================

#=======================================================================
# Packges
#=======================================================================
#devtools::install_github("AlphaGenes/AlphaPart")
library(AlphaPart)
library(tidyverse)
library(tidyr)
library(patchwork)

APMeanSex <- APVarSex <- APMean <- APVar <- NULL


for(i in seq_len(30)){
  cat("Replicate", i, "\n")
  #=======================================================================
  # Reading and organizing Scenario 1
  #=======================================================================
  dataS1 <-  readRDS(paste0("./Results/Scenario1_replicate",i,".rds")) %>%
    dplyr::mutate(across(generation:mother, as.numeric)) %>%
    dplyr::mutate(across(c("sex", "type"), as.factor)) %>%
    arrange(generation, ind) %>%
    select(ind, father, mother, sex, type, generation, tbv, pheno) %>%
    dplyr::mutate(generation = generation - 20) %>%
    droplevels()
  #=======================================================================
  # AlphaPart
  #======================================================================= 
  dataS1 <- as.data.frame(dataS1)
  dataS1$path <- droplevels(dataS1$sex:dataS1$type)
  levels(dataS1$path)[1] <- "F"
  
  #-----------------------------------------------------------------------
  # Sex
  #-----------------------------------------------------------------------
  m1Sex <- AlphaPart(dataS1, colId = "ind", colFid = "father", 
                     colMid = "mother", colBV = "tbv", colPath = "sex",
                     scaleEBV = list(center = TRUE, scale = TRUE))
  S1Sex <- summary(m1Sex, by = "generation")
  S1Sex <- S1Sex$tbv
  S1Sex$rep <- i
  S1Sex$Scenario <- "Medium accuracy"
  APMeanSex <- rbind(APMeanSex, S1Sex)
  S1VarSex <- summary(m1Sex, by = "generation", FUN = var, cov = TRUE)
  S1VarSex <- S1VarSex$tbv
  S1VarSex$rep <- i
  S1VarSex$Scenario <- "Medium accuracy"
  APVarSex <- rbind(APVarSex, S1VarSex)
  
  #-----------------------------------------------------------------------
  # Sex:Type
  #-----------------------------------------------------------------------
  m1 <- AlphaPart(dataS1, colId = "ind", colFid = "father", 
                  colMid = "mother", colBV = "tbv", colPath = "path",
                  scaleEBV = list(center = TRUE, scale = TRUE))
  S1 <- summary(m1, by = "generation")
  S1 <- S1$tbv
  S1$rep <- i
  S1$Scenario <- "Medium accuracy"
  APMean <- rbind(APMean, S1)
  
  S1Var <- summary(m1, by = "generation", FUN = var, cov = TRUE)
  S1Var <- S1Var$tbv
  S1Var$rep <- i
  S1Var$Scenario <- "Medium accuracy"
  APVar <- rbind(APVar, S1Var)
  
  
  #=======================================================================
  # Reading and organizing Scenario 2
  #=======================================================================
  dataS2 <- readRDS(paste0("./Results/Scenario2_replicate",i,".rds")) %>%
    dplyr::mutate(across(generation:mother, as.numeric)) %>%
    dplyr::mutate(across(c("sex", "type"), as.factor)) %>%
    arrange(generation, ind) %>%
    select(ind, father, mother, sex, type, generation, tbv) %>%
    dplyr::mutate(generation = generation - 20) %>%
    droplevels()
  
  #=======================================================================
  # AlphaPart
  #======================================================================= 
  dataS2 <- as.data.frame(dataS2)
  dataS2$path <- droplevels(dataS2$sex:dataS2$type)
  levels(dataS2$path)[1] <- "F"
  
  #-----------------------------------------------------------------------
  # Sex
  #-----------------------------------------------------------------------
  m1Sex <- AlphaPart(dataS2, colId = "ind", colFid = "father", 
                     colMid = "mother", colBV = "tbv", colPath = "sex",
                     scaleEBV = list(center = TRUE, scale = TRUE))
  S2Sex <- summary(m1Sex, by = "generation")
  S2Sex <- S2Sex$tbv
  S2Sex$rep <- i
  S2Sex$Scenario <- "High accuracy"
  APMeanSex <- rbind(APMeanSex, S2Sex)
  
  S2VarSex <- summary(m1Sex, by = "generation", FUN = var, cov = TRUE)
  S2VarSex <- S2VarSex$tbv
  S2VarSex$rep <- i
  S2VarSex$Scenario <- "High accuracy"
  APVarSex <- rbind(APVarSex, S2VarSex)
  
  #-----------------------------------------------------------------------
  # Sex:Type
  #-----------------------------------------------------------------------
  m1 <- AlphaPart(dataS2, colId = "ind", colFid = "father", 
                  colMid = "mother", colBV = "tbv", colPath = "path",
                  scaleEBV = list(center = TRUE, scale = TRUE))
  S2 <- summary(m1, by = "generation")
  S2 <- S2$tbv
  S2$rep <- i
  S2$Scenario <- "High accuracy"
  APMean <- rbind(APMean, S2)
  
  S2Var <- summary(m1, by = "generation", FUN = var, cov = TRUE)
  S2Var <- S2Var$tbv
  S2Var$rep <- i
  S2Var$Scenario <- "High accuracy"
  APVar <- rbind(APVar, S2Var)
}


########################################################################
# Sex
########################################################################
#-----------------------------------------------------------------------
# Path Sex - Mean
#-----------------------------------------------------------------------
APMeanSexCI <- APMeanSex %>%
  group_by(generation, Scenario) %>%
  pivot_longer(cols = Sum:`M`,
               names_to = "type",
               values_to = "Median") %>%
  dplyr::mutate(across(type, factor)) %>%
  group_by(generation, Scenario, type) %>%
  summarise(
    Lower = quantile(Median, probs = 0.025),
    Est = quantile(Median, probs = 0.5),
    Upper = quantile(Median, probs = 0.975)
  )

APMeanSexCI$Scenario <- factor(APMeanSexCI$Scenario, levels = c("Medium accuracy", "High accuracy"))
levels(APMeanSexCI$Scenario)
#-----------------------------------------------------------------------
# Plot
#-----------------------------------------------------------------------
# Genetic Mean
g1 <- APMeanSexCI %>%
  ggplot(aes(y = Est, x = generation, groups = type),
         size = 0.1) +
  facet_wrap(~Scenario)+
  geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = type), alpha = 0.2) + 
  scale_fill_manual(values= c("red", "blue", "black")) +
  geom_line(aes(colour = type), alpha = 0.8) +
  scale_color_manual(values= c("red", "blue", "black")) +
  geom_vline(xintercept = 0, linetype = 2, alpha = 0.3) +
  ylab("Genetic Mean") +
  xlab("Generation") +
  xlim(-20,23)+
  labs(fill = "Selection path:", colour  = "Selection path:" ) +
  theme_bw(base_size = 15) + 
  theme(legend.position = "none")

g2 <- g1 + theme(legend.position = "none") +
  guides(
    color = guide_legend(override.aes = list(linetype = c(0, 1),
                                             shape = c(16, NA),
                                             color = "black")))

textplot <- tibble(
  Est = c(APMeanSexCI[APMeanSexCI$generation==20 & APMeanSexCI$Scenario=="Medium accuracy" & APMeanSexCI$type=="Sum",]$Est,
          APMeanSexCI[APMeanSexCI$generation==20 & APMeanSexCI$Scenario=="High accuracy" & APMeanSexCI$type=="Sum",]$Est,
          APMeanSexCI[APMeanSexCI$generation==20 & APMeanSexCI$Scenario=="Medium accuracy"  & APMeanSexCI$type=="F",]$Est,
          APMeanSexCI[APMeanSexCI$generation==20 & APMeanSexCI$Scenario=="High accuracy" & APMeanSexCI$type=="F",]$Est,
          APMeanSexCI[APMeanSexCI$generation==20 & APMeanSexCI$Scenario=="Medium accuracy"& APMeanSexCI$type=="M",]$Est,
          APMeanSexCI[APMeanSexCI$generation==20 & APMeanSexCI$Scenario=="High accuracy" & APMeanSexCI$type=="M",]$Est),
  Scenario = as.factor(rep(c("Medium accuracy","High accuracy"),3)),
  generation = 22,
  label = c(rep("Sum",2), rep("F",2), rep("M",2)),
  type = c(rep("Sum",2), rep("F",2), rep("M",2)),
  colour = c(rep("black",2), rep("red",2), rep("blue",2)),
  size = 3.6
)

g3 <- g2 + geom_text(
  data = textplot,
  label = textplot$label,
  colour = textplot$colour,
  size = textplot$size
)
g3


#=======================================================================
# Path: Sex - Variance
#=======================================================================
APVarSexCI <- APVarSex %>%
  group_by(generation, Scenario) %>%
  pivot_longer(cols = Sum:`FM`,
               names_to = "type",
               values_to = "Median") %>%
  dplyr::mutate(across(type, factor)) %>%
  group_by(generation, Scenario, type) %>%
  summarise(
    Lower = quantile(Median, probs = 0.025),
    Est = quantile(Median, probs = 0.5),
    Upper = quantile(Median, probs = 0.975)
  )
levels(APVarSexCI$type)[2] <- "F:M"
APVarSexCI$Scenario <- factor(APVarSexCI$Scenario, levels = c("Medium accuracy", "High accuracy"))
levels(APVarSexCI$Scenario)

#-----------------------------------------------------------------------
# Plot
#-----------------------------------------------------------------------
# Genetic Var
g4 <- APVarSexCI %>%
  ggplot(aes(y = Est, x = generation, groups = type),
         size = 0.1) +
  facet_wrap(~Scenario)+
  geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = type), alpha = 0.2) + 
  scale_fill_manual(values= c("red", "#4b7000ff", "blue", "black")) +
  geom_line(aes(colour = type), alpha = 0.8) +
  scale_color_manual(values= c("red", "#4b7000ff", "blue", "black")) +
  geom_vline(xintercept = 0, linetype = 2, alpha = 0.3) +
  ylab("Genetic Variance") +
  xlab("Generation") +
  xlim(-20,23)+
  labs(fill = "Selection path:", colour  = "Selection path:" ) +
  theme_bw(base_size = 15) + 
  theme(legend.position = "none")

g5 <- g4 + theme(legend.position = "none") +
  guides(
    color = guide_legend(override.aes = list(linetype = c(0, 1),
                                             shape = c(16, NA),
                                             color = "black")))

textplot <- tibble(
  Est = c(APVarSexCI[APVarSexCI$generation==20 & APVarSexCI$Scenario=="Medium accuracy" & APVarSexCI$type=="Sum",]$Est,
          APVarSexCI[APVarSexCI$generation==20 & APVarSexCI$Scenario=="High accuracy" & APVarSexCI$type=="Sum",]$Est,
          APVarSexCI[APVarSexCI$generation==20 & APVarSexCI$Scenario=="Medium accuracy"  & APVarSexCI$type=="F",]$Est,
          APVarSexCI[APVarSexCI$generation==20 & APVarSexCI$Scenario=="High accuracy" & APVarSexCI$type=="F",]$Est,
          APVarSexCI[APVarSexCI$generation==20 & APVarSexCI$Scenario=="Medium accuracy"  & APVarSexCI$type=="F:M",]$Est,
          APVarSexCI[APVarSexCI$generation==20 & APVarSexCI$Scenario=="High accuracy" & APVarSexCI$type=="F:M",]$Est,
          APVarSexCI[APVarSexCI$generation==20 & APVarSexCI$Scenario=="Medium accuracy"& APVarSexCI$type=="M",]$Est,
          APVarSexCI[APVarSexCI$generation==20 & APVarSexCI$Scenario=="High accuracy" & APVarSexCI$type=="M",]$Est),
  Scenario = as.factor(rep(c("Medium accuracy","High accuracy"),4)),
  generation = 22,
  label = c(rep("Sum",2), rep("F",2), rep("F:M",2), rep("M",2)),
  type = label,
  colour = c(rep("black",2), rep("red",2), rep("#4b7000ff",2), rep("blue",2)),
  size = 3.6
)

g6 <- g5 + geom_text(
  data = textplot,
  label = textplot$label,
  colour = textplot$colour,
  size = textplot$size
)
g6

g3/g6

ggsave("./Analysis/Figures/AlphaPartSex30Reps.pdf",width = 8, height = 8)

#=======================================================================
# Path: Sex/Selection - Mean
#=======================================================================
APMeanCI <- APMean %>%
  group_by(generation, Scenario) %>%
  pivot_longer(cols = Sum:`M:Selected`,
               names_to = "type",
               values_to = "Median") %>%
  dplyr::mutate(across(type, factor)) %>%
  group_by(generation, Scenario, type) %>%
  summarise(
    Lower = quantile(Median, probs = 0.025),
    Est = quantile(Median, probs = 0.5),
    Upper = quantile(Median, probs = 0.975)
  )
levels(APMeanCI$type)[2:3] <- c("M(N)", "M(S)") 
APMeanCI$Scenario <- factor(APMeanCI$Scenario, levels = c("Medium accuracy", "High accuracy"))
levels(APMeanCI$Scenario)

#-----------------------------------------------------------------------
# Plot
#-----------------------------------------------------------------------
# Genetic Mean
g1 <- APMeanCI %>%
  ggplot(aes(y = Est, x = generation, groups = type),
         size = 0.1) +
  facet_wrap(~Scenario)+
  geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = type), alpha = 0.2) + 
  scale_fill_manual(values= c("red", "#9597a1", "blue", "black")) +
  geom_line(aes(colour = type), alpha = 0.8) +
  scale_color_manual(values= c("red", "#9597a1", "blue", "black")) +
  geom_vline(xintercept = 0, linetype = 2, alpha = 0.3) +
  ylab("Genetic Mean") +
  xlab("Generation") +
  xlim(-20,26.5)+
  labs(fill = "Selection path:", colour  = "Selection path:" ) +
  theme_bw(base_size = 15) + 
  theme(legend.position = "none")

g2 <- g1 + theme(legend.position = "none") +
  guides(
    color = guide_legend(override.aes = list(linetype = c(0, 1),
                                             shape = c(16, NA),
                                             color = "black")))

textplot <- tibble(
  Est = c(APMeanCI[APMeanCI$generation==20 & APMeanCI$Scenario=="Medium accuracy" & APMeanCI$type=="Sum",]$Est,
          APMeanCI[APMeanCI$generation==20 & APMeanCI$Scenario=="High accuracy" & APMeanCI$type=="Sum",]$Est,
          APMeanCI[APMeanCI$generation==20 & APMeanCI$Scenario=="Medium accuracy"  & APMeanCI$type=="F",]$Est,
          APMeanCI[APMeanCI$generation==20 & APMeanCI$Scenario=="High accuracy" & APMeanCI$type=="F",]$Est,
          APMeanCI[APMeanCI$generation==20 & APMeanCI$Scenario=="Medium accuracy"& APMeanCI$type=="M(N)",]$Est,
          APMeanCI[APMeanCI$generation==20 & APMeanCI$Scenario=="High accuracy" & APMeanCI$type=="M(N)",]$Est,
          APMeanCI[APMeanCI$generation==20 & APMeanCI$Scenario=="Medium accuracy"& APMeanCI$type=="M(S)",]$Est,
          APMeanCI[APMeanCI$generation==20 & APMeanCI$Scenario=="High accuracy" & APMeanCI$type=="M(S)",]$Est),
  Scenario = as.factor(rep(c("Medium accuracy","High accuracy"),4)),
  generation = 24.5,
  type = c(rep("Sum",2), rep("F",2), rep("M(N)",2), rep("M(S)",2)),
  colour = c(rep("black",2), rep("red",2), rep("#9597a1",2), rep("blue",2)),
  size = 3.6
)

g3 <- g2 + geom_text(
  data = textplot,
  label = textplot$type,
  colour = textplot$colour,
  size = textplot$size
)
g3

#=======================================================================
# Path: Sex/Selection - Var
#=======================================================================
APVar$cor <- APVar$`FM:Selected`/sqrt(APVar$`F`*APVar$`M:Selected`)
APVarCI <- APVar %>%
  group_by(generation, Scenario) %>%
  pivot_longer(cols = Sum:`M:Non-SelectedM:Selected`,
               names_to = "type",
               values_to = "Median") %>%
  dplyr::mutate(across(type, factor)) %>%
  group_by(generation, Scenario, type) %>%
  summarise(
    Lower = quantile(Median, probs = 0.025),
    Est = quantile(Median, probs = 0.5),
    Upper = quantile(Median, probs = 0.975)
  )
levels(APVarCI$type)[2:6] <- c("F:M(N)", "F:M(S)", "M(N)", "M(N):M(S)", "M(S)")
APVarCI$Scenario <- factor(APVarCI$Scenario, levels = c("Medium accuracy", "High accuracy"))
levels(APVarCI$Scenario)

#-----------------------------------------------------------------------
# Plot
#-----------------------------------------------------------------------
# Genetic Var
g4 <- APVarCI %>%
  ggplot(aes(y = Est, x = generation, groups = type),
         size = 0.1) +
  facet_wrap(~Scenario)+
  geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = type), alpha = 0.2) + 
  scale_fill_manual(values= c("red", "#823175", "#4b7000ff", "#9597a1",
                              "#c4852b", "blue", "black")) +
  geom_line(aes(colour = type), alpha = 0.8) +
  scale_color_manual(values= c("red", "#823175", "#4b7000ff", "#9597a1",
                               "#c4852b", "blue", "black")) +
  geom_vline(xintercept = 0, linetype = 2, alpha = 0.3) +
  ylab("Genetic Variance") +
  xlab("Generation") +
  xlim(-20,26.5)+
  labs(fill = "Selection path:", colour  = "Selection path:" ) +
  theme_bw(base_size = 15) + 
  theme(legend.position = "none")

g5 <- g4 + theme(legend.position = "none") +
  guides(
    color = guide_legend(override.aes = list(linetype = c(0, 1),
                                             shape = c(16, NA),
                                             color = "black")))

textplot <- tibble(
  Est = c(APVarCI[APVarCI$generation==20 & APVarCI$Scenario=="Medium accuracy" & APVarCI$type=="Sum",]$Est,
          APVarCI[APVarCI$generation==20 & APVarCI$Scenario=="High accuracy" & APVarCI$type=="Sum",]$Est,
          APVarCI[APVarCI$generation==20 & APVarCI$Scenario=="Medium accuracy"  & APVarCI$type=="F",]$Est,
          APVarCI[APVarCI$generation==20 & APVarCI$Scenario=="High accuracy" & APVarCI$type=="F",]$Est,
          APVarCI[APVarCI$generation==20 & APVarCI$Scenario=="Medium accuracy"& APVarCI$type=="M(N)",]$Est,
          APVarCI[APVarCI$generation==20 & APVarCI$Scenario=="High accuracy" & APVarCI$type=="M(N)",]$Est,
          APVarCI[APVarCI$generation==20 & APVarCI$Scenario=="Medium accuracy"& APVarCI$type=="M(S)",]$Est,
          APVarCI[APVarCI$generation==20 & APVarCI$Scenario=="High accuracy" & APVarCI$type=="M(S)",]$Est,
          APVarCI[APVarCI$generation==20 & APVarCI$Scenario=="Medium accuracy"& APVarCI$type=="F:M(N)",]$Est,
          APVarCI[APVarCI$generation==20 & APVarCI$Scenario=="High accuracy" & APVarCI$type=="F:M(N)",]$Est,
          APVarCI[APVarCI$generation==20 & APVarCI$Scenario=="Medium accuracy"& APVarCI$type=="F:M(S)",]$Est,
          APVarCI[APVarCI$generation==20 & APVarCI$Scenario=="High accuracy" & APVarCI$type=="F:M(S)",]$Est,
          APVarCI[APVarCI$generation==20 & APVarCI$Scenario=="Medium accuracy"& APVarCI$type=="M(N):M(S)",]$Est,
          APVarCI[APVarCI$generation==20 & APVarCI$Scenario=="High accuracy" & APVarCI$type=="M(N):M(S)",]$Est),
  Scenario = as.factor(rep(c("Medium accuracy","High accuracy"),7)),
  generation = 24.5,
  type = c(rep("Sum",2),  rep("F",2), rep("M(N)",2),  rep("M(S)",2), 
           rep("F:M(N)",2), rep("F:M(S)",2), rep("M(N):M(S)",2)),
  colour = c(rep("black",2), rep("red",2), rep("#9597a1",2), rep("blue",2), 
             rep("#823175",2), rep("#4b7000ff",2), rep("#c4852b",2)),
  size = 3.6
)

g6 <- g5 + geom_text(
  data = textplot,
  label = textplot$type,
  colour = textplot$colour,
  size = textplot$size
)
g6

g3/g6
ggsave("./Analysis/Figures/AlphaPartSexSel30Reps.pdf",width = 8, height = 8)



APVarCI <- APVar %>%
  group_by(generation, Scenario) %>%
  pivot_longer(cols = cor,
               names_to = "type",
               values_to = "Median") %>%
  dplyr::mutate(across(type, factor)) %>%
  group_by(generation, Scenario, type) %>%
  summarise(
    Lower = quantile(Median, probs = 0.025),
    Est = quantile(Median, probs = 0.5),
    Upper = quantile(Median, probs = 0.975)
  )
levels(APVarCI$type) <- c("F:M(S)")
APVarCI$Scenario <- factor(APVarCI$Scenario, levels = c("Medium accuracy", "High accuracy"))
levels(APVarCI$Scenario)

APVarCI %>%
  filter(type == "F:M(S)") %>%
  ggplot(aes(y = Est, x = generation, groups = type),
         size = 0.1) +
  facet_wrap(~Scenario)+
  geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "#4b7000ff",
              alpha = 0.2) + 
  geom_line(colour = "#4b7000ff", alpha = 0.8) +
  geom_vline(xintercept = 0, linetype = 2, alpha = 0.3) +
  ylab(expression(rho['F, M(S)'])) +
  xlab("Generation") +
  xlim(-20,20)+
  labs(fill = "Selection path:", colour  = "Selection path:" ) +
  theme_bw(base_size = 15) + 
  theme(legend.position = "none")
ggsave("./Analysis/Figures/AlphaPartCor30Reps.pdf", width = 8, height = 4)

