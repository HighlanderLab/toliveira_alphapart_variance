########################################################################
## Project: toliveira_alphapart_variance
## Script purpose: Medium accuracy - accounting for inbreeding
## Date: 2022-01-18
## Author: Thiago de Paula Oliveira
########################################################################
#=======================================================================
# Loading packages and functions
#=======================================================================
library(DescTools)
source("./Analysis/CommonFunctions.R")
#=======================================================================
# Structuring the data
#=======================================================================
ebvPath <- "./Analysis/ModelOutput/renf90.inb"
pedPath <- "./Analysis/ModelOutput/renadd01.ped"
dataPath <- "./Results/Scenario1_replicate1.rds"
dataInb <- wrapperDataInb(dataPath, pedPath, ebvPath)

g1 <- dataInb %>%
  ggplot(aes(y=Inbreeding, x = generation)) +
  geom_ribbon(aes(ymin = Lower,
                  ymax = Upper), alpha = 0.2) +
  geom_line() +
  ylim(0, 0.58)+
  geom_vline(xintercept = 0, linetype = 2, alpha = 0.2) +
  xlab("Generation") +
  ylab("Average Inbreeding") +
  ggtitle("Medium accuracy") +
  theme_bw()

#-----------------------------------------------------------------------
# TBV selection
#-----------------------------------------------------------------------
ebvPath <- "./Analysis/ModelOutputTBV/renf90.inb"
pedPath <- "./Analysis/ModelOutputTBV/renadd01.ped"
dataPath <- "./Results/Scenario2_replicate1.rds"
dataInbTBV <- wrapperDataInb(dataPath, pedPath, ebvPath)

g2 <- dataInbTBV %>%
  ggplot(aes(y=Inbreeding, x = generation)) +
  geom_ribbon(aes(ymin = Lower,
                  ymax = Upper), alpha = 0.2) +
  ylim(0, 0.58)+
  geom_line() +
  xlab("Generation") +
  ylab("Average Inbreeding") +
  ggtitle("High accuracy") +
  geom_vline(xintercept = 0, linetype = 2, alpha = 0.2) +
  theme_bw()
g2
g3 <- g1+g2
ggsave("./Analysis/Figures/Inbreeding.pdf", width = 7, height = 4)

rm(list = ls())
#=======================================================================
# Loading packages and functions
#=======================================================================
source("./Analysis/CommonFunctions.R")

#=======================================================================
# Organizing gibbs1f90 database
#=======================================================================
ebvPath <- "./Analysis/ModelOutput/all_solutions"
pedPath <- "./Analysis/ModelOutput/renadd01.ped"
dataPath <- "./Results/Scenario1_replicate1.rds"
dataEBV <- wrapperData(dataPath, pedPath, ebvPath)

#=======================================================================
# test
#=======================================================================
summMean <- dataEBV$data %>%
  group_by(generation, path) %>%
  summarise(meanTBV = mean(tbv),
            across(s1:s1500, ~ mean(.x, na.rm = TRUE)))
summMean$meanEBV <- rowMeans(summMean[,dataEBV$colNames]) 

# Concordance
ccc <- CCC(summMean$meanTBV, summMean$meanEBV)
ccc$rho.c

# Pearson Cor.
ccc$rho.c[1]/ccc$C.b

# C_b
ccc$C.b

#plot
p1 <- summMean %>%
  dplyr::mutate(path = recode(path, 
                              `F`="F", 
                              `M:Non-Selected`= "M(NS)",
                              `M:Selected` = "M(S)")) %>%
  ggplot(aes(x = meanTBV, y = meanEBV)) +
  facet_wrap(~path)+
  geom_point() +
  xlab("Average TBV per generation") +
  ylab("Average EBV per generation") +
  geom_abline(intercept = 0, slope = 1) +
  ggtitle("Medium accuracy") +
  theme_bw()


summVar <- dataEBV$data %>%
  dplyr::mutate(path = recode(path, 
                              `F`="F", 
                              `M:Non-Selected`= "M(NS)",
                              `M:Selected` = "M(S)")) %>%
  group_by(generation, path) %>%
  summarise(varTBV = var(tbv),
            across(s1:s1500, ~ var(.x, na.rm = TRUE)))
summVar$varEBV <- rowMeans(summVar[,dataEBV$colNames]) 


# Concordance
ccc <- CCC(summVar$varTBV, summVar$varEBV)
ccc$rho.c

# Pearson Cor.
ccc$rho.c[1]/ccc$C.b

# C_b
ccc$C.b

# plot
p2 <- summVar %>%
  ggplot(aes(x = varTBV, y = varEBV)) +
  facet_wrap(~path)+
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  xlab("Average TBV per generation") +
  ylab("Average EBV per generation") +
  geom_abline(intercept = 0, slope = 1) +
  ggtitle("Medium accuracy") +
  theme_bw()

p1/p2
ggsave("./Analysis/Figures/EBV_TBV_ScenarioPheno.pdf", 
       width = 9, height = 6)
#=======================================================================
# AlphaPart
#=======================================================================
m1 <- AlphaPart(dataEBV$data, colId = "ind", colFid = "father", 
                colMid = "mother", colPath = "path",
                colBV = dataEBV$colNames, 
                scaleEBV = list(center = TRUE, scale = TRUE))

m2 <- AlphaPart(dataEBV$data, colId = "ind", colFid = "father", 
                colMid = "mother", colPath = "path",
                colBV = "tbv", 
                scaleEBV = list(center = TRUE, scale = TRUE))
saveRDS(m1, "./Analysis/Results/AlphaPart_EBV_ValidationPheno.rds")
saveRDS(m2, "./Analysis/Results/AlphaPart_TBV_ValidationPheno.rds")

rm(dataEBV)
gc()

m2 <- readRDS("./Analysis/Results/AlphaPart_TBV_ValidationPheno.rds")

m2$tbv %>%
  select(ind:path, tbv:`tbv_M:Selected`) %>%
  dplyr::mutate(path = recode(path, 
                              `F`="F", 
                              `M:Non-Selected`= "M(NS)",
                              `M:Selected` = "M(S)")) %>%
  group_by(generation, path) %>%
  summarise(
    Female = mean(tbv_F),
    `M(NS)` = mean(`tbv_M:Non-Selected`),
    `M(S)` = mean(`tbv_M:Selected`)
  ) %>%
  pivot_longer(cols = Female:`M(S)`,
               names_to = "type",
               values_to = "Value") %>%
  dplyr::mutate(across(type, factor)) %>%
  ggplot(aes(y = Value, x = generation, linetype = path, colour = path)) +
  facet_wrap(~type) +
  geom_line() +
  scale_color_manual(values= c("red", "#9597a1", "blue")) +
  labs(linetype = "Path:", colour = "Path:") +
  xlab("Generation") +
  ylab("Contribution to Genetic Mean") +
  theme_bw() +
  theme(
    legend.position = "top"
  )
ggsave("./Analysis/Figures/ContbyPathPheno.pdf", width = 10, height = 4)


m2$tbv %>%
  filter(generation == "39") %>%
  dplyr::mutate(path = recode(path, 
                              `F`="F", 
                              `M:Non-Selected`= "M(NS)",
                              `M:Selected` = "M(S)")) %>%
  select(ind:path, tbv:`tbv_M:Selected`) %>%
  pivot_longer(cols = `tbv_w`:`tbv_M:Selected`,
               names_to = "Contribution",
               values_to = "Value") %>%
  dplyr::mutate(
    Contribution = recode(
      Contribution,tbv_w = "Mendelian sampling terms",
      `tbv_M:Selected` = "Male partition",
      tbv_F = "Female partition")) %>%
  filter(Contribution %in% c("Mendelian sampling terms",
                             "Male partition",
                             "Female partition")) %>%
  ggplot(
    aes( x = Value, fill = Contribution)
  ) +
  facet_wrap(~path) +
  geom_histogram(aes(y = ..density..),
                 bins = 40, position="identity", alpha=0.5,
                 color="black", ) +
  scale_fill_manual(values= c("red", "blue", "gray")) +
  xlab("Genetic gain") +
  ggtitle("Medium accuracy") +
  theme_bw(base_size = 14) +
  xlim(-2.5, 17)+
  labs(fill = "Contribution of:") +
  theme(
    legend.position = "top"
  )

ggsave("./Analysis/Figures/lastGenSelPheno.pdf", width = 8, height = 4)

#=======================================================================
# Summary - Genetic Mean
#=======================================================================
SMean <- wrapperSummary(m1 = m1, m2 = m2, ncols = 1500,
                        levels=c("F", "M(N)", "M(S)", "Sum"))
save(SMean, file = "./Analysis/Results/SMean_TBV_ValidationPheno.rds")

#=======================================================================
# Summary - Genetic Variance
#=======================================================================
SVar <- wrapperSummary(m1 = m1, m2 = m2, ncols = 1500,
                       FUN = var, cov = TRUE,
                       levels=c("F", "F:M(N)", "F:M(S)", "M(N)", 
                                "M(N):M(S)", "M(S)", "Sum"))
save(SVar, file = "./Analysis/Results/SVar_TBV_ValidationPheno.rds")
#load("./Analysis/Results/SMean_TBV_ValidationPheno.rds")
#load("./Analysis/Results/SVar_TBV_ValidationPheno.rds")

#=======================================================================
# Plot
#=======================================================================
# Mean
g1 <- SMean$SummaryEBV %>%
  ggplot(aes(y = Median, x = generation),
         size = 0.1) +
  facet_wrap(~path)+
  geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = path), alpha = 0.2) + 
  scale_fill_manual(values= c("red", "#9597a1", "blue", "black")) +
  geom_line(aes(colour = path), alpha = 0.8) +
  geom_line(data = SMean$SummaryTrue, 
            aes(y = Median, x = generation, colour = path), 
            linetype = 2, alpha = 0.8) +
  scale_color_manual(values= c("red", "#9597a1", "blue", "black")) +
  geom_vline(xintercept = 0, linetype = 2, alpha = 0.3) +
  ylab("Genetic Mean") +
  xlab("Generation") +
  labs(fill = "Selection path:", colour  = "Selection path:" ) +
  ylim(-1, 23) +
  theme_bw(base_size = 15) + 
  ggtitle("Medium accuracy")+
  theme(legend.position = "none")
g1

ggsave("./Analysis/Figures/gibbsGenMeanValidationPheno.pdf", 
       plot = g1, width = 10, height = 7)

g1 <- SMean$mstEBV %>%
  ggplot(aes(y = Median, x = generation),
         size = 0.1) +
  facet_wrap(~path)+
  geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = path), alpha = 0.2) + 
  scale_fill_manual(values= c("red", "#9597a1", "blue", "black")) +
  geom_line(aes(colour = path), alpha = 0.8) +
  geom_line(data = SMean$mstTrue, 
            aes(y = mstTrue, x = generation, colour = path), 
            linetype = 2, alpha = 0.8) +
  scale_color_manual(values= c("red", "#9597a1", "blue", "black")) +
  geom_vline(xintercept = 0, linetype = 2, alpha = 0.3) +
  ylab("Contribution of Mendelian Sampling Term to Genetic Mean") +
  xlab("Generation") +
  labs(fill = "Selection path:", colour  = "Selection path:" ) +
  theme_bw(base_size = 15) + 
  ylim(-1, 1.8) +
  ggtitle("Medium accuracy")+
  theme(legend.position = "none")
g1
ggsave("./Analysis/Figures/gibbsGenMeanValidationPhenoMST.pdf", 
       plot = g1, width = 10, height = 7)


# Variance
g2 <- SVar$SummaryEBV %>%
  ggplot(aes(y = Median, x = generation),
         size = 0.1) +
  facet_wrap(~path)+
  geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = path), alpha = 0.2) + 
  scale_fill_manual(values= c("red", "#823175","#4b7000ff","#9597a1", "#c4852b", "blue", "black")) +
  geom_line(aes(colour = path), alpha = 0.8) +
  geom_line(data = SVar$SummaryTrue, aes(y = Median, x = generation, colour = path), linetype = 2, alpha = 0.8) +
  scale_color_manual(values= c("red", "#823175","#4b7000ff","#9597a1", "#c4852b", "blue", "black")) +
  geom_vline(xintercept = 0, linetype = 2, alpha = 0.3) +
  ylab("Genetic Variance") +
  xlab("Generation") +
  labs(fill = "Selection path:", colour  = "Selection path:" ) +
  ylim(-1, 1.2) +
  ggtitle("Medium accuracy")+
  theme_bw(base_size = 15) + 
  theme(legend.position = "none")
g2
ggsave("./Analysis/Figures/gibbsGenVarValidationPheno.pdf", 
       plot = g2, width = 10, height = 7)

