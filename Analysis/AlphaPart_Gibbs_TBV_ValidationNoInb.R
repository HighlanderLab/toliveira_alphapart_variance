########################################################################
## Project: toliveira_alphapart_variance
## Script purpose: analyse using AlphaPart
## Date: 2022-01-18
## Author: Thiago de Paula Oliveira
########################################################################
#=======================================================================
# Loading packages and functions
#=======================================================================
source("./Analysis/CommonFunctions.R")

#=======================================================================
# Organizing gibbs1f90 database
#=======================================================================
ebvPath <- "./Analysis/ModelOutputTBVNoInb/all_solutions"
pedPath <- "./Analysis/ModelOutputTBVNoInb/renadd01.ped"
dataPath <- "./Results/Scenario2_replicate1.rds"
dataEBV <- wrapperData(dataPath, pedPath, ebvPath)

#=======================================================================
# test
#=======================================================================
summMean <- dataEBV$data %>%
  group_by(generation, path) %>%
  summarise(meanTBV = mean(tbv),
            across(s1:s1500, ~ mean(.x, na.rm = TRUE)))
summMean$meanEBV <- rowMeans(summMean[,dataEBV$colNames]) 


p1 <- summMean %>%
  dplyr::mutate(path = recode(path, 
                              `F`="F", 
                              `M:Non-Selected`= "M(NS)",
                              `M:Selected` = "M(S)")) %>%
  ggplot(aes(x = meanTBV, y = meanEBV)) +
  facet_wrap(~path)+
  geom_point() +
  geom_abline(intercept = 0, slope = 1)


summVar <- dataEBV$data %>%
  group_by(generation, path) %>%
  summarise(varTBV = var(tbv),
            across(s1:s1500, ~ var(.x, na.rm = TRUE)))
summVar$varEBV <- rowMeans(summVar[,dataEBV$colNames]) 


p2 <- summVar %>%
  dplyr::mutate(path = recode(path, 
                              `F`="F", 
                              `M:Non-Selected`= "M(NS)",
                              `M:Selected` = "M(S)")) %>%
  ggplot(aes(x = varTBV, y = varEBV)) +
  facet_wrap(~path)+
  geom_point() +
  geom_abline(intercept = 0, slope = 1)

p1/p2

ggsave("./Analysis/Results/EBV_TBV_ScenarioTBVNoInb.pdf", width = 9, height = 6)
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
saveRDS(m1, "./Analysis/Results/AlphaPart_EBV_ValidationTBVNoInb.rds")
saveRDS(m2, "./Analysis/Results/AlphaPart_TBV_ValidationTBVNoInb.rds")


#=======================================================================
# Summary - Genetic Mean
#=======================================================================
SMean <- wrapperSummary(m1 = m1, m2 = m2, ncols = 1500,
                        levels=c("F", "M(N)", "M(S)", "Sum"))
save(SMean, file = "./Analysis/Results/SMean_TBV_ValidationTBVNoInb.rds")

#=======================================================================
# Summary - Genetic Variance
#=======================================================================
SVar <- wrapperSummary(m1 = m1, m2 = m2, ncols = 1500,
                       FUN = var, cov = TRUE,
                       levels=c("F", "F:M(N)", "F:M(S)", "M(N)", 
                                "M(N):M(S)", "M(S)", "Sum"))
save(SVar, file = "./Analysis/Results/SVar_TBV_ValidationTBVNoInb.rds")
#load("./Analysis/Results/SMean_TBV_ValidationTBVNoInb.rds")
#load("./Analysis/Results/SVar_TBV_ValidationTBVNoInb.rds")

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
  theme_bw(base_size = 15) +   ggtitle("High accuracy")+
  theme(legend.position = "none")
g1

ggsave("./Analysis/Figures/gibbsGenMeanValidationTBVNoInb.pdf", 
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
  ggtitle("High accuracy")+  theme_bw(base_size = 15) + 
  theme(legend.position = "none")
g1
ggsave("./Analysis/Figures/gibbsGenMeanValidationTBVNoInbMST.pdf", 
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
  theme_bw(base_size = 15) +   ggtitle("High accuracy")+
  theme(legend.position = "none")
ggsave("./Analysis/Figures/gibbsGenVarValidationTBVNoInb.pdf", 
       plot = g2, width = 10, height = 7)

