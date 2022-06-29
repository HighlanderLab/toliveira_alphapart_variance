########################################################################
## Project: toliveira_alphapart_variance
## Script purpose: Medium Accuracy - accounting for inbreeding
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
  ggtitle("Medium Accuracy") +
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
  ggtitle("High Accuracy") +
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
tmp <- as.matrix(dataEBV$data[dataEBV$data$father=="0" & 
                                dataEBV$data$mother=="0", -c(1:6,8:11)])
y <- as.matrix(dataEBV$data[, -c(1:6,8:11)])
f <- function(x) {
  sd(x, na.rm = TRUE)
}
scale <- apply(tmp, 2L, f)
y <- y - colMeans(tmp)[1]
y <- y/scale[1]
y <- cbind(y, dataEBV$data[,c(6,9)])

summMean <- y %>%
  group_by(generation, path) %>%
  summarise(meanTBV = mean(tbv),
            across(s1:s1500, ~ mean(.x, na.rm = TRUE)))
summMean$meanEBV <- rowMeans(summMean[,dataEBV$colNames]) 

# Concordance
ccc <- CCC(summMean$meanTBV, summMean$meanEBV)
round(ccc$rho.c,3)

# Pearson Cor.
round(ccc$rho.c[1]/ccc$C.b,3)

# C_b
round(ccc$C.b,3)


#plot
p1 <- summMean %>%
  dplyr::mutate(path = recode(path, 
                              `F`="F", 
                              `M:Non-Selected`= "M(NS)",
                              `M:Selected` = "M(S)")) %>%
  dplyr::mutate(generation = generation-20) %>%
  ggplot(aes(x = meanTBV, y = meanEBV, colour = generation)) +
  facet_wrap(~path)+
  geom_point() +
  xlab("True genetic mean per generation") +
  ylab("Estimated genetic mean per generation") +
  geom_abline(intercept = 0, slope = 1) +
  ggtitle("Medium Accuracy") +
  labs(colour = "Generation:") +
  theme_bw(base_size = 12)+
  theme(legend.position = "top")


summVar <- y %>%
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
round(ccc$rho.c,3)

# Pearson Cor.
round(ccc$rho.c[1]/ccc$C.b,3)

# C_b
round(ccc$C.b,3)

# plot
p2 <- summVar  %>%
  dplyr::mutate(generation = generation-20) %>%
  ggplot(aes(x = varTBV, y = varEBV, colour = generation)) +
  facet_wrap(~path)+
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  xlab("True genetic variance per generation") +
  ylab("Estimated genetic variance per generation") +
  geom_abline(intercept = 0, slope = 1) +
  labs(colour = "Generation:") +
  theme_bw(base_size = 12)+
  theme(legend.position = "none")

p1/p2
ggsave("./Analysis/Figures/EBV_TBV_ScenarioPheno.pdf", 
       width = 9, height = 9)
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

#m1 <- readRDS("./Analysis/Results/AlphaPart_EBV_ValidationPheno.rds")
#m2 <- readRDS("./Analysis/Results/AlphaPart_TBV_ValidationPheno.rds")


rm(dataEBV)
gc()

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
  xlab("Density plot of breeding value partitions") +
  ylab("") +
  ggtitle("Medium Accuracy") +
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
  dplyr::mutate(
    path = factor(path, levels = levels(SMean$SummaryEBV$path)[4:1])
    ) %>%
  ggplot(aes(y = Median, x = generation),
         size = 0.1) +
  facet_wrap(~path)+
  geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = path), alpha = 0.2) + 
  scale_fill_manual(values= c("black", "blue", "#9597a1", "red")) +
  geom_line(aes(colour = path), alpha = 0.8) +
  geom_line(data = SMean$SummaryTrue, 
            aes(y = Median, x = generation, colour = path), 
            linetype = 2, alpha = 0.8) +
  scale_color_manual(values= c("black", "blue", "#9597a1", "red")) +
  geom_vline(xintercept = 0, linetype = 2, alpha = 0.3) +
  ylab("Genetic Mean") +
  xlab("Generation") +
  labs(fill = "Selection path:", colour  = "Selection path:" ) +
  ylim(-1, 23) +
  theme_bw(base_size = 15) + 
  ggtitle("Medium Accuracy")+
  theme(legend.position = "none")
g1

ggsave("./Analysis/Figures/gibbsGenMeanValidationPheno.pdf", 
       plot = g1, width = 6, height = 6)
ggsave("./Analysis/Figures/gibbsGenMeanValidationPhenoPoster.pdf", 
       plot = g1, width = 6, height = 5)


g1 <- SMean$mstEBV %>%
  dplyr::mutate(
    path = factor(path, levels = levels(SMean$mstEBV$path)[4:1])
  ) %>%
  ggplot(aes(y = Median, x = generation),
         size = 0.1) +
  facet_wrap(~path)+
  geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = path), alpha = 0.2) + 
  scale_fill_manual(values= c("black", "blue", "#9597a1", "red")) +
  geom_line(aes(colour = path), alpha = 0.8) +
  geom_line(data = SMean$mstTrue, 
            aes(y = mstTrue, x = generation, colour = path), 
            linetype = 2, alpha = 0.8) +
  scale_color_manual(values= c("black", "blue", "#9597a1", "red")) +
  geom_vline(xintercept = 0, linetype = 2, alpha = 0.3) +
  ylab("Contribution of Mendelian Sampling Term to Genetic Mean") +
  xlab("Generation") +
  labs(fill = "Selection path:", colour  = "Selection path:" ) +
  theme_bw(base_size = 15) + 
  ylim(-1, 1.8) +
  ggtitle("Medium Accuracy")+
  theme(legend.position = "none")
g1
ggsave("./Analysis/Figures/gibbsGenMeanValidationPhenoMST.pdf", 
       plot = g1, width = 6, height = 6)

ggsave("./Analysis/Figures/gibbsGenMeanValidationPhenoMSTPoster.pdf", 
       plot = g1, width = 6, height = 5)


# Variance
g2 <- SVar$SummaryEBV %>%
  dplyr::mutate(
    path = factor(path, levels = levels(SVar$SummaryEBV$path)[c(7,6,4,1,5,3,2)])
  ) %>%
  ggplot(aes(y = Median, x = generation),
         size = 0.1) +
  facet_wrap(~path)+
  geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = path), alpha = 0.2) + 
  scale_fill_manual(values= c("black", "blue","#9597a1", "red", "#c4852b","#4b7000ff","#823175")) +
  geom_line(aes(colour = path), alpha = 0.8) +
  geom_line(data = SVar$SummaryTrue, aes(y = Median, x = generation, colour = path), linetype = 2, alpha = 0.8) +
  scale_color_manual(values= c("black", "blue","#9597a1", "red", "#c4852b","#4b7000ff","#823175")) +
  geom_vline(xintercept = 0, linetype = 2, alpha = 0.3) +
  ylab("Genetic Variance") +
  xlab("Generation") +
  labs(fill = "Selection path:", colour  = "Selection path:" ) +
  ylim(-1, 1.2) +
  ggtitle("Medium Accuracy")+
  theme_bw(base_size = 15) + 
  theme(legend.position = "none")
g2
ggsave("./Analysis/Figures/gibbsGenVarValidationPheno.pdf", 
       plot = g2, width = 8, height = 8)

ggsave("./Analysis/Figures/gibbsGenVarValidationPhenoPoster.pdf", 
       plot = g2, width = 9, height = 7)
