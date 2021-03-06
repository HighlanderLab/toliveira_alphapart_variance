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
ebvPath <- "./Analysis/ModelOutputTBV/all_solutions"
pedPath <- "./Analysis/ModelOutputTBV/renadd01.ped"
dataPath <- "./Results/Scenario2_replicate1.rds"
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
                              `M:Selected` = "M(S)"))  %>%
  dplyr::mutate(generation = generation-20) %>%
  ggplot(aes(x = meanTBV, y = meanEBV, colour = generation)) +
  facet_wrap(~path)+
  geom_point() +
  xlab("True genetic mean per generation") +
  ylab("Estimated genetic mean per generation") +
  geom_abline(intercept = 0, slope = 1) +
  ggtitle("High accuracy") +
  labs(colour = "Generation:") +
  theme_bw(base_size = 12)+
  theme(legend.position = "top")


summVar <- y %>%
  group_by(generation, path) %>%
  summarise(varTBV = var(tbv),
            across(s1:s1500, ~ var(.x, na.rm = TRUE)))
summVar$varEBV <- rowMeans(summVar[,dataEBV$colNames]) 

# Concordance
ccc <- CCC(summVar$varTBV, summVar$varEBV)
round(ccc$rho.c,2)

# Pearson Cor.
round(ccc$rho.c[1]/ccc$C.b,2)

# C_b
round(ccc$C.b,2)

# plot
p2 <- summVar  %>%
  dplyr::mutate(generation = generation-20) %>%
  dplyr::mutate(path = recode(path, 
                              `F`="F", 
                              `M:Non-Selected`= "M(NS)",
                              `M:Selected` = "M(S)")) %>%
  ggplot(aes(x = varTBV, y = varEBV, colour = generation)) +
  facet_wrap(~path)+
  geom_point() +
  xlab("True genetic variance per generation") +
  ylab("Estimated genetic variance per generation") +
  geom_abline(intercept = 0, slope = 1) +
  labs(colour = "Generation:") +
  theme_bw(base_size = 12)+
  theme(legend.position = "none")

p1/p2

ggsave("./Analysis/Figures/EBV_TBV_ScenarioTBV.pdf", width = 9, height = 9)
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
saveRDS(m1, "./Analysis/Results/AlphaPart_EBV_ValidationTBV.rds")
saveRDS(m2, "./Analysis/Results/AlphaPart_TBV_ValidationTBV.rds")


m1 <- readRDS("./Analysis/Results/AlphaPart_EBV_ValidationTBV.rds")
m2 <- readRDS("./Analysis/Results/AlphaPart_TBV_ValidationTBV.rds")

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
ggsave("./Analysis/Figures/ContbyPathTBV.pdf", width = 10, height = 4)


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
  ggtitle("High accuracy") +
  theme_bw(base_size = 14) +  
  xlim(-2.5, 17) +
  labs(fill = "Contribution of:") +
  theme(
    legend.position = "top"
  )
ggsave("./Analysis/Figures/lastGenSelTBV.pdf", width = 8, height = 4)

#=======================================================================
# Summary - Genetic Mean
#=======================================================================
SMean <- wrapperSummary(m1 = m1, m2 = m2, ncols = 1500,
                        levels=c("F", "M(N)", "M(S)", "Sum"))
save(SMean, file = "./Analysis/Results/SMean_TBV_ValidationTBV.rds")

#=======================================================================
# Summary - Genetic Variance
#=======================================================================
SVar <- wrapperSummary(m1 = m1, m2 = m2, ncols = 1500,
                       FUN = var, cov = TRUE,
                       levels=c("F", "F:M(N)", "F:M(S)", "M(N)", 
                                "M(N):M(S)", "M(S)", "Sum"))
save(SVar, file = "./Analysis/Results/SVar_TBV_ValidationTBV.rds")
#load("./Analysis/Results/SMean_TBV_ValidationTBV.rds")
#load("./Analysis/Results/SVar_TBV_ValidationTBV.rds")

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
  ggtitle("High accuracy")+
  theme_bw(base_size = 15) + 
  theme(legend.position = "none")
g1

ggsave("./Analysis/Figures/gibbsGenMeanValidationTBV.pdf", 
       plot = g1, width = 6, height = 6)

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
  ggtitle("High accuracy")+
  ylim(-1, 1.8) +
  theme(legend.position = "none")
g1
ggsave("./Analysis/Figures/gibbsGenMeanValidationTBVMST.pdf", 
       plot = g1, width = 6, height = 6)


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
  ggtitle("High accuracy")+
  ylim(-1, 1.2) +
  theme_bw(base_size = 15) + 
  theme(legend.position = "none")
g2
ggsave("./Analysis/Figures/gibbsGenVarValidationTBV.pdf", 
       plot = g2, width = 8, height = 8)

