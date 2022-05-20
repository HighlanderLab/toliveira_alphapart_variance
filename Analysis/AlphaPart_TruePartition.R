########################################################################
## Project: jobsteter_alphapart_variance
## Script purpose: Dairy Cattle BP Analysis
## Date: 2021-11-30
## Author: Thiago de Paula Oliveira
########################################################################
# Set working directory
rm(list=ls())
#=======================================================================

#=======================================================================
# Packges
#=======================================================================
#devtools::install_github("AlphaGenes/AlphaPart")
library(AlphaPart)
library(tidyverse)
library(tidyr)
library(patchwork)
#=======================================================================
# Reading and organizing Scenario 1
#=======================================================================
dataS1 <-  readRDS("./Results/Scenario1_replicate1.rds") %>%
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
tmpScale <- dataS1 %>% 
  filter(father == 0 & mother == 0) %>%
  summarise(
    mean = mean(tbv),
    sd = sd(tbv)
  )
dataS1$tbv <- (dataS1$tbv - tmpScale$mean)/tmpScale$sd
dataS1$path <- droplevels(dataS1$sex:dataS1$type)
levels(dataS1$path)[1] <- "F"

#-----------------------------------------------------------------------
# Sex
#-----------------------------------------------------------------------
m1Sex <- AlphaPart(dataS1, colId = "ind", colFid = "father", 
                   colMid = "mother", colBV = "tbv", colPath = "sex")
S1Sex <- summary(m1Sex, by = "generation")

S1VarSex <- summary(m1Sex, by = "generation", FUN = var, cov = TRUE)

#-----------------------------------------------------------------------
# Sex:Type
#-----------------------------------------------------------------------
m1 <- AlphaPart(dataS1, colId = "ind", colFid = "father", 
                colMid = "mother", colBV = "tbv", colPath = "path")
S1 <- summary(m1, by = "generation")

S1Var <- summary(m1, by = "generation", FUN = var, cov = TRUE)

#-----------------------------------------------------------------------
# Correlation
#-----------------------------------------------------------------------
m1$tbv$group <- m1$tbv$sex:m1$tbv$type

m1$tbv %>%
  group_by(generation) %>%
  filter(type == "Selected") %>%
  ggplot() +
  geom_point(aes(y = tbv_w, x = generation), shape =2, size = 1.5, alpha = 0.4)+
  geom_point(aes(y = `tbv_M:Selected`, x = generation, color = tbv_F-`tbv_M:Selected`, 
                 size = tbv_F-`tbv_M:Selected`), alpha = 0.6) +
  labs(color = "F - M(S):", size = "F - M(S):") +
  geom_vline(xintercept = 0, alpha = 0.3, linetype = 2) +
  scale_y_continuous(name = 'M(S) contribution',
                     sec.axis = sec_axis(~., name = 'Mendelian Sampling')) +
  xlab("Generation") +
  theme_bw(base_size = 14) +
  theme(legend.position = "none")
#ggsave("./Analysis/Figures/correlation_pheno.pdf",width = 8, height = 6)

#=======================================================================
# Reading and organizing Scenario 2
#=======================================================================
dataS2 <- readRDS("./Results/Scenario2_replicate1.rds") %>%
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
tmpScale <- dataS2 %>% 
  filter(father == 0 & mother == 0) %>%
  summarise(
    mean = mean(tbv),
    sd = sd(tbv)
  )
dataS2$tbv <- (dataS2$tbv - tmpScale$mean)/tmpScale$sd
dataS2$path <- droplevels(dataS2$sex:dataS2$type)
levels(dataS2$path)[1] <- "F"

#-----------------------------------------------------------------------
# Sex
#-----------------------------------------------------------------------
m1Sex <- AlphaPart(dataS2, colId = "ind", colFid = "father", 
                   colMid = "mother", colBV = "tbv", colPath = "sex")
S2Sex <- summary(m1Sex, by = "generation")

S2VarSex <- summary(m1Sex, by = "generation", FUN = var, cov = TRUE)

#-----------------------------------------------------------------------
# Sex:Type
#-----------------------------------------------------------------------
m1 <- AlphaPart(dataS2, colId = "ind", colFid = "father", 
                colMid = "mother", colBV = "tbv", colPath = "path")
S2 <- summary(m1, by = "generation")

S2Var <- summary(m1, by = "generation", FUN = var, cov = TRUE)

m1$tbv %>%
  group_by(generation) %>%
  filter(type == "Selected") %>%
  ggplot() +
  geom_point(aes(y = tbv_w, x = generation), shape =2, size = 1.5, alpha = 0.4)+
  geom_point(aes(y = `tbv_M:Selected`, x = generation, color = tbv_F-`tbv_M:Selected`, 
                 size = tbv_F-`tbv_M:Selected`), alpha = 0.6) +
  labs(color = "F - M(S):", size = "F - M(S):") +
  geom_vline(xintercept = 0, alpha = 0.3, linetype = 2) +
  scale_y_continuous(name = 'M(S) contribution',
                     sec.axis = sec_axis(~., name = 'Mendelian Sampling')) +
  xlab("Generation") +
  theme_bw(base_size = 14) +
  theme(legend.position = "none")

#ggsave("./Analysis/Figures/escorrelation_tbv.pdf",width = 8, height = 6)


########################################################################
# Sex
########################################################################
#-----------------------------------------------------------------------
# Plot Both - Mean
#-----------------------------------------------------------------------
tmp1 <- S1Sex$tbv
tmp1$Selection <- "Medium accuracy"

tmp2 <- S2Sex$tbv
tmp2$Selection <- "High accuracy"

Summ <- rbind(tmp1,tmp2)
Summ$Selection <- factor(Summ$Selection, levels = c("Medium accuracy", "High accuracy"))
levels(Summ$Selection)

rm(tmp1,tmp2)

g1 <- Summ %>%
  mutate(across(Selection, factor, levels=c("Medium accuracy","High accuracy"))) %>%
  ggplot(aes(y = Sum, x = generation),
         size = 0.1) +
  facet_wrap(~Selection)+
  geom_line(aes(colour = "Sum")) +
  geom_line(aes(y = `F`, x = generation, colour = "F"),
            alpha = 0.8) +
  geom_line(aes(y = `M`, x = generation, colour = "M"),
            alpha = 0.8) +
  scale_color_manual(values= c("red", "blue", "black")) +
  geom_vline(xintercept = 0, linetype = 2, alpha = 0.3) +
  ylab("Genetic Mean") +
  xlab("Generation") +
  labs(colour = "Path:") +
  theme_bw(base_size = 18) + 
  theme(legend.position = "none")
g2 <- g1 + theme(legend.position = "none")
g2

#-----------------------------------------------------------------------
# Plot Both - Var
#-----------------------------------------------------------------------
tmp1 <- S1VarSex$tbv
tmp1$Selection <- "Medium accuracy"

tmp2 <- S2VarSex$tbv
tmp2$Selection <- "High accuracy"


Summ <- rbind(tmp1,tmp2)
Summ$Selection <- factor(Summ$Selection, levels = c("Medium accuracy", "High accuracy"))
levels(Summ$Selection)
rm(tmp1,tmp2)

g4 <- Summ %>%
  ggplot(aes(y = Sum, x = generation)) +
  facet_wrap(~Selection)+
  geom_line(aes(colour = "Sum")) +
  geom_line(aes(y = `F`, x = generation, colour = "F"),
            alpha = 0.8) +
  geom_line(aes(y = `FM`, x = generation, colour = "F:M"),
            size = 0.2, alpha = 0.8) +
  scale_color_manual(values= c("red", "orange", "blue", "black")) +
  geom_line(aes(y = `M`, x = generation, colour = "M"),
            alpha = 0.8) +
  geom_vline(xintercept = 0, linetype = 2, alpha = 0.3) +
  ylab("Genetic Variance") +
  labs(colour = "Path:") +
  xlab("Generation") +
  theme_bw(base_size = 18)
g5 <- g4 + theme(legend.position = "none") 

g5

g2/g5
ggsave("./Analysis/Figures/AlphaPartSexP.pdf",width = 8, height = 8)

########################################################################
# Sex:Type
########################################################################
#-----------------------------------------------------------------------
# Plot Both - Mean
#-----------------------------------------------------------------------
tmp1 <- S1$tbv
tmp1$Selection <- "Medium accuracy"

tmp2 <- S2$tbv
tmp2$Selection <- "High accuracy"

Summ <- rbind(tmp1,tmp2)
Summ$Selection <- factor(Summ$Selection, levels = c("Medium accuracy", "High accuracy"))
levels(Summ$Selection)
saveRDS(Summ, "./Analysis/Results/SummaryMean.rds")
rm(tmp1,tmp2)

g1 <- Summ %>%
  ggplot(aes(y = Sum, x = generation),
         size = 0.1) +
  facet_wrap(~Selection)+
  geom_line(aes(colour = "Sum")) +
  geom_line(aes(y = `F`, x = generation, colour = "F"),
            alpha = 0.8) +
  geom_line(aes(y = `M:Selected`, x = generation, colour = "M(S)"),
            alpha = 0.8) +
  geom_line(aes(y = `M:Non-Selected`, x = generation, colour = "M(N)"),
            alpha = 0.8) +
  scale_color_manual(values= c("red", "#9597a1", "blue", "black")) +
  geom_vline(xintercept = 0, linetype = 2, alpha = 0.3) +
  ylab("Genetic Mean") +
  xlab("Generation") +
  labs(colour = "Path:") +
  theme_bw(base_size = 18) + 
  theme(legend.position = "none")

g2 <- g1 + theme(legend.position = "none") 
g2

#ggsave("Figur./Analysis/Figures/esAlphaPartMean.pdf",width = 14, height = 7)
#-----------------------------------------------------------------------
# Plot Both - Var
#-----------------------------------------------------------------------
tmp1 <- S1Var$tbv
tmp1$Selection <- "Medium accuracy"

tmp2 <- S2Var$tbv
tmp2$Selection <- "High accuracy"


Summ <- rbind(tmp1,tmp2)
Summ$Selection <- factor(Summ$Selection, levels = c("Medium accuracy", "High accuracy"))
levels(Summ$Selection)
saveRDS(Summ, "./Analysis/Results/SummaryVar.rds")
rm(tmp1,tmp2)

Summ %>%
  group_by(generation, Selection) %>%
  summarise(
    cor = `FM:Selected` / sqrt(`F` * `M:Selected`)
  ) %>%
  ggplot(aes(y = cor, x = generation)) +
  facet_wrap(~Selection) +
  geom_line() +
  geom_hline(yintercept = 0, linetype = 2, alpha = 0.2)+
  geom_vline(xintercept = 0, linetype = 2, alpha = 0.2) +
  ylab("Correlation between the contribution of F and M(S) to genetic variance") +
  xlab("Generation") +
  theme_bw(base_size = 14)
ggsave("./Analysis/Figures/AlphaPartCor.pdf",width = 8, height = 7)

g4 <- Summ %>%
  ggplot(aes(y = Sum, x = generation)) +
  facet_wrap(~Selection)+
  geom_line(aes(colour = "Sum")) +
  geom_line(aes(y = `F`, x = generation, colour = "F"),
            alpha = 0.8) +
  geom_line(aes(y = `FM:Selected`, x = generation, colour = "F:M(S)"),
            size =0.5, alpha =0.8) +
  geom_line(aes(y = `FM:Non-Selected`, x = generation, colour = "F:M(N)"),
            size =0.5, alpha =0.6) +
  geom_line(aes(y = `M:Non-SelectedM:Selected`, x = generation, 
                colour = "M(S):M(N)"),
            size =0.5, alpha =0.6) +
  geom_line(aes(y = `M:Selected`, x = generation, colour = "M(S)"),
            alpha = 0.8) +
  geom_line(aes(y = `M:Non-Selected`, x = generation, colour = "M(N)"),
            size =0.5, alpha =0.8) +
  scale_color_manual(values= c("red", "#823175","#4b7000ff", "#9597a1", "blue", "#c4852b", "black")) +
  geom_vline(xintercept = 0, linetype = 2, alpha = 0.3) +
  ylab("Genetic Variance") +
  xlab("Generation") +
  labs(colour = "Path:") +
  theme_bw(base_size = 18)
g5 <- g4 + theme(legend.position = "none") 
g5

g2/g5
ggsave("./Analysis/Figures/AlphaPartSexSel.pdf",width = 8, height = 8)

