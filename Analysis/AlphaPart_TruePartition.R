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
library(ggridges)
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
# Groups
#-----------------------------------------------------------------------
m1$tbv %>%
  ggplot(aes(y = as.factor(generation), x= tbv_F)) +
  geom_density_ridges(
    aes(fill = "F - Non-Selected", linetype = "F - Non-Selected"),
    alpha = .4, point_alpha = 1, rel_min_height = 0.01
  ) +
  geom_density_ridges(
    aes(y = as.factor(generation), x= `tbv_M:Non-Selected`, fill = "M - Non-Selected",
        linetype = "M - Non-Selected"),
    alpha = .4, point_alpha = 1, rel_min_height = 0.01
  ) +
  geom_density_ridges(
    aes(y = as.factor(generation), x= `tbv_M:Selected`, fill = "M - Selected",
        linetype = "M - Selected"),
    alpha = .4, point_alpha = 1, rel_min_height = 0.01
  ) +
  geom_density_ridges(
    aes(y = as.factor(generation), x= `tbv`,
        fill = "Sum", linetype = "Sum"),
    alpha = .4, point_alpha = 1, rel_min_height = 0.01
  ) +
  ylab("Generation") +
  xlab("Density plot of breeding value partitions") +
  labs(fill = "Path:", linetype = "Path:") +
  theme_bw(base_size = 20) +
  xlim(c(-5,25))+
  ggtitle("Medium accuracy") +
  theme(
    legend.position = "top"
  ) 
ggsave("./Analysis/Figures/Contribution_pheno.pdf",width = 11, height = 11)

#-----------------------------------------------------------------------
# MST and PA
#-----------------------------------------------------------------------
m1$tbv %>%
  ggplot(aes(y = as.factor(generation), x= tbv_pa)) +
  geom_density_ridges(
    aes(fill = type, linetype = "PA"),
    alpha = .4, point_alpha = 1, rel_min_height = 0.01
  ) +
  geom_density_ridges(
    aes(y = as.factor(generation), x= tbv_w, fill = type,
        linetype = "MST"),
    alpha = .4, point_alpha = 1, rel_min_height = 0.01
  ) +
  facet_wrap(~sex) +
  ylab("Year") +
  xlab("Parent Average") +
  ggtitle("Medium accuracy") +
  labs(fill = "Selection paths:", linetype = "Variable:") +
  xlim(-3,23) +
  theme_bw(base_size = 16) +
  theme(
    legend.position = "top"
  ) 
ggsave("./Analysis/Figures/MST_pheno.pdf",width = 9, height = 9)

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
  theme(legend.position = "top")
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


#-----------------------------------------------------------------------
# Groups
#-----------------------------------------------------------------------
m1$tbv %>%
  ggplot(aes(y = as.factor(generation), x= tbv_F)) +
  geom_density_ridges(
    aes(fill = "F - Non-Selected", linetype = "F - Non-Selected"),
    alpha = .4, point_alpha = 1, rel_min_height = 0.01
  ) +
  geom_density_ridges(
    aes(y = as.factor(generation), x= `tbv_M:Non-Selected`, fill = "M - Non-Selected",
        linetype = "M - Non-Selected"),
    alpha = .4, point_alpha = 1, rel_min_height = 0.01
  ) +
  geom_density_ridges(
    aes(y = as.factor(generation), x= `tbv_M:Selected`, fill = "M - Selected",
        linetype = "M - Selected"),
    alpha = .4, point_alpha = 1, rel_min_height = 0.01
  ) +
  geom_density_ridges(
    aes(y = as.factor(generation), x= `tbv`,
        fill = "Sum", linetype = "Sum"),
    alpha = .4, point_alpha = 1, rel_min_height = 0.01
  ) +
  ylab("Generation") +
  xlab("Density plot of breeding value partitions") +
  labs(fill = "Path:", linetype = "Path:") +
  ggtitle("High accuracy") +
  theme_bw(base_size = 20) +
  xlim(c(-5,25))+
  theme(
    legend.position = "top"
  ) 
ggsave("./Analysis/Figures/Contribution_tbv.pdf",width = 11, height = 11)


#-----------------------------------------------------------------------
# MST and PA
#-----------------------------------------------------------------------
m1$tbv %>%
  ggplot(aes(y = as.factor(generation), x= tbv_pa)) +
  geom_density_ridges(
    aes(fill = type, linetype = "PA"),
    alpha = .4, point_alpha = 1, rel_min_height = 0.01
  ) +
  geom_density_ridges(
    aes(y = as.factor(generation), x= tbv_w, fill = type,
        linetype = "MST"),
    alpha = .4, point_alpha = 1, rel_min_height = 0.01
  ) +
  facet_wrap(~sex) +
  ylab("Year") +
  xlab("Parent Average") +
  ggtitle("High accuracy") +
  labs(fill = "Selection paths:", linetype = "Variable:") +
  xlim(-3,23) +
  theme_bw(base_size = 16) +
  theme(
    legend.position = "top"
  )
ggsave("./Analysis/Figures/MST_tbv.pdf",width = 9, height = 9)

#-----------------------------------------------------------------------
# M(S) contribution and MST
#-----------------------------------------------------------------------
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
  theme(legend.position = "top")

#ggsave("./Analysis/Figures/escorrelation_tbv.pdf",width = 8, height = 6)

m1$tbv %>%
  filter(generation >0) %>%
  ggplot(aes(y = tbv_F, x = `tbv_M:Selected`, colour = path)) +
  geom_point() +
  facet_wrap(~generation) +
  ylab("Contribution from Female") +
  xlab("Contribution from Male") +
  theme_bw() +
  theme(
    legend.position = "top"
  )
ggsave("./Analysis/Figures/negative_cor_generation20.pdf",width = 15, height = 12)
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
  scale_linetype_manual(
    values = c("solid", "longdash", "dashed", "dotted"))+
  geom_line() +
  geom_line(aes(y = `F`, x = generation),
            colour = "red", alpha = 0.8) +
  geom_line(aes(y = `M`, x = generation),
            colour = "blue", alpha = 0.8) +
  geom_vline(xintercept = 0, linetype = 2, alpha = 0.3) +
  ylab("Genetic Mean") +
  xlab("Generation") +
  xlim(-20,23)+
  theme_bw(base_size = 18) + 
  theme(legend.position = "none")
g2 <- g1 + theme(legend.position = "top") +
  guides(
    color = guide_legend(override.aes = list(linetype = c(0, 1),
                                             shape = c(16, NA),
                                             color = "black")))
textplot <- tibble(
  Sum = rep(c(mean(Summ[Summ$generation==20 & Summ$Selection=="Medium accuracy",][,3]),
              mean(Summ[Summ$generation==20 & Summ$Selection=="High accuracy",][,3]),
              mean(Summ[Summ$generation==20 & Summ$Selection=="Medium accuracy",][,5]),
              mean(Summ[Summ$generation==20 & Summ$Selection=="High accuracy",][,5]),
              mean(Summ[Summ$generation==20 & Summ$Selection=="Medium accuracy",][,4]),
              mean(Summ[Summ$generation==20 & Summ$Selection=="High accuracy",][,4])),2),
  Selection = as.factor(rep(rep(c("Medium accuracy","High accuracy"),3),2)),
  Scenario = c(rep("IS", 6), rep("WFS", 6)),
  generation = 22,
  label = rep(c(rep("Sum",2), rep("M",2), rep("F",2)),2),
  colour = rep(c(rep("black",2), rep("blue",2), rep("red",2)),2),
  size = 3.6
)

g3 <- g2 + geom_text(
  data = textplot,
  label = textplot$label,
  colour = textplot$colour,
  size = textplot$size
)
g3
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
  geom_line() +
  geom_line(aes(y = `F`, x = generation),
            colour = "red", alpha = 0.8) +
  geom_line(aes(y = `FM`, x = generation),
            colour = "#4b7000ff", size = 0.2, alpha = 0.6) +
  geom_line(aes(y = `M`, x = generation),
            colour = "blue", alpha = 0.8) +
  geom_vline(xintercept = 0, linetype = 2, alpha = 0.3) +
  ylab("Genetic Variance") +
  xlim(-20,23)+
  xlab("Generation") +
  theme_bw(base_size = 18)
g5 <- g4 + theme(legend.position = "top") +
  guides(
    color = guide_legend(override.aes = list(linetype = c(0, 1),
                                             shape = c(16, NA),
                                             color = "black")))

textplot <- tibble(
  Sum = rep(c(mean(Summ[Summ$generation==20 & Summ$Selection=="Medium accuracy",][,3]),
              mean(Summ[Summ$generation==20 & Summ$Selection=="High accuracy",][,3]),
              mean(Summ[Summ$generation==20 & Summ$Selection=="Medium accuracy",][,5]),
              mean(Summ[Summ$generation==20 & Summ$Selection=="High accuracy",][,5]) - 0.04,
              mean(Summ[Summ$generation==20 & Summ$Selection=="Medium accuracy",][,4]),
              mean(Summ[Summ$generation==20 & Summ$Selection=="High accuracy",][,4])+0.02,
              mean(Summ[Summ$generation==20 & Summ$Selection=="Medium accuracy",][,6]),
              mean(Summ[Summ$generation==20 & Summ$Selection=="High accuracy",][,6])),2),
  Selection = as.factor(rep(rep(c("Medium accuracy","High accuracy"),4),2)),
  Scenario = c(rep("IS", 8), rep("WFS", 8)),
  generation = 22,
  label = rep(c(rep("Sum",2), rep("M",2), rep("F",2), rep("F:M",2)),2),
  colour = rep(c(rep("black",2), rep("blue",2), rep("red",2), 
                 rep("#4b7000ff",2)),2),
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
  scale_linetype_manual(
    values = c("solid", "longdash", "dashed", "dotted"))+
  geom_line() +
  geom_line(aes(y = `F`, x = generation),
            colour = "red", alpha = 0.8) +
  geom_line(aes(y = `M:Selected`, x = generation),
            colour = "blue", alpha = 0.8) +
  geom_line(aes(y = `M:Non-Selected`, x = generation),
            colour = "#9597a1", alpha = 0.8) +
  geom_vline(xintercept = 0, linetype = 2, alpha = 0.3) +
  ylab("Genetic Mean") +
  xlab("Generation") +
  xlim(-20,32)+
  theme_bw(base_size = 18) + 
  theme(legend.position = "none")

g2 <- g1 + theme(legend.position = "top") +
  guides(
    color = guide_legend(override.aes = list(linetype = c(0, 1),
                                             shape = c(16, NA),
                                             color = "black")))
g2
textplot <- tibble(
  Sum = rep(c(mean(Summ[Summ$generation==20 & Summ$Selection=="Medium accuracy",][,3]),
              mean(Summ[Summ$generation==20 & Summ$Selection=="High accuracy",][,3]),
              mean(Summ[Summ$generation==20 & Summ$Selection=="Medium accuracy",][,5]),
              mean(Summ[Summ$generation==20 & Summ$Selection=="High accuracy",][,5]),
              mean(Summ[Summ$generation==20 & Summ$Selection=="Medium accuracy",][,4]),
              mean(Summ[Summ$generation==20 & Summ$Selection=="High accuracy",][,4]),
              mean(Summ[Summ$generation==20 & Summ$Selection=="Medium accuracy",][,6]),
              mean(Summ[Summ$generation==20 & Summ$Selection=="High accuracy",][,6])),2),
  Selection = as.factor(rep(rep(c("Medium accuracy","High accuracy"),4),2)),
  Scenario = c(rep("IS", 8), rep("WFS", 8)),
  generation = 23.5,
  label = rep(c(rep("Sum",2), rep("M(N)",2), rep("F",2), rep("M(S)",2)),2),
  colour = rep(c(rep("black",2), rep("#9597a1",2), rep("red",2), 
                 rep("blue",2)),2),
  size = 3.6
)

g3 <- g2 + geom_text(
  data = textplot,
  label = textplot$label,
  colour = textplot$colour,
  size = textplot$size
)
g3

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
  ylab("Correlation between the F and M(S) breeding values partition") +
  xlab("Generation") +
  theme_bw(base_size = 14)
ggsave("./Analysis/Figures/AlphaPartCor.pdf",width = 8, height = 7)

g4 <- Summ %>%
  ggplot(aes(y = Sum, x = generation)) +
  facet_wrap(~Selection)+
  geom_line() +
  geom_line(aes(y = `F`, x = generation),
            colour = "red", alpha = 0.8) +
  geom_line(aes(y = `FM:Selected`, x = generation),
            colour = "#4b7000ff", size =0.5, alpha =0.8) +
  geom_line(aes(y = `FM:Non-Selected`, x = generation),
            colour = "#823175", size =0.5, alpha =0.6) +
  geom_line(aes(y = `M:Non-SelectedM:Selected`, x = generation),
            colour = "#c4852b", size =0.5, alpha =0.6) +
  geom_line(aes(y = `M:Selected`, x = generation),
            colour = "blue", alpha = 0.8) +
  geom_line(aes(y = `M:Non-Selected`, x = generation),
            colour = "#9597a1", size =0.5, alpha =0.8) +
  geom_vline(xintercept = 0, linetype = 2, alpha = 0.3) +
  ylab("Genetic Variance") +
  xlim(-20,32)+
  xlab("Generation") +
  theme_bw(base_size = 18)
g5 <- g4 + theme(legend.position = "top") +
  guides(
    color = guide_legend(override.aes = list(linetype = c(0, 1),
                                             shape = c(16, NA),
                                             color = "black")))

textplot <- tibble(
  Sum = rep(c(mean(Summ[Summ$generation==20 & Summ$Selection=="Medium accuracy",][,3]),  #Sum
              mean(Summ[Summ$generation==20 & Summ$Selection=="High accuracy",][,3]),
              mean(Summ[Summ$generation==20 & Summ$Selection=="Medium accuracy",][,5]), # M(MS)
              mean(Summ[Summ$generation==20 & Summ$Selection=="High accuracy",][,5]+0.05),
              mean(Summ[Summ$generation==20 & Summ$Selection=="Medium accuracy",][,6]-0.03), # M(S)
              mean(Summ[Summ$generation==20 & Summ$Selection=="High accuracy",][,6]-0.01),
              mean(Summ[Summ$generation==20 & Summ$Selection=="Medium accuracy",][,4]+0.02), # F
              mean(Summ[Summ$generation==20 & Summ$Selection=="High accuracy",][,4]+0.02),
              mean(Summ[Summ$generation==20 & Summ$Selection=="Medium accuracy",][,7]+0.07), # F:M(N)
              mean(Summ[Summ$generation==20 & Summ$Selection=="High accuracy",][,7]-0.07),
              mean(Summ[Summ$generation==20 & Summ$Selection=="Medium accuracy",][,8]-0.04), # F:M(S)
              mean(Summ[Summ$generation==20 & Summ$Selection=="High accuracy",][,8]-0.05),
              mean(Summ[Summ$generation==20 & Summ$Selection=="Medium accuracy",][,9]), #M(N):M(S)
              mean(Summ[Summ$generation==20 & Summ$Selection=="High accuracy",][,9])),2),
  Selection = as.factor(rep(rep(c("Medium accuracy","High accuracy"),7),2)),
  Scenario = c(rep("IS", 7*2), rep("WFS", 7*2)),
  generation = 27,
  label = rep(c(rep("Sum",2),     rep("M(N)",2),  rep("M(S)",2), rep("F",2), 
                rep("F:M(N)",2), rep("F:M(S)",2), rep("M(S):M(N)",2)),2),
  colour = rep(c(rep("black",2), rep("#9597a1",2), rep("blue",2), 
                 rep("red",2),   rep("#823175",2), rep("#4b7000ff",2), 
                 rep("#c4852b",2)),2),
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
ggsave("./Analysis/Figures/AlphaPartSexSel.pdf",width = 8, height = 8)
