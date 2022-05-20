########################################################################
## Project: jobsteter_alphapart_variance
## Script purpose: Dairy Cattle BP Analysis
## Date: 2021-11-30
## Author: Thiago de Paula Oliveira
########################################################################
# Set working directory
rm(list=ls())

# SET OPTIONS ---------------------------------------
cat("SETTING OPTIONS... \n\n", sep = "")
options(encoding = "UTF-8") # sets string encoding to UTF-8

#=======================================================================

#=======================================================================
# Packges
#=======================================================================
#devtools::install_github("AlphaGenes/AlphaPart")
library(tidyverse)
library(tidyr)

#=======================================================================
# Reading and organizing Scenario 1
#=======================================================================
dataS1 <- readRDS("./Results/Scenario1Var_replicate1.rds") %>%
  dplyr::mutate(generation = generation - 20) %>%
  dplyr::mutate(Selection = "Phenotypic") %>%
  dplyr::mutate(VarA_scaled = VarA/(0.3)) %>%
  dplyr::mutate(across(c("sex", "type"), as.factor)) %>%
  droplevels()

dataS2 <- readRDS("./Results/Scenario2Var_replicate1.rds") %>%
  dplyr::mutate(generation = generation - 20) %>%
  dplyr::mutate(Selection = "TBV") %>%
  dplyr::mutate(VarA_scaled = VarA/(0.3)) %>%
  dplyr::mutate(across(c("sex", "type"), as.factor)) %>%
  droplevels()

data <- rbind(dataS1, dataS2)
levels(data$type) <- c("Selection Canditates", "Used as Parents")
data <- unique(data)

data %>%
  ggplot(aes(y = MeanGeno, x = generation)) +
  facet_grid(type~Selection) +
  geom_line(aes(colour = sex)) +
  geom_vline(xintercept = 0, linetype=2, alpha=0.3)+
  theme_bw(base_size = 20)+
  ylab("Genetic Mean") +
  xlab("Generation") +
  labs(linetype = "Progeny:", colour = "Sex:") +
  theme(legend.position = "top")
ggsave("./Analysis/Figures/MeanA.pdf", width = 9, height = 7)


dataSum <- data %>%
  group_by(Selection, generation, type) %>%
  summarise(
    SumM = mean(MeanGeno),
    SumA = mean(VarA_scaled),
    SumAG = mean(GenicVA)
  )
dataSum[dataSum$type == "Used as Parents", ]$SumA <- NA

data %>%
  ggplot(aes(y = VarA_scaled, x = generation)) +
  facet_grid(type~Selection) +
  geom_line(data = dataSum, aes(y = SumA, x = generation),
            size = 0.4, alpha = 0.6) +
  geom_line(aes(colour = sex), alpha = 0.8, size =0.5) +
  scale_colour_manual(values = c("red","blue"))+
  scale_linetype_manual(values = c(3,1,2)) +
  ylab("Genetic Variance") +
  xlab("Generation") +
  labs(linetype = "Progeny:", colour = "Sex:") +
  geom_vline(xintercept = 0, linetype=2, alpha=0.3)+
  theme_bw(base_size = 20) +
  theme(legend.position = "top")
ggsave("./Analysis/Figures/VarA.pdf", width = 9, height = 7)
