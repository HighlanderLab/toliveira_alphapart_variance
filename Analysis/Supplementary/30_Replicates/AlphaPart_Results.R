########################################################################
## Project: toliveira_alphapart_variance
## Script purpose: Summarise the results 
## Date: 2022-03-21
## Author: Thiago de Paula Oliveira
########################################################################
#=======================================================================
# Packages
#=======================================================================
library(tidyverse)


########################################################################
# Phenotype Selection
########################################################################
#=======================================================================
# Genetic Mean
#=======================================================================
data <-  NULL
for(i in 1:30){
  cat("rep", i, "\n")
  load(paste0(
    "./Analysis/Supplementary/30_Replicates/Results/SMean_TBV_ValidationPheno", 
    i,".rds"))
  tmp <- inner_join(SMean[[1]], SMean[[2]], 
                     by = c("generation" = "generation", "path" = "path"))
  tmp <- tmp %>%
    dplyr::mutate(diff = Median.x - Median.y) %>%
    select(generation, path, diff) %>%
    dplyr::mutate(rep = i)
  data <- rbind(data, tmp)
}


data %>%
  group_by(generation, path) %>%
  summarise(
    diffMean = mean(diff),
    diffLower = diffMean - qnorm(0.99) * sd(diff),
    diffUpper = diffMean + qnorm(0.99) * sd(diff)
  ) %>%
  ggplot(aes(y = diffMean, x = generation)) +
  facet_wrap(~path) +
  geom_ribbon(aes(ymin = diffLower, ymax = diffUpper), alpha = 0.2) +
  geom_line() +
  ylab("Difference between estimated and true partition for genetic mean") +
  xlab("Generation")+
  ylim(-3.5, 1.5) +
  ggtitle("Medium accuracy") +
  geom_hline(yintercept = 0, linetype = 2, alpha = 0.4, colour = "blue") +
  geom_vline(xintercept = 0, linetype = 2, alpha = 0.2) +
  theme_bw(base_size = 14)
ggsave("./Analysis/Supplementary/30_Replicates/Figures/diffMeanPheno.pdf", 
       width = 7, height = 7)

#=======================================================================
# Mendelian sampling term
#=======================================================================
data <-  NULL
for(i in 1:30){
  cat("rep", i, "\n")
  load(paste0(
    "./Analysis/Supplementary/30_Replicates/Results/SMean_TBV_ValidationPheno", 
    i,".rds"))
  tmp <- inner_join(SMean[[4]], SMean[[3]], 
                    by = c("generation" = "generation", "path" = "path"))
  tmp <- tmp %>%
    dplyr::mutate(diff = Median - mstTrue) %>%
    select(generation, path, diff) %>%
    dplyr::mutate(rep = i)
  data <- rbind(data, tmp)
}


data %>%
  group_by(generation, path) %>%
  summarise(
    diffMean = mean(diff),
    diffLower = diffMean - qnorm(0.99) * sd(diff),
    diffUpper = diffMean + qnorm(0.99) * sd(diff)
  ) %>%
  ggplot(aes(y = diffMean, x = generation)) +
  facet_wrap(~path) +
  geom_ribbon(aes(ymin = diffLower, ymax = diffUpper), alpha = 0.2) +
  geom_line() +
  ylim(-1.1, 0.6) +
  ggtitle("Medium accuracy") +
  ylab("Difference between estimated and true partition for MST") +
  xlab("Generation")+
  geom_hline(yintercept = 0, linetype = 2, alpha = 0.4, colour = "blue") +
  geom_vline(xintercept = 0, linetype = 2, alpha = 0.2) +
  theme_bw(base_size = 14)
ggsave("./Analysis/Supplementary/30_Replicates/Figures/diffMeanMSTPheno.pdf", 
       width = 7, height = 7)

#=======================================================================
# Genetic Variance
#=======================================================================
data <-  NULL
for(i in 1:30){
  cat("rep", i, "\n")
  load(paste0(
    "./Analysis/Supplementary/30_Replicates/Results/SVar_TBV_ValidationPheno", 
    i,".rds"))
  tmp <- inner_join(SVar[[1]], SVar[[2]], 
                    by = c("generation" = "generation", "path" = "path"))
  tmp <- tmp %>%
    dplyr::mutate(diff = Median.x - Median.y) %>%
    select(generation, path, diff) %>%
    dplyr::mutate(rep = i)
  data <- rbind(data, tmp)
}


data %>%
  group_by(generation, path) %>%
  summarise(
    diffMean = mean(diff),
    diffLower = diffMean - qnorm(0.99) * sd(diff),
    diffUpper = diffMean + qnorm(0.99) * sd(diff)
  ) %>%
  ggplot(aes(y = diffMean, x = generation)) +
  facet_wrap(~path) +
  geom_ribbon(aes(ymin = diffLower, ymax = diffUpper), alpha = 0.2,
              colour = NA) +
  geom_line() +
  ylab("Difference between estimated and true partition for genetic variance") +
  xlab("Generation")+
  ylim(-0.25, 0.25) +
  ggtitle("Medium accuracy") +
  geom_vline(xintercept = 0, linetype = 2, alpha = 0.2) +
  geom_hline(yintercept = 0, linetype = 2, alpha = 0.4, colour = "blue") +
  theme_bw(base_size = 14)
ggsave("./Analysis/Supplementary/30_Replicates/Figures/diffVarPheno.pdf", 
       width = 7, height = 7)

########################################################################
# TBV Selection
########################################################################
#=======================================================================
# Genetic Mean
#=======================================================================
data <-  NULL
for(i in 1:30){
  cat("rep", i, "\n")
  load(paste0(
    "./Analysis/Supplementary/30_Replicates/Results/SMean_TBV_ValidationTBV", 
    i,".rds"))
  tmp <- inner_join(SMean[[1]], SMean[[2]], 
                    by = c("generation" = "generation", "path" = "path"))
  tmp <- tmp %>%
    dplyr::mutate(diff = Median.x - Median.y) %>%
    select(generation, path, diff) %>%
    dplyr::mutate(rep = i)
  data <- rbind(data, tmp)
}


data %>%
  group_by(generation, path) %>%
  summarise(
    diffMean = mean(diff),
    diffLower = diffMean - qnorm(0.99) * sd(diff),
    diffUpper = diffMean + qnorm(0.99) * sd(diff)
  ) %>%
  ggplot(aes(y = diffMean, x = generation)) +
  facet_wrap(~path) +
  geom_ribbon(aes(ymin = diffLower, ymax = diffUpper), alpha = 0.2) +
  geom_line() +
  ylab("Difference between estimated and true partition for genetic mean") +
  xlab("Generation")+
  geom_hline(yintercept = 0, linetype = 2, alpha = 0.4, colour = "blue") +
  ggtitle("High accuracy") +
  geom_vline(xintercept = 0, linetype = 2, alpha = 0.2) +  
  ylim(-3.5, 1.5) +
  theme_bw(base_size = 14)
ggsave("./Analysis/Supplementary/30_Replicates/Figures/diffMeanTBV.pdf", 
       width = 7, height = 7)

#=======================================================================
# Mendelian sampling term
#=======================================================================
data <-  NULL
for(i in 1:30){
  cat("rep", i, "\n")
  load(paste0(
    "./Analysis/Supplementary/30_Replicates/Results/SMean_TBV_ValidationTBV", 
    i,".rds"))
  tmp <- inner_join(SMean[[4]], SMean[[3]], 
                    by = c("generation" = "generation", "path" = "path"))
  tmp <- tmp %>%
    dplyr::mutate(diff = Median - mstTrue) %>%
    select(generation, path, diff) %>%
    dplyr::mutate(rep = i)
  data <- rbind(data, tmp)
}


data %>%
  group_by(generation, path) %>%
  summarise(
    diffMean = mean(diff),
    diffLower = diffMean - qnorm(0.99) * sd(diff),
    diffUpper = diffMean + qnorm(0.99) * sd(diff)
  ) %>%
  ggplot(aes(y = diffMean, x = generation)) +
  facet_wrap(~path) +
  geom_ribbon(aes(ymin = diffLower, ymax = diffUpper), alpha = 0.2) +
  geom_line() +
  ylab("Difference between estimated and true partition for MST") +
  xlab("Generation")+
  geom_hline(yintercept = 0, linetype = 2, alpha = 0.4, colour = "blue") +
  ggtitle("High accuracy") +
  geom_vline(xintercept = 0, linetype = 2, alpha = 0.2) +
  ylim(-1.1, 0.6) +
  theme_bw(base_size = 14)
ggsave("./Analysis/Supplementary/30_Replicates/Figures/diffMeanMSTTBV.pdf", 
       width = 7, height = 7)

#=======================================================================
# Genetic Variance
#=======================================================================
data <-  NULL
for(i in 1:30){
  cat("rep", i, "\n")
  load(paste0(
    "./Analysis/Supplementary/30_Replicates/Results/SVar_TBV_ValidationTBV", 
    i,".rds"))
  tmp <- inner_join(SVar[[1]], SVar[[2]], 
                    by = c("generation" = "generation", "path" = "path"))
  tmp <- tmp %>%
    dplyr::mutate(diff = Median.x - Median.y) %>%
    select(generation, path, diff) %>%
    dplyr::mutate(rep = i)
  data <- rbind(data, tmp)
}


data %>%
  group_by(generation, path) %>%
  summarise(
    diffMean = mean(diff),
    diffLower = diffMean - qnorm(0.99) * sd(diff),
    diffUpper = diffMean + qnorm(0.99) * sd(diff)
  ) %>%
  ggplot(aes(y = diffMean, x = generation)) +
  facet_wrap(~path) +
  geom_ribbon(aes(ymin = diffLower, ymax = diffUpper), alpha = 0.2,
              colour = NA) +
  
  geom_line() +
  ylab("Difference between estimated and true partition for genetic variance") +
  ggtitle("High accuracy") +
  xlab("Generation")+
  geom_hline(yintercept = 0, linetype = 2, alpha = 0.4, colour = "blue") +
  geom_vline(xintercept = 0, linetype = 2, alpha = 0.2) +
  ylim(-0.25, 0.25) +
  theme_bw(base_size = 14)
ggsave("./Analysis/Supplementary/30_Replicates/Figures/diffVarTBV.pdf", 
       width = 7, height = 7)
