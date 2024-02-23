################################################################################
### Objective 1. The effects of metabolite identity on bat diet preferences. ###
################################################################################

rm(list=ls())

library(tidyverse)
library(dplyr)
library(viridis)
library(ggplot2)
library(MetBrewer)
library(glmmTMB) 
library(psych)
library(effects)
library(ggpubr)
library(emmeans)
library(parameters)
library(car)
library(performance)

### Get data ### 
setwd("/Users/marianagelambi/Desktop/bat-metabolomics/objective1")
data <- read.csv("data_objective1.csv")
head(data)

# Replace values in the 'food_consumed' column that are less than 0.0199 with 0. These values were not 'true consumption.'
data$food_consumed[data$food_consumed < 0.0199] <- 0
# Calculate the average food consumption for control treatment (C, concentration 0) for each bat
control_avg <- data %>%
  filter(compound == "C") %>%
  group_by(bat) %>%
  summarize(avg_food_consumed = mean(food_consumed))
# Merge the average control values back into the original dataframe 
data <- data %>%
  left_join(control_avg, by = "bat") %>%
  mutate(ratio= food_consumed/avg_food_consumed)
# Clean data from bats that did not eat or eat to little from the control,  
data <- data %>%
  filter(!(ratio > 25 | is.na(ratio) | compound == 'C')) # Keep rows where the ratio is not greater than 25, is not missing, and the compound is not equal to 'C'.
# Subset in different concentrations 
c01_data <- data %>%
  filter(concentration == "0.1") # 40 obs
c2_data <- data %>%
  filter(concentration == "2")  # 40 obs
c3_data <- data %>%
  filter(concentration == "3")  # 37 obs

hist(c01_data$ratio)
shapiro.test(c01_data$ratio)
hist(c2_data$ratio)
shapiro.test(c2_data$ratio)
hist(c3_data$ratio)
shapiro.test(c3_data$ratio)

### INDIVIDUAL EFFECT OF SECONDARY METABOLITES CONCENTRATION ON BAT PREFERENCE ###

# Subset in different concentrations 
c01_data <- data %>%
  filter(concentration == "0.1") # 40 obs
c2_data <- data %>%
  filter(concentration == "2")  # 40 obs
c3_data <- data %>%
  filter(concentration == "3")  # 37 obs

############
### 0.1% ###
############
piperine_0.1 <- c01_data %>%
  filter(compound == "C1") # 10 obs
ta_0.1 <- c01_data %>%
  filter(compound == "C2") # 10 obs
eugenol_0.1 <- c01_data %>%
  filter(compound == "C3") # 10 obs
phytol_0.1 <- c01_data %>%
  filter(compound == "C4") # 10 obs

### Wilcoxon signed rank test with continuity correction
wilcoxtest_piperine_0.1 <- wilcox.test(piperine_0.1$ratio, mu = 1, exact = FALSE, conf.int = TRUE)
wilcoxtest_piperine_0.1
ci_low_piperine_0.1 <- 0.8343988
ci_high_piperine_0.1 <- 1.2093825

wilcoxtest_ta_0.1 <- wilcox.test(ta_0.1$ratio, mu = 1, exact = FALSE, conf.int = TRUE)
wilcoxtest_ta_0.1
ci_low_ta_0.1 <- 0.788126 
ci_high_ta_0.1 <- 1.328504

wilcoxtest_eugenol_0.1 <- wilcox.test(eugenol_0.1$ratio, mu = 1, exact = FALSE, conf.int = TRUE)
wilcoxtest_eugenol_0.1
ci_low_eugenol_0.1 <- 0.4908292
ci_high_eugenol_0.1 <- 1.0397079

wilcoxtest_phytol_0.1 <- wilcox.test(phytol_0.1$ratio, mu = 1, exact = FALSE, conf.int = TRUE)
wilcoxtest_phytol_0.1
ci_low_phytol_0.1 <- 0.4422667
ci_high_phytol_0.1 <- 0.9746503

### GLMMs
glmm_0.1 <- glmmTMB(ratio ~ compound + (1|date) + (1|bat), data = c01_data)
summary(glmm_0.1)
Anova(glmm_0.1)
parameters(glmm_0.1)
plot(allEffects(glmm_0.1))
shapiro.test(resid(glmm_0.1)) # p-value = 0.9599
hist(resid(glmm_0.1)) # looks good
summary(allEffects(glmm_0.1))
r2(glmm_0.1)
glmm_0.1_emmeans <-emmeans(glmm_0.1,~compound, type="response")
pairwise_comp_0.1 <- pairs(glmm_0.1_emmeans)
pairwise_comp_0.1
write.csv(pairwise_comp_0.1 ,
          "pairwise_comp_0.1.csv")

### Graph 0.1% ~ compounds
graph_0.1 <- ggplot(c01_data, aes(x = compound, y = ratio)) +
  theme_test(base_size = 12) + 
  ggtitle("(A) 0.1%")+ 
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(legend.position = "null") + 
  geom_violin(fill = "light gray", color = "white", alpha = 0.6) + 
  geom_point(size = 4.5, aes(color = "black", shape = compound))  +
  geom_point(size = 3, alpha = 0.8, aes(color = compound, shape = compound))  +
  annotate("segment", x = -Inf, xend = Inf, y = 1, yend = 1, color = "black", linetype = "dashed") + 
  xlab (" ") + 
  ylab ("Ratio of consumption") + 
  scale_y_continuous(breaks = c(0, 1, 2),
                     limits = c(0, 2)) +
  scale_x_discrete(labels=c("C" = "Control", "C1" = "Piperine", "C2" = "Tannic\nacid", "C3" = "Eugenol", "C4" = "Phytol")) +
  scale_shape_manual(values = c("C1" = "circle", "C2" = "square", "C3" = "triangle", "C4" = "diamond"),
                     name = " ", labels = c("Piperine", "Tannic acid", "Eugenol", "Phytol")) + 
  scale_color_manual(values = c("C1" = "#3B528BFF", "C2" = "#21908CFF", "C3" = "#5DC863FF", "C4" = "#FDE725FF"),
                     name = " ", labels = c("Piperine", "Tannic acid", "Eugenol", "Phytol")) + 
  #piperine
  geom_point(aes(x = 1, y = wilcoxtest_piperine_0.1$estimate), color = "black", size = 4) + 
  geom_errorbar(aes(x = 1, ymin = ci_low_piperine_0.1, ymax = ci_high_piperine_0.1),
                color = "black", width = 0.2, position = position_dodge(width = 0.2)) + 
  geom_text(data = NULL, aes(x = 1, y = 2, label = "P = 0.476"), size = 2.5) + 
  #tannic acid
  geom_point(aes(x = 2, y = wilcoxtest_ta_0.1$estimate), color = "black", size = 4) + 
  geom_errorbar(aes(x = 2, ymin = ci_low_ta_0.1, ymax = ci_high_ta_0.1),
                color = "black", width = 0.2, position = position_dodge(width = 0.2)) + 
  geom_text(data = NULL, aes(x = 2, y = 2, label = "P = 0.760"), size = 2.5) + 
  #eugenol
  geom_point(aes(x = 3, y = wilcoxtest_eugenol_0.1$estimate), color = "black", size = 4) + 
  geom_errorbar(aes(x = 3, ymin = ci_low_eugenol_0.1, ymax = ci_high_eugenol_0.1),
                color = "black", width = 0.2, position = position_dodge(width = 0.2)) + 
  geom_text(data = NULL, aes(x = 3, y = 2, label = "P = 0.126"), size = 2.5) + 
  #phytol
  geom_point(aes(x = 4, y = wilcoxtest_phytol_0.1$estimate), color = "black", size = 4) + 
  geom_errorbar(aes(x = 4, ymin = ci_low_phytol_0.1, ymax = ci_high_phytol_0.1),
                color = "black", width = 0.2, position = position_dodge(width = 0.2)) + 
  geom_text(x = 4, y = 2, label = "P = 0.019", size = 2.5)
    

graph_0.1


############
### 2% ###
############
piperine_2 <- c2_data %>%
  filter(compound == "C1") # 10 obs
ta_2 <- c2_data %>%
  filter(compound == "C2") # 10 obs
eugenol_2 <- c2_data %>%
  filter(compound == "C3") # 10 obs
phytol_2 <- c2_data %>%
  filter(compound == "C4") # 10 obs

### Wilcoxon signed rank test with continuity correction
wilcoxtest_piperine_2 <- wilcox.test(piperine_2$ratio, mu = 1, exact = FALSE, conf.int = TRUE)
wilcoxtest_piperine_2
ci_low_piperine_2 <- 0.4538228
ci_high_piperine_2 <- 1.0466989

wilcoxtest_ta_2 <- wilcox.test(ta_2$ratio, mu = 1, exact = FALSE, conf.int = TRUE)
wilcoxtest_ta_2
ci_low_ta_2 <- 0.4543688
ci_high_ta_2 <- 1.0431951

wilcoxtest_eugenol_2 <- wilcox.test(eugenol_2$ratio, mu = 1, exact = FALSE, conf.int = TRUE)
wilcoxtest_eugenol_2
ci_low_eugenol_2 <- 0.1148112
ci_high_eugenol_2 <- 0.6183521

wilcoxtest_phytol_2 <- wilcox.test(phytol_2$ratio, mu = 1, exact = FALSE, conf.int = TRUE)
wilcoxtest_phytol_2
ci_low_phytol_2 <- 0.3336493
ci_high_phytol_2 <- 0.9982877

### GLMMs
glmm_2 <- glmmTMB(ratio ~ compound + (1|date) + (1|bat), data = c2_data)
summary(glmm_2)
Anova(glmm_2)
parameters(glmm_2)
plot(allEffects(glmm_2))
shapiro.test(resid(glmm_2)) # p-value = 0.8079
hist(resid(glmm_2)) # looks good
summary(allEffects(glmm_2))
r2(glmm_2)
glmm_2_emmeans <-emmeans(glmm_2,~compound, type="response")
pairwise_comp_2 <- pairs(glmm_2_emmeans)

write.csv(pairwise_comp_2 ,
          "pairwise_comp_2.csv")
### Graph 2% ~ compounds
graph_2 <- ggplot(c2_data, aes(x = compound, y = ratio)) +
  theme_test(base_size = 12) + 
  ggtitle("(B) 2%")+ 
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(legend.position = "null") + 
  geom_violin(fill = "light gray", color = "white", alpha = 0.6) + 
  geom_point(size = 4.5, aes(color = "black", shape = compound))  +
  geom_point(size = 3, alpha = 0.8, aes(color = compound, shape = compound))  +
  annotate("segment", x = -Inf, xend = Inf, y = 1, yend = 1, color = "black", linetype = "dashed") + 
  xlab (" ") + 
  ylab (" ") + 
  scale_y_continuous(breaks = c(0, 1, 2),
                     limits = c(0, 2)) +
  scale_x_discrete(labels=c("C" = "Control", "C1" = "Piperine", "C2" = "Tannic\nacid", "C3" = "Eugenol", "C4" = "Phytol")) +
  scale_shape_manual(values = c("C1" = "circle", "C2" = "square", "C3" = "triangle", "C4" = "diamond"),
                     name = " ", labels = c("Piperine", "Tannic acid", "Eugenol", "Phytol")) + 
  scale_color_manual(values = c("C1" = "#3B528BFF", "C2" = "#21908CFF", "C3" = "#5DC863FF", "C4" = "#FDE725FF"),
                     name = " ", labels = c("Piperine", "Tannic acid", "Eugenol", "Phytol")) + 
  #piperine
  geom_point(aes(x = 1, y = wilcoxtest_piperine_2$estimate), color = "black", size = 4) + 
  geom_errorbar(aes(x = 1, ymin = ci_low_piperine_2, ymax = ci_high_piperine_2),
                color = "black", width = 0.2, position = position_dodge(width = 0.2)) + 
  geom_text(data = NULL, aes(x = 1, y = 2, label = "P = 0.185"), size = 2.5) + 
  #tannic acid
  geom_point(aes(x = 2, y = wilcoxtest_ta_2$estimate), color = "black", size = 4) + 
  geom_errorbar(aes(x = 2, ymin = ci_low_ta_2, ymax = ci_high_ta_2),
                color = "black", width = 0.2, position = position_dodge(width = 0.2)) + 
  geom_text(data = NULL, aes(x = 2, y = 2, label = "P = 0.185"), size = 2.5) + 
  #eugenol
  geom_point(aes(x = 3, y = wilcoxtest_eugenol_2$estimate), color = "black", size = 4) + 
  geom_errorbar(aes(x = 3, ymin = ci_low_eugenol_2, ymax = ci_high_eugenol_2),
                color = "black", width = 0.2, position = position_dodge(width = 0.2)) + 
  geom_text(data = NULL, aes(x = 3, y = 2, label = "P = 0.006"), size = 2.5) + 
  #phytol
  geom_point(aes(x = 4, y = wilcoxtest_phytol_2$estimate), color = "black", size = 4) + 
  geom_errorbar(aes(x = 4, ymin = ci_low_phytol_2, ymax = ci_high_phytol_2),
                color = "black", width = 0.2, position = position_dodge(width = 0.2)) + 
  geom_text(data = NULL, aes(x = 4, y = 2, label = "P = 0.041"), size = 2.5)


graph_2

###########
### 3% ###
###########

piperine_3 <- c3_data %>%
  filter(compound == "C1") # 9 obs
ta_3 <- c3_data %>%
  filter(compound == "C2") # 9 obs
eugenol_3 <- c3_data %>%
  filter(compound == "C3") # 10 obs
phytol_3 <- c3_data %>%
  filter(compound == "C4") # 9 obs

### Wilcoxon signed rank test with continuity correction
wilcoxtest_piperine_3 <- wilcox.test(piperine_3$ratio, mu = 1, exact = FALSE, conf.int = TRUE)
wilcoxtest_piperine_3
ci_low_piperine_3 <- 0.7637845 
ci_high_piperine_3 <- 1.4615681

wilcoxtest_ta_3 <- wilcox.test(ta_3$ratio, mu = 1, exact = FALSE, conf.int = TRUE)
wilcoxtest_ta_3
ci_low_ta_3 <- 0.7296888 
ci_high_ta_3 <- 1.4483672

wilcoxtest_eugenol_3 <- wilcox.test(eugenol_3$ratio, mu = 1, exact = FALSE, conf.int = TRUE)
wilcoxtest_eugenol_3
ci_low_eugenol_3 <- 0.02969686 
ci_high_eugenol_3 <- 0.90703955

wilcoxtest_phytol_3 <- wilcox.test(phytol_3$ratio, mu = 1, exact = FALSE, conf.int = TRUE)
wilcoxtest_phytol_3
ci_low_phytol_3 <-  0.01293217 
ci_high_phytol_3 <- 0.93842839

### GLMMs
glmm_3 <- glmmTMB(ratio ~ compound + (1|date) + (1|bat), data = c3_data)
summary(glmm_3)
Anova(glmm_3)
parameters(glmm_3)
plot(allEffects(glmm_3))
shapiro.test(resid(glmm_3)) # p-value = 0.9599
hist(resid(glmm_3)) # looks good
summary(allEffects(glmm_3))
r2(glmm_3)
glmm_3_emmeans <-emmeans(glmm_3,~compound, type="response")
pairwise_comp_3 <- pairs(glmm_3_emmeans)
pairwise_comp_3
write.csv(pairwise_comp_3 ,
          "pairwise_comp_3.csv")

### Graph 3% ~ compounds
graph_3 <- ggplot(c3_data, aes(x = compound, y = ratio)) +
  theme_test(base_size = 12) + 
  ggtitle("(C) 3%")+ 
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(legend.position = "null") + 
  geom_violin(fill = "light gray", color = "white", alpha = 0.6) + 
  geom_point(size = 4.5, aes(color = "black", shape = compound))  +
  geom_point(size = 3, alpha = 0.8, aes(color = compound, shape = compound))  +
  annotate("segment", x = -Inf, xend = Inf, y = 1, yend = 1, color = "black", linetype = "dashed") + 
  xlab (" ") + 
  ylab (" ") + 
  scale_y_continuous(breaks = c(0, 1, 2),
                     limits = c(0, 2)) +
  scale_x_discrete(labels=c("C" = "Control", "C1" = "Piperine", "C2" = "Tannic\nacid", "C3" = "Eugenol", "C4" = "Phytol")) +
  scale_shape_manual(values = c("C1" = "circle", "C2" = "square", "C3" = "triangle", "C4" = "diamond"),
                     name = " ", labels = c("Piperine", "Tannic acid", "Eugenol", "Phytol")) + 
  scale_color_manual(values = c("C1" = "#3B528BFF", "C2" = "#21908CFF", "C3" = "#5DC863FF", "C4" = "#FDE725FF"),
                     name = " ", labels = c("Piperine", "Tannic acid", "Eugenol", "Phytol")) + 
  #piperine
  geom_point(aes(x = 1, y = wilcoxtest_piperine_3$estimate), color = "black", size = 4) + 
  geom_errorbar(aes(x = 1, ymin = ci_low_piperine_3, ymax = ci_high_piperine_3),
                color = "black", width = 0.2, position = position_dodge(width = 0.2)) + 
  geom_text(data = NULL, aes(x = 1, y = 2, label = "P = 0.554"), size = 2.5) + 
  #tannic acid
  geom_point(aes(x = 2, y = wilcoxtest_ta_3$estimate), color = "black", size = 4) + 
  geom_errorbar(aes(x = 2, ymin = ci_low_ta_3, ymax = ci_high_ta_3),
                color = "black", width = 0.2, position = position_dodge(width = 0.2)) + 
  geom_text(data = NULL, aes(x = 2, y = 2, label = "P = 0.636"), size = 2.5) + 
  #eugenol
  geom_point(aes(x = 3, y = wilcoxtest_eugenol_3$estimate), color = "black", size = 4) + 
  geom_errorbar(aes(x = 3, ymin = ci_low_eugenol_3, ymax = ci_high_eugenol_3),
                color = "black", width = 0.2, position = position_dodge(width = 0.2)) + 
  geom_text(data = NULL, aes(x = 3, y = 2, label = "P = 0.019"), size = 2.5) + 
  #phytol
  geom_point(aes(x = 4, y = wilcoxtest_phytol_3$estimate), color = "black", size = 4) + 
  geom_errorbar(aes(x = 4, ymin = ci_low_phytol_3, ymax = ci_high_phytol_3),
                color = "black", width = 0.2, position = position_dodge(width = 0.2)) + 
  geom_text(data = NULL, aes(x = 4, y = 2, label = "P = 0.033"), size = 2.5)

graph_3

preference_graph <- ggarrange(graph_0.1,
                              graph_2,
                              graph_3,
                              ncol = 3, nrow = 1)

preference_graph
ggsave(file="preference_graph.jpg", 
       plot= preference_graph,
       width=8,height=4,units="in",dpi=600)

####################################
####################################
####################################
####################################
