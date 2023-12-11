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

# Wilcoxon signed-rank test -- non-parametric alternative to the one sample t-test
wilcox_test_01 <- wilcox.test(c01_data$ratio, mu = 1, exact = FALSE, conf.int = TRUE)
wilcox_test_01

wilcox_test_2 <- wilcox.test(c2_data$ratio, mu = 1, exact = FALSE, conf.int = TRUE)
wilcox_test_2

wilcox_test_3 <- wilcox.test(c3_data$ratio, mu = 1, exact = FALSE, conf.int = TRUE)
wilcox_test_3

### Graphs ###
mean_value_01 <- 0.9521435
ci_low_01 <- 0.8512433
ci_high_01 <- 1.0051647

graph_0.1 <- ggplot(c01_data, aes(y = ratio)) +
  theme_test(base_size = 15) + 
  ggtitle("(A) 0.1%")+ 
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(legend.position = "bottom",
        axis.text.x = element_text(size = 16),  # Adjust x-axis label size
        legend.text = element_text(size = 13))  + # Adjust legend text size 
  geom_violin(aes(x = 0, y = ratio), fill = "light gray", color = "white", alpha = 0.2) + 
  geom_jitter(aes(x = 0, color = compound, shape = compound), size = 3, alpha = 0.8, width = 0.05) +
  annotate("segment", x = -Inf, xend = Inf, y = 1, yend = 1, color = "black", linetype = "dashed") + 
  xlab (" ") + 
  ylab ("Ratio of consumption") + 
  geom_text(data = NULL, aes(x = 0, y = 2, label = "P-value = 0.075"), size = 5) +
  scale_x_continuous(breaks = NULL) +  # Remove x-axis tick marks
  scale_shape_manual(values = c("C1" = "circle", "C2" = "square", "C3" = "triangle", "C4" = "diamond"),
                     name = " ", labels = c("Piperine", "Tannic acid", "Eugenol", "Phytol")) + 
  scale_color_manual(values = c("C1" = "#3B528BFF", "C2" = "#21908CFF", "C3" = "#5DC863FF", "C4" = "#FDE725FF"),
                      name = " ", labels = c("Piperine", "Tannic acid", "Eugenol", "Phytol")) +
  scale_y_continuous(breaks = c(0, 1, 2),
                     limits = c(0, 2)) +  # Set the limits of the Y-axis
  # Add mean point with confidence intervals
  geom_point(aes(x = 0, y = mean_value_01), color = "black", size = 4) +
  geom_errorbar(aes(x = 0, ymin = ci_low_01, ymax = ci_high_01), 
                color = "black", 
                width = 0.05, 
                position = position_nudge(x = 0))

graph_0.1


#
mean_value_2 <- 0.6445069 
ci_low_2 <- 0.4963644
ci_high_2 <- 0.8331745

graph_2 <- ggplot(c2_data, aes(y = ratio)) +
  theme_test(base_size = 15) + 
  ggtitle("(B) 2%")+ 
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(legend.position = "bottom",
        axis.text.x = element_text(size = 16),  # Adjust x-axis label size
        legend.text = element_text(size = 13))  + # Adjust legend text size 
  geom_violin(aes(x = 0, y = ratio), fill = "light gray", color = "white", alpha = 0.2) + 
  geom_jitter(aes(x = 0, color = compound, shape = compound), size = 3, alpha = 0.8, width = 0.05) +
  annotate("segment", x = -Inf, xend = Inf, y = 1, yend = 1, color = "black", linetype = "dashed") + 
  xlab (" ") + 
  ylab (" ") + 
  geom_text(data = NULL, aes(x = 0, y = 2, label = "P-value < 0.001"), size = 5) +
  scale_x_continuous(breaks = NULL) +  # Remove x-axis tick marks
  scale_shape_manual(values = c("C1" = "circle", "C2" = "square", "C3" = "triangle", "C4" = "diamond"),
                     name = " ", labels = c("Piperine", "Tannic acid", "Eugenol", "Phytol")) + 
  scale_color_manual(values = c("C1" = "#3B528BFF", "C2" = "#21908CFF", "C3" = "#5DC863FF", "C4" = "#FDE725FF"),
                     name = " ", labels = c("Piperine", "Tannic acid", "Eugenol", "Phytol")) +
  scale_y_continuous(breaks = c(0, 1, 2),
                     limits = c(0, 2)) +  # Set the limits of the Y-axis
  # Add mean point with confidence intervals
  geom_point(aes(x = 0, y = mean_value_2), color = "black", size = 4) +
  geom_errorbar(aes(x = 0, ymin = ci_low_2, ymax = ci_high_2), 
                color = "black", 
                width = 0.05, 
                position = position_nudge(x = 0))

graph_2
#
mean_value_3 <- 0.7707695 
ci_low_3 <- 0.5626209
ci_high_3 <- 0.9602930

graph_3 <- ggplot(c3_data, aes(y = ratio)) +
  theme_test(base_size = 15) + 
  ggtitle("(C) 3%")+ 
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(legend.position = "bottom",
        axis.text.x = element_text(size = 16),  # Adjust x-axis label size
        legend.text = element_text(size = 13))  + # Adjust legend text size 
  geom_violin(aes(x = 0, y = ratio), fill = "light gray", color = "white", alpha = 0.2) + 
  geom_jitter(aes(x = 0, color = compound, shape = compound), size = 3, alpha = 0.8, width = 0.05) +
  annotate("segment", x = -Inf, xend = Inf, y = 1, yend = 1, color = "black", linetype = "dashed") + 
  xlab (" ") + 
  ylab (" ") + 
  geom_text(data = NULL, aes(x = 0, y = 2, label = "P-value = 0.021"), size = 5) +
  scale_x_continuous(breaks = NULL) +  # Remove x-axis tick marks
  scale_shape_manual(values = c("C1" = "circle", "C2" = "square", "C3" = "triangle", "C4" = "diamond"),
                     name = " ", labels = c("Piperine", "Tannic acid", "Eugenol", "Phytol")) + 
  scale_color_manual(values = c("C1" = "#3B528BFF", "C2" = "#21908CFF", "C3" = "#5DC863FF", "C4" = "#FDE725FF"),
                     name = " ", labels = c("Piperine", "Tannic acid", "Eugenol", "Phytol")) +
  scale_y_continuous(breaks = c(0, 1, 2),
                     limits = c(0, 2)) +  # Set the limits of the Y-axis
  # Add mean point with confidence intervals
  geom_point(aes(x = 0, y = mean_value_3), color = "black", size = 4) +
  geom_errorbar(aes(x = 0, ymin = ci_low_3, ymax = ci_high_3), 
                color = "black", 
                width = 0.05, 
                position = position_nudge(x = 0))

graph_3

#
preference_obj1 <- ggarrange(graph_0.1,
                             graph_2,
                             graph_3,
                             ncol = 3, nrow = 1,
                             align = "hv", 
                             common.legend = TRUE,
                             legend = "bottom")
preference_obj1 
ggsave(file="preference_obj1 .jpg", 
       plot=preference_obj1,
       width=8,height=4,units="in",dpi=600)

####################################
####################################
####################################
####################################



