#################################################################################################
### Objective 2. The effects of secondary metabolites consumption on the bat fecal metabolome ###
#################################################################################################

rm(list=ls())

library(tidyverse)
library(dplyr)
library(vegan)
library(viridis)
library(ggplot2)
library(patchwork)
library(MetBrewer)
library(glmmTMB) 
library(randomForest)
library(varSelRF)
library(VSURF)
library(Boruta)
library(psych)
library(effects)
library(AER)
library(ggpubr)
library(emmeans)
library(classyfireR)

### Get data ### 
setwd("/Users/marianagelambi/Desktop/bat-metabolomics/objective2")
data_final <- read.csv("data_final_metaMS_dec2_normalizedinsresponse.csv")
data_final$concentration <- as.factor(data_final$concentration)
data_final <- data_final %>%
  select(where(~ any(. != 0)))
head(data_final)

### Subset data in different concentrations (all compounds in 0.1, 2, and 3)
concentration_0.1 <- data_final %>% filter(batID %in% c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10"))
concentration_2 <- data_final %>% filter(batID %in% c("11", "12", "13", "14", "15", "16", "17", "18", "19", "20"))
concentration_3 <- data_final %>% filter(batID %in% c("21", "22", "23", "24", "25", "26", "27", "28", "29", "30"))

# eliminate variable with all zeros 
concentration_0.1_no0 <- concentration_0.1 %>%
  select(where(~ any(. != 0)))
concentration_2_no0 <- concentration_2 %>%
  select(where(~ any(. != 0)))
concentration_3_no0 <- concentration_3 %>%
  select(where(~ any(. != 0))) # eliminate variable with all zeros 

################################
### RANDOM FOREST AND BORUTA ###
################################

### 0.1%
x1 <- concentration_0.1_no0[, 2:211] #x1 should only have the compound data, no other factors
c1 <- as.factor(concentration_0.1_no0$treatment)
randfor <- randomForest(x1, c1, importance=TRUE, proximity=TRUE, oob.prox=TRUE, ntree=2000)
plot(randfor)
importance_allcompounds <- importance(randfor)
varImpPlot(randfor)
boruta <- Boruta(x1, c1, pValue=0.05, maxRuns=1000)
plot(boruta)
plotImpHistory(boruta)
boruta_table <- getImpRfZ(x1, c1)
boruta_table
boruta_sigfeatures <- getSelectedAttributes(boruta)
boruta_sigfeatures

### 2%
x1 <- concentration_2_no0[ , 2:151] #x1 should only have the compound data, no other factors
c1 <- as.factor(concentration_2_no0$treatment)
randfor <- randomForest(x1, c1, importance=TRUE, proximity=TRUE, oob.prox=TRUE, ntree=2000)
plot(randfor)
importance_allcompounds <- importance(randfor)
varImpPlot(randfor)
boruta <- Boruta(x1, c1, pValue=0.05, maxRuns=1000)
plot(boruta)
plotImpHistory(boruta)
boruta_table <- getImpRfZ(x1, c1)
boruta_table
boruta_sigfeatures <- getSelectedAttributes(boruta)
boruta_sigfeatures

### 3%
x1 <- concentration_3_no0[ , 2:200] #x1 should only have the compound data, no other factors
c1 <- as.factor(concentration_3_no0$treatment)
randfor <- randomForest(x1, c1, importance=TRUE, proximity=TRUE, oob.prox=TRUE, ntree=2000)
plot(randfor)
importance_allcompounds <- importance(randfor)
varImpPlot(randfor)
boruta <- Boruta(x1, c1, pValue=0.05, maxRuns=1000)
plot(boruta)
plotImpHistory(boruta)
boruta_table <- getImpRfZ(x1, c1)
boruta_table
boruta_sigfeatures <- getSelectedAttributes(boruta)
boruta_sigfeatures

#######################################
### GENERALIZED LINEAR MIXED MODELS ###
#######################################

### 0.1%

# 226
concentration_0.1_no0$Unknown.226_gamma <- concentration_0.1_no0$Unknown.226 + 0.0001
concentration_0.1_no0$Unknown.226_gamma
c226 <- glmmTMB(Unknown.226_gamma ~ treatment + (1|batID), data = concentration_0.1_no0, family=Gamma(link = log))
summary(c226)
plot(allEffects(c226))
c226_emmeans <-emmeans(c226,~ treatment, type="response")
c226_emmeans
c226_emmeans <- as.data.frame(c226_emmeans)
concentration_0.1_no0$predictions_c226 <- predict(c226, concentration_0.1_no0, re.form=NA, type="response")

concentration_226_graph <- ggplot(concentration_0.1_no0, aes(x = treatment, y = predictions_c226, color = treatment, shape = treatment)) +
  theme_test(base_size = 15) +
  scale_color_viridis(option = "D", discrete=TRUE, name = "Compounds") + 
  geom_point (data = concentration_0.1_no0, aes(x = treatment, y = Unknown.226_gamma, color = treatment), size = 2) +
  stat_summary(fun.data = mean_se, color = "black", size = 1) + 
  geom_errorbar(data = c226_emmeans, aes(x = treatment, y = response, ymin = asymp.LCL, ymax = asymp.UCL), width = 0.2, color = "black") +
  ylab ("Eugenol") +
  xlab (" ") +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) + 
  scale_x_discrete(labels=c("C" = "Control", "C1" = "Piperine", "C2" = "Tannic acid", "C3" = "Eugenol", "C4" = "Phytol")) + 
  geom_text(x = 4, y = 0.06, label = "***", size = 8, color = "black")
concentration_226_graph

# 6
concentration_0.1_no0$Unknown.6_gamma <- concentration_0.1_no0$Unknown.6 + 0.0001
c6 <- glmmTMB(Unknown.6_gamma ~ treatment + (1|batID), data = concentration_0.1_no0, family=Gamma(link = log))
summary(c6)
plot(allEffects(c6))
c6_emmeans <-emmeans(c6,~ treatment, type="response")
c6_emmeans
effect_size_c6 <- eff_size(c6_emmeans, sigma= sigma(c6), edf = df.residual(c6)) # get effect sizes 
effect_size_c6
c6_emmeans <- as.data.frame(c6_emmeans)
concentration_0.1_no0$predictions_c6 <- predict(c6, concentration_0.1_no0, re.form=NA, type="response")

concentration_6_graph <- ggplot(concentration_0.1_no0, aes(x = treatment, y = predictions_c6, color = treatment, shape = treatment)) +
  theme_test(base_size = 15) +
  scale_color_viridis(option = "D", discrete=TRUE, name = "Compounds") + 
  geom_point (data = concentration_0.1_no0, aes(x = treatment, y = Unknown.6_gamma, color = treatment), size = 2) +
  stat_summary(fun.data = mean_se, color = "black", size = 1) + 
  geom_errorbar(data = c6_emmeans, aes(x = treatment, y = response, ymin = asymp.LCL, ymax = asymp.UCL), width = 0.2, color = "black") +
  ylab ("A") +
  xlab (" ") +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) + 
  scale_x_discrete(labels=c("C" = "Control", "C1" = "Piperine", "C2" = "Tannic acid", "C3" = "Eugenol", "C4" = "Phytol")) + 
  geom_text(x = 3, y = 0.15, label = "*", size = 8, color = "black")

concentration_6_graph

allcompounds_01_RF <- concentration_6_graph + 
  concentration_226_graph

allcompounds_01_RF

ggsave(file="allcompounds_01_RF.jpg", 
       plot=allcompounds_01_RF,
       width=6,height=3,units="in",dpi=300)

### 2% 

# 353
concentration_2_no0$Unknown.353_gamma <- concentration_2_no0$Unknown.353 + 0.0001
c353 <- glmmTMB(Unknown.353_gamma ~ treatment + (1|batID), data = concentration_2_no0, family=Gamma(link = log))
summary(c353)
plot(allEffects(c353))
c353_emmeans <-emmeans(c353,~ treatment, type="response")
c353_emmeans
effect_size_c353 <- eff_size(c353_emmeans, sigma= sigma(c353), edf = df.residual(c353)) # get effect sizes 
effect_size_c353
c353_emmeans <- as.data.frame(c353_emmeans)
concentration_2_no0$predictions_c353 <- predict(c353, concentration_2_no0, re.form=NA, type="response")

concentration_353_graph <- ggplot(concentration_2_no0, aes(x = treatment, y = predictions_c353, color = treatment, shape = treatment)) +
  theme_test(base_size = 15) +
  scale_color_viridis(option = "D", discrete=TRUE, name = "Compounds") + 
  scale_shape_manual(values = c("C" = "cross", "C1" = "circle", "C2" = "square", "C3" = "triangle", "C4" = "diamond"),
                     name = " ", labels = c("Control", "Piperine", "Tannic acid", "Eugenol", "Phytol")) + 
  stat_summary(fun.data = mean_se, color = "black", size = 1) + 
  geom_point(data = concentration_2_no0, aes(x = treatment, y = Unknown.353_gamma, color = "black"), size = 3) +
  geom_point(data = concentration_2_no0, aes(x = treatment, y = Unknown.353_gamma, color = treatment), size = 2) +
  geom_errorbar(data = c353_emmeans, aes(x = treatment, y = response, ymin = asymp.LCL, ymax = asymp.UCL), width = 0.2, color = "black") +
  ylab ("Eugenol") +
  xlab (" ") +
  theme(legend.position = "none") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
  geom_text(x = 4, y = 0.40, label = "***", size = 8, color = "black")
concentration_353_graph

# 226
concentration_2_no0$Unknown.226_gamma <- concentration_2_no0$Unknown.226 + 0.0001
c226 <- glmmTMB(Unknown.226_gamma ~ treatment + (1|batID), data = concentration_2_no0, family=Gamma(link = log))
summary(c226)
plot(allEffects(c226))
c226_emmeans <-emmeans(c226,~ treatment, type="response")
c226_emmeans
effect_size_c226 <- eff_size(c226_emmeans, sigma= sigma(c226), edf = df.residual(c226)) # get effect sizes 
effect_size_c226
c226_emmeans <- as.data.frame(c226_emmeans)
concentration_2_no0$predictions_c226 <- predict(c226, concentration_2_no0, re.form=NA, type="response")

concentration_226_graph <- ggplot(concentration_2_no0, aes(x = treatment, y = predictions_c226, color = treatment, shape = treatment)) +
  theme_test(base_size = 15) +
  scale_color_viridis(option = "D", discrete=TRUE, name = "Compounds") + 
  geom_point(data = concentration_2_no0, aes(x = treatment, y = Unknown.226_gamma, color = "black"), size = 3) +
  geom_point(data = concentration_2_no0, aes(x = treatment, y = Unknown.226_gamma, color = treatment), size = 2) +
  scale_shape_manual(values = c("C" = "cross", "C1" = "circle", "C2" = "square", "C3" = "triangle", "C4" = "diamond"),
                     name = " ", labels = c("Control", "Piperine", "Tannic acid", "Eugenol", "Phytol")) + 
  stat_summary(fun.data = mean_se, color = "black", size = 1) + 
  geom_errorbar(data = c226_emmeans, aes(x = treatment, y = response, ymin = asymp.LCL, ymax = asymp.UCL), width = 0.2, color = "black") +
  ylab ("Eugenol") +
  xlab (" ") +
  theme(legend.position = "none") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
  geom_text(x = 4, y = 0.35, label = "***", size = 8, color = "black")
concentration_226_graph

# 150
concentration_2_no0$Unknown.150_gamma <- concentration_2_no0$Unknown.150 + 0.0001
c150 <- glmmTMB(Unknown.150_gamma ~ treatment + (1|batID), data = concentration_2_no0, family=Gamma(link = log))
summary(c150)
plot(allEffects(c150))
c150_emmeans <-emmeans(c150,~ treatment, type="response")
c150_emmeans
effect_size_c150 <- eff_size(c150_emmeans, sigma= sigma(c150), edf = df.residual(c150)) # get effect sizes 
effect_size_c150
c150_emmeans <- as.data.frame(c150_emmeans)
concentration_2_no0$predictions_c150 <- predict(c150, concentration_2_no0, re.form=NA, type="response")

concentration_150_graph <- ggplot(concentration_2_no0, aes(x = treatment, y = predictions_c150, color = treatment, shape = treatment)) +
  theme_test(base_size = 15) +
  scale_shape_manual(values = c("C" = "cross", "C1" = "circle", "C2" = "square", "C3" = "triangle", "C4" = "diamond"),
                     name = " ", labels = c("Control", "Piperine", "Tannic acid", "Eugenol", "Phytol")) + 
  scale_color_viridis(option = "D", discrete=TRUE, name = "Compounds") + 
  geom_point (data = concentration_2_no0, aes(x = treatment, y = Unknown.150_gamma, color = "black"), size = 3) +
  geom_point (data = concentration_2_no0, aes(x = treatment, y = Unknown.150_gamma, color = treatment), size = 2) +
  stat_summary(fun.data = mean_se, color = "black", size = 1) + 
  geom_errorbar(data = c150_emmeans, aes(x = treatment, y = response, ymin = asymp.LCL, ymax = asymp.UCL), width = 0.2, color = "black") +
  ylab ("A") +
  xlab (" ") +
  theme(legend.position = "none") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
  geom_text(x = 2, y = 0.03, label = "***", size = 8, color = "black") + 
  geom_text(x = 4, y = 0.165, label = "***", size = 8, color = "black")
  
concentration_150_graph

# 371
concentration_2_no0$Unknown.371_gamma <- concentration_2_no0$Unknown.371 + 0.0001
c371 <- glmmTMB(Unknown.371_gamma ~ treatment + (1|batID), data = concentration_2_no0, family=Gamma(link = log))
summary(c371)
plot(allEffects(c371))
c371_emmeans <-emmeans(c371,~ treatment, type="response")
c371_emmeans
effect_size_c371 <- eff_size(c371_emmeans, sigma= sigma(c371), edf = df.residual(c371)) # get effect sizes 
effect_size_c371
c371_emmeans <- as.data.frame(c371_emmeans)
concentration_2_no0$predictions_c371 <- predict(c371, concentration_2_no0, re.form=NA, type="response")

concentration_371_graph <- ggplot(concentration_2_no0, aes(x = treatment, y = predictions_c371, color = treatment, shape = treatment)) +
  theme_test(base_size = 15) +
  scale_color_viridis(option = "D", discrete=TRUE, name = "Compounds") + 
  scale_shape_manual(values = c("C" = "cross", "C1" = "circle", "C2" = "square", "C3" = "triangle", "C4" = "diamond"),
                     name = " ", labels = c("Control", "Piperine", "Tannic acid", "Eugenol", "Phytol")) + 
  geom_point (data = concentration_2_no0, aes(x = treatment, y = Unknown.371_gamma, color = "black"), size = 3) +
  geom_point (data = concentration_2_no0, aes(x = treatment, y = Unknown.371_gamma, color = treatment), size = 2) +
  stat_summary(fun.data = mean_se, color = "black", size = 1) + 
  geom_errorbar(data = c371_emmeans, aes(x = treatment, y = response, ymin = asymp.LCL, ymax = asymp.UCL), width = 0.2, color = "black") +
  ylab ("B") +
  xlab (" ") +
  theme(legend.position = "none") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
  geom_text(x = 5, y = 0.007, label = "***", size = 8, color = "black") 
concentration_371_graph

# 14
concentration_2_no0$Unknown.14_gamma <- concentration_2_no0$Unknown.14 + 0.0001
c14 <- glmmTMB(Unknown.14_gamma ~ treatment + (1|batID), data = concentration_2_no0, family=Gamma(link = log))
summary(c14)
plot(allEffects(c14))
c14_emmeans <-emmeans(c14,~ treatment, type="response")
c14_emmeans
effect_size_c14 <- eff_size(c14_emmeans, sigma= sigma(c14), edf = df.residual(c14)) # get effect sizes 
effect_size_c14
c14_emmeans <- as.data.frame(c14_emmeans)
concentration_2_no0$predictions_c14 <- predict(c14, concentration_2_no0, re.form=NA, type="response")

concentration_14_graph <- ggplot(concentration_2_no0, aes(x = treatment, y = predictions_c14, color = treatment, shape = treatment)) +
  theme_test(base_size = 15) +
  scale_color_viridis(option = "D", discrete=TRUE, name = "Compounds") + 
  geom_point (data = concentration_2_no0, aes(x = treatment, y = Unknown.14_gamma, color = "black"), size = 3) +
  geom_point (data = concentration_2_no0, aes(x = treatment, y = Unknown.14_gamma, color = treatment), size = 2) +
  scale_shape_manual(values = c("C" = "cross", "C1" = "circle", "C2" = "square", "C3" = "triangle", "C4" = "diamond"),
                     name = " ", labels = c("Control", "Piperine", "Tannic acid", "Eugenol", "Phytol")) + 
  stat_summary(fun.data = mean_se, color = "black", size = 1) + 
  geom_errorbar(data = c14_emmeans, aes(x = treatment, y = response, ymin = asymp.LCL, ymax = asymp.UCL), width = 0.2, color = "black") +
  ylab ("C") +
  xlab (" ") +
  theme(legend.position = "none") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
  geom_text(x = 2, y = 0.06, label = "***", size = 8, color = "black") + 
  geom_text(x = 3, y = 0.02, label = "***", size = 8, color = "black") + 
  geom_text(x = 4, y = 0.012, label = "***", size = 8, color = "black") + 
  geom_text(x = 5, y = 0.02, label = "***", size = 8, color = "black")
concentration_14_graph

# 69
concentration_2_no0$Unknown.69_gamma <- concentration_2_no0$Unknown.69 + 0.0001
c69 <- glmmTMB(Unknown.69_gamma ~ treatment + (1|batID), data = concentration_2_no0, family=Gamma(link = log))
summary(c69)
plot(allEffects(c69))
c69_emmeans <-emmeans(c69,~ treatment, type="response")
c69_emmeans
effect_size_c69 <- eff_size(c69_emmeans, sigma= sigma(c69), edf = df.residual(c69)) # get effect sizes 
effect_size_c69
c69_emmeans <- as.data.frame(c69_emmeans)
concentration_2_no0$predictions_c69 <- predict(c69, concentration_2_no0, re.form=NA, type="response")

concentration_69_graph <- ggplot(concentration_2_no0, aes(x = treatment, y = predictions_c69, color = treatment, shape = treatment)) +
  theme_test(base_size = 15) +
  scale_color_viridis(option = "D", discrete=TRUE, name = "Compounds") + 
  geom_point (data = concentration_2_no0, aes(x = treatment, y = Unknown.69_gamma, color = "black"), size = 3) +
  geom_point (data = concentration_2_no0, aes(x = treatment, y = Unknown.69_gamma, color = treatment), size = 2) +
  scale_shape_manual(values = c("C" = "cross", "C1" = "circle", "C2" = "square", "C3" = "triangle", "C4" = "diamond"),
                     name = " ", labels = c("Control", "Piperine", "Tannic acid", "Eugenol", "Phytol")) + 
  stat_summary(fun.data = mean_se, color = "black", size = 1) + 
  stat_summary(fun.data = mean_se, color = "black", size = 1) + 
  geom_errorbar(data = c69_emmeans, aes(x = treatment, y = response, ymin = asymp.LCL, ymax = asymp.UCL), width = 0.2, color = "black") +
  ylab ("D") +
  xlab (" ") +
  theme(legend.position = "none") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
  geom_text(x = 3, y = 0.0275, label = "***", size = 8, color = "black") + 
  ylim(0, 0.035)
concentration_69_graph

# 18
concentration_2_no0$Unknown.18_gamma <- concentration_2_no0$Unknown.18 + 0.0001
c18 <- glmmTMB(Unknown.18_gamma ~ treatment + (1|batID), data = concentration_2_no0, family=Gamma(link = log))
summary(c18)
plot(allEffects(c18))
c18_emmeans <-emmeans(c18,~ treatment, type="response")
c18_emmeans
effect_size_c18 <- eff_size(c18_emmeans, sigma= sigma(c18), edf = df.residual(c18)) # get effect sizes 
effect_size_c18
c18_emmeans <- as.data.frame(c18_emmeans)
concentration_2_no0$predictions_c18 <- predict(c18, concentration_2_no0, re.form=NA, type="response")

concentration_18_graph <- ggplot(concentration_2_no0, aes(x = treatment, y = predictions_c18, color = treatment, shape = treatment)) +
  theme_test(base_size = 15) +
  scale_color_viridis(option = "D", discrete=TRUE, name = "Compounds") + 
  geom_point (data = concentration_2_no0, aes(x = treatment, y = Unknown.18_gamma, color = "black"), size = 3) +
  geom_point (data = concentration_2_no0, aes(x = treatment, y = Unknown.18_gamma, color = treatment), size = 2) +
  scale_shape_manual(values = c("C" = "cross", "C1" = "circle", "C2" = "square", "C3" = "triangle", "C4" = "diamond"),
                     name = " ", labels = c("Control", "Piperine", "Tannic acid", "Eugenol", "Phytol")) + 
  stat_summary(fun.data = mean_se, color = "black", size = 1) + 
  geom_errorbar(data = c18_emmeans, aes(x = treatment, y = response, ymin = asymp.LCL, ymax = asymp.UCL), width = 0.2, color = "black") +
  ylab ("E") +
  xlab (" ") +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) + 
  scale_x_discrete(labels=c("C" = "Control", "C1" = "Piperine", "C2" = "Tannic acid", "C3" = "Eugenol", "C4" = "Phytol")) + 
  geom_text(x = 3, y = 0.015, label = "**", size = 8, color = "black") 
concentration_18_graph

# 370
concentration_2_no0$Unknown.370_gamma <- concentration_2_no0$Unknown.370 + 0.0001
c370 <- glmmTMB(Unknown.370_gamma ~ treatment + (1|batID), data = concentration_2_no0, family=Gamma(link = log))
summary(c370)
plot(allEffects(c370))
c370_emmeans <-emmeans(c370,~ treatment, type="response")
c370_emmeans
effect_size_c370 <- eff_size(c370_emmeans, sigma= sigma(c370), edf = df.residual(c370)) # get effect sizes 
effect_size_c370
c370_emmeans <- as.data.frame(c370_emmeans)
concentration_2_no0$predictions_c370 <- predict(c370, concentration_2_no0, re.form=NA, type="response")

concentration_370_graph <- ggplot(concentration_2_no0, aes(x = treatment, y = predictions_c370, color = treatment, shape = treatment)) +
  theme_test(base_size = 15) +
  scale_color_viridis(option = "D", discrete=TRUE, name = "Compounds") + 
  geom_point (data = concentration_2_no0, aes(x = treatment, y = Unknown.370_gamma, color = "black"), size = 3) +
  geom_point (data = concentration_2_no0, aes(x = treatment, y = Unknown.370_gamma, color = treatment), size = 2) +
  scale_shape_manual(values = c("C" = "cross", "C1" = "circle", "C2" = "square", "C3" = "triangle", "C4" = "diamond"),
                     name = " ", labels = c("Control", "Piperine", "Tannic acid", "Eugenol", "Phytol")) + 
  stat_summary(fun.data = mean_se, color = "black", size = 1) + 
  geom_errorbar(data = c370_emmeans, aes(x = treatment, y = response, ymin = asymp.LCL, ymax = asymp.UCL), width = 0.2, color = "black") +
  ylab ("Phytol") +
  xlab (" ") +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) + 
  scale_x_discrete(labels=c("C" = "Control", "C1" = "Piperine", "C2" = "Tannic acid", "C3" = "Eugenol", "C4" = "Phytol")) + 
  geom_text(x = 5, y = 0.20, label = "***", size = 8, color = "black") 
concentration_370_graph

# 46
concentration_2_no0$Unknown.46_gamma <- concentration_2_no0$Unknown.46 + 0.0001
c46 <- glmmTMB(Unknown.46_gamma ~ treatment + (1|batID), data = concentration_2_no0, family=Gamma(link = log))
summary(c46)
plot(allEffects(c46))
c46_emmeans <-emmeans(c46,~ treatment, type="response")
c46_emmeans
effect_size_c46 <- eff_size(c46_emmeans, sigma= sigma(c46), edf = df.residual(c46)) # get effect sizes 
effect_size_c46
c46_emmeans <- as.data.frame(c46_emmeans)
concentration_2_no0$predictions_c46 <- predict(c46, concentration_2_no0, re.form=NA, type="response")

concentration_46_graph <- ggplot(concentration_2_no0, aes(x = treatment, y = predictions_c46, color = treatment, shape = treatment)) +
  theme_test(base_size = 15) +
  scale_color_viridis(option = "D", discrete=TRUE, name = "Compounds") + 
  geom_point (data = concentration_2_no0, aes(x = treatment, y = Unknown.46_gamma, color = "black"), size = 3) +
  geom_point (data = concentration_2_no0, aes(x = treatment, y = Unknown.46_gamma, color = treatment), size = 2) +
  scale_shape_manual(values = c("C" = "cross", "C1" = "circle", "C2" = "square", "C3" = "triangle", "C4" = "diamond"),
                     name = " ", labels = c("Control", "Piperine", "Tannic acid", "Eugenol", "Phytol")) + 
  stat_summary(fun.data = mean_se, color = "black", size = 1) + 
  geom_errorbar(data = c46_emmeans, aes(x = treatment, y = response, ymin = asymp.LCL, ymax = asymp.UCL), width = 0.2, color = "black") +
  ylab ("F") +
  xlab (" ") +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) + 
  scale_x_discrete(labels=c("C" = "Control", "C1" = "Piperine", "C2" = "Tannic acid", "C3" = "Eugenol", "C4" = "Phytol")) + 
  geom_text(x = 5, y = 0.025, label = "***", size = 8, color = "black") 

concentration_46_graph

allcompounds_2_RF <- concentration_353_graph +
  concentration_226_graph +
  concentration_150_graph +
  concentration_371_graph +
  concentration_14_graph +
  concentration_69_graph +
  concentration_18_graph +
  concentration_370_graph +
  concentration_46_graph
allcompounds_2_RF

ggsave(file="allcompounds_2_RF.jpg", 
       plot=allcompounds_2_RF,
       width=9,height=9,units="in",dpi=300)

### 3% 

#356
concentration_3_no0$Unknown.356_gamma <- concentration_3_no0$Unknown.356 + 0.0001
c356 <- glmmTMB(Unknown.356_gamma ~ treatment + (1|batID), data = concentration_3_no0, family=Gamma(link = log))
summary(c356)
plot(allEffects(c356))
c356_emmeans <-emmeans(c356,~ treatment, type="response")
c356_emmeans
effect_size_c356 <- eff_size(c356_emmeans, sigma= sigma(c356), edf = df.residual(c356)) # get effect sizes 
effect_size_c356
c356_emmeans <- as.data.frame(c356_emmeans)
concentration_3_no0$predictions_c356 <- predict(c356, concentration_3_no0, re.form=NA, type="response")

concentration_356_graph <- ggplot(concentration_3_no0, aes(x = treatment, y = predictions_c356, color = treatment, shape = treatment)) +
  theme_test(base_size = 15) +
  scale_color_viridis(option = "D", discrete=TRUE, name = "Compounds") + 
  scale_shape_manual(values = c("C" = "cross", "C1" = "circle", "C2" = "square", "C3" = "triangle", "C4" = "diamond"),
                     name = " ", labels = c("Control", "Piperine", "Tannic acid", "Eugenol", "Phytol")) + 
  geom_point (data = concentration_3_no0, aes(x = treatment, y = Unknown.356_gamma, color = "black"), size = 3) +
  geom_point (data = concentration_3_no0, aes(x = treatment, y = Unknown.356_gamma, color = treatment), size = 2) +
  stat_summary(fun.data = mean_se, color = "black", size = 1) + 
  geom_errorbar(data = c356_emmeans, aes(x = treatment, y = response, ymin = asymp.LCL, ymax = asymp.UCL), width = 0.2, color = "black") +
  ylab ("A") +
  xlab ("") +
  theme(legend.position = "none") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
  geom_text(x = 4, y = 0.021, label = "***", size = 8, color = "black") 
concentration_356_graph

#353
concentration_3_no0$Unknown.353_gamma <- concentration_3_no0$Unknown.353 + 0.0001
c353 <- glmmTMB(Unknown.353_gamma ~ treatment + (1|batID), data = concentration_3_no0, family=Gamma(link = log))
summary(c353)
plot(allEffects(c353))
c353_emmeans <-emmeans(c353,~ treatment, type="response")
c353_emmeans
effect_size_c353 <- eff_size(c353_emmeans, sigma= sigma(c353), edf = df.residual(c353)) # get effect sizes 
effect_size_c353
c353_emmeans <- as.data.frame(c353_emmeans)
concentration_3_no0$predictions_c353 <- predict(c353, concentration_3_no0, re.form=NA, type="response")

concentration_353_graph <- ggplot(concentration_3_no0, aes(x = treatment, y = predictions_c353, color = treatment, shape = treatment)) +
  theme_test(base_size = 15) +
  scale_color_viridis(option = "D", discrete=TRUE, name = "Compounds") + 
  scale_shape_manual(values = c("C" = "cross", "C1" = "circle", "C2" = "square", "C3" = "triangle", "C4" = "diamond"),
                     name = " ", labels = c("Control", "Piperine", "Tannic acid", "Eugenol", "Phytol")) + 
  geom_point (data = concentration_3_no0, aes(x = treatment, y = Unknown.353_gamma, color = "black"), size = 3) +
  geom_point (data = concentration_3_no0, aes(x = treatment, y = Unknown.353_gamma, color = treatment), size = 2) +
  stat_summary(fun.data = mean_se, color = "black", size = 1) + 
  geom_errorbar(data = c353_emmeans, aes(x = treatment, y = response, ymin = asymp.LCL, ymax = asymp.UCL), width = 0.2, color = "black") +
  ylab ("Eugenol") +
  xlab ("") +
  theme(legend.position = "none") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
  geom_text(x = 4, y = 0.55, label = "***", size = 8, color = "black") 
concentration_353_graph

#226
concentration_3_no0$Unknown.226_gamma <- concentration_3_no0$Unknown.226 + 0.0001
c226 <- glmmTMB(Unknown.226_gamma ~ treatment + (1|batID), data = concentration_3_no0, family=Gamma(link = log))
summary(c226)
plot(allEffects(c226))
c226_emmeans <-emmeans(c226,~ treatment, type="response")
c226_emmeans
effect_size_c226 <- eff_size(c226_emmeans, sigma= sigma(c226), edf = df.residual(c226)) # get effect sizes 
effect_size_c226
c226_emmeans <- as.data.frame(c226_emmeans)
concentration_3_no0$predictions_c226 <- predict(c226, concentration_3_no0, re.form=NA, type="response")

concentration_226_graph <- ggplot(concentration_3_no0, aes(x = treatment, y = predictions_c226, color = treatment, shape = treatment)) +
  theme_test(base_size = 15) +
  scale_color_viridis(option = "D", discrete=TRUE, name = "Compounds") + 
  scale_shape_manual(values = c("C" = "cross", "C1" = "circle", "C2" = "square", "C3" = "triangle", "C4" = "diamond"),
                     name = " ", labels = c("Control", "Piperine", "Tannic acid", "Eugenol", "Phytol")) + 
  geom_point (data = concentration_3_no0, aes(x = treatment, y = Unknown.226_gamma, color = "black"), size = 3) +
  geom_point (data = concentration_3_no0, aes(x = treatment, y = Unknown.226_gamma, color = treatment), size = 2) +
  stat_summary(fun.data = mean_se, color = "black", size = 1) + 
  geom_errorbar(data = c226_emmeans, aes(x = treatment, y = response, ymin = asymp.LCL, ymax = asymp.UCL), width = 0.2, color = "black") +
  ylab ("Eugenol") +
  xlab ("") +
  theme(legend.position = "none") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
  geom_text(x = 4, y = 3.4, label = "***", size = 8, color = "black") +
  ylim(0, 3.5)
concentration_226_graph

# 354
concentration_3_no0$Unknown.354_gamma <- concentration_3_no0$Unknown.354 + 0.0001
c354 <- glmmTMB(Unknown.354_gamma ~ treatment + (1|batID), data = concentration_3_no0, family=Gamma(link = log))
summary(c354)
plot(allEffects(c354))
c354_emmeans <-emmeans(c354,~ treatment, type="response")
c354_emmeans
effect_size_c354 <- eff_size(c354_emmeans, sigma= sigma(c354), edf = df.residual(c354)) # get effect sizes 
effect_size_c354
c354_emmeans <- as.data.frame(c354_emmeans)
concentration_3_no0$predictions_c354 <- predict(c354, concentration_3_no0, re.form=NA, type="response")

concentration_354_graph <- ggplot(concentration_3_no0, aes(x = treatment, y = predictions_c354, color = treatment, shape = treatment)) +
  theme_test(base_size = 15) +
  scale_color_viridis(option = "D", discrete=TRUE, name = "Compounds") + 
  scale_shape_manual(values = c("C" = "cross", "C1" = "circle", "C2" = "square", "C3" = "triangle", "C4" = "diamond"),
                     name = " ", labels = c("Control", "Piperine", "Tannic acid", "Eugenol", "Phytol")) + 
  geom_point (data = concentration_3_no0, aes(x = treatment, y = Unknown.354_gamma, color = "black"), size = 3) +
  geom_point (data = concentration_3_no0, aes(x = treatment, y = Unknown.354_gamma, color = treatment), size = 2) +
  stat_summary(fun.data = mean_se, color = "black", size = 1) + 
  geom_errorbar(data = c354_emmeans, aes(x = treatment, y = response, ymin = asymp.LCL, ymax = asymp.UCL), width = 0.2, color = "black") +
  ylab ("B") +
  xlab ("") +
  theme(legend.position = "none") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
  geom_text(x = 4, y = 0.0355, label = "***", size = 8, color = "black") 
concentration_354_graph

# 371
concentration_3_no0$Unknown.371_gamma <- concentration_3_no0$Unknown.371 + 0.0001
c371 <- glmmTMB(Unknown.371_gamma ~ treatment + (1|batID), data = concentration_3_no0, family=Gamma(link = log))
summary(c371)
plot(allEffects(c371))
c371_emmeans <-emmeans(c371,~ treatment, type="response")
c371_emmeans
effect_size_c371 <- eff_size(c371_emmeans, sigma= sigma(c371), edf = df.residual(c371)) # get effect sizes 
effect_size_c371
c371_emmeans <- as.data.frame(c371_emmeans)
concentration_3_no0$predictions_c371 <- predict(c371, concentration_3_no0, re.form=NA, type="response")

concentration_371_graph <- ggplot(concentration_3_no0, aes(x = treatment, y = predictions_c371, color = treatment, shape = treatment)) +
  theme_test(base_size = 15) +
  scale_color_viridis(option = "D", discrete=TRUE, name = "Compounds") + 
  scale_shape_manual(values = c("C" = "cross", "C1" = "circle", "C2" = "square", "C3" = "triangle", "C4" = "diamond"),
                     name = " ", labels = c("Control", "Piperine", "Tannic acid", "Eugenol", "Phytol")) + 
  geom_point (data = concentration_3_no0, aes(x = treatment, y = Unknown.371_gamma, color = "black"), size = 3) +
  geom_point (data = concentration_3_no0, aes(x = treatment, y = Unknown.371_gamma, color = treatment), size = 2) +
  stat_summary(fun.data = mean_se, color = "black", size = 1) + 
  geom_errorbar(data = c371_emmeans, aes(x = treatment, y = response, ymin = asymp.LCL, ymax = asymp.UCL), width = 0.2, color = "black") +
  ylab ("C") +
  xlab ("") +
  theme(legend.position = "none") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
  geom_text(x = 5, y = 0.03, label = "***", size = 8, color = "black") 
concentration_371_graph

# 372
concentration_3_no0$Unknown.372_gamma <- concentration_3_no0$Unknown.372 + 0.0001
c372 <- glmmTMB(Unknown.372_gamma ~ treatment + (1|batID), data = concentration_3_no0, family=Gamma(link = log))
summary(c372)
plot(allEffects(c372))
c372_emmeans <-emmeans(c372,~ treatment, type="response")
c372_emmeans
effect_size_c372 <- eff_size(c372_emmeans, sigma= sigma(c372), edf = df.residual(c372)) # get effect sizes 
effect_size_c372
c372_emmeans <- as.data.frame(c372_emmeans)
concentration_3_no0$predictions_c372 <- predict(c372, concentration_3_no0, re.form=NA, type="response")

concentration_372_graph <- ggplot(concentration_3_no0, aes(x = treatment, y = predictions_c372, color = treatment, shape = treatment)) +
  theme_test(base_size = 15) +
  scale_color_viridis(option = "D", discrete=TRUE, name = "Compounds") + 
  scale_shape_manual(values = c("C" = "cross", "C1" = "circle", "C2" = "square", "C3" = "triangle", "C4" = "diamond"),
                     name = " ", labels = c("Control", "Piperine", "Tannic acid", "Eugenol", "Phytol")) + 
  geom_point (data = concentration_3_no0, aes(x = treatment, y = Unknown.372_gamma, color = "black"), size = 3) +
  geom_point (data = concentration_3_no0, aes(x = treatment, y = Unknown.372_gamma, color = treatment), size = 2) +
  stat_summary(fun.data = mean_se, color = "black", size = 1) + 
  geom_errorbar(data = c372_emmeans, aes(x = treatment, y = response, ymin = asymp.LCL, ymax = asymp.UCL), width = 0.2, color = "black") +
  ylab ("D") +
  xlab ("") +
  theme(legend.position = "none") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
  geom_text(x = 5, y = 0.015, label = "***", size = 8, color = "black") 
concentration_372_graph

# 334
concentration_3_no0$Unknown.334_gamma <- concentration_3_no0$Unknown.334 + 0.0001
c334 <- glmmTMB(Unknown.334_gamma ~ treatment + (1|batID), data = concentration_3_no0, family=Gamma(link = log))
summary(c334)
plot(allEffects(c334))
c334_emmeans <-emmeans(c334,~ treatment, type="response")
c334_emmeans
effect_size_c334 <- eff_size(c334_emmeans, sigma= sigma(c334), edf = df.residual(c334)) # get effect sizes 
effect_size_c334
c334_emmeans <- as.data.frame(c334_emmeans)
concentration_3_no0$predictions_c334 <- predict(c334, concentration_3_no0, re.form=NA, type="response")

concentration_334_graph <- ggplot(concentration_3_no0, aes(x = treatment, y = predictions_c334, color = treatment, shape = treatment)) +
  theme_test(base_size = 15) +
  scale_color_viridis(option = "D", discrete=TRUE, name = "Compounds") + 
  scale_shape_manual(values = c("C" = "cross", "C1" = "circle", "C2" = "square", "C3" = "triangle", "C4" = "diamond"),
                     name = " ", labels = c("Control", "Piperine", "Tannic acid", "Eugenol", "Phytol")) + 
  geom_point (data = concentration_3_no0, aes(x = treatment, y = Unknown.334_gamma, color = "black"), size = 3) +
  geom_point (data = concentration_3_no0, aes(x = treatment, y = Unknown.334_gamma, color = treatment), size = 2) +
  stat_summary(fun.data = mean_se, color = "black", size = 1) + 
  geom_errorbar(data = c334_emmeans, aes(x = treatment, y = response, ymin = asymp.LCL, ymax = asymp.UCL), width = 0.2, color = "black") +
  ylab ("E") +
  xlab ("") +
  theme(legend.position = "none") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
  geom_text(x = 3, y = 0.023, label = "***", size = 8, color = "black") 
concentration_334_graph

# 335
concentration_3_no0$Unknown.335_gamma <- concentration_3_no0$Unknown.335 + 0.0001
c335 <- glmmTMB(Unknown.335_gamma ~ treatment + (1|batID), data = concentration_3_no0, family=Gamma(link = log))
summary(c335)
plot(allEffects(c335))
c335_emmeans <-emmeans(c335,~ treatment, type="response")
c335_emmeans
effect_size_c335 <- eff_size(c335_emmeans, sigma= sigma(c335), edf = df.residual(c335)) # get effect sizes 
effect_size_c335
c335_emmeans <- as.data.frame(c335_emmeans)
concentration_3_no0$predictions_c335 <- predict(c335, concentration_3_no0, re.form=NA, type="response")

concentration_335_graph <- ggplot(concentration_3_no0, aes(x = treatment, y = predictions_c335, color = treatment, shape = treatment)) +
  theme_test(base_size = 15) +
  scale_color_viridis(option = "D", discrete=TRUE, name = "Compounds") + 
  scale_shape_manual(values = c("C" = "cross", "C1" = "circle", "C2" = "square", "C3" = "triangle", "C4" = "diamond"),
                     name = " ", labels = c("Control", "Piperine", "Tannic acid", "Eugenol", "Phytol")) + 
  geom_point (data = concentration_3_no0, aes(x = treatment, y = Unknown.335_gamma, color = "black"), size = 3) +
  geom_point (data = concentration_3_no0, aes(x = treatment, y = Unknown.335_gamma, color = treatment), size = 2) +
  stat_summary(fun.data = mean_se, color = "black", size = 1) + 
  geom_errorbar(data = c335_emmeans, aes(x = treatment, y = response, ymin = asymp.LCL, ymax = asymp.UCL), width = 0.2, color = "black") +
  ylab ("F") +
  xlab ("") +
  theme(legend.position = "none") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
  geom_text(x = 3, y = 0.0024, label = "***", size = 8, color = "black") 
concentration_335_graph

# 69
concentration_3_no0$Unknown.69_gamma <- concentration_3_no0$Unknown.69 + 0.0001
c69 <- glmmTMB(Unknown.69_gamma ~ treatment + (1|batID), data = concentration_3_no0, family=Gamma(link = log))
summary(c69)
plot(allEffects(c69))
c69_emmeans <-emmeans(c69,~ treatment, type="response")
c69_emmeans
effect_size_c69 <- eff_size(c69_emmeans, sigma= sigma(c69), edf = df.residual(c69)) # get effect sizes 
effect_size_c69
c69_emmeans <- as.data.frame(c69_emmeans)
concentration_3_no0$predictions_c69 <- predict(c69, concentration_3_no0, re.form=NA, type="response")

concentration_69_graph <- ggplot(concentration_3_no0, aes(x = treatment, y = predictions_c69, color = treatment, shape = treatment)) +
  theme_test(base_size = 15) +
  scale_color_viridis(option = "D", discrete=TRUE, name = "Compounds") + 
  scale_shape_manual(values = c("C" = "cross", "C1" = "circle", "C2" = "square", "C3" = "triangle", "C4" = "diamond"),
                     name = " ", labels = c("Control", "Piperine", "Tannic acid", "Eugenol", "Phytol")) + 
  geom_point (data = concentration_3_no0, aes(x = treatment, y = Unknown.69_gamma, color = "black"), size = 3) +
  geom_point (data = concentration_3_no0, aes(x = treatment, y = Unknown.69_gamma, color = treatment), size = 2) +
  stat_summary(fun.data = mean_se, color = "black", size = 1) + 
  geom_errorbar(data = c69_emmeans, aes(x = treatment, y = response, ymin = asymp.LCL, ymax = asymp.UCL), width = 0.2, color = "black") +
  ylab ("G") +
  xlab ("") +
  theme(legend.position = "none") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
  geom_text(x = 3, y = 0.47, label = "***", size = 8, color = "black") +
  ylim(0, 0.55)
concentration_69_graph

# 199
concentration_3_no0$Unknown.199_gamma <- concentration_3_no0$Unknown.199 + 0.0001
c199 <- glmmTMB(Unknown.199_gamma ~ treatment + (1|batID), data = concentration_3_no0, family=Gamma(link = log))
summary(c199)
plot(allEffects(c199))
c199_emmeans <-emmeans(c199,~ treatment, type="response")
c199_emmeans
effect_size_c199 <- eff_size(c199_emmeans, sigma= sigma(c199), edf = df.residual(c199)) # get effect sizes 
effect_size_c199
c199_emmeans <- as.data.frame(c199_emmeans)
concentration_3_no0$predictions_c199 <- predict(c199, concentration_3_no0, re.form=NA, type="response")

concentration_199_graph <- ggplot(concentration_3_no0, aes(x = treatment, y = predictions_c199, color = treatment, shape = treatment)) +
  theme_test(base_size = 15) +
  scale_color_viridis(option = "D", discrete=TRUE, name = "Compounds") + 
  scale_shape_manual(values = c("C" = "cross", "C1" = "circle", "C2" = "square", "C3" = "triangle", "C4" = "diamond"),
                     name = " ", labels = c("Control", "Piperine", "Tannic acid", "Eugenol", "Phytol")) + 
  geom_point (data = concentration_3_no0, aes(x = treatment, y = Unknown.199_gamma, color = "black"), size = 3) +
  geom_point (data = concentration_3_no0, aes(x = treatment, y = Unknown.199_gamma, color = treatment), size = 2) +
  stat_summary(fun.data = mean_se, color = "black", size = 1) + 
  geom_errorbar(data = c199_emmeans, aes(x = treatment, y = response, ymin = asymp.LCL, ymax = asymp.UCL), width = 0.2, color = "black") +
  ylab ("Gallic acid") +
  xlab ("") +
  theme(legend.position = "none") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
  geom_text(x = 3, y = 0.38, label = "***", size = 8, color = "black") +
  geom_text(x = 4, y = 0.05, label = "*", size = 8, color = "black") 
concentration_199_graph

# 7
concentration_3_no0$Unknown.7_gamma <- concentration_3_no0$Unknown.7 + 0.0001
c7 <- glmmTMB(Unknown.7_gamma ~ treatment + (1|batID), data = concentration_3_no0, family=Gamma(link = log))
summary(c7)
plot(allEffects(c7))
c7_emmeans <-emmeans(c7,~ treatment, type="response")
c7_emmeans
effect_size_c7 <- eff_size(c7_emmeans, sigma= sigma(c7), edf = df.residual(c7)) # get effect sizes 
effect_size_c7
c7_emmeans <- as.data.frame(c7_emmeans)
concentration_3_no0$predictions_c7 <- predict(c7, concentration_3_no0, re.form=NA, type="response")

concentration_7_graph <- ggplot(concentration_3_no0, aes(x = treatment, y = predictions_c7, color = treatment, shape = treatment)) +
  theme_test(base_size = 15) +
  scale_color_viridis(option = "D", discrete=TRUE, name = "Compounds") + 
  scale_shape_manual(values = c("C" = "cross", "C1" = "circle", "C2" = "square", "C3" = "triangle", "C4" = "diamond"),
                     name = " ", labels = c("Control", "Piperine", "Tannic acid", "Eugenol", "Phytol")) + 
  geom_point (data = concentration_3_no0, aes(x = treatment, y = Unknown.7_gamma, color = "black"), size = 3) +
  geom_point (data = concentration_3_no0, aes(x = treatment, y = Unknown.7_gamma, color = treatment), size = 2) +
  stat_summary(fun.data = mean_se, color = "black", size = 1) + 
  geom_errorbar(data = c7_emmeans, aes(x = treatment, y = response, ymin = asymp.LCL, ymax = asymp.UCL), width = 0.2, color = "black") +
  ylab ("H") +
  xlab ("") +
  theme(legend.position = "none") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
  geom_text(x = 2, y = 0.7, label = "*", size = 8, color = "black") +
  geom_text(x = 3, y = 0.7, label = "**", size = 8, color = "black") 
concentration_7_graph

# 33
concentration_3_no0$Unknown.33_gamma <- concentration_3_no0$Unknown.33 + 0.0001
c33 <- glmmTMB(Unknown.33_gamma ~ treatment + (1|batID), data = concentration_3_no0, family=Gamma(link = log))
summary(c33)
plot(allEffects(c33))
c33_emmeans <-emmeans(c33,~ treatment, type="response")
c33_emmeans
effect_size_c33 <- eff_size(c33_emmeans, sigma= sigma(c33), edf = df.residual(c33)) # get effect sizes 
effect_size_c33
c33_emmeans <- as.data.frame(c33_emmeans)
concentration_3_no0$predictions_c33 <- predict(c33, concentration_3_no0, re.form=NA, type="response")

concentration_33_graph <- ggplot(concentration_3_no0, aes(x = treatment, y = predictions_c33, color = treatment, shape = treatment)) +
  theme_test(base_size = 15) +
  scale_color_viridis(option = "D", discrete=TRUE, name = "Compounds") + 
  scale_shape_manual(values = c("C" = "cross", "C1" = "circle", "C2" = "square", "C3" = "triangle", "C4" = "diamond"),
                     name = " ", labels = c("Control", "Piperine", "Tannic acid", "Eugenol", "Phytol")) + 
  geom_point (data = concentration_3_no0, aes(x = treatment, y = Unknown.33_gamma, color = "black"), size = 3) +
  geom_point (data = concentration_3_no0, aes(x = treatment, y = Unknown.33_gamma, color = treatment), size = 2) +
  stat_summary(fun.data = mean_se, color = "black", size = 1) + 
  geom_errorbar(data = c33_emmeans, aes(x = treatment, y = response, ymin = asymp.LCL, ymax = asymp.UCL), width = 0.2, color = "black") +
  ylab ("I") +
  xlab ("") +
  theme(legend.position = "none") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
  geom_text(x = 3, y = 0.02, label = "***", size = 8, color = "black") +
  geom_text(x = 5, y = 0.03, label = "*", size = 8, color = "black") 
concentration_33_graph

# 31
concentration_3_no0$Unknown.31_gamma <- concentration_3_no0$Unknown.31 + 0.0001
c31 <- glmmTMB(Unknown.31_gamma ~ treatment + (1|batID), data = concentration_3_no0, family=Gamma(link = log))
summary(c31)
plot(allEffects(c31))
c31_emmeans <-emmeans(c31,~ treatment, type="response")
c31_emmeans
effect_size_c31 <- eff_size(c31_emmeans, sigma= sigma(c31), edf = df.residual(c31)) # get effect sizes 
effect_size_c31
c31_emmeans <- as.data.frame(c31_emmeans)
concentration_3_no0$predictions_c31 <- predict(c31, concentration_3_no0, re.form=NA, type="response")

concentration_31_graph <- ggplot(concentration_3_no0, aes(x = treatment, y = predictions_c31, color = treatment, shape = treatment)) +
  theme_test(base_size = 15) +
  scale_color_viridis(option = "D", discrete=TRUE, name = "Compounds") + 
  scale_shape_manual(values = c("C" = "cross", "C1" = "circle", "C2" = "square", "C3" = "triangle", "C4" = "diamond"),
                     name = " ", labels = c("Control", "Piperine", "Tannic acid", "Eugenol", "Phytol")) + 
  geom_point(data = concentration_3_no0, aes(x = treatment, y = Unknown.31_gamma, color = "black"), size = 3) +
  geom_point(data = concentration_3_no0, aes(x = treatment, y = Unknown.31_gamma, color = treatment), size = 2) +
  stat_summary(fun.data = mean_se, color = "black", size = 1) + 
  geom_errorbar(data = c31_emmeans, aes(x = treatment, y = response, ymin = asymp.LCL, ymax = asymp.UCL), width = 0.2, color = "black") +
  ylab ("J") +
  xlab ("") +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) + 
  scale_x_discrete(labels=c("C" = "Control", "C1" = "Piperine", "C2" = "Tannic acid", "C3" = "Eugenol", "C4" = "Phytol")) + 
  geom_text(x = 3, y = 0.010, label = "***", size = 8, color = "black") +
  geom_text(x = 5, y = 0.025, label = "***", size = 8, color = "black") 
concentration_31_graph


# 94
concentration_3_no0$Unknown.94_gamma <- concentration_3_no0$Unknown.94 + 0.0001
c94 <- glmmTMB(Unknown.94_gamma ~ treatment + (1|batID), data = concentration_3_no0, family=Gamma(link = log))
summary(c94)
plot(allEffects(c94))
c94_emmeans <-emmeans(c94,~ treatment, type="response")
c94_emmeans
effect_size_c94 <- eff_size(c94_emmeans, sigma= sigma(c94), edf = df.residual(c94)) # get effect sizes 
effect_size_c94
c94_emmeans <- as.data.frame(c94_emmeans)
concentration_3_no0$predictions_c94 <- predict(c94, concentration_3_no0, re.form=NA, type="response")

concentration_94_graph <- ggplot(concentration_3_no0, aes(x = treatment, y = predictions_c94, color = treatment, shape = treatment)) +
  theme_test(base_size = 15) +
  scale_color_viridis(option = "D", discrete=TRUE, name = "Compounds") + 
  scale_shape_manual(values = c("C" = "cross", "C1" = "circle", "C2" = "square", "C3" = "triangle", "C4" = "diamond"),
                     name = " ", labels = c("Control", "Piperine", "Tannic acid", "Eugenol", "Phytol")) + 
  geom_point (data = concentration_3_no0, aes(x = treatment, y = Unknown.94_gamma, color = "black"), size = 3) +
  geom_point (data = concentration_3_no0, aes(x = treatment, y = Unknown.94_gamma, color = treatment), size = 2) +
  stat_summary(fun.data = mean_se, color = "black", size = 1) + 
  geom_errorbar(data = c94_emmeans, aes(x = treatment, y = response, ymin = asymp.LCL, ymax = asymp.UCL), width = 0.2, color = "black") +
  ylab ("K") +
  xlab ("") +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) + 
  scale_x_discrete(labels=c("C" = "Control", "C1" = "Piperine", "C2" = "Tannic acid", "C3" = "Eugenol", "C4" = "Phytol")) + 
  geom_text(x = 4, y = 0.32, label = "***", size = 8, color = "black") +
  ylim(0, 0.4)
concentration_94_graph

# 370
concentration_3_no0$Unknown.370_gamma <- concentration_3_no0$Unknown.370 + 0.0001
c370 <- glmmTMB(Unknown.370_gamma ~ treatment + (1|batID), data = concentration_3_no0, family=Gamma(link = log))
summary(c370)
plot(allEffects(c370))
c370_emmeans <-emmeans(c370,~ treatment, type="response")
c370_emmeans
effect_size_c370 <- eff_size(c370_emmeans, sigma= sigma(c370), edf = df.residual(c370)) # get effect sizes 
effect_size_c370
c370_emmeans <- as.data.frame(c370_emmeans)
concentration_3_no0$predictions_c370 <- predict(c370, concentration_3_no0, re.form=NA, type="response")

concentration_370_graph <- ggplot(concentration_3_no0, aes(x = treatment, y = predictions_c370, color = treatment, shape = treatment)) +
  theme_test(base_size = 15) +
  scale_color_viridis(option = "D", discrete=TRUE, name = "Compounds") + 
  scale_shape_manual(values = c("C" = "cross", "C1" = "circle", "C2" = "square", "C3" = "triangle", "C4" = "diamond"),
                     name = " ", labels = c("Control", "Piperine", "Tannic acid", "Eugenol", "Phytol")) + 
  geom_point (data = concentration_3_no0, aes(x = treatment, y = Unknown.370_gamma, color = "black"), size = 3) +
  geom_point (data = concentration_3_no0, aes(x = treatment, y = Unknown.370_gamma, color = treatment), size = 2) +
  stat_summary(fun.data = mean_se, color = "black", size = 1) + 
  geom_errorbar(data = c370_emmeans, aes(x = treatment, y = response, ymin = asymp.LCL, ymax = asymp.UCL), width = 0.2, color = "black") +
  ylab ("Phytol") +
  xlab ("") +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) + 
  scale_x_discrete(labels=c("C" = "Control", "C1" = "Piperine", "C2" = "Tannic acid", "C3" = "Eugenol", "C4" = "Phytol")) + 
  geom_text(x = 5, y = 0.35, label = "***", size = 8, color = "black") 
concentration_370_graph

allcompounds_3_RF <- (concentration_356_graph+
                               concentration_353_graph+
                               concentration_226_graph+
                               concentration_354_graph+
                               concentration_371_graph+
                               concentration_372_graph+
                               concentration_334_graph+
                               concentration_335_graph+
                               concentration_69_graph+
                               concentration_199_graph+
                               concentration_7_graph+
                               concentration_33_graph+
                               concentration_31_graph+
                               concentration_94_graph+
                               concentration_370_graph) + plot_layout(ncol = 3)

ggsave(file="allcompounds_3_RF.jpg", 
       plot=allcompounds_3_RF,
       width=10,height=13,units="in",dpi=300)

##################################################
### Compound classification using classyfireR ####
##################################################

### Subset data per compound identity 

class_compounds <- read.csv("RF_classification_compoundsID.csv")
head(class_compounds)
InChI_Keys <- class_compounds$InChIKey # a vector with the InChIKey
Classification_List <- purrr::map(InChI_Keys, get_classification)
Classification_List # update RF_classification.csv file with classification results

