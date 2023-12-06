#################################################################################################
### Objective 2. The effects of secondary metabolites consumption on the bat fecal metabolome ###
#################################################################################################

rm(list=ls())

library(tidyverse)
library(dplyr)
library(vegan)
library(viridis)
library(ggplot2)
library(ggpubr)

### Get data ### 
data_final <- read.csv("data_final_metaMS_dec2_normalizedinsresponse.csv")
data_final$concentration <- as.factor(data_final$concentration)
data_final <- data_final %>%
  select(where(~ any(. != 0)))
head(data_final)

### Subset data in different concentrations (all compounds in 0.1, 2, and 3)
concentration_0.1 <- data_final %>% filter(batID %in% c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10"))
concentration_2 <- data_final %>% filter(batID %in% c("11", "12", "13", "14", "15", "16", "17", "18", "19", "20"))
concentration_3 <- data_final %>% filter(batID %in% c("21", "22", "23", "24", "25", "26", "27", "28", "29", "30"))

############
### NMDS ###
############

set.seed(1)

### 0.1
concentration_0.1_no0 <- concentration_0.1 %>%
  select(where(~ any(. != 0))) # elimnate variable with all zeros 

mmatrix <- concentration_0.1_no0[, 2:211] # select just the peaks
matrix <- as.matrix(mmatrix) # turn data frame into matrix
matrix
nmds_results <- metaMDS(mmatrix,
                        trymax = 300,
                        distance = "bray", # Specify a bray-curtis distance
                        try = 1000) # Number of iterations

nmds_results # Stress: 0.1556059 
stressplot(nmds_results)
ordiplot(nmds_results, type="text")
plot(nmds_results, type = "t")

scores(nmds_results$points)

data.scores <- as.data.frame(scores(nmds_results$points))
data.scores$treatment <- concentration_0.1_no0$treatment
data.scores$concentration <- concentration_0.1_no0$concentration

# PERMANOVA, 0.1
adonis_concentration_0.1 <- adonis2(matrix ~ data.scores$treatment, perm=9999)
adonis_concentration_0.1
# BETADISP, 0.1 
dist_matrix <- vegdist(matrix, method = "bray")
groups <- concentration_0.1_no0$treatment
dispersal <- betadisper(dist_matrix, groups, type = c("centroid"))
anova(dispersal)
tukey <- TukeyHSD(dispersal)
tukey

nmdsgraph_concentration_0.1 <- ggplot(data.scores, aes(x = MDS1, y = MDS2, colour = treatment)) +
  theme_classic(base_size = 20) + 
  ggtitle("(A) 0.1%")+ 
  theme(plot.title = element_text(hjust = 0.5)) + 
  geom_point(size = 5, alpha = 0.8) +
  ylab ("NMDS2") +
  xlab ("NMDS1") + 
  stat_ellipse(level = 0.95, linewidth = 1, linetype = "dashed") +  
  theme(legend.position = "bottom") + 
  annotate("text", x = 2, y = 1, size = 6, label = paste("Stress = 0.156\nPERMDISP2, P = 0.621\nPERMANOVA, P = 0.501")) +
  scale_color_viridis(option = "D", discrete=TRUE, name = "Metabolites", labels = c("Control", "Piperine", "Tannic acid", "Eugenol", "Phytol")) 
nmdsgraph_concentration_0.1

### 2
concentration_2_no0 <- concentration_2 %>%
  select(where(~ any(. != 0))) # eliminate variable with all zeros 

mmatrix <- concentration_2_no0[ , 2:151] # select just the peaks
matrix <- as.matrix(mmatrix) # turn data frame into matrix
matrix
nmds_results <- metaMDS(mmatrix,
                        trymax = 300,
                        distance = "bray", # Specify a bray-curtis distance
                        try = 1000) # Number of iterations

nmds_results # Stress:     0.2202312 
stressplot(nmds_results)
ordiplot(nmds_results, type="text")
plot(nmds_results, type = "t")
scores(nmds_results$points)
data.scores <- as.data.frame(scores(nmds_results$points))
data.scores$treatment <- concentration_2_no0$treatment
data.scores$concentration <- concentration_2_no0$concentration

# PERMANOVA, 2
adonis_concentration_2 <- adonis2(matrix ~ data.scores$treatment, perm=9999)
adonis_concentration_2
# BETADISP, 2 
dist_matrix <- vegdist(matrix, method = "bray")
groups <- concentration_2_no0$treatment
dispersal <- betadisper(dist_matrix, groups, type = c("centroid"))
anova(dispersal)
tukey <- TukeyHSD(dispersal)
tukey

nmdsgraph_concentration_2 <- ggplot(data.scores, aes(x = MDS1, y = MDS2, colour = treatment)) +
  theme_classic(base_size = 20) + 
  ggtitle("(B) 2%")+ 
  theme(plot.title = element_text(hjust = 0.5)) + 
  geom_point(size = 5, alpha = 0.8) +
  scale_color_viridis(option = "D", discrete=TRUE, name = "Compounds, 2%") + 
  ylab ("NMDS2") +
  xlab ("NMDS1") + 
  stat_ellipse(level = 0.95, linewidth = 1, linetype = "dashed") +   
  annotate("text", x = 1, y = 1, size = 6, label = paste("Stress = 0.220\nPERMDISP2, P = 0.253\nPERMANOVA, P = 0.069")) +
  theme(legend.position = "none") 
nmdsgraph_concentration_2

### 3
concentration_3_no0 <- concentration_3 %>%
  select(where(~ any(. != 0))) # eliminate variable with all zeros 

mmatrix <- concentration_3_no0[ , 2:200] # select just the peaks
matrix <- as.matrix(mmatrix) # turn data frame into matrix
matrix
nmds_results <- metaMDS(mmatrix,
                        trymax = 300,
                        distance = "bray", # Specify a bray-curtis distance
                        try = 1000) # Number of iterations

nmds_results # Stress:     0.1526375 
stressplot(nmds_results)
ordiplot(nmds_results, type="text")
plot(nmds_results, type = "t")

scores(nmds_results$points)

data.scores <- as.data.frame(scores(nmds_results$points))
data.scores$treatment <- concentration_3_no0$treatment
data.scores$concentration <- concentration_3_no0$concentration

# PERMANOVA, 3
adonis_concentration_3 <- adonis2(matrix ~ data.scores$treatment, perm=9999)
adonis_concentration_3
# BETADISP, 3
dist_matrix <- vegdist(matrix, method = "bray")
groups <- concentration_3_no0$treatment
dispersal <- betadisper(dist_matrix, groups, type = c("centroid"))
anova(dispersal)
tukey <- TukeyHSD(dispersal)
tukey

nmdsgraph_concentration_3 <- ggplot(data.scores, aes(x = MDS1, y = MDS2, colour = treatment)) +
  theme_classic(base_size = 20) + 
  ggtitle("(C) 3%")+ 
  theme(plot.title = element_text(hjust = 0.5)) + 
  geom_point(size = 5, alpha = 0.8) +
  ylab ("NMDS2") +
  xlab ("NMDS1") + 
  stat_ellipse(level = 0.95, linewidth = 1) + 
  scale_color_viridis(option = "D", discrete=TRUE, name = "Compounds", labels = c("Control", "Piperine", "Tannic acid", "Eugenol", "Phytol")) +
  annotate("text", x = 1.5, y = 1, size = 6, label = paste("Stress = 0.152\nPERMDISP2, P = 0.171\nPERMANOVA, P = 0.003**")) +
  theme(legend.position = "right")

nmdsgraph_concentration_3

NMDS_graph <- ggarrange(nmdsgraph_concentration_0.1,
                        nmdsgraph_concentration_2,
                        nmdsgraph_concentration_3,
                        ncol = 1, nrow = 3,
                        common.legend = TRUE)
NMDS_graph

ggsave(file="NMDS_graph.jpg", 
       plot= NMDS_graph,
       width=9,height=20,units="in",dpi=300)
