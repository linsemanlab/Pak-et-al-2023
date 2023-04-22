
##Principal Components Analysis (PCA) - Proteomics Data 
##Paper Title: Non-canonical Metabolic and Molecular Effects of Calorie Restriction Are Revealed by Varying Temporal Conditions

######
#Packages Needed 

#Used for data tidying functions
library(dplyr)
library(tidyr)

#Used for PCA 
install.packages("FactoMineR")
library(FactoMineR)

#Used for PCA visualization 
install.packages("factoextra")
library(factoextra)
install.packages("pamr")
library(pamr)
install.packages("ggfortify")
library(ggfortify)
install_github("vqv/ggbiplot", force = TRUE, dependencies=TRUE)
library(ggbiplot)

#######
#PCA Analaysis 

#Creating new dataframe with only timepoint A - removing text columns 
Pro_log_imputed_scaled_8_hr <- Pro_log_scaled %>% filter(Timepoint == "A")
View(Pro_log_imputed_scaled_8_hr)
Pro_log_8_hr <- Pro_log_imputed_scaled_8_hr[ -c(1:5) ]
View(Pro_log_8_hr)

#Creating new dataframe with only timepoint B - removing text columns 
Pro_log_imputed_scaled_12_hr <- Pro_log_scaled %>% filter(Timepoint == "B")
View(Pro_log_imputed_scaled_12_hr)
Pro_log_12_hr <- Pro_log_imputed_scaled_12_hr[ -c(1:5) ]
View(Pro_log_12_hr)

#Creating new dataframe with only timepoint C - removing text columns 
Pro_log_imputed_scaled_16_hr <- Pro_log_scaled %>% filter(Timepoint == "C")
View(Pro_log_imputed_scaled_16_hr)
Pro_log_16_hr <- Pro_log_imputed_scaled_16_hr[ -c(1:5) ]
View(Pro_log_16_hr)

#Creating new dataframe with only timepoint D - removing text columns 
Pro_log_imputed_scaled_24_hr <- Pro_log_scaled %>% filter(Timepoint == "D")
View(Pro_log_imputed_scaled_24_hr)
Pro_log_24_hr <- Pro_log_imputed_scaled_24_hr[ -c(1:5) ]
View(Pro_log_24_hr)

######
#Running the PCA for each timepoint 

#PCA - Timepoint A (8 hours)
pro.res.pca.8 <- PCA(Pro_log_8_hr, scale.unit = FALSE, ncp = 5, graph = FALSE)
print(pro.res.pca.8)
eig.val <- get_eigenvalue(pro.res.pca.8)
eig.val
fviz_eig(pro.res.pca.8, addlabels = TRUE, ylim = c(0, 35))#Scree plot
#Top contributors to each PCA dimension 
fviz_contrib(pro.res.pca.8, choice = "var", axes = 1, top = 50)
fviz_contrib(pro.res.pca.8, choice = "var", axes = 2, top = 50)
fviz_contrib(pro.res.pca.8, choice = "var", axes = 1:2, top = 50)

pro.res.pca.12 <- PCA(Pro_log_12_hr, scale.unit = FALSE, ncp = 5, graph = FALSE)
print(pro.res.pca.12)
eig.val <- get_eigenvalue(pro.res.pca.12)
eig.val
fviz_eig(pro.res.pca.12, addlabels = TRUE, ylim = c(0, 35))#Scree plot
#Top contributors to each PCA dimension 
fviz_contrib(pro.res.pca.12, choice = "var", axes = 1, top = 50)
fviz_contrib(pro.res.pca.12, choice = "var", axes = 2, top = 50)
fviz_contrib(pro.res.pca.12, choice = "var", axes = 1:2, top = 50)

pro.res.pca.16 <- PCA(Pro_log_16_hr, scale.unit = FALSE, ncp = 5, graph = FALSE)
print(pro.res.pca.16)
eig.val <- get_eigenvalue(pro.res.pca.16)
eig.val
fviz_eig(pro.res.pca.16, addlabels = TRUE, ylim = c(0, 35))#Scree plot
#Top contributors to each PCA dimension 
fviz_contrib(pro.res.pca.16, choice = "var", axes = 1, top = 50)
fviz_contrib(pro.res.pca.16, choice = "var", axes = 2, top = 50)
fviz_contrib(pro.res.pca.16, choice = "var", axes = 1:2, top = 50)

pro.res.pca.24 <- PCA(Pro_log_24_hr, scale.unit = FALSE, ncp = 5, graph = FALSE)
print(pro.res.pca.24)
eig.val <- get_eigenvalue(pro.res.pca.24)
eig.val
fviz_eig(pro.res.pca.24, addlabels = TRUE, ylim = c(0, 35))#Scree plot
#Top contributors to each PCA dimension 
fviz_contrib(pro.res.pca.24, choice = "var", axes = 1, top = 50)
fviz_contrib(pro.res.pca.24, choice = "var", axes = 2, top = 50)
fviz_contrib(pro.res.pca.24, choice = "var", axes = 1:2, top = 50)

#####
#Final plots 

#PCA - Timepoint A (8 hours)
PCA_8 <- fviz_pca_ind(pro.res.pca.8,
                      geom.ind = "point", # show points only (but not "text")
                      col.ind = Pro_log_imputed_scaled_8_hr$Group,# color by groups
                      mean.point = FALSE,
                      axes.linetype = NA,
                      palette = c("#282A74", "#72C8F1"),
                      pointsize = 10,
                      addEllipses = TRUE, # Concentration ellipses
                      ellipse.type = "confidence",
                      legend.title = "Group_Timepoint")+
  theme_minimal()+
  scale_shape_manual(values=c(16,18))+
  scale_y_continuous(limits = c(-12.5, 12.5), breaks = seq(-10, 10, by = 5))+ 
  scale_x_continuous(limits = c(-12.5, 12.5), breaks = seq(-10, 10, by = 5))+  
  annotate("text", x = 10, y = 12, label = "8hr Fast", family = "Arial", fontface = "bold", size = 12)+
  theme(axis.ticks = element_line(linewidth = 1.5), axis.ticks.length=unit(.5, "cm"), axis.line = element_blank(), plot.margin = unit(c(1,1,1,1), "cm"),axis.title.y = element_text(margin = margin(t = 0, r = 25, b = 0, l = 0)),axis.title.x = element_text(margin = margin(t = 25, r = 0, b = 0, l = 0)))+
  theme(panel.border = element_rect(fill = "transparent", color = "black", linewidth = 2, linetype = "solid"), panel.grid = element_blank(), text = element_text(family = "Arial"), legend.position="right", legend.justification = "top", legend.box.background = element_rect(color = "black", linewidth = 1),legend.title = element_blank(),axis.text = element_text(size = 40, family = "Arial", color = "black", face = "bold"), axis.title = element_text(size = 40, family = "Arial", face = "bold"), legend.text = element_text(size = 25, face = "bold", family = "Arial"), legend.key.size = unit(2, "cm"), plot.title = element_blank())+
  labs(x = "Component 1 (62.02%)", y = "Component 2 (14.36%)")
PCA_8

#Exporting plot as a PNG 
png(filename = 'Pro_PCA_8hr.png', type = 'cairo', units = 'in',
    width = 14, height = 10, pointsize = 16, res = 600)
PCA_8
dev.off()

#PCA - Timepoint B (12 hours)
PCA_12 <- fviz_pca_ind(pro.res.pca.12,
                       geom.ind = "point", # show points only (but not "text")
                       col.ind = Pro_log_imputed_scaled_12_hr$Group,# color by groups
                       mean.point = FALSE,
                       axes.linetype = NA,
                       palette = c("#282A74", "#72C8F1"),
                       pointsize = 10, 
                       addEllipses = TRUE, # Concentration ellipses
                       ellipse.type = "confidence",
                       legend.title = "Group_Timepoint")+
  theme_minimal()+
  scale_shape_manual(values=c(16,18))+
  scale_y_continuous(limits = c(-12.5, 12.5), breaks = seq(-10, 10, by = 5))+ 
  scale_x_continuous(limits = c(-12.5, 12.5), breaks = seq(-10, 10, by = 5))+ 
  annotate("text", x = 10, y = 12, label = "12hr Fast", family = "Arial", fontface = "bold", size = 12)+
  theme(axis.ticks = element_line(linewidth = 1.5), axis.ticks.length=unit(.5, "cm"), axis.line = element_blank(), plot.margin = unit(c(1,1,1,1), "cm"),axis.title.y = element_text(margin = margin(t = 0, r = 25, b = 0, l = 0)),axis.title.x = element_text(margin = margin(t = 25, r = 0, b = 0, l = 0)))+
  theme(panel.border = element_rect(fill = "transparent", color = "black", linewidth = 2, linetype = "solid"), panel.grid = element_blank(), text = element_text(family = "Arial"), legend.position="right", legend.justification = "top", legend.box.background = element_rect(color = "black", linewidth = 1),legend.title = element_blank(),axis.text = element_text(size = 40, family = "Arial", color = "black", face = "bold"), axis.title = element_text(size = 40, family = "Arial", face = "bold"), legend.text = element_text(size = 25, face = "bold", family = "Arial"), legend.key.size = unit(2, "cm"), plot.title = element_blank())+
  labs(x = "Component 1 (50.50%)", y = "Component 2 (17.16%)")
PCA_12

#Exporting plot as a PNG 
png(filename = 'Pro_PCA_12hr.png', type = 'cairo', units = 'in',
    width = 14, height = 10, pointsize = 16, res = 600)
PCA_12
dev.off()

#PCA - Timepoint C (16 hours)
PCA_16 <- fviz_pca_ind(pro.res.pca.16,
                       geom.ind = "point", # show points only (but not "text")
                       col.ind = Pro_log_imputed_scaled_16_hr$Group,# color by groups
                       mean.point = FALSE,
                       axes.linetype = NA,
                       palette = c("#282A74", "#72C8F1"),
                       pointsize = 10,
                       addEllipses = TRUE, # Concentration ellipses
                       ellipse.type = "confidence",
                       legend.title = "Group_Timepoint")+
  theme_minimal()+
  scale_shape_manual(values=c(16,18))+
  scale_y_continuous(limits = c(-12.5, 12.5), breaks = seq(-10, 10, by = 5))+ 
  scale_x_continuous(limits = c(-12.5, 12.5), breaks = seq(-10, 10, by = 5))+ 
  annotate("text", x = 10, y = 12, label = "16hr Fast", family = "Arial", fontface = "bold", size = 12)+
  theme(axis.ticks = element_line(linewidth = 1.5), axis.ticks.length=unit(.5, "cm"), axis.line = element_blank(), plot.margin = unit(c(1,1,1,1), "cm"),axis.title.y = element_text(margin = margin(t = 0, r = 25, b = 0, l = 0)),axis.title.x = element_text(margin = margin(t = 25, r = 0, b = 0, l = 0)))+
  theme(panel.border = element_rect(fill = "transparent", color = "black", linewidth = 2, linetype = "solid"), panel.grid = element_blank(), text = element_text(family = "Arial"), legend.position="right", legend.justification = "top", legend.box.background = element_rect(color = "black", linewidth = 1),legend.title = element_blank(),axis.text = element_text(size = 40, family = "Arial", color = "black", face = "bold"), axis.title = element_text(size = 40, family = "Arial", face = "bold"), legend.text = element_text(size = 25, face = "bold", family = "Arial"), legend.key.size = unit(2, "cm"), plot.title = element_blank())+
  labs(x = "Component 1 (33.02%)", y = "Component 2 (22.50%)")
PCA_16

#Exporting plot as a PNG 
png(filename = 'Pro_PCA_16hr.png', type = 'cairo', units = 'in',
    width = 14, height = 10, pointsize = 16, res = 600)
PCA_16
dev.off()

#PCA - Timepoint D (24 hours)
PCA_24 <- fviz_pca_ind(pro.res.pca.24,
                       geom.ind = "point", # show points only (but not "text")
                       col.ind = Pro_log_imputed_scaled_24_hr$Group,# color by groups
                       mean.point = FALSE,
                       axes.linetype = NA,
                       palette = c("#282A74", "#72C8F1"),
                       pointsize = 10,
                       addEllipses = TRUE, # Concentration ellipses
                       ellipse.type = "confidence",
                       legend.title = "Group_Timepoint")+
  theme_minimal()+
  scale_shape_manual(values=c(16,18))+
  scale_y_continuous(limits = c(-12.5, 12.5), breaks = seq(-10, 10, by = 5))+ 
  scale_x_continuous(limits = c(-12.5, 12.5), breaks = seq(-10, 10, by = 5))+ 
  annotate("text", x = 10, y = 12, label = "24hr Fast", family = "Arial", fontface = "bold", size = 12)+
  theme(axis.ticks = element_line(linewidth = 1.5), axis.ticks.length=unit(.5, "cm"), axis.line = element_blank(), plot.margin = unit(c(1,1,1,1), "cm"),axis.title.y = element_text(margin = margin(t = 0, r = 25, b = 0, l = 0)),axis.title.x = element_text(margin = margin(t = 25, r = 0, b = 0, l = 0)))+
  theme(panel.border = element_rect(fill = "transparent", color = "black", linewidth = 2, linetype = "solid"), panel.grid = element_blank(), text = element_text(family = "Arial"), legend.position="right", legend.justification = "top", legend.box.background = element_rect(color = "black", linewidth = 1),legend.title = element_blank(),axis.text = element_text(size = 40, family = "Arial", color = "black", face = "bold"), axis.title = element_text(size = 40, family = "Arial", face = "bold"), legend.text = element_text(size = 25, face = "bold", family = "Arial"), legend.key.size = unit(2, "cm"), plot.title = element_blank())+
  labs(x = "Component 1 (51.63%)", y = "Component 2 (14.81%)")
PCA_24

#Exporting plot as a PNG 
png(filename = 'Pro_PCA_24hr.png', type = 'cairo', units = 'in',
    width = 14, height = 10, pointsize = 16, res = 600)
PCA_24
dev.off()
