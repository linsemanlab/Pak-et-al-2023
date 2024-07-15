
##Metabolomics Enrichment Analysis  
##Paper Title: Non-canonical Metabolic and Molecular Effects of Calorie Restriction Are Revealed by Varying Temporal Conditions

######
#Packages needed 
library(dplyr)
library(tidyr)
library(Cairo)
library(fgsea)
library(ggplot2)

#######
#Creating dataframes containing fold change and significance data for each group and timepoint

#CR vs. AL 
metabolite_data_CR_AL <- Met_log2FC_ttest_CR_AL %>%
  filter(Group == "CR:AL") %>%
  select(metabolite = Metabolite, score = Log2_FC)
View(metabolite_data_CR_AL)

score_list_CR_AL <- metabolite_data_CR_AL$score
names(score_list_CR_AL) <- metabolite_data_CR_AL$metabolite

#CR-A vs. AL-A 
metabolite_data_A <- Met_log2FC_ttest %>%
  filter(Group == "CR-A:AL-A") %>%
  select(metabolite = Metabolite, score = Log2_FC)
View(metabolite_data_A)

score_list_A <- metabolite_data_A$score
names(score_list_A) <- metabolite_data_A$metabolite

#CR-B vs. AL-B
metabolite_data_B <- Met_log2FC_ttest %>%
  filter(Group == "CR-B:AL-B") %>%
  select(metabolite = Metabolite, score = Log2_FC)
View(metabolite_data_B)

score_list_B <- metabolite_data_B$score
names(score_list_B) <- metabolite_data_B$metabolite

#CR-C vs. AL-C
metabolite_data_C <- Met_log2FC_ttest %>%
  filter(Group == "CR-C:AL-C") %>%
  select(metabolite = Metabolite, score = Log2_FC)
View(metabolite_data_C)

score_list_C <- metabolite_data_C$score
names(score_list_C) <- metabolite_data_C$metabolite

#CR-D vs. AL-D
metabolite_data_D <- Met_log2FC_ttest %>%
  filter(Group == "CR-D:AL-D") %>%
  select(metabolite = Metabolite, score = Log2_FC)
View(metabolite_data_D)

score_list_D <- metabolite_data_D$score
names(score_list_D) <- metabolite_data_D$metabolite

######
#Creating lists of pathway names for each dataframe created above for use in enrichment analysis 

filtered_path_names_CR_AL <- Metabolomics_Pathway_Names[Metabolomics_Pathway_Names$Compound %in% metabolite_data_CR_AL$metabolite, ]
View(filtered_path_names_CR_AL)
pathway_vector_CR_AL <- as.character(filtered_path_names_CR_AL$Pathway)

filtered_path_names_A <- Metabolomics_Pathway_Names[Metabolomics_Pathway_Names$Compound %in% metabolite_data_A$metabolite, ]
View(filtered_path_names_A)
pathway_vector_A <- as.character(filtered_path_names_A$Pathway)

filtered_path_names_B <- Metabolomics_Pathway_Names[Metabolomics_Pathway_Names$Compound %in% metabolite_data_B$metabolite, ]
View(filtered_path_names_B)
pathway_vector_B <- as.character(filtered_path_names_B$Pathway)

filtered_path_names_C <- Metabolomics_Pathway_Names[Metabolomics_Pathway_Names$Compound %in% metabolite_data_C$metabolite, ]
View(filtered_path_names_C)
pathway_vector_C <- as.character(filtered_path_names_C$Pathway)

filtered_path_names_D <- Metabolomics_Pathway_Names[Metabolomics_Pathway_Names$Compound %in% metabolite_data_D$metabolite, ]
View(filtered_path_names_D)
pathway_vector_D <- as.character(filtered_path_names_D$Pathway)

pathway_list_CR_AL <- split(filtered_path_names_CR_AL$Compound, filtered_path_names_CR_AL$Pathway)
pathway_list_A <- split(filtered_path_names_A$Compound, filtered_path_names_A$Pathway)
pathway_list_B <- split(filtered_path_names_B$Compound, filtered_path_names_B$Pathway)
pathway_list_C <- split(filtered_path_names_C$Compound, filtered_path_names_C$Pathway)
pathway_list_D <- split(filtered_path_names_D$Compound, filtered_path_names_D$Pathway)

pathway_vector_CR_AL <- names(pathway_list_CR_AL)
pathway_vector_A <- names(pathway_list_A)
pathway_vector_B <- names(pathway_list_B)
pathway_vector_C <- names(pathway_list_C)
pathway_vector_D <- names(pathway_list_D)

######
#CR vs AL - Metabolite Set Enrichment Analysis 

msea_results_CR_AL <- fgsea(pathways = pathway_list_CR_AL, 
                        stats = score_list_CR_AL, 
                        eps = 0.0)

head(msea_results_CR_AL[order(padj), ])

#Filtered subset for significant pathways
significant_pathways <- subset(msea_results_CR_AL, -log10(padj) > 1.30)

#Add asterisk to pathway labels if significant
msea_results_CR_AL$pathway_label <- msea_results_CR_AL$pathway
msea_results_CR_AL$pathway_label[msea_results_CR_AL$pathway %in% significant_pathways$pathway] <- paste0(msea_results_CR_AL$pathway[msea_results_CR_AL$pathway %in% significant_pathways$pathway], "*")

#Plotting results 
Met_Lollipop_CR_AL <- ggplot(msea_results_CR_AL, aes(x = NES, y = reorder(pathway_label, NES))) +
  geom_segment(aes(x = 0, xend = NES, y = pathway_label, yend = pathway_label)) +
  geom_point(aes(size = -log10(padj), color = ES), stroke = 0.5) +  # Add black outline
  scale_color_gradient(low = "yellow", high = "red") +  # Color gradient for ES
  labs(
    x = "Normalized Enrichment Score (NES)",
    y = "Pathway",
    color = "NES",
    size = "-log10(Adjusted p-value)")+
  theme_minimal() +
  theme(
    text = element_text(family = "Arial"), 
    plot.margin = unit(c(1,1,1,1), "cm"),
    legend.position="right", 
    legend.title = element_text(size = 15, family = "Arial", face = "bold"), 
    axis.title.y = element_text(margin = margin(t = 0, r = 25, b = 0, l = 0)),
    axis.title.x = element_text(margin = margin(t = 25, r = 0, b = 0, l = 0)), 
    axis.text = element_text(size = 20, family = "Arial", color = "black"), 
    axis.title = element_text(size = 25, family = "Arial", face = "bold"), 
    legend.text = element_text(size = 15, family = "Arial"))
Met_Lollipop_CR_AL

#Exporting plot 
png(filename = 'Met_Lollipop_CR_AL.png', type = 'cairo', units = 'in',
    width = 14, height = 10, pointsize = 16, res = 600)
Met_Lollipop_CR_AL
dev.off()

######
#CR-A vs. AL-A - Metabolite Set Enrichment Analysis 

msea_results_A <- fgsea(pathways = pathway_list_A, 
                        stats = score_list_A, 
                        eps = 0.0)

head(msea_results_A[order(padj), ])

#Filtered subset for significant pathways
significant_pathways <- subset(msea_results_A, -log10(padj) > 1.30)

#Add asterisk to pathway labels if significant
msea_results_A$pathway_label <- msea_results_A$pathway
msea_results_A$pathway_label[msea_results_A$pathway %in% significant_pathways$pathway] <- paste0(msea_results_A$pathway[msea_results_A$pathway %in% significant_pathways$pathway], "*")

#Plotting results 
Met_Lollipop_CR_A_AL_A <- ggplot(msea_results_A, aes(x = NES, y = reorder(pathway_label, NES))) +
  geom_segment(aes(x = 0, xend = NES, y = pathway_label, yend = pathway_label)) +
  geom_point(aes(size = -log10(padj), color = ES), stroke = 0.5) +  # Add black outline
  scale_color_gradient(low = "yellow", high = "red") +  # Color gradient for ES
  labs(
    x = "Normalized Enrichment Score (NES)",
    y = "Pathway",
    color = "NES",
    size = "-log10(Adjusted p-value)")+
  theme_minimal() +
  theme(
    text = element_text(family = "Arial"), 
    plot.margin = unit(c(1,1,1,1), "cm"),
    legend.position="right", 
    legend.title = element_text(size = 15, family = "Arial", face = "bold"), 
    axis.title.y = element_text(margin = margin(t = 0, r = 25, b = 0, l = 0)),
    axis.title.x = element_text(margin = margin(t = 25, r = 0, b = 0, l = 0)), 
    axis.text = element_text(size = 20, family = "Arial", color = "black"), 
    axis.title = element_text(size = 25, family = "Arial", face = "bold"), 
    legend.text = element_text(size = 15, family = "Arial"))
Met_Lollipop_CR_A_AL_A

#Exporting plot 
png(filename = 'Met_Lollipop_CR_A_AL_A.png', type = 'cairo', units = 'in',
    width = 14, height = 10, pointsize = 16, res = 600)
Met_Lollipop_CR_A_AL_A
dev.off()

######
#CR-B vs. AL-B - Metabolite Set Enrichment Analysis 

msea_results_B <- fgsea(pathways = pathway_list_B, 
                        stats = score_list_B, 
                        eps = 0.0)

head(msea_results_B[order(pval), ])

#Filtered subset for significant pathways
significant_pathways <- subset(msea_results_B, -log10(padj) > 1.30)

#Add asterisk to pathway labels if significant
msea_results_B$pathway_label <- msea_results_B$pathway
msea_results_B$pathway_label[msea_results_B$pathway %in% significant_pathways$pathway] <- paste0(msea_results_B$pathway[msea_results_B$pathway %in% significant_pathways$pathway], "*")

#Plotting results 
Met_Lollipop_CR_B_AL_B <- ggplot(msea_results_B, aes(x = NES, y = reorder(pathway_label, NES))) +
  geom_segment(aes(x = 0, xend = NES, y = pathway_label, yend = pathway_label)) +
  geom_point(aes(size = -log10(padj), color = ES), stroke = 0.5) +  # Add black outline
  scale_color_gradient(low = "yellow", high = "red") +  # Color gradient for ES
  labs(
    x = "Normalized Enrichment Score (NES)",
    y = "Pathway",
    color = "NES",
    size = "-log10(Adjusted p-value)")+
  theme_minimal() +
  theme(
    text = element_text(family = "Arial"), 
    plot.margin = unit(c(1,1,1,1), "cm"),
    legend.position="right", 
    legend.title = element_text(size = 15, family = "Arial", face = "bold"), 
    axis.title.y = element_text(margin = margin(t = 0, r = 25, b = 0, l = 0)),
    axis.title.x = element_text(margin = margin(t = 25, r = 0, b = 0, l = 0)), 
    axis.text = element_text(size = 20, family = "Arial", color = "black"), 
    axis.title = element_text(size = 25, family = "Arial", face = "bold"), 
    legend.text = element_text(size = 15, family = "Arial"))
Met_Lollipop_CR_B_AL_B

#Exporting plot 
png(filename = 'Met_Lollipop_CR_B_AL_B.png', type = 'cairo', units = 'in',
    width = 14, height = 10, pointsize = 16, res = 600)
Met_Lollipop_CR_B_AL_B
dev.off()

######
#CR-C vs. AL-C - Metabolite Set Enrichment Analysis 

msea_results_C <- fgsea(pathways = pathway_list_C, 
                        stats = score_list_C, 
                        eps = 0.0)

head(msea_results_C[order(pval), ])

#Filtered subset for significant pathways
significant_pathways <- subset(msea_results_C, -log10(padj) > 1.30)

#Add asterisk to pathway labels if significant
msea_results_C$pathway_label <- msea_results_C$pathway
msea_results_C$pathway_label[msea_results_C$pathway %in% significant_pathways$pathway] <- paste0(msea_results_C$pathway[msea_results_C$pathway %in% significant_pathways$pathway], "*")

#Plotting results 
Met_Lollipop_CR_C_AL_C <- ggplot(msea_results_C, aes(x = NES, y = reorder(pathway_label, NES))) +
  geom_segment(aes(x = 0, xend = NES, y = pathway_label, yend = pathway_label)) +
  geom_point(aes(size = -log10(padj), color = ES), stroke = 0.5) +  # Add black outline
  scale_color_gradient(low = "yellow", high = "red") +  # Color gradient for ES
  labs(
    x = "Normalized Enrichment Score (NES)",
    y = "Pathway",
    color = "NES",
    size = "-log10(Adjusted p-value)")+
  theme_minimal() +
  theme(
    text = element_text(family = "Arial"), 
    plot.margin = unit(c(1,1,1,1), "cm"),
    legend.position="right", 
    legend.title = element_text(size = 15, family = "Arial", face = "bold"), 
    axis.title.y = element_text(margin = margin(t = 0, r = 25, b = 0, l = 0)),
    axis.title.x = element_text(margin = margin(t = 25, r = 0, b = 0, l = 0)), 
    axis.text = element_text(size = 20, family = "Arial", color = "black"), 
    axis.title = element_text(size = 25, family = "Arial", face = "bold"), 
    legend.text = element_text(size = 15, family = "Arial"))
Met_Lollipop_CR_C_AL_C

#Exporting plot 
png(filename = 'Met_Lollipop_CR_C_AL_C.png', type = 'cairo', units = 'in',
    width = 14, height = 10, pointsize = 16, res = 600)
Met_Lollipop_CR_C_AL_C
dev.off()

######
#CR-D vs. AL-D - Metabolite Set Enrichment Analysis 

msea_results_D <- fgsea(pathways = pathway_list_D, 
                        stats = score_list_D, 
                        eps = 0.0)

head(msea_results_D[order(padj), ])

#Filtered subset for significant pathways
significant_pathways <- subset(msea_results_D, -log10(padj) > 1.30)

#Add asterisk to pathway labels if significant
msea_results_D$pathway_label <- msea_results_D$pathway
msea_results_D$pathway_label[msea_results_D$pathway %in% significant_pathways$pathway] <- paste0(msea_results_D$pathway[msea_results_D$pathway %in% significant_pathways$pathway], "*")

#Plotting results 
Met_Lollipop_CR_D_AL_D <- ggplot(msea_results_D, aes(x = NES, y = reorder(pathway_label, NES))) +
  geom_segment(aes(x = 0, xend = NES, y = pathway_label, yend = pathway_label)) +
  geom_point(aes(size = -log10(padj), color = ES), stroke = 0.5) +  # Add black outline
  scale_color_gradient(low = "yellow", high = "red") +  # Color gradient for ES
  labs(
    x = "Normalized Enrichment Score (NES)",
    y = "Pathway",
    color = "NES",
    size = "-log10(Adjusted p-value)")+
  theme_minimal() +
  theme(
    text = element_text(family = "Arial"), 
    plot.margin = unit(c(1,1,1,1), "cm"),
    legend.position="right", 
    legend.title = element_text(size = 15, family = "Arial", face = "bold"), 
    axis.title.y = element_text(margin = margin(t = 0, r = 25, b = 0, l = 0)),
    axis.title.x = element_text(margin = margin(t = 25, r = 0, b = 0, l = 0)), 
    axis.text = element_text(size = 20, family = "Arial", color = "black"), 
    axis.title = element_text(size = 25, family = "Arial", face = "bold"), 
    legend.text = element_text(size = 15, family = "Arial"))
Met_Lollipop_CR_D_AL_D

#Exporting plot 
png(filename = 'Met_Lollipop_CR_D_AL_D.png', type = 'cairo', units = 'in',
    width = 14, height = 10, pointsize = 16, res = 600)
Met_Lollipop_CR_D_AL_D
dev.off()
