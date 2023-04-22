
##T-Tests, Fold-Changes and volcano plot for metabolomics dataset 
##Paper Title: Non-canonical Metabolic and Molecular Effects of Calorie Restriction Are Revealed by Varying Temporal Conditions

######

#Used to generate tables in R for viewing
library(knitr)
install.packages("kableExtra")
library(kableExtra)

#Function used to generate tables using knitr and kableExtra 
knitr_table <- function(x) {
  x %>% 
    knitr::kable(format = "html", digits = Inf, 
                 format.args = list(big.mark = ",")) %>%
    kableExtra::kable_styling(font_size = 15)
}

#Used to make sure labels don't overlap in volcano plot 
library(ggrepel)
library(dplyr)
library(tidyr)
library(purrr)

#Used to format volcano plot 
install.packages("extrafont")
library(extrafont)
font_import()
fonts()
library(scales)
library(ggplot2)

######
#Fold-change analysis

#Calculating ratios between group means - CR-A/AL-A
met_df_1 <- Metabolomics_Raw_Data %>% 
  filter(Timepoint == "A") %>%
  group_by(Group) %>%
  summarise_if(is.numeric, mean)
met_result_1 <- as.data.frame(met_df_1[2,-1]/met_df_1[1,-1])
met_result_1$Group <- "CR-A:AL-A"

#Calculating ratios between group means - CR-B/AL-B
met_df_2 <- Metabolomics_Raw_Data %>% 
  filter(Timepoint == "B") %>%
  group_by(Group) %>%
  summarise_if(is.numeric, mean)
met_result_2 <- as.data.frame(met_df_2[2,-1]/met_df_2[1,-1])
met_result_2$Group <- "CR-B:AL-B"

#Calculating ratios between group means - CR-C/AL-C
met_df_3 <- Metabolomics_Raw_Data %>% 
  filter(Timepoint == "C") %>%
  group_by(Group) %>%
  summarise_if(is.numeric, mean)
met_result_3 <- as.data.frame(met_df_3[2,-1]/met_df_3[1,-1])
met_result_3$Group <- "CR-C:AL-C"

#Calculating ratios between group means - CR-D/AL-D
met_df_4 <- Metabolomics_Raw_Data %>% 
  filter(Timepoint == "D") %>%
  group_by(Group) %>%
  summarise_if(is.numeric, mean)
met_result_4 <- as.data.frame(met_df_4[2,-1]/met_df_4[1,-1])
met_result_4$Group <- "CR-D:AL-D"

#Merging individual results for each group comparison into one dataframe 
met_results <- rbind(met_result_1, met_result_2, met_result_3, met_result_4)
met_results_final <- met_results %>%
  dplyr::select(Group, everything())
View(met_results_final)

#Creating new dataframe - takes the logarithm base 2
log2_met <- met_results_final
log2_met[2:171] <- log2(log2_met[2:171])
View(log2_met)

#Turns the new dataframe into a long form dataframe 
log2_met_long <- log2_met %>%
  pivot_longer(-Group, names_to = "Metabolite", values_to = "Log2_FC")
View(log2_met_long)

######
#T-tests

#Filters on timepoint A for CR-A/AL-A comparison
met_df_5 <- Metabolomics_Raw_Data %>%
  filter(Timepoint == "A")
View(met_df_5)

#Runs t-tests between CR-A and AL-A 
met_res_A <- met_df_5 %>%
  select_if(is.numeric)%>%
  map_df(~ broom::tidy(t.test(.~Group, data=met_df_5, var.equal=FALSE)), .id = 'var')
met_res_A$Group <- "CR-A:AL-A"

#Filters on timepoint B for CR-B/AL-B comparison
met_df_6 <- Metabolomics_Raw_Data %>%
  filter(Timepoint == "B")

#Runs t-tests between CR-B and AL-B
met_res_B <- met_df_6 %>%
  select_if(is.numeric)%>%
  map_df(~ broom::tidy(t.test(.~Group, data=met_df_6, var.equal=FALSE)), .id = 'var')
met_res_B$Group <- "CR-B:AL-B"

#Filters on timepoint C for CR-C/AL-C comparison
met_df_7 <- Metabolomics_Raw_Data %>%
  filter(Timepoint == "C")

#Runs t-tests between CR-C and AL-C
met_res_C <- met_df_7 %>%
  select_if(is.numeric)%>%
  map_df(~ broom::tidy(t.test(.~Group, data=met_df_7, var.equal=FALSE)), .id = 'var')
met_res_C$Group <- "CR-C:AL-C"

#Filters on timepoint D for CR-D/AL-D comparison
met_df_8 <- Metabolomics_Raw_Data %>%
  filter(Timepoint == "D")

#Runs t-tests between CR-D and AL-D
met_res_D <- met_df_8 %>%
  select_if(is.numeric)%>%
  map_df(~ broom::tidy(t.test(.~Group, data=met_df_8, var.equal=FALSE)), .id = 'var')
met_res_D$Group <- "CR-D:AL-D"

#Merging individual results for each group comparison into one dataframe 
met_results_t_test <- rbind(met_res_A, met_res_B, met_res_C, met_res_D)
View(met_results_t_test)
colnames(met_results_t_test)[which(names(met_results_t_test) == "var")] <- "Metabolite"

#Merges t-test data and fold-change data 
Met_log2FC_ttest <- merge(log2_met_long,met_results_t_test,by=c("Group","Metabolite"))
View(Met_log2FC_ttest)

#View new data frame as knitr table
Met_log2FC_ttest%>% 
  knitr_table()

######
#Creating a volcano plot to display t-test/fold-change data 

#Filtering results by values that increased, decreased or were unchanged - Creates an "Abundance" variable for use in plot 
Met_log2FC_ttest <- Met_log2FC_ttest %>% 
  mutate(Abundance = case_when(Log2_FC >= log(2) & p.value <= 0.05 ~ "Increased",
                     Log2_FC <= -log(2) & p.value <= 0.05 ~ "Decreased",
                     TRUE ~ "Unchanged"))

#knitr table showing Abundance 
Met_log2FC_ttest %>% dplyr::count(Abundance) %>% knitr_table()

#Filtering results by p-value - Creates a "Significance" variable for use in plot 
Met_log2FC_ttest <- Met_log2FC_ttest %>% 
  mutate(
    Significance = case_when(
      abs(Log2_FC) >= log(2) & p.value <= 0.05 & p.value > 0.01 ~ "p-value < 0.05", 
      abs(Log2_FC) >= log(2) & p.value <= 0.01 & p.value > 0.001 ~ "p-value < 0.01",
      abs(Log2_FC) >= log(2) & p.value <= 0.001 ~ "p-value < 0.001", 
      TRUE ~ "Unchanged"))

#knitr table showing Abundance and Significance 
Met_log2FC_ttest %>% dplyr::count(Abundance, Significance) %>% knitr_table()

#Creating new dataframe with all significant metabolites (<0.05 - increased) 
met_significant_genes_pos <- Met_log2FC_ttest %>% 
  filter(Significance == 'p-value < 0.05' & Log2_FC < 0)
View(met_significant_genes_pos)

#Creating new dataframe with all significant metabolites (<0.05 - decreased) 
met_significant_genes_neg <- Met_log2FC_ttest %>% 
  filter(Significance == 'p-value < 0.05' & Log2_FC < 0)
View(met_significant_genes_neg)

#Creating lines for plot 
lines <- data.frame(
  Group = unique(Met_log2FC_ttest$Group),
  line = c(0,0,0,0))

#Creating labels for plot 
labels <- data.frame(
  Group = unique(Met_log2FC_ttest$Group),
  label = c("8hr Fast","12hr Fast","16hr Fast","24hr Fast"))

#Renaming metabolites with long names (don't fit well on the plot)
met_significant_genes_pos[4, 2] <- "acyl-C4"
met_significant_genes_neg[3, 2] <- "a-Linolenic acid"
met_significant_genes_neg[4, 2] <- "acyl-C18:2"
met_significant_genes_neg[7, 2] <- "Arachidonic acid"
met_significant_genes_neg[8, 2] <- "Dodecanoic acid"
met_significant_genes_neg[12, 2] <- "acyl-C16"
met_significant_genes_neg[14, 2] <- "acyl-C12"
met_significant_genes_neg[15, 2] <- "acyl-C18:1"
met_significant_genes_neg[12, 2] <- "O-dodecanoyl-carnitine (acyl-C12)"
met_significant_genes_neg[16, 2] <- "acyl-C14"
met_significant_genes_neg[17, 2] <- "Octadecenoic acid"
met_significant_genes_neg[18, 2] <- "Tetradecanoic acid"
met_significant_genes_neg[23, 2] <- "acetyl-carnitine"
met_significant_genes_neg[26, 2] <- "Arachidonic acid"
met_significant_genes_neg[32, 2] <- "acyl-C18:1"
met_significant_genes_neg[33, 2] <- "Octadecenoic acid"
met_significant_genes_neg[34, 2] <- "Pentanoate"
met_significant_genes_neg[20, 2] <- "7,10,13,16,19-Docosapentaenoic acid"
met_significant_genes_neg[13, 2] <- "Linoleic acid"
met_significant_genes_neg[21, 2] <- "8,11,14-Eicosatrienoic acid"

#Volcano plot 
met_VP <- ggplot(data=Met_log2FC_ttest, aes(x=Log2_FC, y=-log10(p.value))) +
  geom_point(aes(color = Abundance), size = 8) +
  facet_wrap(~Group, scales='free_y') + 
  geom_hline(data = lines, linewidth = 2.5, aes(yintercept = line))+
  geom_vline(data = lines, linewidth = 2.5, aes(xintercept = line))+
  theme_minimal()+
  xlab(expression("log"[2]*"Fold-Change")) + 
  ylab(expression("-log"[10]*"p-value")) + 
  ylim(0, 7) + 
  scale_x_continuous(
  labels = label_number(accuracy = 1))+
  geom_label(data = labels, aes(label = label), x = 4.3, y = 0.4, label.size = NA, size=20, fontface = "bold", family = "Arial")+
  scale_color_manual(values=c("#6E91CB", "#B12625", "#BFBEBE")) +
  theme(plot.margin = unit(c(1,1,1,1), "cm"),axis.title.y = element_text(margin = margin(t = 0, r = 25, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 25, r = 0, b = 0, l = 0)), strip.text.x = element_blank())+
  theme(axis.ticks = element_line(linewidth = 3), axis.ticks.length=unit(1, "cm"), text = element_text(family = "Arial"), 
        panel.spacing = unit(4, "lines"), legend.position="none",axis.text = element_text(size = 60, color = "black", face="bold"), 
        axis.title = element_text(size = 100, color = "black", face="bold"), plot.title = element_blank())+
  geom_label_repel(data = met_significant_genes_pos, max.overlaps = Inf, min.segment.length = 0.5, 
                   point.padding = 0.75, 
                   box.padding = 0.75, 
                   seed = 42, 
                   fill = NA, 
                   mapping = aes(Log2_FC, -log10(p.value), label = Metabolite),
                   size = 8.5,
                   family = "Arial", 
                   color = "#B12625", 
                   xlim = c(0.25, NA),
                   ylim = c(0,7),
                   point.size = 1, 
                   label.size = NA,
                   label.padding = unit(0.5, "lines")) +
  geom_label_repel(data = met_significant_genes_neg, max.overlaps = Inf, min.segment.length = 0.5, 
                   point.padding = 0.75, 
                   box.padding = 0.75, 
                   seed = 42, 
                   fill = NA, 
                   family = "Arial", 
                   mapping = aes(Log2_FC, -log10(p.value), label = Metabolite),
                   size = 8.5,
                   color = "#6E91CB", 
                   xlim = c(NA, -0.25),
                   ylim = c(0,7),
                   point.size = 1, 
                   label.size = NA,
                   label.padding = unit(0.5, "lines"))
met_VP

#Exporting volcano plot as a PNG 
png(filename = 'Met_Volcano_Plots_<0.05.png', type = 'cairo', units = 'in',
    width = 56, height = 32.5, pointsize = 16, res = 300)
met_VP
dev.off()

######
#Source data for metabolomics volcano plot 
write.csv(Met_log2FC_ttest,file="Met_Volcano_Plot_Source_Data.csv") 
