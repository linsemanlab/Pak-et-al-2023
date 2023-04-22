
##T-Tests, Fold-Changes and volcano plot for proteomics dataset 
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
pro_df_1 <- Proteomics_Raw_Data_Spectral_Count_Cleaner %>% 
  filter(Timepoint == "A") %>%
  group_by(Group) %>%
  summarise_if(is.numeric, mean)
pro_result_1 <- as.data.frame(pro_df_1[2,-1]/pro_df_1[1,-1])
pro_result_1$Group <- "CR-A:AL-A"

#Calculating ratios between group means - CR-B/AL-B
pro_df_2 <- Proteomics_Raw_Data_Spectral_Count_Cleaner %>% 
  filter(Timepoint == "B") %>%
  group_by(Group) %>%
  summarise_if(is.numeric, mean)
pro_result_2 <- as.data.frame(pro_df_2[2,-1]/pro_df_2[1,-1])
pro_result_2$Group <- "CR-B:AL-B"

#Calculating ratios between group means - CR-C/AL-C
pro_df_3 <- Proteomics_Raw_Data_Spectral_Count_Cleaner %>% 
  filter(Timepoint == "C") %>%
  group_by(Group) %>%
  summarise_if(is.numeric, mean)
pro_result_3 <- as.data.frame(pro_df_3[2,-1]/pro_df_3[1,-1])
pro_result_3$Group <- "CR-C:AL-C"

#Calculating ratios between group means - CR-D/AL-D
pro_df_4 <- Proteomics_Raw_Data_Spectral_Count_Cleaner %>% 
  filter(Timepoint == "D") %>%
  group_by(Group) %>%
  summarise_if(is.numeric, mean)
pro_result_4 <- as.data.frame(pro_df_4[2,-1]/pro_df_4[1,-1])
pro_result_4$Group <- "CR-D:AL-D"

#Merging individual results for each group comparison into one dataframe 
pro_results <- rbind(pro_result_1, pro_result_2, pro_result_3, pro_result_4)
pro_results_final <- pro_results %>%
  dplyr::select(Group, everything())
View(pro_results_final)

#Creating new dataframe - takes the logarithm base 2
log2_pro <- pro_results_final
log2_pro[2:294] <- log2(log2_pro[2:294])
View(log2_pro)

#Turns the new dataframe into a long form dataframe 
log2_pro_long <- log2_pro %>%
  pivot_longer(-Group, names_to = "Protein", values_to = "Log2_FC")
View(log2_pro_long)

######
#T-tests

#Filters on timepoint A for CR-A/AL-A comparison
pro_df_5 <- Proteomics_Raw_Data_Spectral_Count_Cleaner %>%
  filter(Timepoint == "A")

#Runs t-tests between CR-A and AL-A 
pro_res_A <- pro_df_5 %>%
  select_if(is.numeric)%>%
  map_df(~ broom::tidy(t.test(.~Group, data=pro_df_5, var.equal=FALSE)), .id = 'var')
pro_res_A$Group <- "CR-A:AL-A"

#Filters on timepoint B for CR-B/AL-B comparison
pro_df_6 <- Proteomics_Raw_Data_Spectral_Count_Cleaner %>%
  filter(Timepoint == "B")

#Runs t-tests between CR-B and AL-B
pro_res_B <- pro_df_6 %>%
  select_if(is.numeric)%>%
  map_df(~ broom::tidy(t.test(.~Group, data=pro_df_6, var.equal=FALSE)), .id = 'var')
pro_res_B$Group <- "CR-B:AL-B"

#Filters on timepoint C for CR-C/AL-C comparison
pro_df_7 <- Proteomics_Raw_Data_Spectral_Count_Cleaner %>%
  filter(Timepoint == "C")

#Runs t-tests between CR-C and AL-C
pro_res_C <- pro_df_7 %>%
  select_if(is.numeric)%>%
  map_df(~ broom::tidy(t.test(.~Group, data=pro_df_7, var.equal=FALSE)), .id = 'var')
pro_res_C$Group <- "CR-C:AL-C"

#Filters on timepoint D for CR-D/AL-D comparison
pro_df_8 <- Proteomics_Raw_Data_Spectral_Count_Cleaner %>%
  filter(Timepoint == "D")

#Runs t-tests between CR-D and AL-D
pro_res_D <- pro_df_8 %>%
  select_if(is.numeric)%>%
  map_df(~ broom::tidy(t.test(.~Group, data=pro_df_8, var.equal=FALSE)), .id = 'var')
pro_res_D$Group <- "CR-D:AL-D"

#Merging individual results for each group comparison into one dataframe 
pro_results_t_test <- rbind(pro_res_A, pro_res_B, pro_res_C, pro_res_D)
View(pro_results_t_test)
colnames(pro_results_t_test)[which(names(pro_results_t_test) == "var")] <- "Protein"

#Merges t-test data and fold-change data 
Pro_log2FC_ttest <- merge(log2_pro_long,pro_results_t_test,by=c("Group","Protein"))
View(Pro_log2FC_ttest)

#View new data frame as knitr table
Pro_log2FC_ttest%>% 
  knitr_table()

######
#Creating a volcano plot to display t-test/fold-change data 

#Filtering results by values that increased, decreased or were unchanged - Creates an "Abundance" variable for use in plot 
Pro_log2FC_ttest <- Pro_log2FC_ttest %>% 
  mutate(
    Abundance = case_when(Log2_FC >= log(2) & p.value <= 0.05 ~ "Increased",
                          Log2_FC <= -log(2) & p.value <= 0.05 ~ "Decreased",
                          TRUE ~ "Unchanged"))

#knitr table showing Abundance 
Pro_log2FC_ttest %>% count(Abundance) %>% knitr_table()

#Filtering results by p-value - Creates a "Significance" variable for use in plot 
Pro_log2FC_ttest <- Pro_log2FC_ttest %>% 
  mutate(
    Significance = case_when(
      abs(Log2_FC) >= log(2) & p.value <= 0.05 & p.value > 0.01 ~ "p-value < 0.05", 
      abs(Log2_FC) >= log(2) & p.value <= 0.01 & p.value > 0.001 ~ "p-value < 0.01",
      abs(Log2_FC) >= log(2) & p.value <= 0.001 ~ "p-value < 0.001", 
      TRUE ~ "Unchanged"))

#knitr table showing Abundance and Significance 
Pro_log2FC_ttest %>% count(Abundance, Significance) %>% knitr_table()

#Creating new dataframe with all significant proteins (<0.05 - increased) 
pro_significant_genes_pos <- Pro_log2FC_ttest %>% 
  filter(Significance == 'p-value < 0.05' & Log2_FC > 0)
View(pro_significant_genes_pos)

#Creating new dataframe with all significant proteins (<0.05 - decreased) 
pro_significant_genes_neg <- Pro_log2FC_ttest %>% 
  filter(Significance == 'p-value < 0.05' & Log2_FC < 0)
View(pro_significant_genes_neg)

#Creating lines for plot 
lines <- data.frame(
  Group = unique(Pro_log2FC_ttest$Group),
  line = c(0,0,0,0))

#Creating labels for plot 
labels <- data.frame(
  Group = unique(Pro_log2FC_ttest$Group),
  label = c("8hr Fast","12hr Fast","16hr Fast","24hr Fast"))

#Volcano plot 
pro_VP <- ggplot(data=Pro_log2FC_ttest, aes(x=Log2_FC, y=-log10(p.value))) +
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
  geom_label(data = labels, aes(label = label), x = 2.6, y = 0.4, label.size = NA, size=20, fontface = "bold", family = "Arial")+
  scale_color_manual(values=c("#6E91CB", "#B12625", "#BFBEBE")) +
  theme(plot.margin = unit(c(1,1,1,1), "cm"),axis.title.y = element_text(margin = margin(t = 0, r = 25, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 25, r = 0, b = 0, l = 0)), strip.text.x = element_blank())+
  theme(axis.ticks = element_line(linewidth = 3), axis.ticks.length=unit(1, "cm"), text = element_text(family = "Arial"), 
        panel.spacing = unit(4, "lines"), legend.position="none",axis.text = element_text(size = 60, color = "black", face="bold"),
        axis.title = element_text(size = 100, color = "black", face="bold"), plot.title = element_blank())+
  geom_label_repel(data = pro_significant_genes_pos, max.overlaps = Inf, min.segment.length = 0.5, 
                   point.padding = 0.75, 
                   box.padding = 0.75, 
                   seed = 42, 
                   fill = NA, 
                   mapping = aes(Log2_FC, -log10(p.value), label = Protein),
                   size = 8.5,
                   family = "Arial", 
                   color = "#B12625", 
                   xlim = c(0.25, NA),
                   ylim = c(0,7),
                   point.size = 1, 
                   label.size = NA,
                   label.padding = unit(0.5, "lines")) +
  geom_label_repel(data = pro_significant_genes_neg, max.overlaps = Inf, min.segment.length = 0.5, 
                   point.padding = 0.75, 
                   box.padding = 0.75, 
                   seed = 42, 
                   fill = NA, 
                   family = "Arial", 
                   mapping = aes(Log2_FC, -log10(p.value), label = Protein),
                   size = 8.5,
                   color = "#6E91CB", 
                   xlim = c(NA, -0.25),
                   ylim = c(0,7),
                   point.size = 1, 
                   label.size = NA,
                   label.padding = unit(0.5, "lines"))
pro_VP

#Exporting volcano plot as a PNG 
png(filename = 'Pro_Volcano_Plots_<0.05.png', type = 'cairo', units = 'in',
    width = 56, height = 32.5, pointsize = 16, res = 300)
pro_VP
dev.off()

######

#Source data for proteomics volcano plot 
write.csv(Pro_log2FC_ttest,file="Pro_Volcano_Plot_Source_Data.csv") 
