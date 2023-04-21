
##Data cleaning, normalization (log10 transformation), and scaling (pareto) for proteomics dataset 
##Paper Title: Non-canonical Metabolic and Molecular Effects of Calorie Restriction Are Revealed by Varying Temporal Conditions

######
#R packages used for data cleaning: 

#General use 
library(tidyverse)
library(dplyr)
library(tidyr)
library(broom)

#Used for pareto scaling 
library(IMIFA)

######

#Input data - raw proteomics data
View(Proteomics_Raw_Data_Spectral_Count)

#Creating new dataframe with raw data 
Proteomics_Raw_Data_Spectral_Count_Clean <- Proteomics_Raw_Data_Spectral_Count
View(Proteomics_Raw_Data_Spectral_Count_Clean)

#Filters the dataframe containg raw proteomics data by removing any columns that have more than 25% of their values equal to zero
Proteomics_Raw_Data_Spectral_Count_Cleaner <- Proteomics_Raw_Data_Spectral_Count_Clean[sapply(Proteomics_Raw_Data_Spectral_Count_Clean, function(x) mean(x == 0) <= 0.25)]
View(Proteomics_Raw_Data_Spectral_Count_Cleaner)

#Converts the dataframe to all numeric variables and drops the first columns containing text 
Proteomics_Raw_Data_Spectral_Count_Cleaner_numeric <- Proteomics_Raw_Data_Spectral_Count_Cleaner[5:492]

#Filters the dataframe further by retaining only columns where the standard deviation is greater than 3 
Proteomics_Raw_Data_Spectral_Count_Cleaner_numeric <- Filter(function(x) sd(x) > 3, Proteomics_Raw_Data_Spectral_Count_Cleaner_numeric)

#Adding the text columns back into the new dataframe 
Proteomics_Raw_Data_Spectral_Count_Cleaner <- cbind(Proteomics_Raw_Data_Spectral_Count_Clean[1:4], Proteomics_Raw_Data_Spectral_Count_Cleaner_numeric)
View(Proteomics_Raw_Data_Spectral_Count_Cleaner)

#Adds 1 to all numeric values in the filtered dataframe (zero values cannot be log10 transformed - applying this to all values maintains consistency)
Proteomics_Raw_Data_Spectral_Count_Cleaner[5:297] <- Proteomics_Raw_Data_Spectral_Count_Cleaner[5:297] + 1

#####
#Applying a Log transformation (log10) to the filtered/cleaned data 
Pro_log <- Proteomics_Raw_Data_Spectral_Count_Cleaner
Pro_log[5:297] <- log10(Pro_log[5:297])
View(Pro_log)

######
#Pareto scaling (mean-centered and divided by the square root of the standard deviation of each variable) the log transformed daa 
Pro_log_scaled <- Pro_log
Pro_log_scaled[5:297] <- pareto_scale(Pro_log_scaled[5:297], centering = TRUE)
View(Pro_log_scaled)

#####
#Exporting CSV with final cleaned/normalized data 
write.csv(Pro_log_scaled,file="Pro_log_scaled.csv") 


