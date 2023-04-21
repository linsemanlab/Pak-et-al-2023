
##Data cleaning, normalization (log10 transformation), and scaling (pareto) for metabolomics dataset 
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

#Used for data imputation
library(mice)

#Used to clean column headers
library(janitor)

######
#Input data - raw metabolomics data
View(Metabolomics_Raw_Data)

#Using janitor to clean column header names 
Metabolomics_Raw_Data_Clean <- Metabolomics_Raw_Data %>%
  clean_names()
View(Metabolomics_Raw_Data_Clean)

#Converting zero values to NA for data imputation
Metabolomics_Raw_Data_Clean[Metabolomics_Raw_Data_Clean == 0] <- NA
View(Metabolomics_Raw_Data_Clean)

#Checking to make sure zeroes were correctly converted to NAs 
is.na(Metabolomics_Raw_Data_Clean) %>% table()
sapply(Metabolomics_Raw_Data_Clean, function(y) sum(length(which(is.na(y)))))

#Removing the columns with text prior to imputation (requires a dataframe with only numeric data)
Metabolomics_Raw_Data_Clean <- Metabolomics_Raw_Data_Clean[ -c(1:3) ]
View(Metabolomics_Raw_Data_Clean)

#Using CART based imputation to replace NA values 
Metabolomics_Raw_Data_Imputed <- mice(Metabolomics_Raw_Data_Clean, m=5, maxit = 50, method = 'cart', seed = 245435)
summary(Metabolomics_Raw_Data_Imputed)
Metabolomics_Raw_Data_Imputed <- complete(Metabolomics_Raw_Data_Imputed,1)
View(Metabolomics_Raw_Data_Imputed)

#Exporting CSV with imputed data for re-upload without having to re-do the imputation (computationally heavy - takes some time)
write.csv(Metabolomics_Raw_Data_Imputed,file="Metabolomics_Raw_Data_Imputed.csv") 

#Adding the text columns back into the new dataframe with the imputed data 
Metabolomics_Raw_Data_Imputed_Clean <- cbind(Metabolomics_Raw_Data[1:3], Metabolomics_Raw_Data_Imputed)
View(Metabolomics_Raw_Data_Imputed_Clean)

#Merging the original column headers into the new dataframe with the imputed data 
colnames(Metabolomics_Raw_Data_Imputed_Clean) <- colnames(Metabolomics_Raw_Data)
View(Metabolomics_Raw_Data_Imputed_Clean)

#Checking to make sure data imputation was sucessful 
is.na(Metabolomics_Raw_Data_Imputed_Clean) %>% table()
sapply(Metabolomics_Raw_Data_Imputed_Clean, function(y) sum(length(which(is.na(y)))))

#####
#Applying a Log transformation (log10) to the imputed data 
Met_log <- Metabolomics_Raw_Data_Imputed_Clean
Met_log[4:173] <- log10(Metabolomics_Raw_Data_Imputed_Clean[4:173])
View(Met_log)

######
#Pareto scaling (mean-centered and divided by the square root of the standard deviation of each variable) the log transformed daa 
Met_log_imputed_scaled <- Met_log
Met_log_imputed_scaled[4:173] <- pareto_scale(Met_log_imputed_scaled[4:173], centering = TRUE)
View(Met_log_imputed_scaled)

#####
#Exporting CSV with final cleaned/normalized data 
write.csv(Met_log_imputed_scaled,file="Met_log_imputed_scaled.csv") 

