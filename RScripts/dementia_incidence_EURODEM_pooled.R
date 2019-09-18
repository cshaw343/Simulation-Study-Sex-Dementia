#*******************************************************************************
# Dementia Incidence Rates 
# Source: Gender differences in the incidence of AD and vascular dementia: 
# The EURODEM Studies (Andersen 1999); Table 2
# 
# Life Table birth cohort:  
#*******************************************************************************

#---- Package Loading ----
if (!require("pacman")) 
  install.packages("pacman", repos='http://cran.us.r-project.org')

p_load("tidyverse")

#--- Dementia incidence rates ----
EURODEM_inc_rates <- tibble("Low_Age" = seq(65, 95, by = 5), 
                            "High_Age" = seq(69, 99, by = 5), 
                            "Visit_Age" = High_Age + 1, 
                            "M_Person_Years" = c(2948, 3524, 2751, 1809, 894, 
                                                 344, 344), 
                            "M_All_Dementia" = c(10, 29, 38, 45, 40, 15, 15), 
                            "M_AD" = c(3, 14, 20, 29, 24, 6, 6), 
                            "M_Vasc" = c(4, 9, 11, 10, 8, 5, 5), 
                            "F_Person_Years" = c(3404, 4254, 3777, 2729, 1496, 
                                                 837, 837), 
                            "F_All_Dementia" = c(10, 26, 71, 100, 77, 67, 67), 
                            "F_AD" = c(7, 18, 43, 80, 56, 52, 52), 
                            "F_Vasc" = c(1, 6, 14, 11, 10, 9, 9), 
                            "Total_Person_Years" = M_Person_Years + 
                              F_Person_Years, 
                            "Total_All_Dementia" = M_All_Dementia + 
                              F_All_Dementia, 
                            "Total_AD" = M_AD + F_AD, 
                            "Total_Vasc" = M_Vasc + F_Vasc, 
                            "Total_All_Dementia_1000PY" = 
                              1000*Total_All_Dementia/Total_Person_Years, 
                            "Total_AD_1000PY" = 
                              1000*Total_AD/Total_Person_Years, 
                            "Female_All_Dementia_1000PY" = 
                              1000*F_All_Dementia/F_Person_Years, 
                            "Male_All_Dementia_1000PY" = 
                              1000*M_All_Dementia/M_Person_Years, 
                            "Female_AD_1000PY" = 
                              1000*F_AD/F_Person_Years, 
                            "Male_AD_1000PY" = 
                              1000*M_AD/M_Person_Years) %>%
  mutate("Female_All_Dementia_IR" = F_All_Dementia/F_Person_Years, 
         "Male_All_Dementia_IR" = M_All_Dementia/M_Person_Years, 
         "All_Dementia_IRR_F:M" = Female_All_Dementia_IR/Male_All_Dementia_IR, 
         "Female_AD_IR" = F_AD/F_Person_Years, 
         "Male_AD_IR" = M_AD/M_Person_Years, 
         "AD_IRR_F:M" = Female_AD_IR/Male_AD_IR)
  
  