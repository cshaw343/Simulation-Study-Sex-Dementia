#*******************************************************************************
# Dementia Incidence Rates (2000 - 2013) per 1000 person years
# Source: Supplementary Table 2 and 4, Inequalities in dementia incidence 
# between six racial and ethnic groups over 14 years (Mayeda et al, 2016)
# 
# Life Table birth cohort:  1919 - 1921
#*******************************************************************************

#---- Package Loading ----
if (!require("pacman")) 
  install.packages("pacman", repos='http://cran.us.r-project.org')

p_load("tidyverse")

#---- Source Files ----
source("RScripts/misc_custom_functions.R")

#---- Tables ----
#The original source has incidence rates by race, age, and sex.
#We will start with calibrating to Whites as a first attempt

  #---- Dementia incidence rates ----
  dem_rates_whites <- 
    data_frame("LowAge" = c(64, 70, 75, 80, 85, 90), 
               "HighAge" = c(69, 74, 79, 84, 89, 1000), 
               "Rate" = c(3.45, 8.44, 17.35, 34.51, 61.35, 99.26))
  
  #Divide by 200 to get the probability(diagnosis)/5 years
  dem_incidence <- as.matrix(dem_rates_whites$Rate/200) %>% as.data.frame() %>%
    t()
  
  incidence_data_labels <- tibble("age" = rep("age", ncol(dem_incidence)),
                                  "ages" = seq(70, 95, by = 5)) %>% 
    unite("incidence_labels", sep = "")
  
  colnames(dem_incidence) <- incidence_data_labels$incidence_labels
  
  #---- Cumulative dementia risk ----
  cum_risk_whites <- 
    tibble("age" = seq(70, 100, by = 5), 
           "cum_risk" = c(1.3, 4.7, 10.6, 19.7, 29.8, 37, 40)) %>%
    mutate("risk" = sub_recurse(cum_risk))
  
  dem_risk <- cum_risk_whites$risk %>% t()
  
  risk_data_labels <- tibble("age" = rep("age", ncol(dem_risk)), 
                             "ages" = seq(70, 100, by = 5)) %>% 
    unite("risk_labels", sep = "")
  
  colnames(dem_risk) <- risk_data_labels$risk_labels

