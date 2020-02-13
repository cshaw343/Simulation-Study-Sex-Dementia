#*******************************************************************************
# Dementia Incidence Rates from the Adult Changes in Thought (ACT) Study
# Source: Characterization of Dementia and Alzheimer’s Disease 
# in an Older Population: Updated Incidence and Life Expectancy 
# With and Without Dementia. Tom et al. Am J Public Health 2015
# 
# Life table birth cohort: US 1919-1921
# Note: We calibrate to All Dementia rates in the simulation
#*******************************************************************************

#---- Package Loading ----
if (!require("pacman")) 
  install.packages("pacman", repos='http://cran.us.r-project.org')

p_load("tidyverse")

#--- Dementia incidence rates ----
ACT_inc_rates <- tibble("Low_Age" = seq(65, 90, by = 5), 
                        "High_Age" = seq(69, 94, by = 5), 
                        "Visit_Age" = High_Age + 1, 
                        "M_Person_Years" = c(677, 2367, 2849, 2275, 1091, 
                                             346), 
                        "M_All_Dementia" = c(5, 27, 60, 112, 69, 34), 
                        "M_AD" = c(3, 13, 42, 79, 58, 29), 
                        "M_Non_AD" = c(2, 14, 18, 33, 11, 5), 
                        "F_Person_Years" = c(795, 3158, 4147, 3583, 1929, 
                                             835), 
                        "F_All_Dementia" = c(3, 25, 75, 160, 155, 90), 
                        "F_AD" = c(1, 17, 59, 132, 130, 76), 
                        "F_Non_AD" = c(2, 8, 16, 28, 25, 14), 
                        "Total_Person_Years" = M_Person_Years + 
                          F_Person_Years, 
                        "Total_All_Dementia" = M_All_Dementia + 
                          F_All_Dementia, 
                        "Total_AD" = M_AD + F_AD, 
                        "Total_Non_AD" = M_Non_AD + F_Non_AD, 
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
  mutate("Female_AD_IR" = F_AD/F_Person_Years, 
         "Male_AD_IR" = M_AD/M_Person_Years,
         "AD_IRR_F:M" = Female_AD_IR/Male_AD_IR, 
         "AD_SE_IRR" = sqrt(1/F_AD + 1/M_AD), 
         "AD_95CI_upper" = exp(log(`AD_IRR_F:M`) + 1.96*AD_SE_IRR), 
         "AD_95CI_lower" = exp(log(`AD_IRR_F:M`) - 1.96*AD_SE_IRR), 
         "All_Dem_IRR_F:M" = Female_All_Dementia_1000PY/
           Male_All_Dementia_1000PY, 
         "All_Dem_SE_IRR" = sqrt(1/F_All_Dementia + 1/M_All_Dementia), 
         "All_Dem_95CI_upper" = exp(log(`All_Dem_IRR_F:M`) + 
                                      1.96*All_Dem_SE_IRR), 
         "All_Dem_95CI_lower" = exp(log(`All_Dem_IRR_F:M`) - 
                                      1.96*All_Dem_SE_IRR))

# #---- Plots ----
# ACT_plot_data <- 
#   tibble("Age" = seq(69, 94, by = 5), 
#          "Men" = ACT_inc_rates$Male_All_Dementia_1000PY, 
#          "Women" = ACT_inc_rates$Female_All_Dementia_1000PY) %>%
#   gather(key = "Sex/Gender", value = "value", c("Men", "Women")) %>% 
#   mutate_at("Sex/Gender", as.factor)
# ACT_plot_data$`Sex/Gender` <- 
#   fct_relevel(ACT_plot_data$`Sex/Gender`, "Men", after = 1)
# 
# ggplot(aes(Age, value), data = ACT_plot_data) + 
#   geom_point(aes(colour = `Sex/Gender`)) + theme_minimal() + 
#   geom_line(aes(color = `Sex/Gender`)) + 
#   theme(axis.text = element_text(size = 12)) +
#   #ggtitle("Dementia Incidence Rates by Sex from EURODEM Pooled Analysis ") + 
#   ylab("Dementia Incidence Rate per 1000 Person Years")
