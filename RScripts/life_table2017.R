#***************************************************************
# 2017 Life Table
# Source: National Vital Statistics Reports, Vol. 68, No. 7, 
# June 24, 2019 (pg 49-50)
# Birth cohort: 1919 - 1921
# **Used data with race aggregated**
#***************************************************************

#---- Package Loading and Options ----
if (!require("pacman")) 
  install.packages("pacman", repos='http://cran.us.r-project.org')

p_load("tidyverse", "here")

#---- Source Files ----
source(here("RScripts", "misc_custom_functions.R"))

#---- Hazard Function ----
haz <- function(age, logprobs){
  HZ <- vector(length = length(age))
  for(i in 2:length(HZ)){
    HZ[i] = -(logprobs[i] - logprobs[i - 1])/(age[i] - age[i - 1])
  }
  return(HZ)
}

#---- Life Table Data ----
#"Survivors" represents number surviving out of 100,000 born alive
ages <- seq(from = 0, to = 100, by = 5)

life <- tibble("Birth cohort" = "1919-1921",
               "Age" = ages, 
               "Survivors" = c(100000, 83389, 88129, 87144, 85441, 83146, 80642, 
                               77961, 75114, 72036, 68429, 63947, 58079, 50560, 
                               41090, 29729, 18298, 8683, 2941, 646, 67), 
               "Prob" = Survivors/100000, 
               "logProb" = log(Prob), 
               "CP" = cond_prob(Survivors), 
               "Haz" = haz(age = Age, logprobs = logProb), 
               "Country" = "US", 
               "Sex" = "both") %>% 
  write_csv("Data/US_cohort_table_MF_calcs.csv")

male_life <- tibble("Birth cohort" = "1919-1921",
                    "Age" = ages, 
                    "Survivors" = c(100000, 88505, 87184, 86156, 84440, 82252, 
                                    79890, 77514, 74432, 71244, 67553, 62965, 
                                    56917, 49218, 39668, 28316, 17128, 7920, 
                                    2527, 556, 62), 
                    "Prob" = Survivors/100000,
                    "logProb" = log(Prob),
                    "CP" = cond_prob(Survivors), 
                    "Haz" = haz(age = Age, logprobs = logProb), 
                    "Country" = "US", 
                    "Sex" = "Male") %>% 
  write_csv("Data/US_cohort_table_M_calcs.csv")

female_life <- tibble("Birth cohort" = "1919-1921",
                      "Age" = ages, 
                      "Survivors" = c(100000, 90380, 89186, 88247, 86556, 
                                      84135, 81463, 78713, 75907, 72954, 69452, 
                                      65099, 59438, 52126, 42741, 31344, 19613, 
                                      9515, 3314, 728, 72), 
                      "Prob" = Survivors/100000,
                      "logProb" = log(Prob),
                      "CP" = cond_prob(Survivors), 
                      "Haz" = haz(age = Age, logprobs = logProb), 
                      "Country" = "US", 
                      "Sex" = "Female") %>% 
  write_csv("Data/US_cohort_table_F_calcs.csv")

# #---- Hazard Plots ----
# #Creating plot data
# female_hazards <- female_life %>% dplyr::select(c("Age", "Haz")) %>% 
#   melt(., id.vars = "Age")
# 
# male_hazards <- male_life %>% dplyr::select(c("Age", "Haz")) %>% 
#   melt(., id.vars = "Age")
# 
# haz_ratios <- Hratio %>% cbind(ages) %>% 
#   mutate("Age" = ages) %>% dplyr::select(-ages) %>% 
#   melt(., id.vars = "Age")
# 
# haz_plot_data <- rbind(female_hazards, male_hazards, haz_ratios)
# 
# #Creating plot
# hazard_plot<- ggplot(haz_plot_data, aes(Age, value), color = variable) + 
#   geom_line(data = subset(haz_plot_data, variable == "Haz"), 
#             aes(color = variable), size = 1) + 
#   geom_line(data = subset(haz_plot_data, variable == "Haz"), 
#             aes(color = variable), size = 1, alpha = 0.6) + 
#   geom_line(data = subset(haz_plot_data, variable == "ratio"), 
#             aes(color = variable), size = 1) + ylim(0, 1.25) + 
#   labs(y = "Hazard", 
#        x = "Age", 
#        color = "Hazard") + 
#   theme_minimal()
# 
# #Saving plot output
# #ggsave(filename = "hazard_plot_2014life.jpeg", width = 10, height = 7, 
# #       plot = hazard_plot)


