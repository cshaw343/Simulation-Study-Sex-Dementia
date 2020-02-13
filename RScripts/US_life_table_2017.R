#***************************************************************
# 2017 US Life Table
# Source: National Vital Statistics Reports, Vol. 68, No. 7, 
# June 24, 2019 (pg 49-52)
# Birth cohort: 1919 - 1921
# 
# Note: We calibrate to data for Whites in the simulation
#***************************************************************

#---- Package Loading and Options ----
if (!require("pacman")) 
  install.packages("pacman", repos='http://cran.us.r-project.org')

p_load("tidyverse", "here")

#---- Conditional probabilities function ----
cond_prob <- function(x){
  probs <- vector(length = length(x))
  for(i in 2:length(probs)){
    probs[i] = x[i]/x[i - 1]
  }
  probs[1] <- 1
  return(probs)
}

#---- Hazard Function ----
haz <- function(age, logprobs){
  HZ <- vector(length = length(age))
  for(i in 2:length(HZ)){
    HZ[i] = -(logprobs[i] - logprobs[i - 1])/(age[i] - age[i - 1])
  }
  return(HZ)
}

#---- Life Table Data-- Race-aggregated ----
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

#---- Life Table Data-- White ----
#"Survivors" represents number surviving out of 100,000 born alive
ages <- seq(from = 0, to = 100, by = 5)

white_life <- tibble("Birth cohort" = "1919-1921",
               "Age" = ages, 
               "Survivors" = c(100000, 89771, 88536, 87633, 86159, 84106, 
                               81787, 79277, 76642, 73705, 70250, 65875, 60013, 
                               52411, 42736, 31086, 19149, 9078, 2991, 643, 62), 
               "Prob" = Survivors/100000, 
               "logProb" = log(Prob), 
               "CP" = cond_prob(Survivors), 
               "Haz" = haz(age = Age, logprobs = logProb), 
               "Country" = "US", 
               "Sex" = "both") %>% 
  write_csv("Data/US_cohort_table_white_MF_calcs.csv")

white_male_life <- tibble("Birth cohort" = "1919-1921",
                    "Age" = ages, 
                    "Survivors" = c(100000, 88842, 87530, 86546, 84997, 83061, 
                                    80888, 78441, 75733, 72696, 69107, 64574, 
                                    58498, 50663, 40873, 29205, 17655, 8154, 
                                    2568, 556, 61), 
                    "Prob" = Survivors/100000,
                    "logProb" = log(Prob),
                    "CP" = cond_prob(Survivors), 
                    "Haz" = haz(age = Age, logprobs = logProb), 
                    "Country" = "US", 
                    "Sex" = "Male") %>% 
  write_csv("Data/US_cohort_table_white_M_calcs.csv")

white_female_life <- tibble("Birth cohort" = "1919-1921",
                      "Age" = ages, 
                      "Survivors" = c(100000, 90721, 89564, 88712, 87281, 
                                      85163, 82740, 80206, 77624, 74871, 71547, 
                                      67323, 61704, 54299, 44638, 32777, 20492, 
                                      9909, 3372, 721, 63), 
                      "Prob" = Survivors/100000,
                      "logProb" = log(Prob),
                      "CP" = cond_prob(Survivors), 
                      "Haz" = haz(age = Age, logprobs = logProb), 
                      "Country" = "US", 
                      "Sex" = "Female") %>% 
  write_csv("Data/US_cohort_table_white_F_calcs.csv")


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

# #---- Sanity check ----
# male_life_check = male_life %>% 
#   mutate("cp50" = 
#            Survivors/as.numeric(male_life[which(male_life$Age == 50), 
#                                           "Survivors"]), 
#          "log_cp50" = log(cp50), 
#          "Haz_check" = haz(age = Age, logprobs = log_cp50)) 
# 
# female_life_check = female_life %>% 
#   mutate("cp50" = 
#            Survivors/as.numeric(female_life[which(female_life$Age == 50), 
#                                             "Survivors"]))
