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

#---- Source Files ----
source(here("RScripts", "cond_prob.R"))
source(here("RScripts", "haz.R"))

ages <- seq(from = 0, to = 100, by = 5)

#---- Life Table Data-- Race-aggregated ----
#"Survivors" represents number surviving out of 100,000 born alive
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
               "Sex" = "both") 

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
                    "Sex" = "Male") 

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
                      "Sex" = "Female") 

#---- Life Table Data-- White ----
#"Survivors" represents number surviving out of 100,000 born alive
white_life <- tibble("Birth cohort" = "1919-1921",
                     "Age" = ages, 
                     "Survivors" = c(100000, 89771, 88536, 87633, 86159, 84106, 
                                     81787, 79277, 76642, 73705, 70250, 65875, 
                                     60013, 52411, 42736, 31086, 19149, 9078, 
                                     2991, 643, 62), 
                     "Prob" = Survivors/100000, 
                     "logProb" = log(Prob), 
                     "CP" = cond_prob(Survivors), 
                     "Haz" = haz(age = Age, logprobs = logProb), 
                     "Country" = "US", 
                     "Sex" = "both") 

white_male_life <- tibble("Birth cohort" = "1919-1921",
                          "Age" = ages, 
                          "Survivors" = c(100000, 88842, 87530, 86546, 84997, 
                                          83061, 80888, 78441, 75733, 72696, 
                                          69107, 64574, 58498, 50663, 40873, 
                                          29205, 17655, 8154, 2568, 556, 61), 
                          "Prob" = Survivors/100000,
                          "logProb" = log(Prob),
                          "CP" = cond_prob(Survivors), 
                          "Haz" = haz(age = Age, logprobs = logProb), 
                          "Country" = "US", 
                          "Sex" = "Male") 

white_female_life <- tibble("Birth cohort" = "1919-1921",
                            "Age" = ages, 
                            "Survivors" = c(100000, 90721, 89564, 88712, 87281, 
                                            85163, 82740, 80206, 77624, 74871, 
                                            71547, 67323, 61704, 54299, 44638, 
                                            32777, 20492, 9909, 3372, 721, 63), 
                            "Prob" = Survivors/100000,
                            "logProb" = log(Prob),
                            "CP" = cond_prob(Survivors), 
                            "Haz" = haz(age = Age, logprobs = logProb), 
                            "Country" = "US", 
                            "Sex" = "Female") 

#---- Life Table Data-- Black ----
#"Survivors" represents number surviving out of 100,000 born alive
black_life <- tibble("Birth cohort" = "1919-1921",
                     "Age" = ages, 
                     "Survivors" = c(100000, 86174, 84690, 83180, 79641, 74973, 
                                     70492, 65865, 61244, 56442, 51422, 45803, 
                                     39418, 32738, 25585, 18011, 11376, 5794, 
                                     2317, 689, 129), 
                     "Prob" = Survivors/100000, 
                     "logProb" = log(Prob), 
                     "CP" = cond_prob(Survivors), 
                     "Haz" = haz(age = Age, logprobs = logProb), 
                     "Country" = "US", 
                     "Sex" = "both") 

black_male_life <- tibble("Birth cohort" = "1919-1921",
                          "Age" = ages, 
                          "Survivors" = c(100000, 85195, 83768, 82332, 79057, 
                                          74540, 70344, 65873, 61353, 56589, 
                                          51880, 46581, 40506, 34042, 26923, 
                                          18854, 11615, 5605, 2040, 552, 77), 
                          "Prob" = Survivors/100000,
                          "logProb" = log(Prob),
                          "CP" = cond_prob(Survivors), 
                          "Haz" = haz(age = Age, logprobs = logProb), 
                          "Country" = "US", 
                          "Sex" = "Male") 

black_female_life <- tibble("Birth cohort" = "1919-1921",
                            "Age" = ages, 
                            "Survivors" = c(100000, 87149, 85607, 83954, 80154, 
                                            75359, 70633, 65857, 61130, 56230, 
                                            50780, 44742, 37954, 31044, 24107, 
                                            17216, 11151, 5972, 2579, 818, 179), 
                            "Prob" = Survivors/100000,
                            "logProb" = log(Prob),
                            "CP" = cond_prob(Survivors), 
                            "Haz" = haz(age = Age, logprobs = logProb), 
                            "Country" = "US", 
                            "Sex" = "Female") 

#---- quick calcs ----
#How different are sex diffs in survival between whites and blacks

plot(black_female_life$Haz[12:21]/black_male_life$Haz[12:21], pch = 19, 
     col = "black", ylab = "HR (women:men)", xlab = "Index")
points(white_female_life$Haz[12:21]/white_male_life$Haz[12:21])
