#***************************************************************
# 2014 Life Table
# Source: National Vital Statistics Reports, Vol. 66, No. 4, 
# August 14, 2017 (pg 48-49)
# Birth cohort:  1919 - 1921
#***************************************************************

#---- Package Loading and Options ----
if (!require("pacman")) 
  install.packages("pacman", repos='http://cran.us.r-project.org')

p_load("tidyverse")

#---- Conditional Probabilities Function ----
cond_prob <- function(x){
  probs <- vector(length = length(x))
  for(i in 2:length(probs)){
    probs[i] = x[i]/x[i - 1]
  }
  probs[1] <- 1
  return(probs)
}

#---- Life Table Data ----
#"Survivors" represents number surviving out of 100,000 born alive
life <- tibble("Age" = seq(from = 50, to = 100, by = 5), 
               "Survivors" = c(68429, 63947, 58079, 50560, 41090, 29729, 18298, 
                               8683, 2941, 646, 67)) %>% 
  mutate("CP" = cond_prob(Survivors))

male_life <- tibble("MAge" = seq(from = 50, to = 100, by = 5), 
                    "MSurvivors" = c(67553, 62965, 56917, 49218, 39668, 28316, 
                                    17128, 7920, 2527, 556, 62)) %>% 
  mutate("CP" = cond_prob(MSurvivors))

female_life <- tibble("FAge" = seq(from = 50, to = 100, by = 5), 
                      "FSurvivors" = c(69452, 65099, 59438, 52126, 42741, 31344, 
                                      19613, 9515, 3314, 728, 72)) %>% 
  mutate("CP" = cond_prob(FSurvivors))

#---- Hazard Function Calculations ----

