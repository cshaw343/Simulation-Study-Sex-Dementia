#***************************************************************
# 2014 Life Table
# Source: National Vital Statistics Reports, Vol. 66, No. 4, 
# August 14, 2017 (pg 48-49)
# Birth cohort:  1919 - 1921
#***************************************************************

#"Survivors" represents number surviving out of 100,000 born alive
male_life <- tibble("Age" = seq(from = 50, to = 100, by = 5), 
                    "Survivors" = c(67553, 62965, 56917, 49218, 39668, 28316, 
                                    17128, 7920, 2527, 556, 62))

female_life <- tibble("Age" = seq(from = 50, to = 100, by = 5), 
                      "Survivors" = c(69452, 65099, 59438, 52126, 42741, 31344, 
                                      19613, 9515, 3314, 728, 72))
