#***************************************************************
# Dementia Incidence Rates (2000 - 2013) per 1000 person years
# Source: 
# 
# Life Table birth cohort:  1919 - 1921
#***************************************************************

#The original source has incidence rates by race, age, and sex.
#We will start with calibrating to Whites as a first attempt

dem_rates_whites <- 
  data_frame("Low Age" = c(64, 70, 75, 80, 85, 90), 
             "High Age" = c(69, 74, 79, 84, 89, 1000), 
             "Rate" = c(3.45, 8.44, 17.35, 34.51, 61.35, 99.26))
