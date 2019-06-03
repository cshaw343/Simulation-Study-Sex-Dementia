#*******************************************************************************
# Simulation Senario A
# 25% cumulative incidence of mortality
# Exposure = male gender
# Based on US lifetables, sample is 51% F / 49% M at age 50
# Remaining LE at age 50:  33.2 years for F and 29.6 years for M
# This is a test change
#*******************************************************************************

#---- Number of simulation runs ----
runs = 100

#---- Number of observations ----
#Multiples of 1000
num_obs <- 100000

#---- Baseline age ----
age0 <- 50

#---- Assessment parameters ----
#Time between assessments (measured in years)
int_time <- 5

#Number of assessments
num_tests <- 9

#Resulting visit times
visit_times <- seq(from = 0, to = int_time*num_tests, by = int_time)

#---- Prevalance of exposure (male) ----
#Exposure = male gender
psex <- 0.49

#---- Variances and correlations ----
#Both slope and slope variance has placeholder values because of the covariance
#matrix for slope and intercept... i.e. stop trying to delete these!!!!!

cij_var0 <- 0.2   #Variance of random cognitive intercept
#Need one value for each visit, including baseline
cij_var1 <- c(0.001, #baseline measure (holding place and really doesn't matter)
#              0.001, 0.002, 0.002, 0.003, 0.004, 0.005, 0.011, 0.01775, 0.0195) #Time-dependent variance of random cognitive slope
              0.001, 0.0012, 0.00204, 0.00204, 0.00204, 0.00253161, 0.00253161, 0.00253161, 0.00253161)
cij_cov <- 0.01   #Covariance of random intercept and random slope
cij_var3 <- 1     #Variance of noise for Cij (cognitive function for person i at time j)
cij_r1 <- 0.3     #Correlation between noise terms for Cij; this may need to be adjusted
#cij_var4 <- 0.19     #Variance of measurement error of Cij

#---- Parameters for Cij ----
#Cognitive function for person i at time j
b00 <- 0      #Cognitive intercept for females
b01 <- 0      #Effect of sex on cognitive intercept
b02 <- -0.05  #Effect of age on cognitive intercept; Note: Everyone is the same age so there is no age effect (since baseline centered ages are 0 for everyone)
b03 <- 0      #Effect of U (unmeasured/underlying variable) on cognitive intercept

cij_knots <- seq(55, 90, by = 5) #Specify which ages to place knots

#Need one value for each visit time, including baseline
#First value is cognitive slope, the remaining values are changes in cognitive slopes
#These are: b10a, b10b - b10a, b10c - b10b, etc...
#ie Cognitive slope for females age 50-70, change in cognitive slope for females age 70-85, etc...
#Based on slopes_dem-cut_search.R script (results from 20190202)
cij_slopes <- c(0, 
#                -0.00475552, -0.01820177, -0.02680375, -0.03914714, 
#                -0.03312803, -0.06, -0.03056074, -0.15, -0.40)
                -0.00475552, -0.01630603, -0.00773323, -0.01184527, -0.00968730, -0.09400000,  
                0.00000000,  0.00000000,  0.00000000)

b11 <- 0      #Effect of sex on cognitive slope
b12 <- -0.005 #Effect of age on cognitive slope; Note: Everyone is the same age so there is no age effect
b13 <- -0.05  #Effect of U on cognitive slope (currently age constant)

#---- Parameters for Sij (survival for person i at time j) ----
#Effect of sex on log hazard of death; chosen using calc from euro_life_tables.R 
g1 <- log(c(1.6160, 1.8973, 2.1307, 2.2726, 2.2225, 2.0917, 1.8494, 1.6147, 
            1.3909, 1.2596)) 
g2 <- 0     #Effect of U on log hazard of death
g3 <- 0     #Effect of interaction between sex and U on log hazard of death
g4 <- 0     #Effect of cognitive slope at time j on log hazard of death
g5 <- 0     #Effect of cognitive function at time j on log hazard of death
#Effect of prevalent dementia on log hazard of death
#From Kaiser dataset for now
g6 <- log(c(rep(8, 3), 7.45, 6.29, 4.87, 4.41, 3.84, 3.13, 3.13))      

#---- Baseline hazard of death for unexposed ----
#Computed in lambda_search_euro.R script
#Based on 35x3000 = 105000 observations
lambda <- c(0.004140, 0.004944, 0.004944, 0.009450, 0.012630, 0.021630, 0.003000, 0.003000, 0.003000)
  #c(0.00414, 0.00577, 0.00824, 0.01260, 0.02105, 0.03605, 0.06316, 
            #0.10918, 0.20142) 

#---- Dementia Cut Point ----
#Need one value for each visit time, including baseline
#Based on slopes_dem-cut_search.R script (results from 20190202)
dem_cut <- -5.75
# dem_cuts <- c(-2.98629, 
#               -2.98629, -3.28948, -3.77503, -4.45, -4.975, -5.25, -6.225, 
#               -6.8225, -7.05)
