#*******************************************************************************
# Simulation Senario A
# 25% cumulative incidence of mortality
# Exposure = male gender
# Based on US lifetables, sample is 51% F / 49% M at age 50
# Remaining LE at age 50:  33.2 years for F and 29.6 years for M
# This is a test change
#*******************************************************************************

#---- Number of simulation runs ----
runs = 50

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
#Exposure = female gender
pfemale <- 0.51

#---- Variances and correlations ----
#Both slope and slope variance has placeholder values because of the covariance
#matrix for slope and intercept... i.e. stop trying to delete these!!!!!

cij_var0 <- 0.2   #Variance of random cognitive intercept
#Need one value for each visit, including baseline
#For male matching
#cij_var1 <- c(0.001, #Need this for the noise term for the baseline measure
#              rep(0.001, 3), rep(0.0015, 3), 0.00165, rep(0.004, 2)) #Time-dependent variance of random cognitive slope

#For total matching
cij_var1 <- c(0.001, #Need this for the noise term for the baseline measure
              rep(0.001, 3), 0.009, 0.08, 0.3, 0.950, 2.05, 5) #Time-dependent variance of random cognitive slope
cij_cov <- 0.01   #Covariance of random intercept and random slope
cij_var3 <- 1     #Variance of noise for Cij (cognitive function for person i at time j)
cij_r1 <- 0.3     #Correlation between noise terms for Cij; this may need to be adjusted
#cij_var4 <- 0.19     #Variance of measurement error of Cij

#---- Parameters for Cij ----
#Cognitive function for person i at time j
b00 <- 0      #Cognitive intercept for males
b01 <- 0      #Effect of female on cognitive intercept
b02 <- -0.05  #Effect of age on cognitive intercept; Note: Everyone is the same age so there is no age effect (since baseline centered ages are 0 for everyone)
b03 <- 0      #Effect of U (unmeasured/underlying variable) on cognitive intercept

cij_knots <- seq(55, 90, by = 5) #Specify which ages to place knots

#Need one value for each visit time, including baseline
#First value is cognitive slope, the remaining values are changes in cognitive slopes
#These are: b10a, b10b - b10a, b10c - b10b, etc...
#ie Cognitive slope for females age 50-55, change in cognitive slope for males age 55-60, etc...
#Based on iterative_dem_slope_opt.R script

#For male matching
#cij_slopes <- c(-0.0107143, -0.0397509, -0.0125, 0, 0, 0, -0.055, 0, 0)

#For total matching
cij_slopes <- c(rep(-0.005, 2), -0.003, -0.008, -0.015, -0.015, -0.015, -0.22, 
                -0.35) 
                
b11 <- 0      #Effect of female on cognitive slope
b12 <- -0.005 #Effect of age on cognitive slope; Note: Everyone is the same age so there is no age effect
b13 <- -0.05  #Effect of U on cognitive slope (currently age constant)

#---- Parameters for Sij (survival for person i at time j) ----
#Effect of sex (being female) on log hazard of death; 
#chosen using calc from life_table_calcs.R 
g1 <- log(c(0.5270786, 0.4693230, 0.4400211, 0.4499390, 0.4780844, 0.5407054, 
            0.6193134, 0.7189676, 0.7939100)) 
g2 <- 0     #Effect of U on log hazard of death
g3 <- 0     #Effect of interaction between female and U on log hazard of death
g4 <- 0     #Effect of cognitive slope at time j on log hazard of death
g5 <- 0     #Effect of cognitive function at time j on log hazard of death
#Effect of prevalent dementia on log hazard of death
#From Kaiser dataset for now
g6 <- log(c(rep(1, 10)))      

#---- Baseline hazard of death for unexposed ----
#For male matching
#lambda <- c(0.0073, 0.01175, 0.018, 0.02875, 0.0425, 0.06, 0.089, 0.1275, 
#            0.2175)

#For total matching
lambda <- c(0.0073, 0.0116, 0.018, 0.028, 0.045, 0.07, 0.1, 0.15, 0.27)

#---- Dementia Cut Point ----
#Need one value for each visit time, including baseline
#Based on slopes_dem-cut_search.R script (results from 20190202)
dem_cut <- -6.4
# dem_cuts <- c(-2.98629, 
#               -2.98629, -3.28948, -3.77503, -4.45, -4.975, -5.25, -6.225, 
#               -6.8225, -7.05)
