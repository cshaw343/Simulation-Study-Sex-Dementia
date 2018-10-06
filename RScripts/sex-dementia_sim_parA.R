#*******************************************************************************
# Simulation Senario A
# 25% cumulative incidence of mortality
# Exposure = male gender
# Based on US lifetables, sample is 51% F / 49% M at age 50
# Remaining LE at age 50:  33.2 years for F and 29.6 years for M
# This is a test change
#*******************************************************************************

#---- Number of observations ----
num_obs <- 10

#---- Baseline age ----
age0 <- rep(50, num_obs)

#---- Assessment parameters ----
#Time between assessments (measured in years)
int_time <- 5

#Number of assessments
num_tests <- 10

#---- Prevalance of exposure (male) ----
#Exposure = male gender
psex <- 0.49

#---- Variances and correlations ----
cij_var0 <- 0.2   #Variance of random cognitive intercept
cij_var1 <- 0.005 #Variance of random cognitive slope
cij_cov <- 0.01   #Covariance of random intercept and random slope
cij_var3 <- 1     #Variance of noise for Cij (cognitive function for person i at time j)
cij_r1 <- 0.3     #Correlation between noise terms for Cij; this may need to be adjusted
#var4 <- 0.19     #Variance of measurement error of Cij

#---- Parameters for Cij ----
#Cognitive function for person i at time j
#Experimenting with knots placed at every age
cij_b00 <- 0      #Cognitive intercept for females
cij_b01 <- 0      #Effect of sex on cognitive intercept
cij_b02 <- -0.05  #Effect of age on cognitive intercept; Note: Everyone is the same age so there is no age effect (since baseline centered ages are 0 for everyone)
cij_b03 <- 0      #Effect of U (unmeasured/underlying variable) on cognitive intercept

#First value is cognitive slopes, the remaining values are changes in cognitive slopes
cij_slopes <- c(rep(-0.01, 3), -0.015, -0.075, -0.0825, -0.155, -0.274, -1, 0)

#This is replaced by the vector of slopes above
# b10a <- 0     #Cognitive slope for females age 50-70
# b10delta_b <- 0   #Change in cognitive slope for females age 70-85
# b10delta_c <- 0   #Change in cognitive slope for females age 85+
# b10b <- -0.15 #Cognitive slope for females age 70-85
# b10c <- -0.4  #Cognitive slope for females age 85+

cij_b11 <- 0      #Effect of sex on cognitive slope
cij_b12 <- -0.005 #Effect of age on cognitive slope; Note: Everyone is the same age so there is no age effect
cij_b13 <- -0.05  #Effect of U on cognitive slope

#---- Parameters for Fij ----
#Functional ability of person i at time j




#---- Parameters for Sij (survival for person i at time j) ----
#Effect of sex on log hazard of death; chosen using calc from life_table2014.R 
g1 <- log(c(1.081, 1.087, 1.11, 1.107, 1.087, 1.087, 1.072, 1.066, 1.083, 
            0.999, 0.948)) 
g2 <- 0       #Effect of age at time j on log hazard of death (exp(0.095) = 1.10)
g3 <- 0       #Effect of U on log hazard of death
g4 <- 0       #Interaction effect of sex and age on log hazard of death
g5 <- 0       #Effect of cognitive slope at time j on log hazard of death
g6 <- 0       #Effect of cognitive function at time j on log hazard of death

#---- Baseline hazard of death for unexposed ----
#Computed in lambda_search.R script
#Based on 35x3000 = 105000 observations
lambda <- c(0.0134, 0.0187, 0.0266, 0.0392, 0.0622, 0.0946, 0.1439, 0.2072, 
            0.2903, 0.4571)

#---- Dementia Cut Point ----
#Chosen using demcut_search.R
dem_cut <- -4.75
