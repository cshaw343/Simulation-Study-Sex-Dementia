#***************************************************************
# Simulation Senario A
# 25% cumulative incidence of mortality
# Exposure = male gender
# Based on US lifetables, sample is 51% F / 49% M at age 50
# Remaining LE at age 50:  33.2 years for F and 29.6 years for M
#***************************************************************

#---- Number of observations ----
num_obs <- 3000

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
var0 <- 0.2   #Variance of random cognitive intercept
var1 <- 0.005 #Variance of random cognitive slope
cov <- 0.01   #Covariance of random intercept and random slope
var3 <- 0.2   #Variance of noise for Cij (cognitive function for person i at time j)
r1 <- 0       #Correlation between noise terms for Cij
var4 <- 0.19  #Variance of measurement error of Cij

#---- Parameters for Cij ----
#Knots placed at ages 70 and 85
b00 <- 0      #Cognitive intercept for females
b01 <- 0      #Effect of sez on cognitive intercept
b02 <- -0.05  #Effect of age on cognitive intercept; Note: Everyone is the same age so there is no age effect
b03 <- 0      #Effect of U (unmeasured/underlying variable) on cognitive intercept
b10a <- 0     #Cognitive slope for females age 50-70
b10b <- -0.15 #Cognitive slope for females age 70-85
b10c <- -0.4  #Cognitive slope for females age 85+
b11 <- 0      #Effect of sex on cognitive slope
b12 <- -0.005 #Effect of age on cognitive slope; Note: Everyone is the same age so there is no age effect
b13 <- -0.05  #Effect of U on cognitive slope

#---- Parameters for Sij (survival for person i at time j) ----
g1 <- 0.47    #Effect of sex on log hazard of death
g2 <- 0.095   #Effect of age at time j on log hazard of death (exp(0.095) = 1.10)
g3 <- 0       #Effect of U on log hazard of death
g4 <- (0.095)*(-0.01) #Interaction effect of sex and age on log hazard of death
g5 <- 0       #Effect of cognitive slope at time j on log hazard of death
g6 <- 0       #Effect of cognitive function at time j on log hazard of death

#---- Baseline hazard of death ----
#Should vary by time j
lambda <- 0.0031

#---- Dementia Cut Point ----
#Initially set to a value
#Choose reasonable value based on find_dem_cut script
dem_cut <- -3
