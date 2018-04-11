#***************************************************************
# Simulation Senario A:  No anticipated bias
# 25% cumulative incidence of mortality
# Exposure = male gender
# Based on US lifetables, sample is 51% F / 49% M at age 50
# Remaining LE at age 50:  33.2 years for F and 29.6 years for M
#***************************************************************

#---- Time between assessments ----
#Time is measured in years
int_time <- 5

#---- Prevalance of exposure ----

pexp <- 0.49

#---- Variances and correlations ----
var0 <- 0.2   #Variance of random cognitive intercept
var1 <- 0.005 #Variance of random cognitive slope
cov <- 0.01   #Covariance of random intercept and random slope
var3 <- 0.2  #Variance of noise for Cij (cognitive function for person i at time j)
r1 <- 0.40    #Correlation between noise terms for Cij
var4 <- 0.19  #Variance of measurement error of Cij

#---- Parameters for Cij (cognitive function for person i at time j) ----
b00 <- 0      #Cognitive intercept for unexposed
b01 <- 0      #Effect of exposure on cognitive intercept
b02 <- -0.05  #Effect of age on cognitive intercept
b03 <- 0      #Effect of U (unmeasured/underlying variable) on cognitive intercept
b10 <- -0.05  #Cognitive slope for unexposed
b11 <- 0      #Effect of exposure on cognitive slope
b12 <- -0.005 #Effect of age on cognitive slope
b13 <- -0.05  #Effect of U on cognitive slope

#---- Parameters for Sij (survival for person i at time j) ----
g1 <- 0.47    #Effect of exposure on log hazard of death
g2 <- 0.095   #Effect of TD (time of diagnosis) age on log hazard of death
g3 <- 0       #Effect of U on log hazard of death
g4 <- (0.095)*(-0.01) #Interaction effect of exposure and age on log hazard of death
g5 <- 0       #Effect of TD cognitive slope on log hazard of death
g6 <- 0       #Effect of TD cognitive function on log hazard of death

#---- Baseline hazard of death ----
lambda <- 0.004315
