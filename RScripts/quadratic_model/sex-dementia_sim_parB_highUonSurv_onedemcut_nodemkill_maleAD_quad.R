#*******************************************************************************
# Simulation Senario B
# 25% cumulative incidence of mortality
# Exposure = male gender
# Based on US lifetables, sample is 51% F / 49% M at age 50
# Remaining LE at age 50:  33.2 years for F and 29.6 years for M
# This is a test change
#*******************************************************************************

#---- Number of simulation runs ----
runs = 1000

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

#---- Prevalance of exposure (female) ----
#Exposure = female gender
pfemale <- 0.51

#---- Variances and correlations ----
#Calculated based on the work in test_quadratic_trajectories.R
cij_var0 <- 0.05       #Variance of random cognitive intercept
cij_var1 <- 0.001      #Variance of random linear term
cij_var2 <- 0.000009   #Variance of random quadratic term (use tiny value b/c calcs give 0)

cij_cov01 <- -0.00009  #Covariance between random intercept and random linear term
cij_cov12 <- 0         #Covariance between random linear and random quadratic term
cij_cov02 <- 0         #Covariance between random intercept and random quadratic term (use tiny value b/c calcs give 0)
                         

cij_var3 <- 1     #Variance of noise for Cij (cognitive function for person i at time j)

#cij_r1 <- 0.3     #Correlation between noise terms for Cij; this may need to be adjusted
#cij_var4 <- 0.19     #Variance of measurement error of Cij

#---- Parameters for Cij ----
#Cognitive function for person i at time j
b00 <- 0      #Cognitive intercept for males 
b01 <- 0      #Effect of female on cognitive intercept
b02 <- 0      #Effect of age on cognitive intercept; Note: Everyone is the same age so there is no age effect (since baseline centered ages are 0 for everyone)
b03 <- -0.5   #Effect of U (unmeasured/underlying variable) on cognitive intercept (taken from Marden et. al. 2017)


b10 <- 0.049   #Cognitive linear term for males                
b11 <- 0       #Effect of female on cognitive linear term
b12 <- 0       #Effect of age on cognitive linear term; Note: Everyone is the same age so there is no age effect
b13 <- 0       #Effect of U on cognitive linear term 

b20 <- -0.003075   #Cognitive quadratic term for males 
b21 <- 0         #Effect of female on cognitive quadratic term
b22 <- 0         #Effect of age on cognitive quadratic term; Note: Everyone is the same age so there is no age effect 
b23 <- 0         #Effect of U on cognitive quadratic term

#---- Parameters for Sij (survival for person i at time j) ----
#Effect of sex (being female) on log hazard of death; 
#chosen using calc from life_table_calcs.R 
g1 <- log(c(0.911228, 0.884474, 0.888991, 0.887634, 0.870331, 0.870331, 
            0.878532, 0.887716, 0.906304))

g2 <- log(3.5)  #Effect of U on log hazard of death
g3 <- 0         #Effect of interaction between male and U on log hazard of death

#---- Baseline hazard of death for unexposed ----
#For male matching
lambda <- c(0.00715524, 0.01244461, 0.02231527, 0.04281595, 0.09159394, 
            0.20048593, 0.49220949, 1.22742358, 3.27047484)

#---- Baseline hazard of random dementia ----
lambda_dj <- 0.004

#---- Dementia Cut Point ----
dem_cut <- -6.5

