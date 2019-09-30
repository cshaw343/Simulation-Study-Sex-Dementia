#*******************************************************************************
# Simulation Senario A
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
cij_var1 <- 0.00009    #Variance of random linear term
cij_var2 <- 0.000005   #Variance of random quadratic term (use tiny value b/c calcs give 0)

cij_cov01 <- -0.00009  #Covariance between random intercept and random linear term
cij_cov12 <- 0         #Covariance between random linear and random quadratic term
cij_cov02 <- 0         #Covariance between random intercept and random quadratic term (use tiny value b/c calcs give 0)
                         

cij_var3 <- 1     #Variance of noise for Cij (cognitive function for person i at time j)

#cij_r1 <- 0.3     #Correlation between noise terms for Cij; this may need to be adjusted
#cij_var4 <- 0.19     #Variance of measurement error of Cij

#---- Parameters for Cij ----
#Cognitive function for person i at time j
b00 <- 0      #Cognitive intercept for males (taken from quad fit to linear model)
b01 <- 0      #Effect of female on cognitive intercept
b02 <- 0      #Effect of age on cognitive intercept; Note: Everyone is the same age so there is no age effect (since baseline centered ages are 0 for everyone)
b03 <- -0.1   #Effect of U (unmeasured/underlying variable) on cognitive intercept (taken from Marden et. al. 2017)


b10 <- 0.08    #Cognitive linear term for males (taken from quad fit to linear model)                
b11 <- 0       #Effect of female on cognitive linear term
b12 <- 0       #Effect of age on cognitive linear term; Note: Everyone is the same age so there is no age effect
b13 <- 0       #Effect of U on cognitive linear term 

b20 <- -0.00425  #Cognitive quadratic term for males (taken from quad fit to linear model)
b21 <- 0         #Effect of female on cognitive quadratic term
b22 <- 0         #Effect of age on cognitive quadratic term; Note: Everyone is the same age so there is no age effect 
b23 <- 0         #Effect of U on cognitive quadratic term

#---- Parameters for Sij (survival for person i at time j) ----
#Effect of sex (being female) on log hazard of death; 
#chosen using calc from life_table_calcs.R 
g1 <- log(c(0.905976, 0.907032, 0.901912, 0.913730, 0.927659, 0.923475, 
            0.931010, 0.924748, 0.999272)) 

g2 <- 0     #Effect of U on log hazard of death
g3 <- 0     #Effect of interaction between female and U on log hazard of death

#---- Baseline hazard of death for unexposed ----
#For male matching
lambda <- c(0.0141156, 0.0201450, 0.0291019, 0.0432166, 0.0671710, 0.1016312, 
            0.1546990, 0.2280313, 0.3027693)

# #For total matching
# lambda <- c(0.0073, 0.0116, 0.018, 0.028, 0.045, 0.07, 0.1, 0.15, 0.27)

#---- Baseline hazard of random dementia ----
lambda_dj <- 0.004

#---- Dementia Cut Point ----
#Based on slopes_dem-cut_search.R script (results from 20190202)
dem_cut <- -5.5

