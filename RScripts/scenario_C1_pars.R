#*******************************************************************************
# Simulation Senario C1
# 99.1% cumulative incidence of mortality
# Exposure = female gender
# Based on US lifetables, sample is 51% F / 49% M at age 50
# Remaining LE (median) at age 50:  23.5 years for F and 22.5 years for M
#*******************************************************************************

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
b03 <- -0.1   #Effect of U (unmeasured/underlying variable) on cognitive intercept (taken from Marden et. al. 2017)


b10 <- 0.0425 #Cognitive linear term for males                
b11 <- 0      #Effect of female on cognitive linear term
b12 <- 0      #Effect of age on cognitive linear term; Note: Everyone is the same age so there is no age effect
b13 <- 0      #Effect of U on cognitive linear term 

b20 <- -0.002875 #Cognitive quadratic term for males 
b21 <- 0         #Effect of female on cognitive quadratic term
b22 <- 0         #Effect of age on cognitive quadratic term; Note: Everyone is the same age so there is no age effect 
b23 <- 0         #Effect of U on cognitive quadratic term

#---- Parameters for Sij (survival for person i at time j) ----
#Effect of sex (being female) on log hazard of death; 
#chosen using calc from life_table_calcs.R 
g1 <- log(c(1.116428, 1.053103, 0.992379, 0.933184, 0.834850, 0.723279, 
            0.592259, 0.452043, 0.360606)) 

g2 <- 0       #Effect of U on log hazard of death
g3 <- log(2)  #Effect of interaction between female and U on log hazard of death

#---- Baseline hazard of death for unexposed ----
#For male matching
lambda <- c(0.0108112, 0.0166977, 0.0258697, 0.0421338, 0.0742636, 0.1313945, 
            0.2481264, 0.4932790, 0.8745257)


#---- Baseline hazard of random dementia ----
lambda_dj <- 0.007

#---- Dementia Cut Point ----
#Based on slopes_dem-cut_search.R script (results from 20190202)
dem_cut <- -6.5

