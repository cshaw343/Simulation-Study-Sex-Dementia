#*******************************************************************************
# Simulation Senario C2
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

#Resulting visit times (years from baseline)
visit_times <- seq(from = 0, to = int_time*num_tests, by = int_time)

#---- Prevalance of exposure (female) ----
#Exposure = female gender
pfemale <- 0.51

#---- Variances and correlations ----
ci_var0 <- 0.05       #Variance of random cognitive intercept
ci_var1 <- 0.001      #Variance of random linear term
ci_var2 <- 0.000009   #Variance of random quadratic term 

ci_cov01 <- -0.00009  #Covariance between random intercept and random linear term
ci_cov12 <- 0         #Covariance between random linear and random quadratic term
ci_cov02 <- 0         #Covariance between random intercept and random quadratic term 

cij_var3 <- 1     #Variance of noise for Ci model (cognitive function for person i)

#---- Parameters for Ci ----
#Cognitive function for person i
b00 <- 0      #Cognitive intercept for males 
b01 <- 0      #Effect of female on cognitive intercept
b02 <- 0      #Effect of age on cognitive intercept; Note: Everyone is the same age so there is no age effect (since baseline centered ages are 0 for everyone)
b03 <- -0.5   #Effect of U (unmeasured/underlying variable) on cognitive intercept (loosely based on Marden et. al. 2017)


b10 <- 0.047859  #Cognitive linear term for males                
b11 <- 0         #Effect of female on cognitive linear term
b12 <- 0         #Effect of age on cognitive linear term; Note: Everyone is the same age so there is no age effect
b13 <- 0         #Effect of U on cognitive linear term 

b20 <- -0.0033263  #Cognitive quadratic term for males 
b21 <- 0           #Effect of female on cognitive quadratic term
b22 <- 0           #Effect of age on cognitive quadratic term; Note: Everyone is the same age so there is no age effect 
b23 <- 0           #Effect of U on cognitive quadratic term

#---- Parameters for Sij (survival for person i in age-band j) ----
#Effect of sex (being female) on log hazard of death; 
#Values from survival_opt.R 
g1 <- log(c( 1.775376, 1.435483, 1.162078, 0.930176, 0.689173, 0.479025, 
             0.305908, 0.174879, 0.100314)) 

g2 <- 0         #Effect of U on log hazard of death
g3 <- log(3.5)  #Effect of interaction between male and U on log hazard of death

#---- Baseline hazard of death for unexposed ----
lambda <- c(0.007000, 0.012127, 0.022133, 0.042273, 0.090888, 0.199953, 
            0.491885, 1.315793, 3.378532)

#---- Baseline hazard of random/shock dementia ----
lambda_dj <- 0.007

#---- Dementia Cut Point ----
dem_cut <- -6.5

