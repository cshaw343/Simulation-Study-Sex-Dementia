#Simulation Senario A:  No anticipated bias
#25% cumulative incidence of mortality
#Exposure = male gender


#---- Time between assessments ----
#Time is measured in years
int_time <- 5

#---- Prevalance of exposure ----

pexp <- 0.49

#---- Variances and correlations ----
var0 <- 0.2   #Variance of random cognitive intercept
var1 <- 0.005 #Variance of random cognitive slope
cov <- 0.01   #Covariance of random intercept and random slope
var3 <- 0.70  #Variance of noise for Cij (cognitive function for person i at time j)
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


