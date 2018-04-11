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




