#*******************************************************************************
# Simulation Senario A
# 25% cumulative incidence of mortality
# Exposure = male gender
# Based on US lifetables, sample is 51% F / 49% M at age 50
# Remaining LE at age 50:  33.2 years for F and 29.6 years for M
# This is a test change
#*******************************************************************************
#---- Package Loading and Options ----
if (!require("pacman")) 
  install.packages("pacman", repos='http://cran.us.r-project.org')

p_load("tidyverse", "here")

#---- Number of observations ----
num_obs <- 1000

#---- Baseline age ----
age0 <- rep(50, num_obs)

#---- Assessment parameters ----
#Time between assessments (measured in years)
int_time <- 5

#Number of assessments
num_tests <- 10

#Resulting visit times
visit_times <- seq(from = 0, to = int_time*num_tests, by = int_time)

#---- Prevalance of exposure (male) ----
#Exposure = male gender
psex <- 0.49

#---- Variances and correlations ----
cij_var0 <- 0.2   #Variance of random cognitive intercept
cij_var1 <- c(seq(0.001, 0.0064, len = 7), 0.009, 0.0149, rep(0.025, 2)) #Time-dependent variance of random cognitive slope
cij_cov <- 0.01   #Covariance of random intercept and random slope
cij_var3 <- 1     #Variance of noise for Cij (cognitive function for person i at time j)
cij_r1 <- 0.3     #Correlation between noise terms for Cij; this may need to be adjusted
#cij_var4 <- 0.19     #Variance of measurement error of Cij

#---- Parameters for Cij ----
#Read in the parameter table from the slopes_dem-cut_search.R script
search_results <-
  read_csv(paste0(here("Data", 
                  paste0("best_slopes_cuts_20190124.csv"))))

#Filling in the intermediate results table (just to check)
if(nrow(search_results != 10)){
  row <- nrow(search_results)
  search_results[(row + 1):10, "slope"] <- 0
  search_results[(row + 1):10, "dem_cut"] <- search_results[row, "dem_cut"]
  search_results[, "age"] <- seq(55, 100, by = 5)
}

#Cognitive function for person i at time j
b00 <- 0      #Cognitive intercept for females
b01 <- 0      #Effect of sex on cognitive intercept
b02 <- -0.05  #Effect of age on cognitive intercept; Note: Everyone is the same age so there is no age effect (since baseline centered ages are 0 for everyone)
b03 <- 0      #Effect of U (unmeasured/underlying variable) on cognitive intercept

#First value is cognitive slope, the remaining values are changes in cognitive slopes
#These are: b10a, b10b - b10a, b10c - b10b, etc...
#ie Cognitive slope for females age 50-70, change in cognitive slope for females age 70-85, etc...
cij_knots <- seq(55, 95, by = 5) #Specify which ages to place knots
#test slopes
#cij_slopes <- c(0, 0, 0, 0, -0.15, 0, 0, -0.25, 0, 0)

cij_slopes <- head(c(0, search_results$slope), -1)

b11 <- 0      #Effect of sex on cognitive slope
b12 <- -0.005 #Effect of age on cognitive slope; Note: Everyone is the same age so there is no age effect
b13 <- -0.05  #Effect of U on cognitive slope (currently age constant)

#---- Parameters for Sij (survival for person i at time j) ----
#Effect of sex on log hazard of death; chosen using calc from euro_life_tables.R 
g1 <- log(c(1.62, 1.90, 2.13, 2.27, 2.22, 2.09, 1.85, 1.61, 1.39, 1.26, 1.18)) 
g2 <- 0       #Effect of age at time j on log hazard of death (exp(0.095) = 1.10)
g3 <- log(2)  #Effect of U on log hazard of death
g4 <- 0       #Interaction effect of sex and age on log hazard of death
g5 <- 0       #Effect of cognitive slope at time j on log hazard of death
g6 <- 0       #Effect of cognitive function at time j on log hazard of death

#---- Baseline hazard of death for unexposed ----
#Computed in lambda_search_euro.R script
#Based on 35x3000 = 105000 observations
#lambda <- c(0.00414, 0.00577, 0.00824, 0.01260, 0.02105, 0.03605, 0.06316, 
            #0.10918, 0.20142, 0.33345)

#Test US hazards
lambda <- c(0.0211, 0.0227, 0.0243, 0.0262, 0.0280, 0.0299, 0.0327, 0.0359, 
            0.0387, 0.0420)

#---- Dementia Cut Point ----
#Based on slopes_dem-cut_search.R script

#dem_cuts <- head(c(-2.5, best_slopes_cuts$dem_cut), -1)
#test dem_cuts
#dem_cuts <- rep(-3, 10)
dem_cuts <- head(c(search_results$dem_cut[1], search_results$dem_cut), -1)

