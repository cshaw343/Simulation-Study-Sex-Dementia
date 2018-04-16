#---- Package loading, options, seed ----
if (!require("pacman")) 
  install.packages("pacman", repos='http://cran.us.r-project.org')

p_load("tidyverse")

#Using standard notation (as opposed to scientific), rounded to three 
#decimal places
options(scipen = 999)
options(digits = 3)

set.seed(10789)

#---- Specify the parameter file ----
source("sex-demensia_sim_parA.R")

#---- Generating assessment timepoint data ----
visit_times <- seq(from = int_time, to = int_time*num_tests, by = int_time)

#---- Generating age data ----
#Creating ages at each timepoint j
ages <- as_tibble(matrix(NA, nrow = num_obs, ncol = length(age_varnames)))
for(i in 1:length(age_varnames)){
  if(i == 1){
    ages[, i] = age0
  } else ages[, i] = ages[, (i-1)] + int_time
}

#Age labels at each assessment timepoint
age_varnames <- vector(length = num_tests)
for(i in 1:num_tests){
  age_varnames[i] = paste("age", i, sep = "")
}
age_varnames <- append("age0", age_varnames)

#Assigning column labels
colnames(ages) <- age_varnames

#---- Generating centered age data ----
#Creating centered ages at each timepoint j
c_ages <- as_tibble(ages - mean(age0))

#Centered age labels at each assessment timepoint
agec_varnames <- vector(length = length(age_varnames))
for(i in 1:length(agec_varnames)){
  agec_varnames[i] = paste(age_varnames[i], "_c50", sep = "")
}

#Assigning column labels
colnames(c_ages) <- agec_varnames




#---- Generating variable names at assessment timepoints ----

#Generating alpha labels at each assessment timepoint
alpha_varnames <- vector(length = num_tests)
for(i in 1:num_tests){
  alpha_varnames[i] = paste("alpha", i, sep = "")
}

#---- Model for Cognitive Function ----
Cij <- function(t){
  knots = c(0, 20, 35)
  if(t >= 0 & t < 20){
    
  }
}

#---- Generate Covariance Matrix for random slope and intercept terms ----
slope_int_cov <- matrix(c(var0, cov, cov, var1), nrow = 2, byrow = TRUE)

#---- The simulation function ----
sex_dem_sim <- function(){
  #Generating the data
  #Generating IDs, sex, U
  obs <- tibble("id" = seq(from = 1, to = num_obs, by = 1), 
                "sex" = rbinom(num_obs, size = 1, prob = psex), 
                "U" = rnorm(num_obs, mean = 0, sd = 1))
  
  #Creating "true" and "measured" cognitive function at each study assessment
  #measured = true + error
  
  #
  
  
}


sex_dem_sim()
