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

#---- Generating variable names at assessment timepoints ----
#Generating age labels at each assessment timepoint
age_varnames <- vector(length = num_tests)
for(i in 1:num_tests){
  age_varnames[i] = paste("age", i, sep = "")
}
age_varnames <- append("age0", age_varnames)

#Generating centered age labels at each assessment timepoint
agec_varnames <- vector(length = length(age_varnames))
for(i in 1:length(agec_varnames)){
  agec_varnames[i] = paste(age_varnames[i], "_c50", sep = "")
}

#---- The simulation function ----
sex_dem_sim <- function(){
  #Generating the data
  #Generating IDs, sex, U
  obs <- tibble("id" = seq(from = 1, to = num_obs, by = 1), 
                "sex" = rbinom(num_obs, size = 1, prob = psex), 
                "U" = rnorm(num_obs, mean = 0, sd = 1))
  #Creating ages at each visit
  ages <- as_tibble(matrix(NA, nrow = num_obs, ncol = length(age_varnames)))
  for(i in 1:length(age_varnames)){
    if(i == 1){
      ages[, i] = age0
    } else ages[, i] = ages[, (i-1)] + int_time
  }
  colnames(ages) <- age_varnames
  
  #Creating centered ages
  c_ages <- as_tibble(ages - mean(age0))
  colnames(c_ages) <- agec_varnames
}

sex_dem_sim()
