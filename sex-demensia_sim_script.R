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

age_varnames <- vector(length = num_tests)
for(i in 1:num_tests){
  age_varnames[i] = paste("age", i, sep = "")
}
age_varnames <- append("age0", age_varnames)

#---- The simulation function ----
sex_dem_sim <- function(){
  #Generating the data
  #Generating IDs, sex, U
  obs <- tibble("id" = seq(from = 1, to = num_obs, by = 1), 
                "sex" = rbinom(num_obs, size = 1, prob = psex), 
                "U" = rnorm(num_obs, mean = 0, sd = 1))
  #Creating ages at each visit
  ages <- matrix(NA, nrow = num_obs, ncol = num_tests + 1)
  #Creating centered ages
  
  
  
}

sex_dem_sim()
