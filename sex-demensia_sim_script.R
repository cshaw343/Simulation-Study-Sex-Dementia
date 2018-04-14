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

#---- Generating assessment timepoints ---- 
test_times <- seq(from = 0, to = int_time*num_tests, by = int_time)

#---- The simulation function ----
sex_dem_sim <- function(){
  #Generating the data
  #Generating IDs, baseline age, sex, U
  #Creating centered ages
  obs <- tibble("id" = seq(from = 1, to = num_obs, by = 1), 
                "age0" = age0,
                "sex" = rbinom(num_obs, size = 1, prob = psex), 
                "U" = rnorm(num_obs, mean = 0, sd = 1)) %>% 
    mutate("age_c" = age0 - mean(age0))
  return(obs)
}

sex_dem_sim()
