#---- Package loading, options ----
if (!require("pacman")) 
  install.packages("pacman", repos='http://cran.us.r-project.org')

p_load("tidyverse", "here")

options(scipen = 999) #Standard Notation
options(digits = 6)   #Round to 6 decimal places
options(warn = -1)    #Suppress warnings

#---- Source files ----
source(here("RScripts", "sex-dementia_sim_parA.R"))
source(here("RScripts", "variable_names.R"))
source(here("RScripts", "sex-dementia_sim_data_gen.R"))

#---- Plug in newly optimized data ----
cij_slopes <- opt_cij_slopes
cij_var1 <- opt_cij_var1
lambda <- opt_base_haz

#---- Generate the data ----
num_obs = 200000
obs <- data_gen(num_obs) %>% as.data.frame()

#---- Compute incidence rates ----
male_sim_inc_rates <- matrix(ncol = 9, nrow = 1)
colnames(male_sim_inc_rates) <- variable_names$interval_ages[1:num_tests]
rownames(male_sim_inc_rates) <- c("")

for(slot in 1:num_tests){
  if(slot == 1){
    dem_last_wave <- paste0("dem", (slot - 1))
    dem_this_wave <- paste0("dem", (slot - 1), "-", slot)
    death_last_wave <- paste0("death", (slot - 1))
    death_this_wave <- paste0("death", (slot - 1), "-", slot)
    contributed <- paste0("contributed", (slot - 1), "-", slot)
  } else {
    dem_last_wave <- paste0("dem", (slot - 2), "-", (slot - 1))
    dem_this_wave <- paste0("dem", (slot - 1), "-", slot)
    death_last_wave <- paste0("death", (slot - 2), "-", (slot - 1))
    death_this_wave <- paste0("death", (slot - 1), "-", slot)
    contributed <- paste0("contributed", (slot - 1), "-", slot)
  }
  PY_data <- obs %>% 
    dplyr::select(death_last_wave, death_this_wave, 
                  dem_last_wave, dem_this_wave, contributed) %>% 
    filter(!! as.name(death_last_wave) == 0 & 
             !! as.name(dem_last_wave) == 0) 
  
  male_sim_inc_rates[1, slot] = round(1000*(sum(PY_data[, dem_this_wave], 
                                       na.rm = TRUE)/
                                     sum(PY_data[, contributed])), 3)
}

#---- Compute survival data ----
#---- Data by sex ----
male_data <- obs %>% filter(female == 0)

#---- Cohort size ----
num_obs <- nrow(obs)
num_males <- nrow(male_data)

#---- Survival probabilities ----
p_alive <- obs %>% 
  dplyr::select(variable_names$deathij_varnames[1:num_tests]) %>% 
  map_dbl(~sum(. == 0, na.rm = TRUE))/num_obs

#---- Survival probabilities by sex ----
p_alive_males <- male_data %>%
  dplyr::select(variable_names$deathij_varnames[1:num_tests]) %>% 
  map_dbl(~sum(. == 0, na.rm = TRUE))/num_males

