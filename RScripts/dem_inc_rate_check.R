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
obs <- data_gen()

#---- Compute incidence rates ----
sim_inc_rates <- matrix(ncol = 9, nrow = 1)
colnames(sim_inc_rates) <- na.omit(variable_names$interval_ages)
rownames(sim_inc_rates) <- c("")

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
  
  sim_inc_rates[1, slot] = round(1000*(sum(PY_data[, dem_this_wave], 
                                       na.rm = TRUE)/
                                     sum(PY_data[, contributed])), 3)
}

#---- Compute survival data ----
#---- Data by sex ----
female_data <- obs %>% filter(sex == 0)

#---- Cohort size ----
num_obs <- nrow(obs)
num_females <- nrow(female_data)

#---- Survival probabilities ----
p_alive <- obs %>% 
  dplyr::select(na.omit(variable_names$deathij_varnames)) %>% 
  map_dbl(~sum(. == 0, na.rm = TRUE))/num_obs

#---- Survival probabilities by sex ----
p_alive_females <- female_data %>%
  dplyr::select(na.omit(variable_names$deathij_varnames)) %>% 
  map_dbl(~sum(. == 0, na.rm = TRUE))/num_females

