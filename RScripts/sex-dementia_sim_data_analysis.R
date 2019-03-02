#---- Package loading, options ----
if (!require("pacman")) 
  install.packages("pacman", repos='http://cran.us.r-project.org')

p_load("tidyverse", "here")

options(scipen = 999) #Standard Notation
options(digits = 6)   #Round to 6 decimal places
options(warn = -1)    #Suppress warnings

#---- Source files ----
source(here("RScripts", "variable_names.R"))
source(here("RScripts", "sex-dementia_sim_data_gen.R"))

#---- Simulation function ----
sex_dem_sim <- function(){
  
  #---- Generate the data ----
  data <- data_gen()
  
  #---- Dementia Incidence Rates by Sex ----
  sim_rates_by_sex <- matrix(ncol = 9, nrow = 2)
  colnames(sim_rates_by_sex) <- na.omit(variable_names$interval_ages)
  rownames(sim_rates_by_sex) <- c("Female", "Male")
  
  female_data <- data %>% filter(sex == 0)
  male_data <- data %>% filter(sex == 1)
  
  #Computing female incidence rates
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
    PY_data <- female_data %>% 
      dplyr::select(death_last_wave, death_this_wave, 
                    dem_last_wave, dem_this_wave, contributed) %>% 
      filter(!! as.name(death_last_wave) == 0 & 
               !! as.name(dem_last_wave) == 0) 
    
    sim_rates_by_sex[1, slot] = round(1000*(sum(PY_data[, dem_this_wave], 
                                                na.rm = TRUE)/
                                              sum(PY_data[, contributed])), 3)
  }
  
  #Computing male incidence rates
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
    PY_data <- male_data %>% 
      dplyr::select(death_last_wave, death_this_wave, 
                    dem_last_wave, dem_this_wave, contributed) %>% 
      filter(!! as.name(death_last_wave) == 0 & 
               !! as.name(dem_last_wave) == 0) 
    
    sim_rates_by_sex[2, slot] = round(1000*sum(PY_data[, dem_this_wave], 
                                               na.rm = TRUE)/
                                        sum(PY_data[, contributed]), 3)
  }
  
  IRRs <- unlist(sim_rates_by_sex[1, ]/sim_rates_by_sex[2, ]) %>% t()
  logIRRs <- log(IRRs) %>% as.data.frame()
  
  #---- Return ----
  return("logIRRs" = logIRRs)
}





  