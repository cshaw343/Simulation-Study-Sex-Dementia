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
  
  #---- Data by sex ----
  female_data <- data %>% filter(sex == 0)
  male_data <- data %>% filter(sex == 1)
  
  #---- Cohort size ----
  num_obs <- nrow(data)
  num_females <- nrow(female_data)
  num_males <- nrow(male_data)
  
  #---- Survival probabilities ----
  p_alive <- data %>% 
    dplyr::select(na.omit(variable_names$deathij_varnames)) %>% 
    map_dbl(~sum(. == 0, na.rm = TRUE))/num_obs
  
  #---- Survival probabilities by sex ----
  p_alive_females <- female_data %>% 
    dplyr::select(na.omit(variable_names$deathij_varnames)) %>% 
    map_dbl(~sum(. == 0, na.rm = TRUE))/num_females
  
  p_alive_males <- male_data %>% 
    dplyr::select(na.omit(variable_names$deathij_varnames)) %>% 
    map_dbl(~sum(. == 0, na.rm = TRUE))/num_males
  
  #---- Number at risk by sex ----
  at_risk_by_sex <- matrix(ncol = 9, nrow = 2)
  colnames(at_risk_by_sex) <- seq(50, 90, by = 5)
  rownames(at_risk_by_sex) <- c("Female", "Male")
  
  #Computing female at_risk
  for(slot in 1:num_tests){
    if(slot == 1){
      dem_last_wave <- paste0("dem", (slot - 1))
      death_last_wave <- paste0("death", (slot - 1))
    } else {
      dem_last_wave <- paste0("dem", (slot - 2), "-", (slot - 1))
      death_last_wave <- paste0("death", (slot - 2), "-", (slot - 1))
    }
    at_risk <- female_data %>% 
      dplyr::select(death_last_wave, dem_last_wave) %>% 
      filter(!! as.name(death_last_wave) == 0 & 
               !! as.name(dem_last_wave) == 0) 
    
    at_risk_by_sex[1, slot] = nrow(at_risk)
  }
  
  #Computing male at_risk
  for(slot in 1:num_tests){
    if(slot == 1){
      dem_last_wave <- paste0("dem", (slot - 1))
      death_last_wave <- paste0("death", (slot - 1))
    } else {
      dem_last_wave <- paste0("dem", (slot - 2), "-", (slot - 1))
      death_last_wave <- paste0("death", (slot - 2), "-", (slot - 1))
    }
    at_risk <- male_data %>% 
      dplyr::select(death_last_wave, dem_last_wave) %>% 
      filter(!! as.name(death_last_wave) == 0 & 
               !! as.name(dem_last_wave) == 0) 
    
    at_risk_by_sex[2, slot] = nrow(at_risk)
  }
  
  #---- Dementia incidence rates ----
  sim_rates <- matrix(ncol = 9, nrow = 1)
  colnames(sim_rates) <- na.omit(variable_names$interval_ages)
  
  #Computing incidence rates
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
    PY_data <- data %>% 
      dplyr::select(death_last_wave, death_this_wave, 
                    dem_last_wave, dem_this_wave, contributed) %>% 
      filter(!! as.name(death_last_wave) == 0 & 
               !! as.name(dem_last_wave) == 0) 
    
    sim_rates[1, slot] = round(1000*(sum(PY_data[, dem_this_wave], 
                                                na.rm = TRUE)/
                                              sum(PY_data[, contributed])), 3)
  }
  
  #---- Dementia incidence cases, rates, and PY by sex ----
  sim_rates_by_sex <- matrix(ncol = 9, nrow = 2)
  colnames(sim_rates_by_sex) <- na.omit(variable_names$interval_ages)
  rownames(sim_rates_by_sex) <- c("Female", "Male")
  
  inc_cases_by_sex <- matrix(ncol = 9, nrow = 2)
  colnames(inc_cases_by_sex) <- na.omit(variable_names$interval_ages)
  rownames(inc_cases_by_sex) <- c("Female", "Male")
  
  PY_by_sex <- matrix(ncol = 9, nrow = 2)
  colnames(PY_by_sex) <- na.omit(variable_names$interval_ages)
  rownames(PY_by_sex) <- c("Female", "Male")
  
  #Computing female incidence cases, rates, PY
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
    
    inc_cases_by_sex[1, slot] = sum(PY_data[, dem_this_wave], 
                                    na.rm = TRUE)
    
    PY_by_sex[1, slot] = sum(PY_data[, contributed])
    
    sim_rates_by_sex[1, slot] = round(1000*(sum(PY_data[, dem_this_wave], 
                                                na.rm = TRUE)/
                                              sum(PY_data[, contributed])), 3)
  }
  
  #Computing male incidence cases and rates
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
    
    inc_cases_by_sex[2, slot] = sum(PY_data[, dem_this_wave], 
                                    na.rm = TRUE)
    
    PY_by_sex[2, slot] = sum(PY_data[, contributed])
    
    sim_rates_by_sex[2, slot] = round(1000*sum(PY_data[, dem_this_wave], 
                                               na.rm = TRUE)/
                                        sum(PY_data[, contributed]), 3)
  }
  
  IRRs <- unlist(sim_rates_by_sex[1, ]/sim_rates_by_sex[2, ]) %>% t()
  logIRRs <- log(IRRs) %>% as.data.frame()
  
  #---- Total dementia cases (incident + prevalent) ----
  dem_cases_female <- female_data %>% 
    dplyr::select(variable_names$dem_varnames[-1]) %>% 
    colSums(na.rm = TRUE) 
  
  dem_cases_male <- male_data %>% 
    dplyr::select(variable_names$dem_varnames[-1]) %>% 
    colSums(na.rm = TRUE)
  
  dem_cases_by_sex <- rbind(dem_cases_female, dem_cases_male) %>% 
    as.data.frame() %>%
    set_colnames(na.omit(variable_names$interval_ages)) %>% 
    set_rownames(c("Female", "Male"))
  
  #---- Probability of dementia ----
  #Prevalence of dementia by sex
  dem_prev_females <- female_data %>% 
    dplyr::select(variable_names$dem_varnames[-1]) %>% 
    colMeans(na.rm = TRUE) 
  
  dem_prev_males <- male_data %>% 
    dplyr::select(variable_names$dem_varnames[-1]) %>% 
    colMeans(na.rm = TRUE) 
  
  dem_prob_by_sex <- rbind(dem_prev_females, dem_prev_males) %>% 
    as.data.frame() %>%
    set_colnames(na.omit(variable_names$interval_ages)) %>% 
    set_rownames(c("Female", "Male"))
  
  #---- Return ----
  return(list("num_obs" = num_obs, "num_females" = num_females, 
              "num_males" = num_males, "p_alive" = p_alive, 
              "p_alive_females" = p_alive_females, 
              "p_alive_males" = p_alive_males,
              "at_risk_by_sex" = at_risk_by_sex, 
              "inc_cases_by_sex" = inc_cases_by_sex, 
              "dem_cases_by_sex" = dem_cases_by_sex, "PY_by_sex" = PY_by_sex,
              "dem_prob_by_sex" = dem_prob_by_sex, "sim_rates" = sim_rates, 
              "logIRRs" = logIRRs))
}





  