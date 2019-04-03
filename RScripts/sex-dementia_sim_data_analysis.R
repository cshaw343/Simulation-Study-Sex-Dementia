#---- Package loading, options ----
if (!require("pacman")) 
  install.packages("pacman", repos='http://cran.us.r-project.org')

p_load("tidyverse", "here", "survival")

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
  num_females <- data %>% filter(female == 1) %>% nrow()
  num_males <- data %>% filter(female == 0) %>% nrow()
  
  #---- Survival probabilities ----
  p_alive <- data %>% 
    dplyr::select(na.omit(variable_names$deathij_varnames)) %>% 
    map_dbl(~sum(. == 0, na.rm = TRUE))/num_obs
  
  #---- Survival probabilities by sex ----
  p_alive_females <- data %>% filter(female == 1) %>%
    dplyr::select(na.omit(variable_names$deathij_varnames)) %>% 
    map_dbl(~sum(. == 0, na.rm = TRUE))/num_females
  
  p_alive_males <- data %>% filter(female == 0) %>%
    dplyr::select(na.omit(variable_names$deathij_varnames)) %>% 
    map_dbl(~sum(. == 0, na.rm = TRUE))/num_males
  
  #---- Mortality logHR (F:M) ----
  simulated_mortality_logHRs <- vector(length = num_tests)
  
  for(i in 1:length(simulated_mortality_logHRs)){
    survtime_name <- paste0("survtime", i - 1, "-", i)
    death_indicator_name <- paste0("death", i - 1, "-", i)
    cox_model <- coxph(Surv(data[, survtime_name], 
                            data[, death_indicator_name]) ~ data$female)
    simulated_mortality_logHRs[i] <- cox_model$coefficients
  }
  
  #---- Number at risk for dementia by sex + distribution of U ----
  at_risk_females <- vector(length = num_tests) 
  at_risk_males <- vector(length = num_tests) 
  
  mean_U_at_risk_females <- vector(length = num_tests) 
  mean_U_at_risk_males <- vector(length = num_tests) 
  
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
      dplyr::select(death_last_wave, dem_last_wave, U) %>% 
      filter(!! as.name(death_last_wave) == 0 & !! as.name(dem_last_wave) == 0) 
    
    at_risk_females[slot] = nrow(at_risk)
    mean_U_at_risk_females[slot] <- mean(at_risk$U)
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
      dplyr::select(death_last_wave, dem_last_wave, U) %>% 
      filter(!! as.name(death_last_wave) == 0 & 
               !! as.name(dem_last_wave) == 0) 
    
    at_risk_males[slot] = nrow(at_risk)
    mean_U_at_risk_males[slot] <- mean(at_risk$U)
  }
  
  #---- Dementia incidence rates ----
  sim_rates <- vector(length = num_tests) 
  
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
    
    sim_rates[slot] = round(1000*(sum(PY_data[, dem_this_wave], na.rm = TRUE)/
                                    sum(PY_data[, contributed])), 3)
  }
  
  #---- Dementia incidence cases, rates, and PY by sex ----
  sim_rates_females <- vector(length = num_tests) 
  sim_rates_males <- vector(length = num_tests) 
  
  inc_cases_females <- vector(length = num_tests) 
  inc_cases_males <- vector(length = num_tests) 
  
  PY_females <- vector(length = num_tests) 
  PY_males <- vector(length = num_tests) 
  
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
    
    inc_cases_females[slot] = sum(PY_data[, dem_this_wave], na.rm = TRUE)
    
    PY_females[slot] = sum(PY_data[, contributed])
    
    sim_rates_females[slot] = round(1000*(sum(PY_data[, dem_this_wave], 
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
    
    inc_cases_males[slot] = sum(PY_data[, dem_this_wave], na.rm = TRUE)
    
    PY_males[slot] = sum(PY_data[, contributed])
    
    sim_rates_males[slot] = round(1000*sum(PY_data[, dem_this_wave], 
                                           na.rm = TRUE)/
                                    sum(PY_data[, contributed]), 3)
  }
  
  IRRs <- sim_rates_females/sim_rates_males
  logIRRs <- log(IRRs) 
  
  # #---- Dementia log(IRRs) Poisson Regression ----
  # modeled_dementia_logIRRs <- vector(length = num_tests)
  # poisson_reg_data <- rbind(cbind(rep("1", length(inc_cases_females)),
  #                                 na.omit(variable_names$interval_ages),
  #                                 inc_cases_females, PY_females), 
  #                           cbind(rep("0", length(inc_cases_females)), 
  #                                 na.omit(variable_names$interval_ages),
  #                                 inc_cases_males, PY_males)) %>% 
  #   as.tibble() %>% 
  #   set_colnames(c("female", "age_interval", "inc_cases", "PY")) %>%
  #   mutate_at(c("female", "inc_cases", "PY"), as.numeric) %>%
  #   mutate("logPY" = log(PY))
  # 
  # for(i in 1:length(modeled_dementia_logIRRs)){
  #   model_data <- tibble("female" = c(1, 0), 
  #                        "py_years" = c(PY_females[i], PY_males[i]), 
  #                        "cases" = c(inc_cases_females[i], inc_cases_males[i]), 
  #                        "logPY" = log(c(PY_females[i], PY_males[i])))
  #   
  #   poisson_model_rate <- glm(cases ~ female + offset(logPY), 
  #                             family = poisson(link = "log"), data = model_data)
  #   modeled_dementia_logIRRs[i] <- poisson_model_rate$coefficients["female"]
  # }
  # 
  # #With ages
  # poisson_model_rate <- 
  #   glm(inc_cases ~ female + age_interval + offset(logPY), 
  #       family = poisson(link = "log"), data = poisson_reg_data)
  # 
  # 
  
  #---- Dementia logHR (F:M) ----
  simulated_dementia_logHRs <- vector(length = num_tests) 
  
  #Need to feed data to this model but without creating new datasets
  
  for(i in 1:length(simulated_dementia_logHRs)){
    PY_contributed_name <- paste0("contributed", i - 1, "-", i)
    dem_indicator_name <- paste0("dem", i - 1, "-", i)
    cox_model <- coxph(Surv(data[, PY_contributed_name], 
                            data[, dem_indicator_name]) ~ data$female)
    simulated_dementia_logHRs[i] <- cox_model$coefficients
  }
  
  #---- Total dementia cases (incident + prevalent) ----
  dem_cases_female <- female_data %>% 
    dplyr::select(variable_names$dem_varnames[-1]) %>% 
    colSums(na.rm = TRUE) 
  
  dem_cases_male <- male_data %>% 
    dplyr::select(variable_names$dem_varnames[-1]) %>% 
    colSums(na.rm = TRUE)
  
  #---- Probability of dementia ----
  #Prevalence of dementia by sex
  prop_dem_females <- female_data %>% 
    dplyr::select(variable_names$dem_varnames[-1]) %>% 
    colMeans(na.rm = TRUE) 
  
  prop_dem_males <- male_data %>% 
    dplyr::select(variable_names$dem_varnames[-1]) %>% 
    colMeans(na.rm = TRUE) 
  
  #---- Vector to return ----
  results_vec <- c(num_obs, num_females, num_males, p_alive, p_alive_females, 
                   p_alive_males, simulated_mortality_logHRs, at_risk_females, 
                   at_risk_males, sim_rates, sim_rates_females, sim_rates_males, 
                   inc_cases_females, inc_cases_males, PY_females, PY_males, 
                   logIRRs, simulated_dementia_logHRs, dem_cases_female, 
                   dem_cases_male, prop_dem_females, prop_dem_males, 
                   mean_U_at_risk_females, mean_U_at_risk_males)
  
  names(results_vec) <- c("num_obs_baseline", "num_females_baseline", 
                          "num_males_baseline", 
                          na.omit(variable_names$p_alive_varnames), 
                          na.omit(variable_names$p_alive_females_varnames), 
                          na.omit(variable_names$p_alive_males_varnames), 
                          na.omit(variable_names$mortality_logHR_varnames),
                          na.omit(variable_names$at_risk_females_varnames), 
                          na.omit(variable_names$at_risk_males_varnames), 
                          na.omit(variable_names$dem_inc_rate_varnames), 
                          na.omit(variable_names$dem_inc_rate_females_varnames), 
                          na.omit(variable_names$dem_inc_rate_males_varnames), 
                          na.omit(variable_names$inc_cases_females_varnames), 
                          na.omit(variable_names$inc_cases_males_varnames), 
                          na.omit(variable_names$PY_females_varnames), 
                          na.omit(variable_names$PY_males_varnames), 
                          na.omit(variable_names$logIRR_varnames), 
                          na.omit(variable_names$dementia_logHR_varnames), 
                          na.omit(variable_names$dem_cases_females_varnames), 
                          na.omit(variable_names$dem_cases_males_varnames), 
                          na.omit(variable_names$prop_dem_females_varnames), 
                          na.omit(variable_names$prop_dem_males_varnames), 
                          na.omit(
                            variable_names$mean_U_at_risk_females_varnames), 
                          na.omit(variable_names$mean_U_at_risk_males_varnames))
  
  #---- Return ----
  return(results_vec)
}





  