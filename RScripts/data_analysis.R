#---- Package loading, options ----
if (!require("pacman")) 
  install.packages("pacman", repos='http://cran.us.r-project.org')

p_load("tidyverse", "here", "survival")

options(scipen = 999) #Standard Notation
options(digits = 6)   #Round to 6 decimal places
options(warn = -1)    #Suppress warnings

#---- Simulation function ----
sex_dem_sim <- function(num_obs){
  
  #---- Generate the data ----
  data <- data_gen(num_obs)
  
  #---- Data by sex ----
  female_data <- data %>% filter(female == 1)
  male_data <- data %>% filter(female == 0)
  
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
  
  #---- Conditional survival probabilities ----
  cp_alive <- vector(length = length(na.omit(variable_names$deathij_varnames)))
  alive <- data %>% 
    dplyr::select(variable_names$deathij_varnames[1:num_tests])
  
  for(i in 1:length(cp_alive)){
    if(i == 1){
      cp_alive[i] = 
        sum(alive[, variable_names$deathij_varnames[i]] == 0)/nrow(alive)
    } else{
      alive <- alive[alive[, variable_names$deathij_varnames[i - 1]] == 0, ]
      cp_alive[i] = 
        sum(alive[, variable_names$deathij_varnames[i]] == 0)/nrow(alive)
    }
  }
  
  #---- Conditional survival probabilities by sex ----
  cp_alive_males <- 
    vector(length = length(na.omit(variable_names$deathij_varnames)))
  alive_males <- male_data %>% 
    dplyr::select(variable_names$deathij_varnames[1:num_tests])
  
  for(i in 1:length(cp_alive_males)){
    if(i == 1){
      cp_alive_males[i] = 
        sum(alive_males[, variable_names$deathij_varnames[i]] == 0)/
        nrow(alive_males)
    } else{
      alive_males <- 
        alive_males[
          alive_males[, variable_names$deathij_varnames[i - 1]] == 0, ]
      cp_alive_males[i] = 
        sum(alive_males[, variable_names$deathij_varnames[i]] == 0)/
        nrow(alive_males)
    }
  }
  
  cp_alive_females <- 
    vector(length = length(na.omit(variable_names$deathij_varnames)))
  alive_females <- female_data %>% 
    dplyr::select(variable_names$deathij_varnames[1:num_tests])
  
  for(i in 1:length(cp_alive_females)){
    if(i == 1){
      cp_alive_females[i] = 
        sum(alive_females[, variable_names$deathij_varnames[i]] == 0)/
        nrow(alive_females)
    } else{
      alive_females <- 
        alive_females[
          alive_females[, variable_names$deathij_varnames[i - 1]] == 0, ]
      cp_alive_females[i] = 
        sum(alive_females[, variable_names$deathij_varnames[i]] == 0)/
        nrow(alive_females)
    }
  }
  
  #---- Mortality logHR (F:M) ----
  survival_data <- data[, variable_names$Sij_varnames[1:num_tests]]
  survival_data[survival_data == 5] <- 4.99999
  data[, variable_names$Sij_varnames[1:num_tests]] <- survival_data
  
  simulated_mortality_logHRs <- vector(length = num_tests)
  
  for(i in 1:length(simulated_mortality_logHRs)){
    survtime_name <- paste0("survtime", i - 1, "-", i)
    death_indicator_name <- paste0("death", i - 1, "-", i)
    cox_model <- coxph(Surv(data[, survtime_name], 
                            data[, death_indicator_name]) ~ data$female)
    simulated_mortality_logHRs[i] <- cox_model$coefficients
  }
  
  #---- Total dementia incidence rates (5-year bands) ----
  inc_cases <- colSums(data[, variable_names$dem_varnames])
  PY <- colSums(data[, variable_names$contributed_varnames[1:9]])
  sim_rates <- inc_cases[-1]/PY*1000
  
  #---- Dementia incidence cases, rates, and PY by sex (5-year bands) ----
  inc_cases_females <- colSums(female_data[, variable_names$dem_varnames])
  inc_cases_males <- colSums(male_data[, variable_names$dem_varnames])
  
  PY_females <- colSums(female_data[, variable_names$contributed_varnames[1:9]]) 
  PY_males <- colSums(male_data[, variable_names$contributed_varnames[1:9]])
  
  sim_rates_females <- inc_cases_females[-1]/PY_females*1000
  sim_rates_males <- inc_cases_males[-1]/PY_males*1000
  
  #---- Dementia logIRRs (5-year bands) ----
  IRRs <- sim_rates_females/sim_rates_males

  simulated_dementia_logIRRs_data <- matrix(nrow = num_tests, ncol = 5)
  #Point estimate
  simulated_dementia_logIRRs_data[, 1] <- log(IRRs)
  #SE calculation for log(IRR)
  simulated_dementia_logIRRs_data[, 2] <-
    sqrt(1/inc_cases_females[-1] + 1/inc_cases_males[-1])
  #95% CI Lower Bound
  simulated_dementia_logIRRs_data[, 3] <-
    simulated_dementia_logIRRs_data[, 1] -
    1.96*simulated_dementia_logIRRs_data[, 2]
  #95% CI Upper Bound
  simulated_dementia_logIRRs_data[, 4] <-
    simulated_dementia_logIRRs_data[, 1] +
    1.96*simulated_dementia_logIRRs_data[, 2]
  #95% CI Coverage
  simulated_dementia_logIRRs_data[, 5] <-
    (simulated_dementia_logIRRs_data[, 3] < 0 &
    simulated_dementia_logIRRs_data[, 4] > 0)*1

  simulated_dementia_logIRRs_data <-
    as.vector(t(simulated_dementia_logIRRs_data))
  
  #---- Dementia incidence cases, rates, and PY by sex (1-year bands) ----
  inc_cases_females_1year <- 
    colSums(female_data[, variable_names_1year$dem_varnames])
  inc_cases_males_1year <- 
    colSums(male_data[, variable_names_1year$dem_varnames])
  
  PY_females_1year <- 
    colSums(female_data[, variable_names_1year$contributed_varnames]) 
  PY_males_1year <- 
    colSums(male_data[, variable_names_1year$contributed_varnames])
  
  sim_rates_females_1year <- inc_cases_females_1year/PY_females_1year*1000
  sim_rates_males_1year <- inc_cases_males_1year/PY_males_1year*1000

  #---- Dementia logIRRs (1-year bands) ----
  IRRs_1year <- sim_rates_females_1year/sim_rates_males_1year
  
  simulated_dementia_logIRRs_data_1yr <- matrix(nrow = num_tests*5, ncol = 5)
  #Point estimate
  simulated_dementia_logIRRs_data_1yr[, 1] <- log(IRRs_1year)
  #SE calculation for log(IRR)
  simulated_dementia_logIRRs_data_1yr[, 2] <-
    sqrt(1/inc_cases_females_1year + 1/inc_cases_males_1year)
  #95% CI Lower Bound
  simulated_dementia_logIRRs_data_1yr[, 3] <-
    simulated_dementia_logIRRs_data_1yr[, 1] -
    1.96*simulated_dementia_logIRRs_data_1yr[, 2]
  #95% CI Upper Bound
  simulated_dementia_logIRRs_data_1yr[, 4] <-
    simulated_dementia_logIRRs_data_1yr[, 1] +
    1.96*simulated_dementia_logIRRs_data_1yr[, 2]
  #95% CI Coverage
  simulated_dementia_logIRRs_data_1yr[, 5] <-
    (simulated_dementia_logIRRs_data_1yr[, 3] < 0 &
       simulated_dementia_logIRRs_data_1yr[, 4] > 0)*1

  simulated_dementia_logIRRs_data_1yr <-
    as.vector(t(simulated_dementia_logIRRs_data_1yr))
  
  #---- Dementia types ----
  demented <- data %>% filter(dem == 1) %>% 
    dplyr::select("female", "dem", "dem_Ci", "dem_random", 
                  variable_names$dem_varnames[-1]) %>% 
    mutate("dem_both" = ifelse(dem_Ci == 1 & dem_random == 1, 1, 0))
  
  n_demented_by_age <- demented %>% 
    dplyr::select(variable_names$dem_varnames[-1]) %>% colSums() 
  
  demented_women <- demented %>% filter(female == 1)
  n_demented_by_age_W <- demented_women %>% 
    dplyr::select(variable_names$dem_varnames[-1]) %>% colSums() 
  
  demented_men <- demented %>% filter(female == 0) 
  n_demented_by_age_M <- demented_men %>% 
    dplyr::select(variable_names$dem_varnames[-1]) %>% colSums() 
  
  #---- Random dementia people ----
  #Overall
  n_dem_random <- demented %>% filter(dem_both != 1) %>% 
    summarise_at("dem_random", sum) %>% as.numeric()
  prop_dem_random <- n_dem_random/nrow(demented)
  
  n_dem_random_W <- demented_women %>% filter(dem_both != 1) %>% 
    summarise_at("dem_random", sum) %>% as.numeric()
  prop_dem_random_W <- n_dem_random_W/nrow(demented_women)
  
  n_dem_random_M <- demented_men %>% filter(dem_both != 1) %>% 
    summarise_at("dem_random", sum) %>% as.numeric()
  prop_dem_random_M <- n_dem_random/nrow(demented)
  
  #Age-band specific
  n_dem_random_by_age <- demented %>% 
    filter(dem_random == 1 & dem_both != 1) %>% 
    dplyr::select(variable_names$dem_varnames[-1]) %>% colSums() 
  prop_dem_random_by_age <- n_dem_random_by_age/n_demented_by_age
  
  n_dem_random_by_age_W <- demented_women %>% 
    filter(dem_random == 1 & dem_both != 1) %>% 
    dplyr::select(variable_names$dem_varnames[-1]) %>% colSums() 
  prop_dem_random_by_age_W <- n_dem_random_by_age_W/n_demented_by_age_W
  
  n_dem_random_by_age_M <- demented_men %>% 
    filter(dem_random == 1 & dem_both != 1) %>% 
    dplyr::select(variable_names$dem_varnames[-1]) %>% colSums() 
  prop_dem_random_by_age_M <- n_dem_random_by_age_M/n_demented_by_age_M
  
  #---- Ci dementia people ----
  #Overall
  n_dem_Ci <- demented %>% filter(dem_both != 1) %>% 
    summarise_at("dem_Ci", sum) %>% as.numeric()
  prop_dem_Ci <- n_dem_Ci/nrow(demented)
  
  n_dem_Ci_W <- demented_women %>% filter(dem_both != 1) %>% 
    summarise_at("dem_Ci", sum) %>% as.numeric()
  prop_dem_Ci_W <- n_dem_Ci_W/nrow(demented_women)
  
  n_dem_Ci_M <- demented_men %>% filter(dem_both != 1) %>% 
    summarise_at("dem_Ci", sum) %>% as.numeric()
  prop_dem_Ci_M <- n_dem_Ci_M/nrow(demented_men)
  
  #Age-band specific
  n_dem_Ci_by_age <- demented %>% filter(dem_Ci == 1 & dem_both != 1) %>% 
    dplyr::select(variable_names$dem_varnames[-1]) %>% colSums() 
  prop_dem_Ci_by_age <- n_dem_Ci_by_age/n_demented_by_age
  
  n_dem_Ci_by_age_W <- demented_women %>% 
    filter(dem_Ci == 1 & dem_both != 1) %>% 
    dplyr::select(variable_names$dem_varnames[-1]) %>% colSums() 
  prop_dem_Ci_by_age_W <- n_dem_Ci_by_age_W/n_demented_by_age_W
  
  n_dem_Ci_by_age_M <- demented_men %>% 
    filter(dem_Ci == 1 & dem_both != 1) %>% 
    dplyr::select(variable_names$dem_varnames[-1]) %>% colSums() 
  prop_dem_Ci_by_age_M <- n_dem_Ci_by_age_M/n_demented_by_age_M
  
  #---- Both types dementia people ----
  dem_both <- demented %>% filter(dem_both == 1)
  dem_both_women <- dem_both %>% filter(female == 1)
  dem_both_men <- dem_both %>% filter(female == 0)
  
  #Overall
  prop_dem_both <- nrow(dem_both)/nrow(demented)
  prop_dem_both_W <- nrow(dem_both_women)/nrow(demented_women)
  prop_dem_both_M <-  nrow(dem_both_men)/nrow(demented_men)
  
  #Age-band specific
  n_dem_both_by_age <- dem_both %>% 
    dplyr::select(variable_names$dem_varnames[-1]) %>% colSums() 
  prop_dem_both_by_age <- n_dem_both_by_age/n_demented_by_age
  
  n_dem_both_by_age_W <- dem_both_women %>% 
    dplyr::select(variable_names$dem_varnames[-1]) %>% colSums() 
  prop_dem_both_by_age_W <- n_dem_both_by_age_W/n_demented_by_age_W
  
  n_dem_both_by_age_M <- dem_both_men %>% 
    dplyr::select(variable_names$dem_varnames[-1]) %>% colSums() 
  prop_dem_both_by_age_M <- n_dem_both_by_age_M/n_demented_by_age_M
  
  #---- Vector to return ----
  results_vec <- c(num_obs, num_females, num_males, p_alive, p_alive_females, 
                   p_alive_males, cp_alive, cp_alive_females, cp_alive_males,
                   simulated_mortality_logHRs, sim_rates, sim_rates_females, 
                   sim_rates_males, inc_cases_females[-1], inc_cases_males[-1], 
                   PY_females, PY_males, inc_cases_females_1year, 
                   inc_cases_males_1year, PY_females_1year, PY_males_1year,
                   simulated_dementia_logIRRs_data, 
                   simulated_dementia_logIRRs_data_1yr,
                   prop_dem_random, prop_dem_random_W, prop_dem_random_M, 
                   prop_dem_random_by_age, prop_dem_random_by_age_W, 
                   prop_dem_random_by_age_M,prop_dem_Ci, prop_dem_Ci_W, 
                   prop_dem_Ci_M, prop_dem_Ci_by_age, prop_dem_Ci_by_age_W, 
                   prop_dem_Ci_by_age_M, prop_dem_both, prop_dem_both_W, 
                   prop_dem_both_M, prop_dem_both_by_age, 
                   prop_dem_both_by_age_W, prop_dem_both_by_age_M)
  
  names(results_vec) <- 
    c("num_obs_baseline", "num_females_baseline", 
      "num_males_baseline", 
      variable_names$p_alive_varnames[1:num_tests], 
      variable_names$p_alive_females_varnames[1:num_tests], 
      variable_names$p_alive_males_varnames[1:num_tests],
      variable_names$cp_alive_varnames[1:num_tests], 
      variable_names$cp_alive_females_varnames[1:num_tests], 
      variable_names$cp_alive_males_varnames[1:num_tests],
      variable_names$mortality_logHR_varnames[1:num_tests],
      variable_names$dem_inc_rate_varnames[1:num_tests], 
      variable_names$dem_inc_rate_females_varnames[1:num_tests], 
      variable_names$dem_inc_rate_males_varnames[1:num_tests], 
      variable_names$inc_cases_females_varnames[1:num_tests], 
      variable_names$inc_cases_males_varnames[1:num_tests], 
      variable_names$PY_females_varnames[1:num_tests], 
      variable_names$PY_males_varnames[1:num_tests],
      variable_names_1year$inc_cases_females_varnames, 
      variable_names_1year$inc_cases_males_varnames,
      variable_names_1year$PY_females_varnames, 
      variable_names_1year$PY_males_varnames, 
      t(
        cbind(
          na.omit(variable_names$logIRR_varnames),
          na.omit(variable_names$logIRR_SE_varnames),
          na.omit(variable_names$logIRR_95CI_Lower_varnames),
          na.omit(variable_names$logIRR_95CI_Upper_varnames),
          na.omit(variable_names$logIRR_95CI_Coverage_varnames))) %>%
        as.vector(),
      t(
        cbind(
          variable_names_1year$logIRR_varnames,
          variable_names_1year$logIRR_SE_varnames,
          variable_names_1year$logIRR_95CI_Lower_varnames,
          variable_names_1year$logIRR_95CI_Upper_varnames,
          variable_names_1year$logIRR_95CI_Coverage_varnames)) %>%
        as.vector(),
      "prop_dem_random", "prop_dem_random_W", "prop_dem_random_M", 
      variable_names$prop_dem_random_by_age[1:num_tests], 
      variable_names$prop_dem_random_W_by_age[1:num_tests], 
      variable_names$prop_dem_random_M_by_age[1:num_tests],
      "prop_dem_Ci", "prop_dem_Ci_W", "prop_dem_Ci_M", 
      variable_names$prop_dem_Ci_by_age[1:num_tests], 
      variable_names$prop_dem_Ci_W_by_age[1:num_tests], 
      variable_names$prop_dem_Ci_M_by_age[1:num_tests], 
      "prop_dem_both", "prop_dem_both_W", "prop_dem_both_M", 
      variable_names$prop_dem_both_by_age[1:num_tests], 
      variable_names$prop_dem_both_W_by_age[1:num_tests], 
      variable_names$prop_dem_both_M_by_age[1:num_tests]
    )
  
  #---- Return ----
  return(results_vec)
}






  