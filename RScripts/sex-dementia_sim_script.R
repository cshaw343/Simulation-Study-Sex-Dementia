#---- Package loading, options ----
if (!require("pacman")) 
  install.packages("pacman", repos='http://cran.us.r-project.org')

p_load("tidyverse", "MASS", "reshape", "magrittr")

options(scipen = 999) #Standard Notation
options(digits = 6)   #Round to 6 decimal places
options(warn = -1)    #Suppress warnings

#---- Source files ----
source("RScripts/cognitive_function_model.R")
source("RScripts/survival_times.R")
source("RScripts/dementia_onset.R")

#---- The simulation function ----
sex_dem_sim <- function(){
  #---- Generating IDs, sex, U ----
  obs <- tibble("id" = seq(from = 1, to = num_obs, by = 1), 
                "sex" = rbinom(num_obs, size = 1, prob = psex), 
                "U" = rnorm(num_obs, mean = 0, sd = 1))
  
  #---- Generating age data ----
  #Creating ages at each timepoint j
  ages <- as_tibble(matrix(NA, nrow = num_obs, 
                           ncol = length(nrow(variable_names))))
  for(j in 1:nrow(variable_names)){
    if(j == 1){
      ages[, j] = age0 #Creates column of baseline ages
    } else ages[, j] = ages[, (j-1)] + int_time #Creates ages at following timepoints
  }
  colnames(ages) <- variable_names$age_varnames
  
  #---- Generating centered age data ----
  #Creating centered ages at each timepoint j
  c_ages <- as_tibble(ages - mean(age0)) %>% 
    set_colnames(., variable_names$agec_varnames)
  
  #---- Generating "true" cognitive function Cij ----
  #Refer to Manuscript/manuscript_equations.pdf for equation
  
  #Generating random terms for slope and intercept
    #Covariance matrices for random slope and intercept terms
    cij_slope_int_cov <- lapply(1:(num_tests + 1), 
                                function(x) matrix(NA, nrow = 2, ncol = 2))
    for(i in 1:(num_tests + 1)){
      cij_slope_int_cov[[i]] <- matrix(c(cij_var0, cij_cov, 
                                         cij_cov, cij_var1[i]), 
                                       nrow = 2, byrow = TRUE)
    }
    
    #Generate random terms for each individual
    cij_slope_int_noise <- matrix(nrow = num_obs, ncol = 2*(num_tests + 1)) 
    for(i in 1:(num_tests + 1)){
      noise <- mvrnorm(n = num_obs, mu = rep(0, 2), 
                       Sigma = cij_slope_int_cov[[i]]) 
      cij_slope_int_noise[ , (2*i - 1):(2*i)] <- noise
    }
    
    noise_names <- vector(length = 2*(num_tests + 1))
    for(i in 1:(num_tests + 1)){
      noise_names[(2*i - 1):(2*i)] <- c(paste0("z0_", (i - 1), "i"), 
                                        paste0("z1_", (i - 1), "i"))
    }
    
    cij_slope_int_noise %<>% as.data.frame() %>% set_colnames(noise_names)
    
  #Generating noise term (unexplained variance in Cij) for each visit
    #Creating AR(1) correlation matrix
    num_visits = num_tests + 1
    powers <- abs(outer(1:(num_visits), 1:(num_visits), "-")) #Exponents for autoregressive terms
    corr <- sqrt(cij_var3)*(cij_r1^powers)                    #Correlation matrix
    S <- diag(rep(sqrt(cij_var3)), nrow(corr))                #Diagonal matrix of SDs
    cij_cov_mat <- S%*%corr%*%S                               #Covariance matrix
  
    #Generating noise terms
    eps <- as_tibble(mvrnorm(n = num_obs, 
                            mu = rep(0, num_visits), Sigma = cij_cov_mat)) %>%
      set_colnames(., variable_names$eps_varnames)
  
  #First data-binding
  
  #Calculating Cij for each individual
  #Store Cij values and slope values for each assessment
  compute_Cij <- cog_func(cij_knots, cij_slopes, obs)
  Cij <- as.data.frame(compute_Cij$Cij) %>% 
    set_colnames(., variable_names$Cij_varnames)
  cij_slopeij <- as.data.frame(compute_Cij$slopes) %>% 
    #remove the last variable name because there are only 10 intervals
    set_colnames(., head(variable_names$cij_slopeij_varnames, -1)) 
  obs %<>% bind_cols(., Cij, cij_slopeij)
  
  #---- Generate survival time for each person ----
  #Refer to Manuscript/manuscript_equations.pdf for equation
  
  #---- Generating uniform random variables per interval for Sij ----
  rij <- as_tibble(replicate(num_tests, 
                              runif(num_obs, min = 0, max = 1))) %>%
    #remove the last variable name because there are only 10 intervals
    set_colnames(head(variable_names$rij_varnames, -1))
  obs %<>% bind_cols(., rij)
  
  #---- Calculating Sij for each individual ----
  #Store Sij values
  Sij <- as.data.frame(survival(obs, lambda)) %>% 
    set_colnames(head(variable_names$Sij_varnames, -1))
  obs %<>% bind_cols(., Sij)
  
  #---- Survival censoring matrix ----
  censor <- Sij/Sij
  
  #---- Calculating death data for each individual ----
  #Indicator of 1 means the individual died in that interval
  #NAs mean the individual died in a prior interval
  deathij <- (Sij < int_time)*1 
  
  deathij %<>% as.data.frame() %>% 
    set_colnames(head(variable_names$deathij_varnames, -1))
  
  obs %<>% bind_cols(., deathij) %>% 
    mutate("death0" = 0) %>%
    mutate("study_death" = rowSums(deathij, na.rm = TRUE)) #Study death indicator
  
  #Compute overall survival times
  survtime <- vector(length = num_obs)
  survtime[which(obs$study_death == 0)] = num_tests*int_time
  for(i in 1:length(survtime)){
    if(survtime[i] == 0){
      death_int <- min(which(deathij[i, ] == 1)) 
      survtime[i] = int_time*(death_int - 1) + Sij[i, death_int]
    } 
  }
  obs %<>% mutate("survtime" = survtime, 
                  "age_death" = age0 + survtime) #Age at death
  
  # #---- Standardize Cij values ----
  # std_Cij <- obs %>% dplyr::select(variable_names$Cij_varnames) %>%
  #   map_df(~(. - mean(., na.rm = TRUE))/sd(., na.rm = TRUE)) %>%
  #   set_colnames(variable_names$std_Cij_varnames)
  # 
  # obs %<>% bind_cols(., std_Cij)
  
  #---- Create a competing risk outcome ----
  dem_cuts_mat <- matrix(dem_cuts, nrow = nrow(obs), ncol = length(dem_cuts), 
                         byrow = TRUE)
  
  demij <- (Cij < dem_cuts_mat)*1 
  
  for(i in 1:nrow(demij)){
    dem_int <- min(which(demij[i, ] == 1))
    if(is.finite(dem_int)){
      demij[i, dem_int:ncol(demij)] = 1 #Changes dementia indicator to 1 after dementia diagnosis
    }
  }
  
  demij %<>% as.data.frame() %>% 
    set_colnames(variable_names$dem_varnames)
  
  obs %<>% cbind(., demij) #%>% filter(`dem0` == 0)
  
  #---- Censor Cij, Sij, and demij based on death data ----
  filtered_deathij <- obs %>% #Need filtered death because we got rid of people with dementia at baseline
    dplyr::select(head(variable_names$deathij_varnames, -1))
  for(i in 1:nrow(obs)){
    if(obs[i, "study_death"] == 1){
      death_slot <- (min(which(filtered_deathij[i, ] == 1)))
      death_varname <- variable_names$deathij_varnames[death_slot]
      death_int <- str_remove(death_varname, "death")
      end_int <- str_sub(death_int, str_length(death_int)) 
      
      first_C <- paste0("Ci", end_int)
      first_dem <- paste0("dem", death_int)
      
      #Things to censor
      Cs <- 
        variable_names$Cij_varnames[
          c(which(variable_names$Cij_varnames == first_C):
              length(variable_names$Cij_varnames))]
      dems <- 
        variable_names$dem_varnames[
          c(which(variable_names$dem_varnames == first_dem):
              length(variable_names$dem_varnames))]
      
      if(death_int != "8-9"){
        end_int_num <- end_int %>% as.numeric()
        next_end_num <- end_int_num + 1
        first_S <- paste0("survtime", end_int_num, "-", next_end_num)
        
        Ss <- 
          variable_names$Sij_varnames[
            c(which(variable_names$Sij_varnames == first_S):
                (length(variable_names$Sij_varnames) - 1))]
        #Censoring
        obs[i, c(Ss)] <- NA
      }
      #Censoring
      obs[i, c(Cs, dems)] <- NA
      demij[i, dems] <- NA
    }
  }
  
  filtered_dem <- obs %>% dplyr::select(variable_names$dem_varnames)
  dem_wave <- vector(length = nrow(filtered_dem))
  for(i in 1:length(dem_wave)){
    if(is.finite(min(which(filtered_dem[i, ] == 1)))){
      dem_wave[i] <- (min(which(filtered_dem[i, ] == 1)) - 1)
    } else {
      dem_wave[i] = NA
    }
  }
  
  obs %<>%
    mutate("dem_wave" = dem_wave) %>%
    mutate("dem" = (1 - is.na(dem_wave)), #Dementia diagnosis indicator
           "timetodem" = dem_onset(., dem_cuts),    #Time to dementia diagnosis
           "ageatdem" = age0 + timetodem, #Age at dementia diagnosis
           "dem_death" =                  #Dementia status at death
             case_when(dem == 1 & timetodem <= survtime ~ 1,
                       study_death == 1 &
                         (dem == 0 | (dem == 1 & timetodem > survtime)) ~
                         2)) %>%
    mutate_at("dem_death", funs(replace(., is.na(.), 0))) %>%
    mutate("timetodem_death" = if_else(dem == 1, pmin(timetodem, survtime),
                                       survtime),
           "ageatdem_death" = age0 + timetodem_death,
           "dem_alive" = case_when(dem_death == 1 ~ 1,
                                   TRUE ~ 0))
  
  #---- Contributed time ----
  contributed_time <- matrix(nrow = nrow(obs), ncol = num_tests)
  for(i in 1:nrow(contributed_time)){
    last_full_slot <- floor(obs[i, "timetodem_death"]/5)
    contributed_time[i, 1:last_full_slot] <- 5
    if(last_full_slot != 9){
      partial_slot <- last_full_slot + 1
      contributed_time[i, partial_slot] <- obs[i, "timetodem_death"]%%5
    }
  }
  
  contributed_time %<>% as.data.frame() %>% 
    set_colnames(head(variable_names$contributed_varnames, -1))
  
  obs %<>% cbind(., contributed_time)
  
  #---- Bind all data into one dataframe ----
  obs %<>% bind_cols(ages, c_ages)
  
  #---- Dementia Incidence Rates by Sex ----
  sim_rates_by_sex <- matrix(ncol = 9, nrow = 2)
  colnames(sim_rates_by_sex) <- head(variable_names$interval_ages, -1)
  rownames(sim_rates_by_sex) <- c("Female", "Male")
  
  female_data <- obs %>% filter(sex == 0)
  male_data <- obs %>% filter(sex == 1)
  
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
  
  return("logIRRs" = logIRRs)
}



