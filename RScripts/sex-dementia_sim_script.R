#---- Package loading, options, seed ----
if (!require("pacman")) 
  install.packages("pacman", repos='http://cran.us.r-project.org')

p_load("tidyverse", "MASS", "reshape", "magrittr")

options(scipen = 999) #Standard Notation
options(digits = 3)   #Round to 3 decimal places
options(warn = -1)    #Suppress warnings

set.seed(10789)

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
  obs %<>% bind_cols(., ages)
  
  #---- Generating centered age data ----
  #Creating centered ages at each timepoint j
  c_ages <- as_tibble(ages - mean(age0)) %>% 
    set_colnames(., variable_names$agec_varnames)
  obs %<>% bind_cols(., c_ages)
  
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
    for(i in 1:(num_tests + 1)){
      cij_slope_int_noise <- as_tibble(mvrnorm(n = num_obs, mu = rep(0, 2), 
                                               Sigma = 
                                                 cij_slope_int_cov[[i]])) %>% 
      set_colnames(., c(paste0("z0_", (i - 1), "i"), 
                        paste0("z1_", (i - 1), "i")))
      obs %<>% bind_cols(., cij_slope_int_noise)
    }
    
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
    obs %<>% bind_cols(., eps)
  
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
  
  #---- Calculating death data for each individual ----
  #Compute death indicator for each interval
  #Change death indicators to 1 after death
  deathij <- as.tibble((Sij < int_time)*1) %>% 
    set_colnames(head(variable_names$deathij_varnames, -1))
  
  for(i in 1:nrow(deathij)){
    death <- min(which(deathij[i, ] == 1))
    if(is.finite(death)){
      deathij[i, death:ncol(deathij)] = 1 
    }
  }
  obs %<>% bind_cols(., deathij) %>% 
    mutate("study_death" = (rowSums(deathij) > 0)*1) #Study death indicator
  
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
  
  #---- Censor Cij based on death data ----
  for(i in 1:num_obs){
    if(obs[i, "study_death"] == 1){
      death_int <- (min(which(deathij[i, ] == 1)))
      Cs <- c(variable_names$Cij_varnames[(death_int + 1):nrow(variable_names)])
      obs[i, Cs] <- NA
    }
  }
  
  # #---- Standardize Cij values ----
  # std_Cij <- obs %>% dplyr::select(variable_names$Cij_varnames) %>%
  #   map_df(~(. - mean(., na.rm = TRUE))/sd(., na.rm = TRUE)) %>%
  #   set_colnames(variable_names$std_Cij_varnames)
  # 
  # obs %<>% bind_cols(., std_Cij)
  
  #---- Create a competing risk outcome ----
  dem_cuts_mat <- matrix(dem_cuts, nrow = nrow(obs), ncol = length(dem_cuts), 
                         byrow = TRUE)
  
  demij <- obs %>% dplyr::select(variable_names$Cij_varnames) %>% 
    set_colnames(variable_names$dem_varnames)
  demij <- (demij < dem_cuts_mat)*1
  
  dem_wave <- vector(length = nrow(obs))  #Wave at which dementia was diagnosed
  for(i in 1:nrow(demij)){
    dem_time <- min(which(demij[i, ] == 1))
    if(is.finite(dem_time)){
      demij[i, dem_time:ncol(demij)] = 1  #Changes dementia indicators to 1 after initial diagnosis
      dem_wave[i] = dem_time - 1          #Fills in wave of dementia diagnosis
    } else{
      dem_wave[i] = NA
    }
  }
  
  obs %<>% cbind(., demij) %>%
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
  
  #---- Censor demij based on death data ----
  for(i in 1:num_obs){
    if(obs[i, "study_death"] == 1){
      death_int <- (min(which(deathij[i, ] == 1)))
      dems <- c(variable_names$dem_varnames[(death_int + 1):
                                              nrow(variable_names)])
      obs[i, dems] <- NA
    }
  }
  

  return("obs" = obs)
}



