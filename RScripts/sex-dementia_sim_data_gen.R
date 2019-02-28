#---- Package loading, options ----
if (!require("pacman")) 
  install.packages("pacman", repos='http://cran.us.r-project.org')

p_load("tidyverse", "MASS", "reshape", "magrittr")

options(scipen = 999) #Standard Notation
options(digits = 6)   #Round to 6 decimal places
options(warn = -1)    #Suppress warnings

#---- Source files ----
source("RScripts/variable_names.R")
source("RScripts/cognitive_function_model.R")
source("RScripts/survival_times.R")
source("RScripts/dementia_onset.R")

#---- The simulation function ----
data_gen <- function(){
  #---- Create a blank dataset ----
  obs <- matrix(NA, nrow = num_obs, ncol = length(column_names)) %>% 
    as.data.frame() %>% set_colnames(column_names)
  
  #---- Generating IDs, sex, U ----
  obs$id <- seq(from = 1, to = num_obs, by = 1)
  obs$sex <- rbinom(num_obs, size = 1, prob = psex)
  obs$U <- rnorm(num_obs, mean = 0, sd = 1)
  
  #---- Generating age data ----
  #Creating ages at each timepoint j
  for(j in 1:length(variable_names$age_varnames)){
    if(j == 1){
      obs[, variable_names$age_varnames[j]] = age0 #Creates column of baseline ages
    } else 
      obs[, variable_names$age_varnames[j]] = 
        obs[, variable_names$age_varnames[j - 1]] + int_time #Creates ages at following timepoints
  }
  
  #---- Generating centered age data ----
  #Creating baseline-mean-centered ages at each timepoint j
  obs[, variable_names$agec_varnames] <- 
    obs[, variable_names$age_varnames] - mean(age0)
  
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
    noise <- mvrnorm(n = num_obs, mu = rep(0, 2), 
                     Sigma = cij_slope_int_cov[[i]]) 
    obs[, c(paste0("z0_", (i - 1), "i"), paste0("z1_", (i - 1), "i"))] <- noise
  }
  
  #Generating noise term (unexplained variance in Cij) for each visit
  #Creating AR(1) correlation matrix
  num_visits = num_tests + 1
  powers <- abs(outer(1:(num_visits), 1:(num_visits), "-")) #Exponents for autoregressive terms
  corr <- sqrt(cij_var3)*(cij_r1^powers)                    #Correlation matrix
  S <- diag(rep(sqrt(cij_var3)), nrow(corr))                #Diagonal matrix of SDs
  cij_cov_mat <- S%*%corr%*%S                               #Covariance matrix
  
  #Generating noise terms
  obs[, variable_names$eps_varnames] <- 
    mvrnorm(n = num_obs, mu = rep(0, num_visits), Sigma = cij_cov_mat)
  
  #Calculating Cij for each individual
  #Store Cij values and slope values for each assessment
  compute_Cij <- cog_func(cij_knots, cij_slopes, obs)
  obs[, variable_names$Cij_varnames] <- compute_Cij$Cij
  obs[, na.omit(variable_names$cij_slopeij_varnames)] <- compute_Cij$slopes
  
  #---- Generate survival time for each person ----
  #Refer to Manuscript/manuscript_equations.pdf for equation
  
  #---- Generating uniform random variables per interval for Sij ----
  obs[, na.omit(variable_names$rij_varnames)]<- 
    replicate(num_tests, runif(num_obs, min = 0, max = 1))
  
  #---- Calculating Sij for each individual ----
  #Store Sij values
  Sij <- as.data.frame(survival(obs, lambda, rij, c_ages, cij_slopeij, Cij)) %>% 
    set_colnames(head(variable_names$Sij_varnames, -1))
  
  #---- Survival censoring matrix ----
  censor <- Sij/Sij 
  censor %<>% cbind(1, .)
  
  #---- Calculating death data for each individual ----
  #Indicator of 1 means the individual died in that interval
  #NAs mean the individual died in a prior interval
  deathij <- (Sij < int_time)*1 
  
  deathij %<>% as.data.frame() %>% 
    set_colnames(head(variable_names$deathij_varnames, -1)) %>% 
    mutate("death0" = 0) %>% dplyr::select("death0", everything()) %>%
    mutate("study_death" = rowSums(deathij, na.rm = TRUE)) #Study death indicator
  
  #Compute overall survival times
  survtime <- vector(length = num_obs)
  survtime[which(deathij$study_death == 0)] = num_tests*int_time
  for(i in 1:length(survtime)){
    if(survtime[i] == 0){
      death_int <- min(which(deathij[i, ] == 1)) 
      survtime[i] = int_time*(death_int - 1) + Sij[i, (death_int - 1)]
    } 
  }
  
  survtime %<>% as.data.frame() %>% 
    mutate("age_death" = age0 + survtime) %>% #Age at death
    set_colnames(c("survtime", "age_death"))
  
  #---- Censor Cij data ----
  Cij <- Cij*censor
  
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
  
  #---- Censor Cij data ----
  demij <- demij*censor
  
  dem_wave <- vector(length = nrow(obs))
  for(i in 1:length(dem_wave)){
    if(is.finite(min(which(demij[i, ] == 1)))){
      dem_wave[i] <- (min(which(demij[i, ] == 1)) - 1)
    } else {
      dem_wave[i] = NA
    }
  }
  
  demij %<>%
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

  return("data" = obs)
}
  