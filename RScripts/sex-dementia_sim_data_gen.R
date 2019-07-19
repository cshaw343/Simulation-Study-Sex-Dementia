#---- Package loading, options ----
if (!require("pacman")) 
  install.packages("pacman", repos='http://cran.us.r-project.org')

p_load("tidyverse", "MASS", "reshape", "magrittr", "here")

options(scipen = 999) #Standard Notation
options(digits = 6)   #Round to 6 decimal places
options(warn = -1)    #Suppress warnings

#---- Source files ----
source(here("RScripts", "variable_names.R"))
source(here("RScripts", "cognitive_function_model.R"))
source(here("RScripts", "survival_times.R"))
source(here("RScripts", "dementia_onset.R"))
source(here("RScripts", "compare_survtime_timetodem.R"))

#---- The small batch data generation function ----
data_gen <- function(num_obs){
  #---- Create a blank dataset ----
  obs <- matrix(NA, nrow = num_obs, ncol = length(column_names)) %>% 
    as.data.frame() %>% set_colnames(column_names)
  
  #---- Generating IDs, female, U ----
  obs$id <- seq(from = 1, to = num_obs, by = 1)
  obs$female <- rbinom(num_obs, size = 1, prob = pfemale)
  obs$U <- rnorm(num_obs, mean = 0, sd = 1)
  
  #---- Generating age data ----
  #Creating ages at each timepoint j
  ages = matrix(seq(50, 95, by = 5), nrow = 1)
  ones = matrix(1, nrow = num_obs, ncol = 1)
  obs[, variable_names$age_varnames] <- ones %*% ages
  
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
    if(i == 1){
      obs[, c("z0_i", paste0("z1_", (i - 1), "i"))] <- noise
    } else{
      obs[, paste0("z1_", (i - 1), "i")] <- noise[, 2]
    }
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
  obs[, variable_names$cij_slopeij_varnames[1:num_tests]] <- compute_Cij$slopes
  
  #---- Create a competing risk outcome ----
  dem_cuts_mat <- matrix(dem_cut, nrow = nrow(obs), 
                         ncol = length(variable_names$Cij_varnames), 
                         byrow = TRUE)
  
  obs[, variable_names$dem_varnames] <- 
    (obs[, variable_names$Cij_varnames] <= dem_cuts_mat)*1
  
  obs %<>% filter(dem0 == 0)
  
  #---- Generate survival time for each person ----
  #Refer to Manuscript/manuscript_equations.pdf for equation
  
  #---- Generating uniform random variables per interval for Sij ----
  obs[, variable_names$rij_varnames[1:num_tests]]<- 
    replicate(num_tests, runif(num_obs, min = 0, max = 1))
  
  #---- Transpose the matrix for subsequent calculations ----
  obs = t(obs)
  
  #---- Dementia indicators ----
  for(i in 1:ncol(obs)){
    dem_int <- min(which(obs[variable_names$dem_varnames, i] == 1))
    if(is.finite(dem_int)){
      obs["dem_wave", i] <- (dem_int - 1)
      first_censor <- min(which(is.na(obs[variable_names$dem_varnames, i])))
      if(dem_int < 10 & is.finite(first_censor)){
        obs[variable_names$dem_varnames[dem_int:(first_censor - 1)], i] <- 1 #Changes dementia indicator to 1 after dementia diagnosis
      }
      if(dem_int < 10 & !is.finite(first_censor)){
        last_1 <- length(variable_names$dem_varnames)
        obs[variable_names$dem_varnames[dem_int:last_1], i] <- 1 #Changes dementia indicator to 1 after dementia diagnosis
      }
    } else {
      obs["dem_wave", i] = NA
    }
  }
  
  #---- Calculating Sij for each individual ----
  #Store Sij values and survival time
  survival_data <- survival(obs)
  obs[variable_names$Sij_varnames[1:num_tests], ] <- survival_data$Sij
  obs["survtime", ] <- survival_data$survtimes
  
  #---- Calculating death data for each individual ----
  #Indicator of 1 means the individual died in that interval
  #NAs mean the individual died in a prior interval
  obs["death0", ] <- 0
  obs[variable_names$deathij_varnames[1:num_tests], ] <- 
    (obs[variable_names$Sij_varnames[1:num_tests], ] < int_time)*1 
  
  obs["study_death", ] <- 
    colSums(obs[variable_names$deathij_varnames[1:num_tests], ], na.rm = TRUE) #Study death indicator

  obs["age_death", ] <- obs["age0", ] + obs["survtime", ]
  
  # #---- Standardize Cij values ----
  # std_Cij <- obs %>% dplyr::select(variable_names$Cij_varnames) %>%
  #   map_df(~(. - mean(., na.rm = TRUE))/sd(., na.rm = TRUE)) %>%
  #   set_colnames(variable_names$std_Cij_varnames)
  # 
  # obs %<>% bind_cols(., std_Cij)
  
  #---- Survival censoring matrix ----
  censor <- (obs[variable_names$Sij_varnames[1:num_tests], ] == 5)*1
  censor[censor == 0] <- NA
  censor %<>% rbind(1, .)
  
  shifted_censor <- rbind(1, censor[1:(nrow(censor) - 1), ])
  
  #---- Censor Cij, slope, and dem data ----
  obs[variable_names$Cij_varnames, ] <- 
    obs[variable_names$Cij_varnames, ]*censor
  
  # obs[variable_names$cij_slopeij_varnames[1:9], ] <- 
  #   obs[variable_names$cij_slopeij_varnames[1:9], ]*shifted_censor[1:9, ]
  
  obs[variable_names$dem_varnames, ] <- 
    obs[variable_names$dem_varnames, ]*shifted_censor
  
  #---- Dementia calcs ----
  obs["dem", ] <- (1 - is.na(obs["dem_wave", ])) #Dementia diagnosis indicator
  obs["timetodem", ] <- dem_onset(obs, dem_cut) #Time to dementia diagnosis
  obs <- compare_survtime_timetodem(obs)
  obs["ageatdem", ] <- obs["age0", ] + obs["timetodem", ] #Age at dementia diagnosis
  
  #Dementia status at death
  for(i in 1:ncol(obs)){
    if(obs["dem", i] == 1 & obs["timetodem", i] <= obs["survtime", i]){
      obs["dem_death", i] <- 1
    } else if(obs["study_death", i] == 1 & 
              (obs["dem", i] == 0 | (obs["dem", i] == 1 & 
               obs["timetodem", i] > obs["survtime", i]))){
      obs["dem_death", i] <- 2
    } else {
      obs["dem_death", i] <- 0
    }
  }
  
  #Time to dem_death
  for(i in 1:ncol(obs)){
    if(obs["dem", i] == 0){
      obs["timetodem_death", i] <- obs["survtime", i]
    } else {
      obs["timetodem_death", i] <- min(obs["timetodem", i], obs["survtime", i])
    }
  }
  
  obs["ageatdem_death", ] <- obs["age0", ] + obs["timetodem_death", ]
  obs["dem_alive", obs["dem_death", ] == 1] <- 1
  obs["dem_alive", is.na(obs["dem_alive", ])] <- 0
  
  #---- Last Cij value ----
  for(i in 1:ncol(obs)){
    int_start = floor(obs["survtime", i]/5)
    if(int_start == 9){
      next
    }
    int_remain = obs["survtime", i] %% 5
    last_intercept <- paste0("Ci", int_start)
    last_slope <- paste0("cij_slope",int_start , "-", int_start + 1)
    obs["last_Cij", i] = obs[last_intercept, i] + obs[last_slope, i]*int_remain
  }
  
  #---- Contributed time ----
  for(i in 1:ncol(obs)){
    #5-year bands
    last_full_slot <- floor(obs["timetodem_death", i]/5)
    full_slots <- variable_names$contributed_varnames[1:last_full_slot]
    obs[full_slots, i] <- 5
    if(last_full_slot != 9){
      partial_slot <- last_full_slot + 1
      obs[variable_names$contributed_varnames[partial_slot], i] <- 
        obs["timetodem_death", i]%%5
    }
  }
  
  #---- Contributed time (1-year bands) ----
  for(i in 1:ncol(obs)){
    for(j in 1:num_tests){
      contributed_var <- variable_names$contributed_varnames[j]
      contributed_vars_1year_block <-
        variable_names_1year$contributed_varnames[(5*(j-1) + 1):(5*j)]
      
      if(is.na(obs[contributed_var, i])){
        break
      } else if(obs[contributed_var, i] == int_time){
        obs[contributed_vars_1year_block, i] <- 1
      } else {
        last_full_slot <- floor(obs[contributed_var, i])
        if(last_full_slot == 0){
          obs[contributed_vars_1year_block[1], i] <- obs[contributed_var, i]
        } else {
          full_slots <-
            contributed_vars_1year_block[1:last_full_slot]
          obs[full_slots, i] <- 1
          partial_slot <- last_full_slot + 1
          obs[contributed_vars_1year_block[partial_slot], i] <-
            (obs[contributed_var, i] - last_full_slot)
        }
      }
    }
  }
  
  #---- Dementia indicators (1-year bands) ----
  for(i in 1:ncol(obs)){
    for(j in 2:(num_tests + 1)){
      dem_var <- variable_names$dem_varnames[j]
      dem_vars_1year_block <-
        variable_names_1year$dem_varnames[(5*(j-2) + 1):(5*(j-1))]
      contributed_vars_1year_block <-
        variable_names_1year$contributed_varnames[(5*(j-2) + 1):(5*(j-1))]
      if(is.na(obs[dem_var, i])){
        break
      } else if(obs[dem_var, i] == 0){
        obs[dem_vars_1year_block, i] <- 0
      } else{
        obs[dem_vars_1year_block, i] <-
          (obs[contributed_vars_1year_block, i] < 1)*1
      }
    }
  }
  
  #---- Values to return ----
  return(t(obs))
}

