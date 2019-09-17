#---- Package loading, options ----
if (!require("pacman")) 
  install.packages("pacman", repos='http://cran.us.r-project.org')

p_load("tidyverse", "MASS", "reshape", "magrittr", "here")

options(scipen = 999) #Standard Notation
options(digits = 6)   #Round to 6 decimal places
options(warn = -1)    #Suppress warnings

#---- Source files ----
source(here("RScripts", "linear_model", "variable_names.R"))
source(here("RScripts", "linear_model", "create_ages.R"))
source(here("RScripts", "linear_model", "gen_random_effects.R"))
source(here("RScripts", "linear_model", "cognitive_function_model.R"))
source(here("RScripts", "linear_model", "survival_times.R"))
source(here("RScripts", "linear_model", "dementia_onset.R"))
source(here("RScripts", "linear_model", "compare_survtime_timetodem.R"))

#---- The small batch data generation function ----
small_batch_gen <- function(num_obs){
  #---- Create a blank dataset ----
  obs <- matrix(NA, nrow = num_obs, ncol = length(column_names))
  colnames(obs) <- column_names
  
  #---- Generating IDs, female, U ----
  obs[, "id"] <- seq(from = 1, to = num_obs, by = 1)
  obs[, "female"] <- rbinom(num_obs, size = 1, prob = pfemale)
  obs[, "U"] <- rnorm(num_obs, mean = 0, sd = 1)
  
  #---- Generating age data ----
  #Creating ages at each timepoint j
  ages = matrix(seq(50, 95, by = 5), nrow = 1)
  obs[, variable_names$age_varnames] <- create_ages(ages, num_obs)
  
  #---- Generating centered age data ----
  #Creating baseline-mean-centered ages at each timepoint j
  obs[, variable_names$agec_varnames] <- 
    obs[, variable_names$age_varnames] - mean(age0)
  
  #---- Generating "true" cognitive function Cij ----
  #Refer to Manuscript/manuscript_equations.pdf for equation
  
  #Generating random terms, slope, intercept, and overall variation
  obs <- gen_random_effects(obs, cij_var0, cij_var1, cij_cov, cij_var3)
  
  #Calculating Cij for each individual
  obs[, variable_names$Cij_varnames] <- cog_func(cij_knots, cij_slopes, obs)
  
  #Check for dementia at baseline
  obs <- obs[obs[, "Ci0"] > dem_cut, ]
  
  #Calculate observed slopes
  for(i in 1:9){
    obs[, variable_names$cij_slopeij_varnames[i]] <- 
      (obs[, variable_names$Cij_varnames[i + 1]] - 
         obs[, variable_names$Cij_varnames[i]])/int_time
  }
  
  #---- Generating uniform random variables per interval for Sij ----
  obs[, variable_names$rij_varnames[1:num_tests]] <- 
    replicate(num_tests, runif(nrow(obs), min = 0, max = 1))
  
  obs[, "death0"] <- 0
  
  #---- Transpose the matrix for subsequent calculations ----
  obs = t(obs)

  #---- Calculating Sij and survival times for each individual ----
  obs[variable_names$Sij_varnames[1:num_tests], ] <- survival(obs)
  obs["survtime", ] <- colSums(obs[variable_names$Sij_varnames[1:num_tests], ], 
                               na.rm = TRUE)
  
  #---- Calculating death data for each individual ----
  #Indicator of 1 means the individual died in that interval
  #NAs mean the individual died in a prior interval
  obs[variable_names$deathij_varnames[1:num_tests], ] <- 
    (obs[variable_names$Sij_varnames[1:num_tests], ] < int_time)*1 
  
  obs["study_death", ] <- 
    colSums(obs[variable_names$deathij_varnames[1:num_tests], ], na.rm = TRUE) #Study death indicator
  
  obs["age_death", ] <- obs["age0", ] + obs["survtime", ]
  
  #---- Dementia indicators ----
  max_dem <- length(variable_names$dem_varnames)
  for(i in 1:ncol(obs)){
    below_dem <- min(which(obs[variable_names$Cij_varnames, i] < dem_cut))
    if(is.finite(below_dem)){
      obs[variable_names$dem_varnames[1:(below_dem - 1)], i] <- 0
      obs[variable_names$dem_varnames[below_dem:max_dem], i] <- 1
      obs["dem_wave", i] <- (below_dem - 1)
      obs["dem", i] <- 1
    } else{
      obs[variable_names$dem_varnames[1:max_dem], i] <- 0
      obs["dem", i] <- 0
    }
  }
  
  #---- Dementia calculations ----
  obs["timetodem", ] <- dem_onset(obs, dem_cut)
  obs <- compare_survtime_timetodem(obs)
  obs["ageatdem", ] <- obs["age0", ] + obs["timetodem", ] #Age at dementia diagnosis
  
  #---- Censor Cij and dem data ----
  obs <- survival_censor(obs)
  
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
  
  #Dementia status at death
  for(i in 1:ncol(obs)){
    if(obs["dem", i] == 1 & obs["timetodem", i] <= obs["survtime", i]){
      obs["dem_death", i] <- 1
    } else if(obs["study_death", i] == 1 &
              (obs["dem", i] == 0 | (obs["dem", i] == 1 &
                                     obs["timetodem", i] >
                                     obs["survtime", i]))){
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
  return(as.data.frame(t(obs)))
}

# data_gen <- function(num_obs){
#   batch_size = 1000
#   
#   if(num_obs%%batch_size != 0){
#     stop(paste("Number of runs must be a multiple of ", batch_size, "."))
#   }
#   
#   data <- replicate(num_obs/1000, small_batch_gen(batch_size))
#   data <- apply(data, 1, function(x) t(x))
#   data <- matrix(unlist(data), ncol = length(column_names), byrow = FALSE) %>%
#     as.data.frame() %>% set_colnames(column_names)
#   data[, 1] <- seq(from = 1, to = nrow(data), by = 1)
#   
#   return(data)
# }

