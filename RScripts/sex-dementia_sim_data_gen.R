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

#---- The small batch data generation function ----
small_batch_gen <- function(small_batch_n){
  #---- Create a blank dataset ----
  obs <- matrix(NA, nrow = small_batch_n, ncol = length(column_names)) %>% 
    as.data.frame() %>% set_colnames(column_names)
  
  #---- Generating IDs, sex, U ----
  obs$id <- seq(from = 1, to = small_batch_n, by = 1)
  obs$sex <- rbinom(small_batch_n, size = 1, prob = psex)
  obs$U <- rnorm(small_batch_n, mean = 0, sd = 1)
  
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
    noise <- mvrnorm(n = small_batch_n, mu = rep(0, 2), 
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
    mvrnorm(n = small_batch_n, mu = rep(0, num_visits), Sigma = cij_cov_mat)
  
  #Calculating Cij for each individual
  #Store Cij values and slope values for each assessment
  compute_Cij <- cog_func(cij_knots, cij_slopes, obs)
  obs[, variable_names$Cij_varnames] <- compute_Cij$Cij
  obs[, na.omit(variable_names$cij_slopeij_varnames)] <- compute_Cij$slopes
  
  #---- Generate survival time for each person ----
  #Refer to Manuscript/manuscript_equations.pdf for equation
  
  #---- Generating uniform random variables per interval for Sij ----
  obs[, na.omit(variable_names$rij_varnames)]<- 
    replicate(num_tests, runif(small_batch_n, min = 0, max = 1))
  
  #---- Calculating Sij for each individual ----
  #Store Sij values and survival time
  survival_data <- survival(obs)
  obs[, na.omit(variable_names$Sij_varnames)] <- survival_data$Sij
  obs[, "survtime"] <- survival_data$survtimes
  
  #---- Calculating death data for each individual ----
  #Indicator of 1 means the individual died in that interval
  #NAs mean the individual died in a prior interval
  obs[, "death0"] <- 0
  obs[, na.omit(variable_names$deathij_varnames)] <- 
    (obs[, na.omit(variable_names$Sij_varnames)] < int_time)*1 
  
  obs[, "study_death"] <- 
    rowSums(obs[, na.omit(variable_names$deathij_varnames)], na.rm = TRUE) #Study death indicator

  obs[, "age_death"] <- obs[, "age0"] + obs[, "survtime"]
  
  # #---- Standardize Cij values ----
  # std_Cij <- obs %>% dplyr::select(variable_names$Cij_varnames) %>%
  #   map_df(~(. - mean(., na.rm = TRUE))/sd(., na.rm = TRUE)) %>%
  #   set_colnames(variable_names$std_Cij_varnames)
  # 
  # obs %<>% bind_cols(., std_Cij)
  
  #---- Create a competing risk outcome ----
  dem_cuts_mat <- matrix(dem_cuts, nrow = nrow(obs), ncol = length(dem_cuts), 
                         byrow = TRUE)
  
  obs[, variable_names$dem_varnames] <- 
    (obs[, variable_names$Cij_varnames] < dem_cuts_mat)*1
  
  obs %<>% filter(dem0 == 0)
  
  #---- Survival censoring matrix ----
  censor <- (obs[, na.omit(variable_names$Sij_varnames)] == 5)*1
  censor[censor == 0] <- NA
  censor %<>% cbind(1, .)
  
  #---- Censor Cij and dem data ----
  obs[, variable_names$Cij_varnames] <- 
    obs[, variable_names$Cij_varnames]*censor
  
  obs[, variable_names$dem_varnames] <- 
    obs[, variable_names$dem_varnames]*censor
  
  #---- Dementia indicators ----
  for(i in 1:nrow(obs)){
    dem_int <- min(which(obs[i, variable_names$dem_varnames] == 1))
    if(is.finite(dem_int)){
      obs[i, "dem_wave"] <- (dem_int - 1)
      first_censor <- min(which(is.na(obs[i, variable_names$dem_varnames])))
      if(dem_int < 9 & is.finite(first_censor)){
        obs[i, variable_names$dem_varnames[dem_int:(first_censor - 1)]] <- 1 #Changes dementia indicator to 1 after dementia diagnosis
      }
      if(dem_int < 9 & !is.finite(first_censor)){
        last_1 <- length(variable_names$dem_varnames)
        obs[i, variable_names$dem_varnames[dem_int:last_1]] <- 1 #Changes dementia indicator to 1 after dementia diagnosis
      }
    } else {
      obs[i, "dem_wave"] = NA
    }
  }
  
  #---- Dementia calcs ----
  obs[, "dem"] <- (1 - is.na(obs[, "dem_wave"])) #Dementia diagnosis indicator
  obs[, "timetodem"] <- dem_onset(obs, dem_cuts) #Time to dementia diagnosis
  obs[, "ageatdem"] <- obs[, "age0"] + obs[, "timetodem"] #Age at dementia diagnosis
  
  #Dementia status at death
  for(i in 1:nrow(obs)){
    if(obs[i, "dem"] == 1 & obs[i, "timetodem"] <= obs[i, "survtime"]){
      obs[i, "dem_death"] <- 1
    } else if(obs[i, "study_death"] == 1 & 
              (obs[i, "dem"] == 0 | (obs[i, "dem"] == 1 & 
               obs[i, "timetodem"] > obs[i, "survtime"]))){
      obs[i, "dem_death"] <- 2
    } else {
      obs[i, "dem_death"] <- 0
    }
  }
  
  #Time to dem_death
  for(i in 1:nrow(obs)){
    if(obs[i, "dem"] == 0){
      obs[i, "timetodem_death"] <- obs[i, "survtime"]
    } else {
      obs[i, "timetodem_death"] <- min(obs[i, "timetodem"], obs[i, "survtime"])
    }
  }
  
  obs[, "ageatdem_death"] <- obs[, "age0"] + obs[, "timetodem_death"]
  obs[obs[, "dem_death"] == 1, "dem_alive"] <- 1
  obs[is.na(obs[, "dem_alive"]), "dem_alive"] <- 0
  
  #---- Contributed time ----
  for(i in 1:nrow(obs)){
    last_full_slot <- floor(obs[i, "timetodem_death"]/5)
    obs[i, na.omit(variable_names$contributed_varnames)[1:last_full_slot]] <- 5
    if(last_full_slot != 9){
      partial_slot <- last_full_slot + 1
      obs[i, na.omit(variable_names$contributed_varnames)[partial_slot]] <- 
        obs[i, "timetodem_death"]%%5
    }
  }
  
  #---- Values to return ----
  return(obs)
}

#---- Data Generation ----
data_gen <- function(){
  small_batch_n <- 1000
  num_reps <- num_obs/small_batch_n
  
  if(num_obs %% small_batch_n != 0){
    stop(paste0("Number of observations must be a multiple of ", 
                small_batch_n, "."))
  }
  
  data <- replicate(num_reps, small_batch_gen(small_batch_n))
  data <- apply(data, 1, function(x) t(x))
  data_mat <- matrix(unlist(data), ncol = length(column_names), byrow = FALSE)
  
  data_mat %<>% as.data.frame() %>% set_colnames(names(data))
  data_mat[, 1] <- seq(from = 1, to = nrow(data_mat), by = 1)
  
  return(data_mat)
}






  