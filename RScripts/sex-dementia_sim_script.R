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

#---- Generating assessment timepoint data ----
visit_times <- seq(from = 0, to = int_time*num_tests, by = int_time)

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
  #Cij = b00 + z0i + bo1*sexi + b02*age0i + b03*Ui + (b10 + z1i + b11*sexi + 
  #b12*age0i + b13*Ui)*timej + epsilonij
  
  #---- Generating random terms for slope and intercept ----
  #Covariance matrix for random slope and intercept terms
  slope_int_cov <- matrix(c(var0, cov, cov, var1), nrow = 2, byrow = TRUE)
  
  #Generate random terms for each individual
  slope_int_noise <- as_tibble(mvrnorm(n = num_obs, mu = rep(0, 2), 
                                       Sigma = slope_int_cov)) %>% 
    set_colnames(., c("z0i", "z1i"))
  obs %<>% bind_cols(., slope_int_noise)
  
  #---- Generating noise term (unexplained variance in Cij) for each visit ----
  #Creating AR(1) correlation matrix
  num_visits = num_tests + 1
  powers <- abs(outer(1:(num_visits), 1:(num_visits), "-")) #Exponents for autoregressive terms
  corr <- sqrt(var3)*(r1^powers)                            #Correlation matrix
  S <- diag(rep(sqrt(var3)), nrow(corr))                    #Diagonal matrix of SDs
  cov_mat <- S%*%corr%*%S                                   #Covariance matrix
  
  #Generating noise terms
  eps <- as_tibble(mvrnorm(n = num_obs, 
                           mu = rep(0, num_visits), Sigma = cov_mat)) %>%
    set_colnames(., variable_names$eps_varnames)
  obs %<>% bind_cols(., eps)
  
  #---- Calculating Cij for each individual ----
  #Store Cij values and slope values for each assessment
  compute_Cij <- cog_func(slopes, obs)
  Cij <- as.data.frame(compute_Cij$Cij) %>% 
    set_colnames(., Cij_varnames)
  slopeij <- as.data.frame(compute_Cij$slopes) %>% 
    set_colnames(., slopeij_varnames)
  obs %<>% bind_cols(., Cij, slopeij)
  
  #---- Calculating mean Cij by sex ----
  mean_Cij <- obs %>% mutate_at("sex", as.factor) %>% group_by(sex) %>% 
    dplyr::select(sex, Cij_varnames) %>% summarise_all(mean) %>%
    set_colnames(mean_Cij_varnames)
  
  #---- Generate survival time for each person ----
  #Individual hazard functions
  #h(tij|x) = lambda*exp(g1*sexi + g2*ageij + g3*Ui + g4*sexi + 
  #g5*slopeij + g6Cij)
  #See Additional notes in README file
  
  #---- Generating uniform random variables per interval for Sij ----
  USij <- as_tibble(replicate(num_tests, 
                              runif(num_obs, min = 0, max = 1))) %>%
    set_colnames(USij_varnames)
  obs %<>% bind_cols(., USij)
  
  #---- Calculating Sij for each individual ----
  #Store Sij values
  Sij <- as.data.frame(survival(obs, lambda)) %>% set_colnames(Sij_varnames)
  obs %<>% bind_cols(., Sij)
  
  #---- Calculating death data for each individual ----
  #Compute death indicator for each interval
  #Change death indicators to 1 after death
  deathij <- as.tibble((Sij < int_time)*1) %>% set_colnames(deathij_varnames)
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
    death_int <- (min(which(deathij[i, ] == 1)) - 1)
    if(is.finite(death_int)){
      Cs <- vector(length = num_tests)
      for(j in death_int:num_tests){
        Cs[j] <- paste("Ci", j, sep = "")
      }
      Cs <- Cs[Cs != "FALSE"]
      obs[i, Cs] <- NA
    }
  }
  
  #---- Create a competing risk outcome ----
  demij <- obs %>% dplyr::select(Cij_varnames) %>% 
    mutate_all(funs((. < dem_cut)*1)) %>% set_colnames(dem_varnames)
  
  dem_wave <- vector(length = num_obs)  #Wave at which dementia was diagnosed
  for(i in 1:nrow(demij)){
    dem_time <- min(which(demij[i, ] == 1))
    if(is.finite(dem_time)){
      demij[i, dem_time:ncol(demij)] = 1  #Changes dementia indicators to 1 after initial diagnosis
      dem_wave[i] = dem_time - 1          #Fills in wave of dementia diagnosis
    } else{
      dem_wave[i] = NA
    }
  }
  
  obs %<>% bind_cols(., demij) %>% 
    mutate("dem_wave" = dem_wave) %>%
    mutate("dem" = (1 - is.na(dem_wave)), #Dementia diagnosis indicator
           "timetodem" = dem_onset(.),    #Time to dementia diagnosis
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
  
  #---- Compute person years ----
  contributed <- (obs$timetodem_death)%%5
  cases_py1000 <- vector(length = num_tests)
  for(j in 1:num_tests){
    last_test = j - 1
    last_wave <- paste("dem", last_test, sep = "")
    this_wave <- paste("dem", j, sep = "")
    dem_data <- demij %>% dplyr::select(c(last_wave, this_wave)) 
    dem_data %<>% cbind(., contributed)
    dem_data %<>% filter(!! as.name(last_wave) == 0) %>%
      mutate("PY" = case_when(!! as.name(this_wave) == 0 ~ 5, 
                              TRUE ~ contributed)) 
    cases_py1000[j] = 1000*
      sum(dem_data[, this_wave], na.rm = TRUE)/sum(dem_data$PY)  
  }
  
  #---- Edit for actual simulation ----
  #Comment out for simulation checks
  #obs <- obs %>% filter(dem_wave != 0)
  
  #Desired simulation return values    
  #return(list("mean_Cij" = mean_Cij))
  
  #Alternative return function for code checking
  return(list("obs" = obs, "mean_Cij" = mean_Cij, "dem_cases" = cases_py1000))
}



