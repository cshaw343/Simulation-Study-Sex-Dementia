#*******************************************************************************
# This is the function that performs the search for appropriate slope values in 
# the model for Cij
#
# This function is attempting to find slopes that will create dementia incidence
# rates that match the Supplemental Table 2 in Inequalities in dementia 
# incidence between six racial and ethnic groups over 14 years 
# (Mayeda et al, 2016)
#*******************************************************************************

#---- Package Loading and Options ----
if (!require("pacman")) 
  install.packages("pacman", repos='http://cran.us.r-project.org')

p_load("tidyverse", "magrittr")

#Suppress warnings
options(warn = -1)

#---- Source Files ----
source("RScripts/dementia_incidence2000-2013.R")
source("RScripts/sex-dementia_sim_parA.R")
source("RScripts/variable_names.R")
source("RScripts/cognitive_function_model.R")
source("RScripts/survival_times.R")
source("RScripts/dementia_onset.R")

#---- Data generation ----
generate_data <- function(num_obs){
  #---- Generating IDs, sex, U ----
  obs <- tibble("id" = seq(from = 1, to = num_obs, by = 1),
                "sex" = rbinom(num_obs, size = 1, prob = psex), 
                "U" = rnorm(num_obs, mean = 0, sd = 1), 
                "age0" = rep(50, num_obs))
  
  #---- Generating age data ----
  #Creating ages at each timepoint j
  ages <- as_tibble(matrix(NA, nrow = num_obs, 
                           ncol = length(nrow(variable_names))))
  for(j in 1:nrow(variable_names)){
    if(j == 1){
      ages[, j] = obs$age0 #Creates column of baseline ages
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
  #Covariance matrix for random slope and intercept terms
  cij_slope_int_cov <- matrix(c(cij_var0, cij_cov, cij_cov, cij_var1), 
                              nrow = 2, byrow = TRUE)
  
  #Generate random terms for each individual
  cij_slope_int_noise <- as_tibble(mvrnorm(n = num_obs, mu = rep(0, 2), 
                                           Sigma = cij_slope_int_cov)) %>% 
    set_colnames(., c("z0i", "z1i"))
  obs %<>% bind_cols(., cij_slope_int_noise)
  
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
  
  return(obs)
}

#---- Function that we're trying to optimize ----
dem_rate_1000py <- function(SLOPES_DEM_CUTS, dem_rate, index){
  obs <- generate_data(num_obs)
  
  #---- Generating assessment timepoint data ----
  visit_times <- seq(from = 0, to = int_time*index, by = int_time)

  #---- Generating Cij data ----
  #Calculating Cij for each individual
  #Store Cij values and slope values for each assessment
  compute_Cij <- cog_func(cij_knots[1:(index - 1)], 
                          SLOPES_DEM_CUTS[1:index], obs)
  Cij <- as.data.frame(compute_Cij$Cij) %>% 
    set_colnames(., variable_names$Cij_varnames[1:(index + 1)])
  cij_slopeij <- as.data.frame(compute_Cij$slopes) %>% 
    #remove the last variable name because there are only 10 intervals
    set_colnames(., 
                 head(variable_names$cij_slopeij_varnames, -1)[1:(index)]) 
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
    set_colnames(head(variable_names$Sij_varnames, -1)[1:index])
  obs %<>% bind_cols(., Sij)
  
  #---- Calculating death data for each individual ----
  #Compute death indicator for each interval
  #Change death indicators to 1 after death
  deathij <- as.tibble((Sij < int_time)*1) %>% 
    set_colnames(head(variable_names$deathij_varnames, -1)[1:index])
  
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
  
  #---- Create a competing risk outcome ----
  demij <- obs %>% dplyr::select(variable_names$Cij_varnames[1:(index + 1)]) 
  dem_cuts <- matrix(SLOPES_DEM_CUTS[11:(11 + index)], 
                     nrow = nrow(demij), ncol = (index + 1), 
                     byrow = TRUE)
  demij <- (demij < dem_cuts)*1 
  demij %<>% as.data.frame() %>% 
    set_colnames(variable_names$dem_varnames[1:(index + 1)]) 
  
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
           "timetodem" = 
             dem_onset(., 
                       dem_cuts = SLOPES_DEM_CUTS[11:length(SLOPES_DEM_CUTS)]), #Time to dementia diagnosis
           "ageatdem" = age0 + timetodem, #Age at dementia diagnosis
           "dem_death" = case_when(dem == 1 & timetodem <= survtime ~ 1, 
                                   study_death == 1 & 
                                     (dem == 0 | 
                                        (dem == 1 & timetodem > survtime)) ~ 
                                     2)) %>%  #Dementia status at death
    mutate_at("dem_death", funs(replace(., is.na(.), 0))) %>%
    mutate("timetodem_death" = if_else(dem == 1, pmin(timetodem, survtime),
                                       survtime),
           "ageatdem_death" = age0 + timetodem_death,
           "dem_alive" = case_when(dem_death == 1 ~ 1,
                                   TRUE ~ 0))
  
  #---- Censor Cij and dementia data based on death data ----
  for(i in 1:num_obs){
    if(obs[i, "study_death"] == 1){
      death_int <- (min(which(deathij[i, ] == 1)))
      Cs <- c(variable_names$Cij_varnames[(death_int + 1):(index + 1)])
      Ds <- c(variable_names$dem_varnames[(death_int + 1):(index + 1)]) 
      obs[i, Cs] <- NA
      obs[i, Ds] <- NA
    }
  }
  
  # #---- Standardize Cij data ----
  # #**We're not sure about this yet**
  # std_Cij <- obs %>% dplyr::select(variable_names$Cij_varnames) %>% 
  #   map_df(~(. - mean(., na.rm = TRUE))/sd(., na.rm = TRUE)) %>%
  #   set_colnames(variable_names$std_Cij_varnames)
  # 
  # obs %<>% bind_cols(., std_Cij)
  
  #---- Compute person years ----
  contributed <- (obs$timetodem_death)%%5
  cases_py1000 <- vector(length = index)
  for(j in 1:index){
    last_test = j - 1
    last_dem_wave <- paste("dem", last_test, sep = "")
    this_dem_wave <- paste("dem", j, sep = "")
    last_death_wave <- paste("death", (last_test - 1), "-", last_test, sep = "")
    this_death_wave <- paste("death", last_test, "-", j, sep = "")
    if(j == 1){
      dem_data <- obs %>% 
        dplyr::select(c(last_dem_wave, this_dem_wave, this_death_wave))
      dem_data %<>% cbind(., contributed)
      dem_data %<>% 
        filter(!! as.name(last_dem_wave) == 0) %>%
        mutate("PY" = case_when(!! as.name(this_dem_wave) == 0 & 
                                  !! as.name(this_death_wave) == 0 ~ 5,
                                TRUE ~ contributed))
    } else {
      dem_data <- obs %>% 
        dplyr::select(c(last_dem_wave, this_dem_wave, last_death_wave, 
                        this_death_wave))
      dem_data %<>% cbind(., contributed)
      dem_data %<>% 
        filter(!! as.name(last_dem_wave) == 0 & 
                 !! as.name(last_death_wave) == 0) %>%
        mutate("PY" = case_when(!! as.name(this_dem_wave) == 0 & 
                                  !! as.name(this_death_wave) == 0 ~ 5,
                                TRUE ~ contributed))
    }
    
    cases_py1000[j] = 1000*
      sum(dem_data[, this_dem_wave], na.rm = TRUE)/sum(dem_data$PY)
  }
  return(abs(cases_py1000[index] - dem_rate))
}

#---- Function that looks for the slopes ----
find_slopes <- function(num_obs, dem_rate){
  best_slopes_dem_cuts <- c(rep(0, 10), rep(-1.5, 11))
  start = 3 #Represents age right before we have rate data (age 70)
  slot = start + 1 #Picks the slot we are trying to optimize
  incident_rate_diffs <- vector(length = length(dem_rate))
  for(i in 1:length(dem_rate)){ 
    if(i == 1){
      opt <- optim(par = best_slopes_dem_cuts, fn = dem_rate_1000py, 
                   dem_rate = dem_rate[[i]], index = slot, 
                   upper = best_slopes_dem_cuts, 
                   lower = c(rep(-0.15, 7), rep(-0.3, 3), rep(-3, 11))) 
      
      best_slopes_dem_cuts[1:(start + i)] <- opt$par[1:(start + i)]
      best_slopes_dem_cuts[11:(start + i + 10)] <- opt$par[11:(start + i + 10)]
      
    } else if((length(best_slopes_dem_cuts) - (start + (i - 1)) > 13)) {
      opt <- optim(par = best_slopes_dem_cuts, fn = dem_rate_1000py, 
                   dem_rate = dem_rate[[i]], index = slot,
                   upper = best_slopes_dem_cuts, 
                   lower = c(best_slopes_dem_cuts[1:(start + (i - 1))], 
                             rep(-0.15, 
                                 (10 - (start + (i - 1)))), 
                             best_slopes_dem_cuts[11:(start + (i - 1) + 10)], 
                             rep(-10, 11 - (start + (i - 1))))) 
      
      best_slopes_dem_cuts[start + i] <- opts$par[10 + start + i]
      best_slopes_dem_cuts[10 + start + i] <- opts$par[10 + start + i]
    } else {
      opt <- optim(par = best_slopes_dem_cuts, fn = dem_rate_1000py, 
                   dem_rate = dem_rate[[i]], index = slot,
                   upper = best_slopes_dem_cuts, 
                   lower = c(best_slopes_dem_cuts[1:(start + (i - 1))], 
                             rep(-0.30, 
                                 (10 - (start + (i - 1)))), 
                             best_slopes_dem_cuts[11:(start + (i - 1) + 10)], 
                             rep(-20, 11 - (start + (i - 1)))))
      
      best_slopes_dem_cuts[start + i] <- opts$par[10 + start + i]
      best_slopes_dem_cuts[10 + start + i] <- opts$par[10 + start + i]
    }
    incident_rate_diffs[i] <- opt$value
    slot = slot + 1
  }
    return(list("best_slopes_dem_cuts" = best_slopes_dem_cuts, 
                "incident_rates" = incident_rates))
  }

#---- Perform the search ----
fingers_crossed <- find_slopes(num_obs = 10, dem_rate = dem_rates_whites$Rate)


