#---- Package loading, options ----
if (!require("pacman")) 
  install.packages("pacman", repos='http://cran.us.r-project.org')

p_load("here", "tidyverse", "MASS", "optimParallel")

#---- Source files ----
source(here("RScripts", "sex-dementia_sim_parA.R"))
source(here("RScripts", "dementia_incidence_EURODEM_pooled.R"))
source(here("RScripts", "euro_life_tables.R"))
source(here("RScripts", "variable_names.R"))
source(here("RScripts", "cognitive_function_model.R"))
source(here("RScripts", "survival_times.R"))
source(here("RScripts", "dementia_onset.R"))
source(here("RScripts", "compare_survtime_timetodem.R"))

#For survival re-optimization
source(here("RScripts", "sex-dementia_sim_data_gen.R"))

#---- Objective Function ----
dem_inc_rate_match <- function(PARAMETER, which_opt,
                               batch_n, timepoint, 
                               dem_inc_data, cp_survival,
                               opt_cij_slopes, opt_cij_var1, opt_base_haz){
  
  #---- Plug in parameters to optimize for ----
  cij_slopes <- opt_cij_slopes
  cij_var1 <- opt_cij_var1
  
  if(which_opt == "slope"){
    cij_slopes[timepoint] <- PARAMETER
  } else {
    cij_var1[timepoint + 1] <- PARAMETER
  }
  
  #---- Create a blank dataset ----
  obs <- matrix(NA, nrow = batch_n, ncol = length(column_names)) %>% 
    as.data.frame() %>% set_colnames(column_names)
  
  #---- Generating IDs, sex, U ----
  obs$id <- seq(from = 1, to = batch_n, by = 1)
  obs$sex <- rbinom(batch_n, size = 1, prob = psex)
  obs$female <- 1 - obs$sex
  obs$U <- rnorm(batch_n, mean = 0, sd = 1)
  
  #---- Generating age data ----
  #Creating ages at each timepoint j
  ages = matrix(seq(50, 95, by = 5), nrow = 1)
  ones = matrix(1, nrow = batch_n, ncol = 1)
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
    noise <- mvrnorm(n = batch_n, mu = rep(0, 2), 
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
    mvrnorm(n = batch_n, mu = rep(0, num_visits), Sigma = cij_cov_mat)
  
  #Calculating Cij for each individual
  #Store Cij values and slope values for each assessment
  compute_Cij <- cog_func(cij_knots, cij_slopes, obs)
  obs[, variable_names$Cij_varnames] <- compute_Cij$Cij
  obs[, na.omit(variable_names$cij_slopeij_varnames)] <- compute_Cij$slopes
  
  #---- Create a competing risk outcome ----
  dem_cuts_mat <- matrix(dem_cut, nrow = nrow(obs), 
                         ncol = length(variable_names$Cij_varnames), 
                         byrow = TRUE)
  
  obs[, variable_names$dem_varnames] <- 
    (obs[, variable_names$Cij_varnames] < dem_cuts_mat)*1
  
  obs %<>% filter(dem0 == 0)
  
  #---- Generate survival time for each person ----
  #Refer to Manuscript/manuscript_equations.pdf for equation
  
  #---- Generating uniform random variables per interval for Sij ----
  obs[, na.omit(variable_names$rij_varnames)]<- 
    replicate(num_tests, runif(batch_n, min = 0, max = 1))
  
  #---- Survival optimization ----
  #Objective Function
  survival_match <- function(LAMBDA, obs, cp_survival){
    lambda <- opt_base_haz
    lambda[timepoint] <- LAMBDA
    survival_data <- survival(obs)
    obs[, na.omit(variable_names$Sij_varnames)] <- survival_data$Sij
    obs[, "survtime"] <- survival_data$survtimes
    
    #---- Calculating death data for each individual ----
    #Indicator of 1 means the individual died in that interval
    #NAs mean the individual died in a prior interval
    obs[, "death0"] <- 0
    obs[, na.omit(variable_names$deathij_varnames)] <- 
      (obs[, na.omit(variable_names$Sij_varnames)] < int_time)*1 
    
    #---- Survival Analysis ----
    female_data <- obs %>% filter(sex == 0)
    
    p_alive_females <- female_data %>%
      dplyr::select(na.omit(variable_names$deathij_varnames)) %>% 
      map_dbl(~sum(. == 0, na.rm = TRUE))/nrow(female_data)
    
    return(abs(p_alive_females[timepoint] - cp_survival[timepoint]))
  }
  
  #Replace lambda with optimized lambda value
  if(timepoint == 1){
    opt_base_haz[timepoint] <- optim(par = opt_base_haz[timepoint], 
                                     fn = survival_match, 
                                     obs = obs, cp_survival = cp_survival, 
                                     upper = opt_base_haz[timepoint], 
                                     lower = 0.8*opt_base_haz[timepoint], 
                                     method = "L-BFGS-B")$par
  } else {
    opt_base_haz[timepoint] <- optim(par = opt_base_haz[timepoint - 1], 
                                     fn = survival_match, 
                                     obs = obs, cp_survival = cp_survival, 
                                     upper = 2.75*opt_base_haz[timepoint - 1], 
                                     lower = opt_base_haz[timepoint - 1], 
                                     method = "L-BFGS-B")$par
  }
  
  
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
  
  #---- Survival censoring matrix ----
  censor <- (obs[, na.omit(variable_names$Sij_varnames)] == 5)*1
  censor[censor == 0] <- NA
  censor %<>% cbind(1, .)
  
  shifted_censor <- cbind(1, censor[, 1:(ncol(censor) - 1)])
  
  #---- Censor Cij and dem data ----
  obs[, variable_names$Cij_varnames] <- 
    obs[, variable_names$Cij_varnames]*censor
  
  obs[, variable_names$dem_varnames] <- 
    obs[, variable_names$dem_varnames]*shifted_censor
  
  #---- Dementia indicators ----
  for(i in 1:nrow(obs)){
    dem_int <- min(which(obs[i, variable_names$dem_varnames] == 1))
    if(is.finite(dem_int)){
      obs[i, "dem_wave"] <- (dem_int - 1)
      first_censor <- min(which(is.na(obs[i, variable_names$dem_varnames])))
      if(dem_int < 10 & is.finite(first_censor)){
        obs[i, variable_names$dem_varnames[dem_int:(first_censor - 1)]] <- 1 #Changes dementia indicator to 1 after dementia diagnosis
      }
      if(dem_int < 10 & !is.finite(first_censor)){
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
  obs <- compare_survtime_timetodem(obs)
  
  #Time to dem_death
  for(i in 1:nrow(obs)){
    if(obs[i, "dem"] == 0){
      obs[i, "timetodem_death"] <- obs[i, "survtime"]
    } else {
      obs[i, "timetodem_death"] <- min(obs[i, "timetodem"], obs[i, "survtime"])
    }
  }
  
  #---- Contributed time ----
  for(i in 1:nrow(obs)){
    #5-year bands
    last_full_slot <- floor(obs[i, "timetodem_death"]/5)
    full_slots <- na.omit(variable_names$contributed_varnames)[1:last_full_slot]
    obs[i, full_slots] <- 5
    if(last_full_slot != 9){
      partial_slot <- last_full_slot + 1
      obs[i, na.omit(variable_names$contributed_varnames)[partial_slot]] <- 
        obs[i, "timetodem_death"]%%5
    }
  }
  
  #---- Dementia Incidence Analysis ----
  if(timepoint == 1){
    dem_last_wave <- paste0("dem", (timepoint - 1))
    dem_this_wave <- paste0("dem", (timepoint - 1), "-", timepoint)
    contributed <- paste0("contributed", (timepoint - 1), "-", timepoint)
  } else {
    dem_last_wave <- paste0("dem", (timepoint - 2), "-", (timepoint - 1))
    dem_this_wave <- paste0("dem", (timepoint - 1), "-", timepoint)
    contributed <- paste0("contributed", (timepoint - 1), "-", timepoint)
  }
  PY_data <- obs %>% 
    dplyr::select(dem_last_wave, dem_this_wave, contributed) %>% 
    filter(!! as.name(dem_last_wave) == 0) 
  
  dem_inc = round(1000*(sum(PY_data[, dem_this_wave], na.rm = TRUE)/
                          sum(PY_data[, contributed], na.rm = TRUE)), 3)
  
  return(abs(dem_inc - dem_inc_data[timepoint]))
}

#---- Pre-allocation ----

opt_cij_slopes <- rep(0, 10)    #Start with 0 slopes
opt_cij_var1 <- rep(0.001, 10)  #Start with tiny variances in slopes
opt_base_haz <- rep(0.00414, 9)   #Testing replacement of values

timepoint = 5

#---- Values to match ----
#The first values are inputs based on climbing to desired inc rate at age 70
dem_inc_data <- c(0.03, 0.25, 1.00, 
                  head(EURODEM_inc_rates$Total_All_Dementia_1000PY, -1))
cp_survival <- head(female_life_netherlands$CP[-1], -1)

#---- Setup cluster ----
#Setting up cluster for parallel optimization
#If using a Mac/Linux system, it's highly recommended to use the type = "FORK"
#option instead and comment out the clusterEvalQ() lines

stopCluster(cluster)
cluster <- makeCluster(0.5*detectCores(), type = "PSOCK")
clusterEvalQ(cl = cluster, {
  if (!require("pacman")) 
    install.packages("pacman", repos='http://cran.us.r-project.org')
  
  p_load("here", "tidyverse", "MASS")
  
  #Source files
  source(here("RScripts", "sex-dementia_sim_parA.R"))
  source(here("RScripts", "dementia_incidence_EURODEM_pooled.R"))
  source(here("RScripts", "euro_life_tables.R"))
  source(here("RScripts", "variable_names.R"))
  source(here("RScripts", "cognitive_function_model.R"))
  source(here("RScripts", "survival_times.R"))
  source(here("RScripts", "dementia_onset.R"))
  source(here("RScripts", "compare_survtime_timetodem.R"))
  source(here("RScripts", "sex-dementia_sim_data_gen.R"))
}) 

clusterExport(cl = cluster, 
              varlist = c("opt_cij_slopes", "opt_cij_var1", "opt_base_haz", 
                          "dem_inc_data", "timepoint", "column_names", 
                          "variable_names"), 
              envir = environment())

#---- Slope Optimization ----
for(time in timepoint:timepoint){
  optim_values <- 
    replicate(7, 
              optimParallel(par = -0.05, 
                            fn = dem_inc_rate_match,
                            which_opt = "slope",
                            batch_n = 20000, timepoint = time, 
                            dem_inc_data = dem_inc_data, 
                            cp_survival = cp_survival, 
                            opt_cij_slopes = opt_cij_slopes, 
                            opt_cij_var1 = opt_cij_var1, 
                            opt_base_haz = opt_base_haz, 
                            lower = -0.07, 
                            upper = -0.041, 
                            method = "L-BFGS-B", 
                            parallel = list(cl = cluster))$par)
  
  avg_optim_values <- mean(optim_values)
  opt_cij_slopes[time] <- avg_optim_values
}

#---- New cluster ----
stopCluster(cluster)
cluster <- makeCluster(0.5*detectCores(), type = "PSOCK")
clusterEvalQ(cl = cluster, {
  if (!require("pacman")) 
    install.packages("pacman", repos='http://cran.us.r-project.org')
  
  p_load("here", "tidyverse", "MASS")
  
  #Source files
  source(here("RScripts", "sex-dementia_sim_parA.R"))
  source(here("RScripts", "dementia_incidence_EURODEM_pooled.R"))
  source(here("RScripts", "euro_life_tables.R"))
  source(here("RScripts", "variable_names.R"))
  source(here("RScripts", "cognitive_function_model.R"))
  source(here("RScripts", "survival_times.R"))
  source(here("RScripts", "dementia_onset.R"))
  source(here("RScripts", "compare_survtime_timetodem.R"))
  source(here("RScripts", "sex-dementia_sim_data_gen.R"))
}) 

clusterExport(cl = cluster, 
              varlist = c("opt_cij_slopes", "opt_cij_var1", "opt_base_haz", 
                          "dem_inc_data", "timepoint", "column_names", 
                          "variable_names"), 
              envir = environment())

#---- Variance Optimization ----
for(time in timepoint:timepoint){
  optim_values <- 
    replicate(7, 
              optimParallel(par = opt_cij_var1[time], 
                            fn = dem_inc_rate_match,
                            which_opt = "variance",
                            batch_n = 20000, timepoint = time, 
                            dem_inc_data = dem_inc_data, 
                            cp_survival = cp_survival, 
                            opt_cij_slopes = opt_cij_slopes, 
                            opt_cij_var1 = opt_cij_var1, 
                            opt_base_haz = opt_base_haz, 
                            lower = opt_cij_var1[time], 
                            upper = 1.5*opt_cij_var1[time], 
                            method = "L-BFGS-B", 
                            parallel = list(cl = cluster))$par)
  
  avg_optim_values <- mean(optim_values)
  opt_cij_var1[(time + 1):length(opt_cij_var1)] <- avg_optim_values
}

#---- Generate data with newly optimized slopes and variances ----
cij_slopes <- opt_cij_slopes
cij_var1 <- opt_cij_var1

#---- Survival Re-optimization ----
#Objective Function
survival_match <- function(LAMBDA, cp_survival, time){
  obs <- data_gen()
  lambda <- opt_base_haz
  lambda[time] <- LAMBDA
  survival_data <- survival(obs)
  obs[, na.omit(variable_names$Sij_varnames)] <- survival_data$Sij
  obs[, "survtime"] <- survival_data$survtimes

  #---- Calculating death data for each individual ----
  #Indicator of 1 means the individual died in that interval
  #NAs mean the individual died in a prior interval
  obs[, "death0"] <- 0
  obs[, na.omit(variable_names$deathij_varnames)] <-
    (obs[, na.omit(variable_names$Sij_varnames)] < int_time)*1

  #---- Survival Analysis ----
  female_data <- obs %>% filter(sex == 0)

  p_alive_females <- female_data %>%
    dplyr::select(na.omit(variable_names$deathij_varnames)) %>%
    map_dbl(~sum(. == 0, na.rm = TRUE))/nrow(female_data)

  return(abs(p_alive_females[time] - cp_survival[time]))
}

#---- New cluster ----
stopCluster(cluster)
cluster <- makeCluster(0.5*detectCores(), type = "PSOCK")
clusterEvalQ(cl = cluster, {
  if (!require("pacman")) 
    install.packages("pacman", repos='http://cran.us.r-project.org')
  
  p_load("here", "tidyverse", "MASS")
  
  #Source files
  source(here("RScripts", "sex-dementia_sim_parA.R"))
  source(here("RScripts", "dementia_incidence_EURODEM_pooled.R"))
  source(here("RScripts", "euro_life_tables.R"))
  source(here("RScripts", "variable_names.R"))
  source(here("RScripts", "cognitive_function_model.R"))
  source(here("RScripts", "survival_times.R"))
  source(here("RScripts", "dementia_onset.R"))
  source(here("RScripts", "compare_survtime_timetodem.R"))
  source(here("RScripts", "sex-dementia_sim_data_gen.R"))
}) 

clusterExport(cl = cluster, 
              varlist = c("opt_cij_slopes", "opt_cij_var1", "opt_base_haz", 
                          "dem_inc_data", "timepoint", "column_names", 
                          "variable_names", "survival_match"), 
              envir = environment())

#---- Baseline hazard optimization ----
for(time in timepoint:timepoint){
  if (time == 1) {
    base_haz <- replicate(10, 
                          optimParallel(#par = 1.25*opt_base_haz[1],
                            par = 0.0035,
                            fn = survival_match,
                            cp_survival = cp_survival,
                            time = time,
                            #upper = 2*opt_base_haz[1],
                            upper = 0.0041,
                            #lower = opt_base_haz[1],
                            lower = 0.002,
                            method = "L-BFGS-B", 
                            parallel = list(cl = cluster))$par)
  } else {
    base_haz <- replicate(10, 
                          optimParallel(#par = 1.25*opt_base_haz[time - 1],
                                par = 0.0052,
                                fn = survival_match,
                                cp_survival = cp_survival,
                                time = time,
                                upper = 0.00893749,
                                #upper = 2.75*opt_base_haz[time - 1],
                                #lower = opt_base_haz[time - 1],
                                lower = 0.0051,
                                method = "L-BFGS-B", 
                                parallel = list(cl = cluster))$par)
  }
  
  #Replace lambda with optimized lambda value
  opt_base_haz[time:length(opt_base_haz)] <- mean(base_haz)
}


