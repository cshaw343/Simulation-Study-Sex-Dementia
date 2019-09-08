#***************************************************************
# Performs a search for 
#***************************************************************

#---- Package Loading and Options ----
if (!require("pacman")) 
  install.packages("pacman", repos='http://cran.us.r-project.org')

p_load("here", "parallel", "optimParallel")

#---- Source scripts ----
source(here(
  "RScripts", "quadratic_model",
  "sex-dementia_sim_parA_onedemcut_nodemkill_maleAD_quad.R"))
source(here("RScripts", "dementia_incidence_EURODEM_pooled.R"))

#---- Function we want to optimize ----
# We literally have to optimize the entire data generation function, 
#but we will remove irrelevant calculations

dem_rates <- function(PARAMETERS, dem_rates){
  cij_var0 <- PARAMETERS[1]    #Variance of random cognitive intercept
  cij_var1 <- PARAMETERS[2]    #Variance of random linear term
  cij_var2 <- PARAMETERS[3]    #Variance of random quadratic term (use tiny value b/c calcs give 0)
  
  cij_cov01 <- PARAMETERS[4]    #Covariance between random intercept and random linear term
  cij_cov12 <- PARAMETERS[5]    #Covariance between random linear and random quadratic term
  cij_cov02 <- PARAMETERS[6]    #Covariance between random intercept and random quadratic term (use tiny value b/c calcs give 0)
  
  b10 <- PARAMETERS[7]   #Cognitive linear term for males (taken from quad fit to linear model)          
  b20 <- PARAMETERS[8]   #Cognitive quadratic term for males (taken from quad fit to linear model)
  
  dem_cut <- PARAMETERS[9]  #dementia cut point
  
  #---- Create a blank dataset ----
  num_obs = 100
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
  
  #Generating random terms for quadratic model
  #Defining the covariance matrix
  quad_coeff_cov <- matrix(c(cij_var0, cij_cov01, cij_cov02, 
                             cij_cov01, cij_var1, cij_cov12, 
                             cij_cov02, cij_cov12, cij_var2), 
                           nrow = 3, byrow = TRUE)
  
  #Generate random terms for each individual
  obs[, c("z_0i", "z_1i", "z_2i")] <- 
    mvrnorm(n = num_obs, mu = rep(0, 3), Sigma = quad_coeff_cov) 
  
  #Generating noise term (unexplained variance in Cij) 
  obs[, "eps_ij"] <- rnorm(n = num_obs, mean = 0, sd = sqrt(cij_var3))
  
  #Calculating quadratic coefficients for each individual
  obs[, c("a0", "a1", "a2")] <- calc_coeff(obs)
  
  #Calculating Cij at each time point for each individual
  obs[, variable_names$Cij_varnames] <- compute_Cij(obs)
  
  #---- Check for dementia at baseline ----
  obs <- obs[obs[, "Ci0"] > dem_cut, ]
  
  #---- Generate survival time for each person ----
  #Refer to Manuscript/manuscript_equations.pdf for equation
  
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
  
  #Time to dem_death
  for(i in 1:ncol(obs)){
    if(obs["dem", i] == 0){
      obs["timetodem_death", i] <- obs["survtime", i]
    } else {
      obs["timetodem_death", i] <- min(obs["timetodem", i], obs["survtime", i])
    }
  }
  
  #---- Censor Cij and dem data ----
  obs <- survival_censor(obs)
  
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
  
  #---- Covert to dataframe ----
  obs <- as.data.frame(t(obs))
  
  #---- Dem inc rate calculations ----
  male_data <- obs %>% filter(female == 0)
  
  male_sim_inc_rates <- matrix(ncol = 6, nrow = 1)
  
  for(slot in 4:num_tests){
      dem_last_wave <- paste0("dem", (slot - 2), "-", (slot - 1))
      dem_this_wave <- paste0("dem", (slot - 1), "-", slot)
      death_last_wave <- paste0("death", (slot - 2), "-", (slot - 1))
      death_this_wave <- paste0("death", (slot - 1), "-", slot)
      contributed <- paste0("contributed", (slot - 1), "-", slot)
  
    PY_data <- male_data %>% 
      dplyr::select(death_last_wave, death_this_wave, 
                    dem_last_wave, dem_this_wave, contributed) %>% 
      filter(!! as.name(death_last_wave) == 0 & 
               !! as.name(dem_last_wave) == 0) 
    
    male_sim_inc_rates[1, (slot - 3)] = 
      round(1000*(sum(PY_data[, dem_this_wave], 
                      na.rm = TRUE)/sum(PY_data[, contributed])), 3)
  }
  
  #---- Values to return ----
  return(sum(abs(male_sim_inc_rates - dem_rates)))

}

#---- Setting up the cluster for optimization ----
cluster <- makeCluster(0.5*detectCores(), type = "PSOCK")

clusterEvalQ(cl = cluster, {
  
  #---- Package Loading and Options ----
  if (!require("pacman")) 
    install.packages("pacman", repos='http://cran.us.r-project.org')
  
  p_load("here", "MASS")
  
  #Suppress warnings
  options(warnings = -1)
  
  #---- Source scripts ----
  source(here(
    "RScripts", "quadratic_model",
    "sex-dementia_sim_parA_onedemcut_nodemkill_maleAD_quad.R"))
  source(here("RScripts", "quadratic_model", "variable_names_quad.R"))
  source(here("RScripts", "quadratic_model", "create_ages.R"))
  source(here("RScripts", "quadratic_model", "calc_coeff.R"))
  source(here("RScripts", "quadratic_model", "compute_Cij.R"))
  source(here("RScripts", "quadratic_model", "last_Cij.R"))
  source(here("RScripts", "quadratic_model", "survival_times_quad.R"))
  source(here("RScripts", "quadratic_model", "survival_censor.R"))
  source(here("RScripts", "quadratic_model", "dementia_onset_quad.R"))
  source(here("RScripts", "quadratic_model", 
              "compare_survtime_timetodem_quad.R"))
})

#---- Setting up and exporting variables ----
maleAD_rates <- head(EURODEM_inc_rates$Male_AD_1000PY, -1)
parameter_start <- c(cij_var0, cij_var1, cij_var2, 
                     cij_cov01, cij_cov12, cij_cov02, 
                     b10, b20, dem_cut)

clusterExport(cl = cluster,
              varlist = c("maleAD_rates", "parameter_start"),
              envir = environment())

#---- Doing the optimization ----
optim_dem <- optimParallel(par = parameter_start, 
                           fn = dem_rates, dem_rates = maleAD_rates,
                           lower = c(0, 0.05, 0.001, -0.3, -0.1, 0.001,
                                     0.05, -0.003, -6.4),
                           upper = c(2, 1, 0.05, -0.01, -0.001, 0.1,
                                     0.05, -0.003, -6.4),
                           method = "L-BFGS-B", 
                           parallel = list(cl = cluster))

opt_par <- optim_dem$par
how_close <- optim_dem$value


