#***************************************************************
# Performs a search for baseline hazards and 
# log(g1) (affect of female sex on survival)
# Reference data found in euro_life_tables.R 
# (using Netherlands data)
#***************************************************************

#---- Package Loading and Options ----
if (!require("pacman")) 
  install.packages("pacman", repos='http://cran.us.r-project.org')

p_load("here", "MASS")

#Suppress warnings
options(warnings = -1)

#---- Specify source files ----
source(here(
  "RScripts", "linear_model",
  "sex-dementia_sim_parC_onedemcut_uniform_timetodem_nodemkill_UonInt.R"))
source(here("RScripts", "life_table_calcs.R"))
source(here("RScripts", "linear_model", "variable_names.R"))
source(here("RScripts", "linear_model", "create_ages.R"))
source(here("RScripts", "linear_model", "gen_random_effects.R"))
source(here("RScripts", "linear_model", "cognitive_function_model.R"))
source(here("RScripts", "linear_model", "survival_times.R"))

#---- Create datasets for the search ----
pre_survival_data_gen <- function(num_obs){
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
  
  return(obs)
}

survival_temp <- function(obs_matrix, lambda, opt_log_g1 = g1){
  #Calculate survival times for each interval
  Sij <- matrix(ncol = ncol(obs_matrix), nrow = (length(visit_times) - 1))
  for(i in 1:ncol(obs_matrix)){
    survtimes <- matrix(NA, nrow = (length(visit_times) - 1), ncol = 1)
    for(j in 1:length(survtimes)){
      r_name <- variable_names$rij_varnames[j]
      
      survtime = -log(obs_matrix[r_name, i])/
        (lambda[j]*exp(opt_log_g1[j]*obs_matrix["female", i] + 
                         g2*obs_matrix["U", i] + 
                         g3*obs_matrix["female", i]*obs_matrix["U", i]))
      
      survtimes[j] <- as.numeric(survtime)
      if(survtime < 5){
        break
      }
    }
    Sij[, i] <- survtimes
  }
  Sij[Sij > 5] <- 5
  return(Sij)
}


#---- Search for Baseline Hazard ----
#Function for finding lambdas
opt_lambdas <- function(sim_data_unexposed, cp50_unexposed){
  
  #Function we are trying to optimize
  survivors <- function(L, obs, pop_size, cp){
    survtime = -log(sim_data_unexposed[, variable_names$rij_varnames[j]])/
      (L*
         exp(g1[j]*sim_data_unexposed[, "female"] + 
               g2*sim_data_unexposed[, "U"] + 
               g3*sim_data_unexposed[, "female"]*sim_data_unexposed[, "U"]))
    
    alive <- (survtime >= 5) * 1
    return(abs((sum(alive)/pop_size) - cp50_unexposed[j]))
  }
  
  #Create vectors to return
  sim_cp50_unexposed <- vector(length = num_tests)
  opt_lambdas <- vector(length = num_tests)
  
  #Begin search
  num_unexposed = nrow(sim_data_unexposed)
  for(j in 1:length(opt_lambdas)){
    if(j == 1){
      opt_lambdas[j:length(opt_lambdas)] = 
        optimise(survivors, interval = c(0, 1), obs = sim_data_unexposed, 
                 pop_size = num_unexposed, cp = cp50_unexposed)$minimum
      Sij <- t(survival_temp(obs_matrix = t(sim_data_unexposed), 
                             lambda = opt_lambdas))
      sim_data_unexposed = cbind(sim_data_unexposed, (Sij[, j] >= 5)*1)
      colnames(sim_data_unexposed)[ncol(sim_data_unexposed)] <- "alive_now"
      sim_cp50_unexposed[j] = 
        sum(sim_data_unexposed[, "alive_now"]/num_unexposed)
    } else {
      sim_data_unexposed <- 
        sim_data_unexposed[sim_data_unexposed[, "alive_now"] == 1, ] 
      opt_lambdas[j:length(opt_lambdas)] = 
        optimise(survivors, 
                 interval = c(opt_lambdas[j - 1], 2*opt_lambdas[j - 1]), 
                 obs = sim_data_unexposed, pop_size = num_unexposed, 
                 cp = cp50_unexposed)$minimum
      Sij <- t(survival_temp(obs_matrix = t(sim_data_unexposed), 
                             lambda = opt_lambdas))
      sim_data_unexposed[, "alive_now"] <- (Sij[, j] >= 5)*1
      sim_cp50_unexposed[j] = 
        sum(sim_data_unexposed[, "alive_now"]/num_unexposed)
    }
  }
  return(list("opt_lambdas" = opt_lambdas, 
              "sim_cp50_unexposed" = sim_cp50_unexposed))
}

#---- Search for effect of "female" on baseline hazard ----
opt_log_g1s <- function(sim_data_exposed, cp50_exposed){
  #Function we are trying to optimize
  survivors <- function(LOG_G1, obs, pop_size, cp, opt_lambdas){
    survtime = -log(sim_data_exposed[, variable_names$rij_varnames[j]])/
      (opt_lambdas[j]*
         exp(LOG_G1*sim_data_exposed[, "female"] + 
               g2*sim_data_exposed[, "U"] + 
               g3*sim_data_exposed[, "female"]*sim_data_exposed[, "U"]))
    
    alive <- (survtime >= 5) * 1
    return(abs((sum(alive)/pop_size) - cp50_exposed[j]))
  }
  
  #Create vectors to return
  sim_cp50_exposed <- vector(length = num_tests)
  opt_log_g1s <- vector(length = num_tests)
  
  #Begin search
  num_exposed = nrow(sim_data_exposed)
  for(j in 1:length(opt_log_g1s)){
    if(j == 1){
      opt_log_g1s[j:length(opt_log_g1s)] = 
        optimise(survivors, interval = c(-1, 0), obs = sim_data_exposed, 
                 pop_size = num_exposed, cp = cp50_exposed, 
                 opt_lambdas = optim_lambda$opt_lambdas)$minimum
      Sij <- t(survival_temp(obs_matrix = t(sim_data_exposed),
                             lambda = optim_lambda$opt_lambdas,
                             opt_log_g1 = opt_log_g1s))
      sim_data_exposed = cbind(sim_data_exposed, (Sij[, j] >= 5)*1)
      colnames(sim_data_exposed)[ncol(sim_data_exposed)] <- "alive_now"
      sim_cp50_exposed[j] = sum(sim_data_exposed[, "alive_now"]/num_exposed)
    } else {
      sim_data_exposed <- 
        sim_data_exposed[sim_data_exposed[, "alive_now"] == 1, ] 
      opt_log_g1s[j:length(opt_log_g1s)] = 
        optimise(survivors, 
                 interval = c(-1, 0), 
                 obs = sim_data_exposed, pop_size = num_exposed, 
                 cp = cp50_exposed, 
                 opt_lambdas = optim_lambda$opt_lambdas)$minimum
      Sij <- t(survival_temp(obs_matrix = t(sim_data_exposed), 
                             lambda = optim_lambda$opt_lambdas,
                             opt_log_g1 = opt_log_g1s))
      sim_data_exposed[, "alive_now"] <- (Sij[, j] >= 5)*1
      sim_cp50_exposed[j] = sum(sim_data_exposed[, "alive_now"]/num_exposed)
    }
  }
  return(list("opt_log_g1s" = opt_log_g1s, 
              "sim_cp50_exposed" = sim_cp50_exposed))
}
  
#---- Doing the optimization ----
sim_data <- pre_survival_data_gen(500000)

cp50_unexposed <- male_life_netherlands$cum_surv_cond50[-1]
cp50_exposed <- female_life_netherlands$cum_surv_cond50[-1]

sim_data_unexposed <- sim_data[sim_data[, "female"] == 0, ]
sim_data_exposed <- sim_data[sim_data[, "female"] == 1, ]

optim_lambda <- opt_lambdas(sim_data_unexposed, cp50_unexposed)
optim_g1s <- opt_log_g1s(sim_data_exposed, cp50_exposed)


