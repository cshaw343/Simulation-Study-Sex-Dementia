#***************************************************************
# Performs a search for baseline hazards and 
# log(g1) (affect of female sex on survival)
# Reference data found in US_life_table_calcs.R 
# (this is the 1919-1921 birth cohort)
#***************************************************************

#---- Package Loading and Options ----
if (!require("pacman")) 
  install.packages("pacman", repos='http://cran.us.r-project.org')

p_load("here", "MASS", "survival", "future.apply")

#Suppress warnings
options(warnings = -1)

#---- Specify source files ----
source(here(
  "RScripts", "quadratic_model",
  "sex-dementia_sim_parB_onedemcut_nodemkill_male_AllDem_quad.R"))
source(here("RScripts", "US_life_table_calcs.R"))
source(here("RScripts", "quadratic_model", "variable_names_quad.R"))
source(here("RScripts", "quadratic_model", "create_ages.R"))
source(here("RScripts", "quadratic_model", "calc_coeff.R"))
source(here("RScripts", "quadratic_model", "compute_Cij.R"))

#---- Create datasets for the search ----
pre_survival_data_gen <- function(num_obs){
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
  
  #Take the negative absolute value of the quadratic noise so that this is drawn
  #from a negative half normal distribution
  obs[, "z_2i"] <- -abs(obs[, "z_2i"])
  
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
  obs[, variable_names$r1ij_varnames[1:num_tests]] <- 
    replicate(num_tests, runif(nrow(obs), min = 0, max = 1))
  
  obs[, "death0"] <- 0
  
  return(obs)
}

survival_temp <- function(obs_matrix, lambda, g1){
  #Calculate survival times for each interval
  Sij <- matrix(ncol = ncol(obs_matrix), nrow = (length(visit_times) - 1))
  for(i in 1:ncol(obs_matrix)){
    survtimes <- matrix(NA, nrow = (length(visit_times) - 1), ncol = 1)
    for(j in 1:length(survtimes)){
      r_name <- variable_names$r1ij_varnames[j]
      
      survtime = -log(obs_matrix[r_name, i])/
        (lambda[j]*exp(g1[j]*obs_matrix["female", i] + 
                         g2*obs_matrix["U", i] + 
                         g3*(1 - obs_matrix["female", i])*obs_matrix["U", i]))
      
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
opt_lambdas <- function(cp_unexposed, lambdas, exp_g1s, start, sim_data){
  
  #---- Function we are trying to optimize ----
  sim_cp_surv <- function(L, obs, cp_unexposed, exp_g1s, j){
    survival_int <- variable_names$Sij_varnames[j]
    dead <- which(is.na(obs[, survival_int]))
    survtime = -log(obs[, variable_names$r1ij_varnames[j]])/
      (L*
         exp(log(exp_g1s[j])*obs[, "female"] + g2*obs[, "U"] + 
               g3*(1 - obs[, "female"])*obs[, "U"]))
    
    survtime[dead] <- NA
    alive_now <- (survtime >= 5) * 1
    
    return(abs(mean(alive_now, na.rm = TRUE) - cp_unexposed[j]))
  }
  
  sim_data[, na.omit(variable_names$Sij_varnames)] <- 
    t(survival_temp(t(sim_data), lambdas, log(exp_g1s)))

  #Begin search
  for(i in 1:length(lambdas)){
    if(i == 1){
      warm_start = 0.011
      
      lambdas[i:length(lambdas)] =
        optim(warm_start, sim_cp_surv, method = "L-BFGS-B",
              lower = 0.001,
              upper = 0.05,
              obs = sim_data,
              cp_unexposed = cp_unexposed, exp_g1s = exp(g1), j = i)$par
      
      sim_data[, na.omit(variable_names$Sij_varnames)] <- 
        t(survival_temp(t(sim_data), lambdas, log(exp_g1s)))
      
    } else if(i >= 2 && i <= 9){
      warm_start = 1.5*lambdas[i - 1]
      
      lambdas[i:length(lambdas)] =
        optim(warm_start, sim_cp_surv, method = "L-BFGS-B",
              lower = lambdas[i - 1],
              upper = 5*lambdas[i - 1],
              obs = sim_data,
              cp_unexposed = cp_unexposed, exp_g1s = exp(g1), j = i)$par
      
      sim_data[, na.omit(variable_names$Sij_varnames)] <- 
        t(survival_temp(t(sim_data), lambdas, log(exp_g1s)))
    } 
  }
  return(lambdas)
}

#---- Search for effect of "female" on baseline hazard ----
opt_exp_g1s <- function(hr, opt_lambdas, exp_g1s, start, sim_data){
  
  #---- Function we are trying to optimize ----
  sim_hr <- function(EXP_G1, obs, hr, opt_lambdas, j){
    survival_int <- variable_names$Sij_varnames[j]
    dead <- which(is.na(obs[, survival_int]))
    survtime = -log(obs[, variable_names$r1ij_varnames[j]])/
      (opt_lambdas[j]*
         exp(log(EXP_G1)*obs[, "female"] + g2*obs[, "U"] + 
               g3*(1 - obs[, "female"])*obs[, "U"]))
    
    survtime[dead] <- NA
    dead_now <- (survtime < 5) * 1
    survtime[survtime == 5] <- 4.99999
    
    #Modeled mortality
    cox_model <- coxph(Surv(survtime, dead_now) ~ obs[, "female"])
    sim_HR <- as.numeric(exp(cox_model$coefficients))
    
    return(abs(sim_HR - hr[j]))
  }
  
  sim_data[, na.omit(variable_names$Sij_varnames)] <- 
    t(survival_temp(t(sim_data), opt_lambdas, log(exp_g1s)))
  
  for(i in start:length(exp_g1s)){
    if(i == 1){
      warm_start = 1
      
      exp_g1s[i:length(exp_g1s)] = 
        optim(warm_start, sim_hr, method = "L-BFGS-B",
              lower = 0.005,
              upper = 2,
              obs = sim_data, hr = hr, 
              opt_lambdas = opt_lambdas, 
              j = i)$par
      
      sim_data[, na.omit(variable_names$Sij_varnames)] <- 
        t(survival_temp(t(sim_data), opt_lambdas, log(exp_g1s)))
    } else {
      warm_start = exp_g1s[i - 1]
      
      exp_g1s[i:length(exp_g1s)] = 
        optim(warm_start, sim_hr, method = "L-BFGS-B",
              lower = 0.005,
              upper = 2,
              obs = sim_data, hr = hr, 
              opt_lambdas = opt_lambdas, 
              j = i)$par
      
      sim_data[, na.omit(variable_names$Sij_varnames)] <- 
        t(survival_temp(t(sim_data), opt_lambdas, log(exp_g1s)))
    }
  }
  return(exp_g1s)
}

#---- Optimization inputs ----
#Do multiple optimizations and average over runs
cp_unexposed <- male_life_US$CP[-1]
hr_wm <- Hratio_US[-1, ]

#---- Lambda optimization ----
lambda_optimization <- function(cp_unexposed, lambdas, exp_g1s, start){
  sim_data <- pre_survival_data_gen(200000)
  sim_data_unexposed <- sim_data[sim_data[, "female"] == 0, ]
  optim_lambda <- opt_lambdas(cp_unexposed, lambdas, exp_g1s, start = 1, 
                              sim_data_unexposed)
  return(optim_lambda)
}

start <- Sys.time()
plan(multiprocess, workers = 0.5*availableCores())
lambda_runs <- 
  future_replicate(10, lambda_optimization(cp_unexposed, lambda, exp(g1), 
                                           start = 1))
opt_lambdas <- rowMeans(lambda_runs)
plan(sequential)
Sys.time() - start

#---- Gamma optimization ----
exp_g1_optimization <- function(hr, opt_lambdas, exp_g1s, start){
  sim_data <- pre_survival_data_gen(100000)
  optim_exp_g1s <- opt_exp_g1s(hr, opt_lambdas, exp_g1s, start = 1, 
                               sim_data)
}

start <- Sys.time()
plan(multiprocess, workers = 0.5*availableCores())
exp_g1_runs <- 
  future_replicate(50, exp_g1_optimization(hr_wm, opt_lambdas, exp(g1), 
                                           start = 1))
opt_exp_g1s <- rowMeans(exp_g1_runs)
plan(sequential)
Sys.time() - start


