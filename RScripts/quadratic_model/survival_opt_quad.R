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
  "sex-dementia_sim_parB_highUonSurv_onedemcut_nodemkill_male_AllDem_quad.R"))
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
opt_lambdas <- function(sim_data_unexposed, cp_unexposed){

  #Function we are trying to optimize
  survivors <- function(L, obs, pop_size, cp){
    survtime = -log(sim_data_unexposed[, variable_names$r1ij_varnames[j]])/
      (L*
         exp(g1[j]*sim_data_unexposed[, "female"] +
               g2*sim_data_unexposed[, "U"] +
               g3*
               (1 - sim_data_unexposed[, "female"])*sim_data_unexposed[, "U"]))

    alive <- (survtime >= 5) * 1
    return(abs((sum(alive)/pop_size) - cp_unexposed[j]))
  }

  #Create vectors to return
  sim_cp_unexposed <- vector(length = num_tests)
  opt_lambdas <- vector(length = num_tests)

  #Begin search
  for(j in 1:length(opt_lambdas)){
    if(j == 1){
      opt_lambdas[j:length(opt_lambdas)] =
        optim(lambda[j], survivors,
              lower = 0.005,
              upper = 0.007,
              obs = sim_data_unexposed,
              pop_size = nrow(sim_data_unexposed),
              cp = cp_unexposed)$par
      Sij <- t(survival_temp(obs_matrix = t(sim_data_unexposed),
                             lambda = opt_lambdas))
      sim_data_unexposed = cbind(sim_data_unexposed, (Sij[, j] >= 5)*1)
      colnames(sim_data_unexposed)[ncol(sim_data_unexposed)] <- "alive_now"
      sim_cp_unexposed[j] =
        sum(sim_data_unexposed[, "alive_now"]/nrow(sim_data_unexposed))
    } else if(j >= 2 && j <= 2){
      opt_lambdas[j:length(opt_lambdas)] =
        optim(1.74*opt_lambdas[j - 1], survivors,
              lower = 1.7325*opt_lambdas[j - 1],
              upper = 1.75*opt_lambdas[j - 1],
              obs = sim_data_unexposed,
              pop_size = nrow(sim_data_unexposed),
              cp = cp_unexposed)$par
      Sij <- t(survival_temp(obs_matrix = t(sim_data_unexposed),
                             lambda = opt_lambdas))
      sim_data_unexposed = cbind(sim_data_unexposed, (Sij[, j] >= 5)*1)
      colnames(sim_data_unexposed)[ncol(sim_data_unexposed)] <- "alive_now"
      sim_cp_unexposed[j] =
        sum(sim_data_unexposed[, "alive_now"]/nrow(sim_data_unexposed))
    } else if(j >= 3 && j <= 3){
      opt_lambdas[j:length(opt_lambdas)] =
        optim(1.83*opt_lambdas[j - 1], survivors,
              lower = 1.825*opt_lambdas[j - 1],
              upper = 1.85*opt_lambdas[j - 1],
              obs = sim_data_unexposed,
              pop_size = nrow(sim_data_unexposed),
              cp = cp_unexposed)$par
      Sij <- t(survival_temp(obs_matrix = t(sim_data_unexposed),
                             lambda = opt_lambdas))
      sim_data_unexposed = cbind(sim_data_unexposed, (Sij[, j] >= 5)*1)
      colnames(sim_data_unexposed)[ncol(sim_data_unexposed)] <- "alive_now"
      sim_cp_unexposed[j] =
        sum(sim_data_unexposed[, "alive_now"]/nrow(sim_data_unexposed))
    } else if(j >= 4 && j <= 4){
      opt_lambdas[j:length(opt_lambdas)] =
        optim(1.95*opt_lambdas[j - 1], survivors,
              lower = 1.91*opt_lambdas[j - 1],
              upper = 2*opt_lambdas[j - 1],
              obs = sim_data_unexposed,
              pop_size = nrow(sim_data_unexposed),
              cp = cp_unexposed)$par
      Sij <- t(survival_temp(obs_matrix = t(sim_data_unexposed),
                             lambda = opt_lambdas))
      sim_data_unexposed = cbind(sim_data_unexposed, (Sij[, j] >= 5)*1)
      colnames(sim_data_unexposed)[ncol(sim_data_unexposed)] <- "alive_now"
      sim_cp_unexposed[j] =
        sum(sim_data_unexposed[, "alive_now"]/nrow(sim_data_unexposed))
    } else if(j >= 5 && j <= 5){
      opt_lambdas[j:length(opt_lambdas)] =
        optim(2.2*opt_lambdas[j - 1], survivors,
              lower = 2.15*opt_lambdas[j - 1],
              upper = 2.25*opt_lambdas[j - 1],
              obs = sim_data_unexposed,
              pop_size = nrow(sim_data_unexposed),
              cp = cp_unexposed)$par
      Sij <- t(survival_temp(obs_matrix = t(sim_data_unexposed),
                             lambda = opt_lambdas))
      sim_data_unexposed = cbind(sim_data_unexposed, (Sij[, j] >= 5)*1)
      colnames(sim_data_unexposed)[ncol(sim_data_unexposed)] <- "alive_now"
      sim_cp_unexposed[j] =
        sum(sim_data_unexposed[, "alive_now"]/nrow(sim_data_unexposed))
    } else if(j >= 6 && j <= 6){
      opt_lambdas[j:length(opt_lambdas)] =
        optim(2.225*opt_lambdas[j - 1], survivors,
              lower = 2.2*opt_lambdas[j - 1],
              upper = 2.25*opt_lambdas[j - 1],
              obs = sim_data_unexposed,
              pop_size = nrow(sim_data_unexposed),
              cp = cp_unexposed)$par
      Sij <- t(survival_temp(obs_matrix = t(sim_data_unexposed),
                             lambda = opt_lambdas))
      sim_data_unexposed = cbind(sim_data_unexposed, (Sij[, j] >= 5)*1)
      colnames(sim_data_unexposed)[ncol(sim_data_unexposed)] <- "alive_now"
      sim_cp_unexposed[j] =
        sum(sim_data_unexposed[, "alive_now"]/nrow(sim_data_unexposed))
    } else if(j >= 7 && j <= 7){
    opt_lambdas[j:length(opt_lambdas)] =
      optim(2.4625*opt_lambdas[j - 1], survivors,
            lower = 2.46*opt_lambdas[j - 1],
            upper = 2.47*opt_lambdas[j - 1],
            obs = sim_data_unexposed,
            pop_size = nrow(sim_data_unexposed),
            cp = cp_unexposed)$par
    Sij <- t(survival_temp(obs_matrix = t(sim_data_unexposed),
                           lambda = opt_lambdas))
    sim_data_unexposed = cbind(sim_data_unexposed, (Sij[, j] >= 5)*1)
    colnames(sim_data_unexposed)[ncol(sim_data_unexposed)] <- "alive_now"
    sim_cp_unexposed[j] =
      sum(sim_data_unexposed[, "alive_now"]/nrow(sim_data_unexposed))
    } else if(j >= 8 && j <= 8){
      opt_lambdas[j:length(opt_lambdas)] =
        optim(2.69*opt_lambdas[j - 1], survivors,
              lower = 2.675*opt_lambdas[j - 1],
              upper = 2.7025*opt_lambdas[j - 1],
              obs = sim_data_unexposed,
              pop_size = nrow(sim_data_unexposed),
              cp = cp_unexposed)$par
      Sij <- t(survival_temp(obs_matrix = t(sim_data_unexposed),
                             lambda = opt_lambdas))
      sim_data_unexposed = cbind(sim_data_unexposed, (Sij[, j] >= 5)*1)
      colnames(sim_data_unexposed)[ncol(sim_data_unexposed)] <- "alive_now"
      sim_cp_unexposed[j] =
        sum(sim_data_unexposed[, "alive_now"]/nrow(sim_data_unexposed))
    } else if(j >= 9 && j <= 9){
      opt_lambdas[j:length(opt_lambdas)] =
        optim(2.575*opt_lambdas[j - 1], survivors,
              lower = 2.55*opt_lambdas[j - 1],
              upper = 2.617*opt_lambdas[j - 1],
              obs = sim_data_unexposed,
              pop_size = nrow(sim_data_unexposed),
              cp = cp_unexposed)$par
      Sij <- t(survival_temp(obs_matrix = t(sim_data_unexposed),
                             lambda = opt_lambdas))
      sim_data_unexposed = cbind(sim_data_unexposed, (Sij[, j] >= 5)*1)
      colnames(sim_data_unexposed)[ncol(sim_data_unexposed)] <- "alive_now"
      sim_cp_unexposed[j] =
        sum(sim_data_unexposed[, "alive_now"]/nrow(sim_data_unexposed))
    }
  }
  return(list("opt_lambdas" = opt_lambdas,
              "sim_cp_unexposed" = sim_cp_unexposed))
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
      warm_start = 0.8775
      
      exp_g1s[i:length(exp_g1s)] = 
        optim(warm_start, sim_hr, method = "L-BFGS-B",
              lower = 0.5,
              upper = 1,
              obs = sim_data, hr = hr, 
              opt_lambdas = opt_lambdas, 
              j = i)$par
    } else {
      warm_start = exp_g1s[i - 1]
      
      exp_g1s[i:length(exp_g1s)] = 
        optim(warm_start, sim_hr, method = "L-BFGS-B",
              lower = 0.5,
              upper = 1,
              obs = sim_data, hr = hr, 
              opt_lambdas = opt_lambdas, 
              j = i)$par
    }
    
    sim_data[, na.omit(variable_names$Sij_varnames)] <- 
      t(survival_temp(t(sim_data), opt_lambdas, log(exp_g1s)))
    dead_now <- (sim_data[, "survtime0-1"] < 5) * 1
    
    cox_model <- coxph(Surv(sim_data[, "survtime0-1"], 
                            dead_now) ~ sim_data[, "female"])
    sim_HR <- as.numeric(exp(cox_model$coefficients))
  }
  return(exp_g1s)
}

#---- Doing the optimization ----
#Do multiple optimizations and average over runs
cp_unexposed <- male_life_US$CP[-1]
hr_wm <- Hratio_US[-1, ]

lambda_optimization <- function(cp_unexposed){
  sim_data <- pre_survival_data_gen(500000)
  sim_data_unexposed <- sim_data[sim_data[, "female"] == 0, ]
  optim_lambda <- opt_lambdas(sim_data_unexposed, cp_unexposed)$opt_lambdas
  
  return(optim_lambda)
}

lambda_runs <- replicate(1, lambda_optimization(cp_unexposed))
opt_lambdas <- colMeans(t(lambda_runs))

exp_g1_optimization <- function(hr, opt_lambdas, exp_g1s, start){
  sim_data <- pre_survival_data_gen(100000)
  optim_exp_g1s <- opt_exp_g1s(hr, opt_lambdas, exp_g1s, start = 1, sim_data)
}

plan(multiprocess, workers = 0.5*availableCores())
exp_g1_runs <- 
  future_replicate(10, exp_g1_optimization(hr_wm, opt_lambdas, exp(g1), 
      start = 1))
opt_exp_g1s <- rowMeans(exp_g1_runs)




