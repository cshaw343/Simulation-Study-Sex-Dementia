#***************************************************************
# Performs a search for baseline hazards 
# Reference data found in euro_life_tables.R 
# (using Netherlands data)
#***************************************************************

#---- Package Loading and Options ----
if (!require("pacman")) 
  install.packages("pacman", repos='http://cran.us.r-project.org')

p_load("tidyverse", "magrittr")

#Suppress warnings
options(warnings = -1)

#---- Specify source files ----
source("RScripts/sex-dementia_sim_parA.R")
source("RScripts/sex-dementia_sim_script.R")
source("RScripts/euro_life_tables.R")
source("RScripts/variable_names.R")

#---- Create datasets for the search ----
data_gen <- function(){
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
  
  #Calculating Cij for each individual
  #Store Cij values and slope values for each assessment
  compute_Cij <- cog_func(cij_knots, cij_slopes, obs)
  Cij <- as.data.frame(compute_Cij$Cij) %>% 
    set_colnames(., variable_names$Cij_varnames)
  cij_slopeij <- as.data.frame(compute_Cij$slopes) %>% 
    #remove the last variable name because there are only 10 intervals
    set_colnames(., head(variable_names$cij_slopeij_varnames, -1)) 
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
    set_colnames(head(variable_names$Sij_varnames, -1))
  obs %<>% bind_cols(., Sij)
  
  return("obs" = obs)
}

#---- Search for Baseline Hazard ----
#Function for finding lambdas
lambdas <- function(sim_data, cp){
  #Function we are trying to optimize
  survivors <- function(L, obs, cp){
    time_left = -log(obs[rij])/(L*exp(g1*obs[, "sex"] + g2*obs[, agec] + 
                                       g3*obs[, "U"] + 
                                        g4*obs[, "sex"]*obs[, agec] + 
                                        g5*obs[, slope] + g6*obs[, C]))
    alive <- (time_left > 5) * 1
    return(abs(mean(alive) - cp))
  }
  #Create dataframes to return
  cp_alive <- vector(length = num_tests)
  lambdas <- vector(length = num_tests)
  Sij <- vector(length = num_tests)
  #Begin search
  for(j in 1:length(lambdas)){
    test_num = j - 1
    rij <- variable_names$rij_varnames[j]
    agec <- variable_names$agec_varnames[j]
    slope <- variable_names$cij_slopeij_varnames[j]
    C <- paste("Ci", test_num, sep = "")
    life_prob = as.double(cp[j + 1, "CP"])
    if(j == 1){
      lambdas[j] = optimise(survivors, interval = c(0, 1), 
                            obs = sim_data, cp = life_prob)$minimum
      Sij <- survival(obs = sim_data, lambda = lambdas)
      alive_now <- (Sij[[j]] > 5)*1
      sim_data %<>% cbind(., "alive" = (Sij[[j]] > 5)*1)
      cp_alive[j] = mean(alive_now)
    } else {
      sim_data %<>% filter(alive == 1)
      lambdas[j] = optimise(survivors, 
                            interval = c(lambdas[j - 1], 2*lambdas[j - 1]), 
                            obs = sim_data, cp = life_prob)$minimum
      Sij <- survival(obs = sim_data, lambda = lambdas)
      alive_now <- (Sij[[j]] > 5)*1
      sim_data %<>% mutate("alive" = alive_now)
      cp_alive[j] = mean(alive_now)
    }
  }
  return(list("lambdas" = lambdas, "cp_alive" = cp_alive))
}

#---- Averaging over baseline hazard searches----
find_lambda <- function(unexposed, life_table){
  simdata <- data_gen() %>% filter(sex == unexposed)
  search <- lambdas(sim_data = simdata, cp = life_table)
  return(search)
}

#---- Check conditional probabilities using baseline hazards ----
#Make sure to rerun parameter file with desired baseline hazards before running 
#actual simulation
lambda_searches <- 
  replicate(35, find_lambda(unexposed = 0, 
                           life_table = female_life_netherlands))

avg_lambdas <- as_tibble(do.call(rbind, lambda_searches["lambdas", ])) %>%
  colMeans()
avg_cps <- as_tibble(do.call(rbind, lambda_searches["cp_alive", ])) %>%
  colMeans()


