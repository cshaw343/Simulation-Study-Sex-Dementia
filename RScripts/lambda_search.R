#***************************************************************
# Performs a search for baseline hazards 
# Reference data found in life_table2014.R
#***************************************************************

#---- Package Loading and Options ----
if (!require("pacman")) 
  install.packages("pacman", repos='http://cran.us.r-project.org')

p_load("tidyverse", "magrittr")

#Suppress warnings
options(warn = -1)

#---- Specify source files ----
source("sex-dementia_sim_parA.R")
source("sex-dementia_sim_script.R")
source("life_table2014.R")

#---- Create datsets for the search ----
data_gen <- function(){
  #---- Generating IDs, sex, U ----
  obs <- tibble("id" = seq(from = 1, to = num_obs, by = 1), 
                "sex" = rbinom(num_obs, size = 1, prob = psex), 
                "U" = rnorm(num_obs, mean = 0, sd = 1))
  
  #---- Generating age data ----
  #Creating ages at each timepoint j
  ages <- as_tibble(matrix(NA, nrow = num_obs, ncol = length(age_varnames)))
  for(j in 1:length(age_varnames)){
    if(j == 1){
      ages[, j] = age0 #Creates column of baseline ages
    } else ages[, j] = ages[, (j-1)] + int_time #Creates ages at following timepoints
  }
  colnames(ages) <- age_varnames
  obs %<>% bind_cols(., ages)
  
  #---- Generating centered age data ----
  #Creating centered ages at each timepoint j
  c_ages <- as_tibble(ages - mean(age0)) %>% set_colnames(., agec_varnames)
  obs %<>% bind_cols(., c_ages)
  
  #---- Generating "true" cognitive function Cij ----
  #Cij = b00 + z0i + bo1*sexi + b02*age0i + b03*Ui + (b10 + z1i + b11*sexi + 
  #b12*age0i + b13*Ui)*timej + epsilonij
  
  #---- Generating random terms for slope and intercept ----
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
    set_colnames(., eps_varnames)
  obs %<>% bind_cols(., eps)
  
  #---- Calculating Cij for each individual ----
  #Store Cij values and slope values for each assessment
  compute_Cij <- cog_func(obs)
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
  
  return("obs" = obs)
}

#---- Search for Baseline Hazard ----
#Function for finding lambdas
lambdas <- function(sim_data, cp){
  #Function we are trying to optimize
  survivors <- function(L, obs, cp){
    time_left = -log(obs[US])/(L*exp(g1*obs["sex"] + g2*obs[agec] + 
                                       g3*obs["U"] + g4*obs["sex"]*obs[agec] + 
                                       g5*obs[slope] + g6*obs[C]))
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
    US <- paste("U", test_num, test_num + 1, sep = "")
    agec <- paste("age", test_num, "_c50", sep = "")
    slope <- paste("slope", test_num, test_num + 1, sep = "")
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
                            interval = c(lambdas[j - 1], 1.75*lambdas[j - 1]), 
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
lambda_searches <- replicate(35, find_lambda(unexposed = 0, 
                                            life_table = female_life[-1, ]))

avg_lambdas <- as_tibble(do.call(rbind, lambda_searches["lambdas", ])) %>%
  colMeans()
avg_cps <- as_tibble(do.call(rbind, lambda_searches["cp_alive", ])) %>%
  colMeans()

  
