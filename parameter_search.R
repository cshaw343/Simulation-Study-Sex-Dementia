#***************************************************************
# Performs a parameter search using a specified function and 
# specified file containing reference data
#***************************************************************

#---- Package Loading and Options ----
if (!require("pacman")) 
  install.packages("pacman", repos='http://cran.us.r-project.org')

p_load("tidyverse", "magrittr")

set.seed(5318)

#---- Specify source files ----
source("sex-dementia_sim_parA.R")
source("sex-dementia_sim_script.R")
source("life_table2014.R")

#---- Search for Baseline Hazard ----
#Function that creates the dataset for the search
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
      ages[, j] = seq(from = 1, to = num_obs, by = 1) #Creates column of ids
    } else if(j == 2){
      ages[, j] = age0 #Creates column of baseline ages
    } else ages[, j] = ages[, (j-1)] + int_time #Creates ages at following timepoints
  }
  colnames(ages) <- age_varnames
  
  #---- Generating centered age data ----
  #Creating centered ages at each timepoint j
  c_ages <- as_tibble(ages - mean(age0)) %>% 
    mutate("id" = seq(from = 1, to = num_obs, by = 1)) #Creates column of ids
  colnames(c_ages) <- agec_varnames
  
  #---- Generating "true" cognitive function Cij ----
  #Cij = b00 + z0i + bo1*sexi + b02*age0i + b03*Ui + (b10 + z1i + b11*sexi + 
  #b12*age0i + b13*Ui)*timej + epsilonij
  
  #---- Generating random terms for slope and intercept ----
  #Generate random terms for each individual
  slope_int_noise <- as_tibble(mvrnorm(n = num_obs, mu = rep(0, 2), 
                                       Sigma = slope_int_cov)) %>% 
    cbind("id" = seq(from = 1, to = num_obs, by = 1), .) #Creates column of ids
  colnames(slope_int_noise) <- c("id", "z0i", "z1i")
  
  #---- Generating noise term (unexplained variance in Cij) for each visit ----
  sd_eps <- sqrt(var3)
  eps <- as_tibble(replicate(num_tests + 1, 
                             rnorm(n = num_obs, mean = 0, sd = sd_eps))) %>% 
    cbind("id" = seq(from = 1, to = num_obs, by = 1), .) #Creates column of ids
  colnames(eps) <- eps_varnames
  
  #---- Creating complete matrix of observation data ----
  obs <- left_join(obs, ages, by = "id") %>% left_join(c_ages, by = "id") %>%
    left_join(slope_int_noise, by = "id") %>% left_join(eps, by = "id")
  
  #---- Calculating Cij for each individual ----
  #Store Cij values
  Cij <- as.data.frame(cog_func(obs)$Cij) %>% 
    cbind("id" = seq(from = 1, to = num_obs, by = 1), .) #Creating column of ids
  colnames(Cij) <- Cij_varnames
  
  #Store slope values per interval per individual
  slopeij <- as.data.frame(cog_func(obs)$slopes) %>% 
    cbind("id" = seq(from = 1, to = num_obs, by = 1), .) #Creating column of ids
  colnames(slopeij) <- slopeij_varnames
  
  #---- Generating uniform random variables per interval for Sij ----
  USij <- as_tibble(replicate(num_tests, 
                              runif(num_obs, min = 0, max = 1))) %>%
    cbind("id" = seq(from = 1, to = num_obs, by = 1), .) #Creating column of ids
  colnames(USij) <- USij_varnames
  
  #---- Merging Cij, slopeij, and USij with observation data ----
  #Used as input for survival function
  obs <- left_join(obs, Cij, by = "id") %>% left_join(slopeij, by = "id") %>% 
    left_join(USij, by = "id")

  return("obs" = obs)
}

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
  live_bysex <-data_frame("male" = rep(NA, num_tests), 
                          "female" = rep(NA, num_tests))
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
      lambdas[j] = optimise(survivors, interval = c(0, 0.005), 
                            obs = sim_data, cp = life_prob)$minimum
      Sij <- survival(obs = sim_data, lambda = lambdas)
      alive_now <- (Sij[[j]] > 5)*1
      sim_data %<>% cbind(., "alive" = (Sij[[j]] > 5)*1)
      cp_alive[j] = mean(alive_now)
      alive_bysex = sim_data %>% group_by(sex) %>% 
        summarise_at("alive", mean)
      live_bysex[j, ] = alive_bysex$alive
    } else {
      sim_data %<>% filter(alive == 1)
      lambdas[j] = optimise(survivors, 
                            interval = c(lambdas[j - 1], 1.025*lambdas[j - 1]), 
                            obs = sim_data, cp = life_prob)$minimum
      Sij <- survival(obs = sim_data, lambda = lambdas)
      alive_now <- (Sij[[j]] > 5)*1
      sim_data %<>% mutate("alive" = alive_now)
      cp_alive[j] = mean(alive_now)
      alive_bysex = sim_data %>% group_by(sex) %>% 
        summarise_at("alive", mean)
      live_bysex[j, ] = alive_bysex$alive
    }
  }
  return(list("lambdas" = lambdas, "cp_alive" = cp_alive, 
              "live_bysex" = live_bysex))
}

#---- Averaging over baseline hazard searches----
find_lambda <- function(life_table){
  simdata <- data_gen()
  search <- lambdas(sim_data = simdata, cp = life_table)
  return(search)
}

#---- Check conditional probabilities using baseline hazards ----
#Make sure to rerun parameter file with desired baseline hazards before running 
#actual simulation
lambda_searches <- replicate(5, find_lambda(life_table = female_life))

avg_lambdas <- as_tibble(do.call(rbind, lambda_searches["lambdas", ])) %>%
  colMeans()
avg_cps <- as_tibble(do.call(rbind, lambda_searches["cp_alive", ])) %>%
  colMeans()
avg_cps_bysex <- as_tibble(matrix(unlist(lambda_searches["live_bysex", ]), 
                           ncol = 5, byrow = FALSE)) %>% rowMeans() %>% 
  matrix(., ncol = 2, byrow = FALSE)
colnames(avg_cps_bysex) <- c("Males", "Females")

#---- Search for dementia cut-point ----

#Looking for a reasonable dementia cut point
#Use 4 simulated datasets and find quantiles of baseline Cij
dem_cut <- replicate(4, sex_dem_sim()) %>% map("Ci0") %>% unlist() %>% 
  quantile(0.05)
