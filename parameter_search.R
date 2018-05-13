#***************************************************************
# Performs a parameter search using a specified function and 
# specified file containing reference data
#***************************************************************

#---- Package Loading and Options ----
if (!require("pacman")) 
  install.packages("pacman", repos='http://cran.us.r-project.org')

p_load("tidyverse", "magrittr")

#---- Specify source files ----
source("sex-dementia_sim_parA.R")
source("sex-dementia_sim_script.R")
source("life_table2014.R")
source("dementia_incidence2000-2013.R")

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

#---- Search for dementia cut-point ----
demcut <- function(sim_data, demrates){
  #Function we are trying to optimize
  dementia <- function(demcut, cogfunc, rate){
    cases = (cogfunc < demcut)*1
    cases_1000py = 1000*sum(cases)/(int_time*length(cases))
    return(abs(cases_1000py - rate))
  }
  #Dataframes to return
  dem_cuts <- vector(length = num_tests + 1)
  dem_1000py <- vector(length = num_tests + 1)
  #Begin search
  for(j in 1:length(dem_cuts)){
    test_num = j - 1
    age <- paste("age", test_num, sep = "")
    C <- paste("Ci", test_num, sep = "")
    index <- which(demrates$LowAge <= sim_data[1, age] & 
                     demrates$HighAge > sim_data[1, age])
    if(length(index) > 0){
      alive_now <- filter(sim_data, !is.na(sim_data[, C])) %>% dplyr::select(C)
      match_dem <- as.double(demrates[index, "Rate"])
      dem_cuts[j] = optimise(dementia, 
                             interval = c(min(alive_now), max(alive_now)), 
                             cogfunc = alive_now, rate = match_dem)$minimum
      
      dem_1000py[j] = 
        1000*sum((alive_now < dem_cuts[j])*1)/(int_time*nrow(alive_now))
    }
  }
  zeros <- which(dem_cuts == 0)
  dem_cuts[zeros] = dem_cuts[max(zeros) + 1]
  return(list("dem_cuts" = dem_cuts, "dem_1000py" = dem_1000py))
}

find_demcut <- function(dem_table){
  simdata <- sex_dem_sim()
  age_Cijs <- simdata$obs %>% 
    dplyr::select(c(dput(age_varnames[-1]), dput(Cij_varnames[-1])))
  search <- demcut(sim_data = age_Cijs, demrates = dem_table)
  return(search)
}

#This occassionally gives a weird error that I haven't figured out yet
#If it breaks during the search, rerun it
demcut_searches <- replicate(5, find_demcut(dem_rates_whites))

dem_cut <- as_tibble(do.call(rbind, demcut_searches["dem_cuts", ])) %>% 
  colMeans() 
avg_matches <- as_tibble(do.call(rbind, demcut_searches["dem_1000py", ])) %>%
  colMeans()


  
  