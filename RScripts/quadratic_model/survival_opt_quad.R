#***************************************************************
# Performs a search for baseline hazards 
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
  "RScripts", "quadratic_model",
  "sex-dementia_sim_parA_onedemcut_nodemkill_maleAD_quad.R"))
source(here("RScripts", "euro_life_tables.R"))
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
  
  return(obs)
}

#---- Search for Baseline Hazard ----
#Function for finding lambdas
lambdas <- function(sim_data, cp50_unexposed){
  
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
                            interval = c(lambdas[j - 1], 2.75*lambdas[j - 1]), 
                            obs = sim_data, cp = life_prob)$minimum
      Sij <- survival(obs = sim_data, lambda = lambdas)
      alive_now <- (Sij[[j]] > 5)*1
      sim_data %<>% mutate("alive" = alive_now)
      cp_alive[j] = mean(alive_now)
    }
  }
  return(list("lambdas" = lambdas, "cp_alive" = cp_alive))
}

#---- Search for effect of "female" on baseline hazard ----
g1s <- function(sim_data, cp50_exposed){
  
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


